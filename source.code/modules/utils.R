suppressMessages(library(reticulate))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Matrix))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(monocle3))
suppressMessages(library(patchwork))

#-------------------------------------------------------------------------------
# Load own modules

use_python("/public/home/glht01/softwares/Python-3.8.5/python", required = TRUE)
source_python("modules/scanpy_methods.py")

#------------------------------------------------------------------------------

#' mergeMultipleSeuratObjects

#' Merge multiple single cell datasets and as a seurat object
#' @param input_folders Input fold of 10X data.
#' @param do.normalize Normalize the count data or not, default: FALSE.
#' @return Merged seurat object.
#' @export

mergeMultipleSeuratObjects <- function(input_folders, do.normalize = FALSE, ...) {
    seurat_data <- purrr::map(input_folders, Read10X)
    add_sample_name_to_cell <- function(x, y) {
        colnames(x) <- paste(y, colnames(x), sep = "_")
        return(x)
    }

    sample_names <- strsplit(input_folders, "/") %>%
        as.data.frame() %>%
        t() %>%
        .[, dim(.)[2] - 1] %>%
        as.character()
    seurat_data <- purrr::map2(seurat_data, sample_names, add_sample_name_to_cell)
    seurat_objects <- purrr::map2(
        seurat_data,
        sample_names, function(x, y) {
              CreateSeuratObject(count = x, project = y, min.cells = 5, ...)
          }
    )

    merged_seurat <- purrr::reduce(
        seurat_objects,
        function(x, y) {
            merge(x, y, do.normalize = do.normalize)
        }
    )
}

#-------------------------------------------------------------------------------
# Extract cellcycle genes

extractCellCycleGenes <- function(SeuratObj, file.path = "./configs/regev_lab_cell_cycle_genes.txt") {
    cc.genes <- readLines(con = file.path)
    s.genes <- cc.genes[1:43]
    s.genes <- intersect(s.genes, rownames(SeuratObj))
    g2m.genes <- cc.genes[44:97]
    g2m.genes <- intersect(g2m.genes, rownames(SeuratObj))
    return(list(s.genes = s.genes, g2m.genes = g2m.genes))
}

#--------------------------------------------------------------------------------
# Averaged the expression data of each sample

CorrBasedOnAvgExprOfSamples <- function(scrnaSeuratObj) {
    SplitObject.lst <- SplitObject(scrnaSeuratObj, split.by = "source_sn")
    sn.corrs <- lapply(SplitObject.lst, function(x) {
        GetAssayData(x) %>%
            as.matrix() %>%
            rowMeans()
    }) %>%
        as.data.frame() %>%
        cor()
    return(sn.corrs)
}

#---------------------------------------------------------------
# Shannon extropy

shannonEntropy <- function(X) {
    nrep.sns <- dim(X)[2]
    X.norm <- X / rowSums(X)
    X.mu <- rowMeans(X.norm)
    shannon.dist <- 1 / nrep.sns * rowSums(X.norm / X.mu * log2(X.norm / X.mu))
    return(shannon.dist)
}

#--------------------------------------------------------------
# Read reference datasets

readScience <- function(science.expr.file, science.markers.file, ...) {
    ref.data <- read.table(science.expr.file, sep = "\t", header = TRUE) %>% .[!duplicated(.[, 1]), ]
    rownames(ref.data) <- ref.data[, 1]
    ref.data <- ref.data[, -1]

    marker.data <- read.table(science.markers.file, sep = ",", header = TRUE) %>% unlist(use.names = FALSE)
    ref.data[marker.data, ]
}

#-------------------------------------------------------------
# Correlation between Our data and Science data

corrBewtweenOursAndScience <- function(ref.data, merged.obj) {
    common.genes <- intersect(rownames(ref.data), rownames(merged.obj))
    average.exp <- AverageExpression(merged.obj, features = common.genes, assays = "RNA")
    ref.data <- ref.data[common.genes, ]
    corr.ref <- apply(average.exp$RNA, 2, function(x) {
        apply(ref.data, 2, function(y) {
            cor(log2(1 + x), y)
        })
    })
    return(corr.ref)
}

#-----------------------------------------------------------------------
#' mapAgesStatus

#' Map age status, including prental, children, and adult.
#' @param obj Seurat object.
#' @param SampleClassify Sample information from global setrings.
#' @return Seurat object.
#' @export

mapAgesStatus <- function(obj, SampleClassify) {
    obj@meta.data$AgeStatus <- "N"
    obj@meta.data$AgeStatus[obj@meta.data$orig.ident %in% SampleClassify[["Prental"]]] <- "Prental"
    obj@meta.data$AgeStatus[obj@meta.data$orig.ident %in% SampleClassify[["Children"]]] <- "Children"
    obj@meta.data$AgeStatus[obj@meta.data$orig.ident %in% SampleClassify[["Adult"]]] <- "Adult"
    obj@meta.data$AgeStatus <- factor(obj@meta.data$AgeStatus, levels = c("Prental", "Children", "Adult"))
    return(obj)
}


#' downSamplSeurat

#' Down-sampling for Seurat obeject.
#' @param obj Seurat object
#' @param cnt Sample size for each ident, default: 2000.
#' @param seed Randon seed, default: 123.
#' @return Subset of seurat object.
#' @export

downSamplSeurat <- function(obj, cnt = 2000, seed = 123, percent = NULL) {
    set.seed(seed)
    cells <- Idents(obj) %>% table()
    sub.cells <- sapply(names(cells), function(xx) {
        sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names()
        cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
        if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
        return(sub.cells)
    }) %>% unlist(use.names = F)
    subset(obj, cells = sub.cells)
}

#' saveAsH5adFile

#' Convert seurat object to H5AD object.
#' @param obj Seurat object.
#' @param prefix Prefix name of output file.
#' @param out.data.dir Output directory, default: ./.
#' @return Converted H5AD file.
#' @export

saveAsH5adFile <- function(obj, prefix, out.data.dir = "./") {
    file.name <- file.path(out.data.dir, prefix)
    SaveH5Seurat(obj, filename = sprintf("%s.h5Seurat", file.name), overwrite = T)
    Convert(sprintf("%s.h5Seurat", file.name), dest = "h5ad", overwrite = T)
    return(sprintf("%s.h5ad", file.name))
}

#' readH5ADFile

#' Read H5AD data format file.
#' @param file.name  Specify the path of H5AD data format file.
#' @return Seurat object.
#' @export

readH5ADFile <- function(file.name) {
    if ((grep("h5ad", file.name) %>% length()) < 1) file.name <- paste0(file.name, ".h5ad")
    prefix <- gsub(".h5ad", "", file.name)
    Convert(file.name, dest = "h5seurat", overwrite = TRUE)
    obj <- LoadH5Seurat(sprintf("%s.h5seurat", prefix))
    return(obj)
}

#' distBetweenDfs

#' Calculate the distances between two data.frames
#' @param df1 Gene expression data frame.
#' @param df2 Gene expression data frame.
#' @return dist.mat
#' @export

distBetweenDfs <- function(df1, df2, method = "bhattacharyya") {
    require(philentropy)
    ovp.genes <- intersect(rownames(df1), rownames(df2))
    df1.sub <- df1[ovp.genes, ]
    df2.sub <- df2[ovp.genes, ]
    dist.mat <- apply(df1.sub, 2, function(x) {
        dist.val <- apply(df2.sub, 2, function(y) {
            x <- x / sum(x)
            y <- y / sum(y)
            dat <- rbind.data.frame(x, y)
            dist.diversity(dat, p = 2) %>% .[method]
        })
    })
    return(dist.mat)
}


#' readSPMergedData

#' Read merged spatial sequencing data.
#' @return obj.st
#' @export

readSPMergedData <- function() {
    obj.sc.st <- readRDS("../transfer/XQ/data_sc_st/merged.obj_sc_st.rds")
    obj.st <- subset(
        obj.sc.st,
        orig.ident == "stT1" |
            orig.ident == "stT2" |
            orig.ident == "stT3" |
            orig.ident == "stT4" |
            orig.ident == "stT5" |
            orig.ident == "stT6" |
            orig.ident == "stT7" |
            orig.ident == "stT8"
    )
    DefaultAssay(obj.st) <- "RNA"
    return(obj.st)
}

upgma <- function(D, method = "average", ...) {
    DD <- as.dist(D)
    hc <- hclust(DD, method = method, ...)
    result <- as.phylo(hc)
    result <- reorder(result, "postorder")
    result
}

rotate_x <- function(data, labels_vec, rot_angle, cex = 1, ...) {
    plt <- barplot(data, col = "steelblue", xaxt = "n", ...)
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1, 1.1), xpd = TRUE, cex = cex)
    return(plt)
}

mappCellNames <- function(obj, src.names, to.names, meta.name, col.name) {
    xx <- obj@meta.data[, meta.name] %>% as.vector()
    for (idx in 1:length(src.names)) {
        xx[xx == src.names[idx]] <- paste0(to.names[idx], ".1")
    }
    xx <- gsub("\\.1", "", xx)
    if (col.name %in% colnames(obj@meta.data)) {
        idxes <- which(obj@meta.data[, meta.name] %in% src.names)
        obj@meta.data[idxes, col.name] <- xx[idxes]
    } else {
        obj@meta.data[, col.name] <- xx
    }
    return(obj)
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return(intersection / union)
}

#' goAnno

#' Gene ontology analysis for DEGs identified by FindAllMarkers.
#' @param de.markers A data frame of DEGs identified by FindAllMarkers.
#' @param topn Top N genes used for gene ontology analysis, default: 100.
#' @return A list of annotated results.
#' @export

goAnno <- function(de.markers, topn = 100) {
    require(org.Hs.eg.db)
    ego.lst <- lapply(de.markers$cluster %>% unique() %>% as.vector(), function(xx) {
        up.genes <- subset(de.markers, cluster == xx)$gene[1:topn]
        up.genes.map <- clusterProfiler::bitr(up.genes, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
        ego.up.1 <- clusterProfiler::enrichGO(
            gene = up.genes.map$ENTREZID %>% unique(),
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        )
    }) %>% `names<-`(de.markers$cluster %>% unique() %>% as.vector())
    return(ego.lst)
}

#' predictCtecOrMtec
# ‘ predict cTEC or mTEC for TEC cells.
#' @param obj.tec A seurat object containing TEC cells.
#' @return
#' @export

predictCtecOrMtec <- function(obj.qry.tec) {
    obj.ref.tec <- readRDS(".obj.ref.tec.rds")
    Idents(obj.ref.tec) <- obj.ref.tec$Anno_level_3
    ovp.genes <- intersect(rownames(obj.qry.tec), rownames(obj.ref.tec))
    obj.ref.tec <- obj.ref.tec[intersect(genes, ovp.genes), ]
    obj.tec <- obj.qry.tec[intersect(genes, ovp.genes), ]
    sample.size <- table(Idents(obj.ref.tec)) %>%
        min() %>%
        {
            . * 0.5
        } %>%
        as.integer()
    lr.model <- trainLRmodel(obj.ref.tec, c("cTEC", "mTEC"), sample.size = sample.size, prefix = "tec", out.figs.dir = out.figs.dir)
    pred.res <- lr.model$predict(obj.tec %>% GetAssayData(.) %>% t()) %>% `names<-`(colnames(obj.tec))
}

#' smoothSpatialData

#' Smooth spatial-spot data.
#' @param obj Seurat object.
#' @param features A vector of features.
#' @return A data.frame of smoothed data.
#' @export

smoothSpatialData <- function(obj, features, method = "ks") {
    res <- lapply(features, function(cell) {
        w <- obj@meta.data[, cell]
        x <- Seurat::GetTissueCoordinates(object = obj@images[[1]])
        w <- Nebulosa:::calculate_density(w, x, method = method, adjust = 1) %>% scale()
    }) %>%
        do.call(cbind, .) %>%
        `colnames<-`(features)
    return(res)
}

#' ovpMarkers

#' Overlap marker genes between query and reference.
#' @param de.markers.qry A data.frame of differential markers identified FindAllMarkers, Query.
#' @param de.markers.ref A data.frame of differential markers identified FindAllMarkers, Ref.
#' @param top.n Top number of genes, default: 50.
#' @return A list of overlap markers.
#' @export

ovpMarkers <- function(de.markers.qry, de.markers.ref, top.n = 50) {
    topn.qry <- de.markers.qry %>%
        group_by(cluster) %>%
        top_n(n = top.n, wt = avg_log2FC)
    topn.ref <- de.markers.ref %>%
        group_by(cluster) %>%
        top_n(n = top.n, wt = avg_log2FC)
    qry.lst <- split(topn.qry$gene, topn.qry$cluster)
    qry.ref <- split(topn.ref$gene, topn.ref$cluster)
    lapply(qry.lst, function(xx) {
        lapply(qry.ref, function(yy) {
            intersect(xx, yy)
        })
    }) -> res.ovp
    return(res.ovp)
}

#' searchCentricSpotForCARD

#' Search centric spot for CARD segmentation results.
#' @param sc.map.lst A list of CARD segmentation results.
#' @param obj.sp.lst A list of spatial-spots results.
#' @return
#' @export

euclidean <- function(a, b) sqrt(sum((a - b)^2))
searchCentricSpotForCARD <- function(sc.map.lst, obj.sp.lst) {
    res <- lapply(names(sc.map.lst), function(sn.name) {
        obj.1 <- sc.map.lst[[sn.name]]
        obj.2 <- obj.sp.lst[[sn.name]]
        coord.df <- GetTissueCoordinates(obj.2)
        card.coord <- as.data.frame(colData(obj.1))

        dist.df <- parallel::mclapply(obj.2$Assign.Centric %>% unique(), function(centric.spot) {
            ref.coord <- coord.df[centric.spot, ]
            qry.coord <- card.coord[, c("x", "y")]
            dist.vec <- apply(qry.coord, 1, function(x.df) euclidean(x.df, ref.coord[1, ]))
            dist.vec
        }, mc.cores = 10) %>%
            do.call(cbind, .) %>%
            `colnames<-`(obj.2$Assign.Centric %>% unique())
        assign.spots <- apply(dist.df, 1, which.min) %>% colnames(dist.df)[.]
        assign.spots
    }) %>% `names<-`(names(sc.map.lst))
    return(res)
}

#' corrBewtweenObjs

#' Correlations between Seurat objects
#' @param obj.qry Seurat object of query.
#' @param obj.ref Seurat object of reference.
#' @param anno.vec A vector of annotation results.
#' @return A matrix of correlations.
#' @export

corrBewtweenObjs <- function(obj.qry, obj.ref) {
    DefaultAssay(obj.qry) <- "RNA"
    DefaultAssay(obj.ref) <- "RNA"
    common.genes <- intersect(rownames(obj.ref), rownames(obj.qry))
    exp.ref <- AverageExpression(obj.ref, features = common.genes, assays = "RNA")$RNA
    exp.qry <- AverageExpression(obj.qry, features = common.genes, assays = "RNA")$RNA

    ref.data <- exp.ref[common.genes, ]
    qry.data <- exp.qry[common.genes, ]

    corr.ref <- apply(ref.data, 2, function(x) {
        apply(qry.data, 2, function(y) {
            cor(x, y)
        })
    })
    return(corr.ref)
}

#' Re-format Fb-cyling data.
#' @param

reformatThymoma <- function(obj.merged.thymoma) {
    stromal <- readRDS(".stromal.cells.rds")
    ovp.cells <- stromal[intersect(names(stromal), colnames(obj.merged.thymoma))]
    rm.cells <- ovp.cells[ovp.cells == "REMOVE"]
    xx <- subset(obj.merged.thymoma, cells = rm.cells %>% names())
    xx.expr <- FetchData(xx, vars = c("FCN1", "MKI67", "TYMS"))
    cycling.cells <- xx.expr[rowSums(xx.expr > 0) > 0, ] %>% rownames()
    ovp.cells[cycling.cells] <- "Fb_cycling"
    obj.merged.thymoma <- updateAnnoRes(ovp.cells, obj.merged.thymoma, "Anno.Level.Fig.1", "Anno.Level.Fig.1")
    obj.merged.thymoma <- subset(obj.merged.thymoma, Anno.Level.Fig.1 != "REMOVE")
    obj.merged.thymoma
}

#' readLoadScience

readLoadScience <- function() {
    obj.science <- LoadH5Seurat("../4.extdata/HTA08.v01.A05.Science_human_figX.h5Seurat")
    obj.science$AgeStatus <- "Prental"
    obj.science$AgeStatus[obj.science$Age %in% c("3m", "6m", "10m", "15m", "30m", "13y")] <- "Children"
    obj.science$AgeStatus[obj.science$Age %in% c("24y", "35y")] <- "Adult"
    return(obj.science)
}

reformatCARD <- function(MapCellCords, xx.obj) {
    meta.info <- xx.obj@meta.data
    if ("SpotName" %in% colnames(MapCellCords)) {
        xx.df <- xx.obj@meta.data[MapCellCords$SpotName, "Assign.Centric"]
        MapCellCords$Assign.Centric <- xx.df
    } else {
        xx.cord <- GetTissueCoordinates(xx.obj)
        xx.centerSPOT <- paste0(xx.cord[, 1], "x", xx.cord[, 2])
        assign.cen <- meta.info[match(MapCellCords$centerSPOT, xx.centerSPOT), "Assign.Centric"]
        spot.name <- rownames(meta.info)[match(MapCellCords$centerSPOT, xx.centerSPOT)]
        he.label <- meta.info[match(MapCellCords$centerSPOT, xx.centerSPOT), "HE.Labels"]
        MapCellCords$SpotName <- spot.name
        MapCellCords$HE.Labels <- he.label
        MapCellCords$Assign.Centric <- assign.cen
    }
    return(MapCellCords)
}

coExistIndex <- function(sc.map.lst) {
    qry.vec <- rep(0, length(ANNO_ENTIRE_IDNET_FIG1)) %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
    df.mat <- lapply(sc.map.lst, function(coord) {
        coord <- colData(coord)
        data.bak <- table(coord$CT, coord$centerSPOT) %>%
            t() %>%
            as.data.frame.matrix()
        names(data.bak)[which(names(data.bak) == "abT_entry")] <- "abT(entry)"
        diff.cells <- setdiff(names(qry.vec), colnames(data.bak))
        if (length(diff.cells) > 0) {
            data.bak[, diff.cells] <- 0
        }
        data.bak <- data.bak[, names(qry.vec)]
        data.bak[data.bak > 0] <- 1
        data.bak
    }) %>% do.call(rbind, .)

    index <- lapply(colnames(df.mat), function(xx) {
        ref <- df.mat[, xx]
        idxes <- lapply(colnames(df.mat), function(yy) {
            qry <- df.mat[, yy]
            tmp <- ref + qry
            idx <- sum(tmp == 2) / min(sum(ref > 0) + 1, sum(qry > 0) + 1)
            return(idx)
        }) %>%
            unlist() %>%
            as.vector() %>%
            `names<-`(colnames(df.mat))
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(colnames(df.mat))
}

statCMCTypes <- function(sc.map.lst, obj.sp.lst, group = "Ours") {
    qry.vec <- rep(0, length(ANNO_ENTIRE_IDNET_FIG1)) %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
    df.mat <- lapply(names(sc.map.lst), function(sn) {
        xx.obj <- obj.sp.lst[[sn]]
        coord <- colData(sc.map.lst[[sn]]) %>% reformatCARD(., xx.obj)
        coord$HE.Labels[coord$HE.Labels == "Medulla_lo"] <- "Cortex"
        coord$HE.Labels[coord$HE.Labels == "Medulla_hi"] <- "Medulla"
        coord$HE.Labels[coord$HE.Labels == "Medulla_centric"] <- "Medulla"
        coord
    }) %>% do.call(rbind, .)

    data.bak <- table(df.mat$CT, df.mat$HE.Labels) %>%
        t() %>%
        as.data.frame.matrix()
    names(data.bak)[which(names(data.bak) == "abT_entry")] <- "abT(entry)"
    he.spots <- lapply(obj.sp.lst, function(xx) {
        xx$HE.Labels[xx$HE.Labels %in% c("Medulla_hi", "Medulla_centric")] <- "Medulla"
        xx$HE.Labels[xx$HE.Labels %in% c("Medulla_lo")] <- "Cortex"
        table(xx$HE.Labels)
    }) %>%
        do.call(rbind, .) %>%
        colSums()
    data.bak <- sweep(data.bak, 2, he.spots, "/")
    diff.cells <- setdiff(names(qry.vec), colnames(data.bak))
    if (length(diff.cells) > 0) {
        data.bak[, diff.cells] <- 0
    }
    data.bak <- data.bak[, names(qry.vec)]
    data.bak
}

quantileNorm <- function(x) {
    Y <- as.matrix(x)
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
    return(Y)
}
