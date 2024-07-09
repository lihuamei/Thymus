library(Seurat)
library(CellTrek)
library(CARD)

require(future)
options(future.globals.maxSize = 500000 * 1024^2)
plan("multiprocess", workers = 20)

#' --------------------------------------------------------------
#' By Seurat

projBySeurat <- function(obj.sc, obj.st.lst) {
    ovp.genes <- intersect(rownames(obj.sc), rownames(obj.st.lst[[1]]))
    obj.sc <- obj.sc[ovp.genes, ]

    obj.sc <- NormalizeData(obj.sc)
    obj.sc <- FindVariableFeatures(obj.sc)
    obj.sc <- ScaleData(obj.sc, vars.to.regress = c("S.Score", "G2M.Score", "orig.ident", "percent.mt", "percent.ribo"))
    obj.ref <- obj.sc

    obj.srtct <- parallel::mclapply(obj.st.lst, function(obj.sp) {
        sn.name <- sp2SampleName(obj.sp@images %>% names())
        require(future)
        options(future.globals.maxSize = 500000 * 1024^2)
        plan("multiprocess", workers = 20)
        obj.sp <- obj.sp[ovp.genes, ]
        obj.sp <- NormalizeData(obj.sp)
        obj.sp <- FindVariableFeatures(obj.sp)
        obj.sp <- ScaleData(obj.sp)
        anchors <- FindTransferAnchors(reference = obj.ref, query = obj.sp)

        predictions <- TransferData(anchorset = anchors, refdata = obj.ref$Anno.Level.Fig.1)
        obj.sp <- AddMetaData(object = obj.sp, metadata = predictions)
    }, mc.cores = 8)

    saveRDS(obj.srtct, file = file.path(out.data.dir, "obj.SrtCT.rds"))
}

#--------------------------------------------------------------
#' By Celltrek

projByCellTrek <- function(obj.sc, obj.st.lst) {
    obj.celltrek <- parallel::mclapply(obj.st.lst, function(obj.sp) {
        require(future)
        options(future.globals.maxSize = 500000 * 1024^2)
        plan("multiprocess", workers = 20)
        brain_traint <- CellTrek::traint(st_data = obj.sp, sc_data = obj.sc, sc_assay = "RNA", cell_names = "Anno.Level.Fig.1")
        brain_celltrek <- CellTrek::celltrek(
            st_sc_int = brain_traint,
            int_assay = "traint",
            sc_data = obj.sc.sub,
            sc_assay = "RNA",
            reduction = "pca",
            intp = T,
            intp_pnt = 5000,
            intp_lin = F,
            nPCs = 30,
            ntree = 1000,
            dist_thresh = 0.55,
            top_spot = 5,
            spot_n = 5,
            repel_r = 20,
            repel_iter = 20,
            keep_model = T
        )$celltrek
        brain_celltrek$Anno.Level.Fig.1 <- factor(brain_celltrek$Anno.Level.Fig.1, levels = names(ANNO_ENTIRE_COLOR_FIG1))
        brain_celltrek
    }, mc.cores = 8)

    saveRDS(obj.celltrek, file = file.path(out.data.dir, "obj.celltrek.rds"))
}

#-----------------------------------------------------------------------
# By CARD

projByCARD <- function(obj.sc, obj.sp.lst, proj = TRUE) {
    require(CARD)
    CARD.obj.lst <- parallel::mclapply(obj.sp.lst, function(xx) {
        sp.counts <- GetAssayData(xx, slot = "count")
        sc.counts <- GetAssayData(obj.sc, slot = "count")
        sc.meta <- obj.sc@meta.data
        sc.meta$Anno.Level.Fig.1 <- as.vector(sc.meta$Anno.Level.Fig.1)
        sc.meta$Anno.Level.Fig.1[sc.meta$Anno.Level.Fig.1 == "abT(entry)"] <- "abT_entry"
        sp.loc <- GetTissueCoordinates(object = xx@images[[1]]) %>% `colnames<-`(c("x", "y"))
        sc.meta$sampleInfo <- "ALL"
        CARD.obj <- createCARDObject(
            sc_count = sc.counts,
            sc_meta = sc.meta,
            spatial_count = sp.counts,
            spatial_location = sp.loc,
            ct.varname = "Anno.Level.Fig.1",
            ct.select = unique(sc.meta$Anno.Level.Fig.1),
            sample.varname = "sampleInfo",
            minCountGene = 100,
            minCountSpot = 5
        )
        CARD.obj <- CARD_deconvolution(CARD_object = CARD.obj)
        CARD.obj <- CARD.imputation(CARD.obj, NumGrids = 2000, ineibor = 10, exclude = NULL)
        if (proj) CARD.obj <- CARD_SCMapping(CARD.obj, shapeSpot = "Square", numCell = 20, ncore = 1)
        CARD.obj
    }, mc.cores = getOption("mc.cores", 8L)) %>% `names<-`(names(obj.sp.lst))

    saveRDS(CARD.obj.lst, file = file.path(out.data.dir, "CARD.obj.lst.rds"))
}
