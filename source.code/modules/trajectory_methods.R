suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(monocle3))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))


#----------------------------------------------------------------
#' runCellphoneDB

#' Cell-cell communication using CellphoneDB.
#' @param obj Seurat object.
#' @param prefix Prefix name of output files, default: out.
#' @param out.figs.dir Output directory, default: ./.
#' @return NULL
#' @export

runCellphoneDB <- function(obj, prefix = "out", out.figs.dir = "./", ...) {
    GetAssayData(obj) %>% as.matrix() -> anno.exprs
    write.table(anno.exprs, ".Seurat_anno.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

    meta.infos <- as.data.frame(Idents(obj))
    write.table(meta.infos, ".meta_infos.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    system("cellphonedb method statistical-analysis .meta_infos.txt .Seurat_anno.txt --threads=20 --counts-data=hgnc_symbol --threshold=0.25 --pvalue=0.05")
    system("cellphonedb plot dot-plot")
    system("cellphonedb plot heatmap-plot .meta_infos.txt")
    system(sprintf("mv out %s", file.path(out.figs.dir, prefix)))
}

# monocle3Analysis

#' Pesudotime analysis using Monocle3
#' @param obj Seurat object.
#' @param root.cells A vector of root cells.
#' @param color.by Render color using pseudotime or others, default: pseudotime.
#' @return A list object, including Seurat, cds, and ggplot2 objects.
#' @export

monocle3Analysis <- function(obj, root.cells, color.by = "pseudotime", ...) {
    obj.cds <- as.cell_data_set(obj)
    obj.cds <- cluster_cells(cds = obj.cds, reduction_method = "UMAP")
    obj.cds <- learn_graph(obj.cds, use_partition = TRUE)
    obj.cds <- order_cells(obj.cds, reduction_method = "UMAP", root_cells = root.cells)

    pt.cell <- plot_cells(
        cds = obj.cds,
        color_cells_by = color.by,
        show_trajectory_graph = TRUE,
        cell_size = 0.5,
        ...
    )
    obj <- AddMetaData(
        object = obj,
        metadata = obj.cds@principal_graph_aux@listData$UMAP$pseudotime,
        col.name = color.by
    )
    gp <- FeaturePlot(obj, color.by, pt.size = 0.5) & scale_color_viridis_c()
    return(list(obj = obj, gp = gp, cds = obj.cds, pt.cell = pt.cell))
}

diffusionMapScanpy <- function(obj) {
    DefaultAssay(obj) <- "RNA"
    obj <- obj %>% FindVariableFeatures()
    obj <- ScaleData(object = obj, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"))
    de.markers <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
    top.degs <- de.markers %>%
        group_by(cluster) %>%
        top_n(n = 50, wt = avg_log2FC)
    obj <- RunPCA(obj, features = top.degs, npcs = 50)
    obj <- saveAsH5adFile(obj, "test", "./")
    import("modules/")
}

fitModelPseud <- function(obj, root.cells, ncores = 40) {
    expr_matrix <- Seurat::GetAssayData(obj, assay = "RNA", slot = "count")
    my_cds <- monocle3::new_cell_data_set(
        expr_matrix,
        cell_metadata = obj@meta.data
    )
    my_cds <- preprocess_cds(my_cds)
    my_cds <- align_cds(my_cds, alignment_group = "orig.ident")
    my_cds <- reduce_dimension(my_cds)
    my_cds <- cluster_cells(cds = my_cds, reduction_method = "UMAP")
    my_cds <- learn_graph(my_cds, use_partition = TRUE)
    my_cds <- order_cells(my_cds, root_cells = root.cells)

    gene_fits <- monocle3::fit_models(my_cds, model_formula_str = "~pseudotime", core = ncores)
    fit_coefs <- coefficient_table(gene_fits)
    fit.res <- fit_coefs %>%
        filter(term == "pseudotime" & q_value <= 0.05) %>%
        select(gene_id, term, q_value, estimate)
    return(fit.res)
}

#' createMonocleObject

#' Create single-cell experiment object for monocle2 analysis.
#' @param obj.seu Seurat object.
#' @return single-cell experiment object.
#' @export

createMonocleObject <- function(obj.seu) {
    sce <- obj.seu
    sample.ann <- sce@meta.data
    gene.ann <- data.frame(
        gene_short_name = rownames(sce@assays$RNA),
        row.names = rownames(sce@assays$RNA)
    )
    pd <- new("AnnotatedDataFrame", data = sample.ann)
    fd <- new("AnnotatedDataFrame", data = gene.ann)
    ct <- as.data.frame(sce@assays$RNA@counts)
    sc.cds <- newCellDataSet(
        as.matrix(ct),
        phenoData = pd,
        featureData = fd,
        expressionFamily = negbinomial.size(),
        lowerDetectionLimit = 1
    )
    return(sc.cds)
}

#' cellChatAnlysis

#' @param obj Seurat object.
#' @param group.by Annotation cell type names, default: Anno.Level.Fig.1.
#' @param db.chat CellChatDB for analysis, default: CellChatDB.human.
#' @return cellchat object.
#' @export

cellChatAnlysis <- function(obj, group.by = "Anno.Level.Fig.1", db.chat = CellChatDB.human) {
    data.input <- GetAssayData(obj, slot = "data")
    meta.data <- obj@meta.data
    cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = group.by)
    cellchat <- addMeta(cellchat, meta = meta.data)
    cellchat <- setIdent(cellchat, ident.use = group.by)
    cellchat@DB <- db.chat

    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 10)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 2)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    return(cellchat)
}
