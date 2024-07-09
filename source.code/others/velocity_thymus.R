require(dplyr)
require(reticulate)

#-----------------------------------------------------------------
# Load own modules

source("modules/utils.R")
source("modules/global_params.R")
source("modules/seurat_methods.R")
source("modules/visualization.R")
source("modules/findMedullaClusters.R")
source("modules/trajectory_methods.R")

#-----------------------------------------------------------------
# Read single-cell RNA-seq data and save as H5AD format file.

obj.merged.sc <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')
obj.xx <- LoadH5Seurat(".merged.obj.UMAP.Scanpy.h5seurat")

out.data.dir <- file.path("../3.results", "2.fine_anno/data")

uamp.new <- obj.xx[["umap"]]@cell.embeddings %>% `colnames<-`(c("umapn_1", "umapn_2"))
uamp.new.sub <- uamp.new[intersect(colnames(obj.merged.sc), rownames(uamp.new)), ]
obj.merged.sc <- obj.merged.sc[, rownames(uamp.new.sub)]
obj.merged.sc[["umapn"]] <- CreateDimReducObject(embeddings = uamp.new.sub, key = "UMAPN_", assay = DefaultAssay(obj.merged.sc))

obj.sub <- subset(obj.merged.sc, (Anno.Level.11 %in% ANNO_TCR_IDENT) & tcr == 1)
obj.sub$CellTypeNew <- obj.sub$Anno.Level.11
saveAsH5adFile(obj.sub, prefix = "obj.T", out.data.dir = "./")

loom.files <- c(
    "F21TM0625G-D10",
    "TM0820G-C1",
    "TM0820G-C4",
    "TM0902G-F10",
    "TM0902G-F6",
    "F13TM0517G-D3",
    "F23TM0528G-A9",
    "TM0820G-C2",
    "TM0820G-C5",
    "TM0902G-F4",
    "TM0923G-B7",
    "F17TM0517G-D5",
    "TM0804G-F5",
    "TM0820G-C3",
    "TM0820G-C6",
    "TM0902G-F5"
) %>% file.path("/public/home/glht01/projects/Thymus/transfer/XQ/thymus_results_11.29.2020//", ., "velocyto", paste0(., ".loom"))

#------------------------------------------------------------------
# Read loom file created using velocity.

scv <- import("scvelo")
sc <- import("scanpy")
np <- import("numpy")
loompy <- import("loompy")

# loompy$combine(loom.files, output_file = file.path('./', 'thymus_combined.loom'))

# adata.vec <- sc$read(loom.files[1], sparse= TRUE, cache = TRUE)
adata.vec <- sc$read(file.path("./", "thymus_combined.loom"), sparse = TRUE, cache = TRUE)
# adata.sc <- sc$read_h5ad(file.path('./', 'obj.blast.h5ad'))
adata.sc <- sc$read_h5ad(file.path("./", "obj.T.h5ad"))

rownames(adata.vec$obs) <- adata.vec$obs %>%
    rownames() %>%
    gsub(":", "_", .) %>%
    gsub("x", "-1", .)
adata.vec$obs$index <- adata.vec$obs %>%
    rownames() %>%
    gsub(":", "_", .) %>%
    gsub("x", "-1", .)
adata.vec$var_names_make_unique()
adata <- scv$utils$merge(adata.sc, adata.vec)

scv$pp$filter_and_normalize(adata, min_shared_counts = 30, n_top_genes = 2000)
# scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
scv$tl$velocity_graph(adata)

cols <- ANNO_TREND_IDENT_COLOR[ANNO_TCR_IDENT]
scv$pl$velocity_embedding_stream(adata, basis = "umap", color = "CellTypeNew", colors = cols)

scv$pl$velocity_graph(adata, threshold = .1, color = "Anno.Level.Fig.1")
scv$tl$latent_time(adata)
scv$pl$scatter(adata, color = "velocity_pseudotime", cmap = "gnuplot")

adata$uns[["neighbors"]][["distances"]] <- adata$obsp[["distances"]]
adata$uns[["neighbors"]][["connectivities"]] <- adata$obsp[["connectivities"]]
scv$tl$paga(adata, groups = "Anno.Level.Fig.1")
scv$pl$paga(adata, basis = "umap", size = 50, alpha = .1, min_edge_width = 2, node_size_scale = 1.5)

topgenes <- adata$var["fit_likelihood"]
topgenes_vals <- topgenes[, 1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing = TRUE)
scv$pl$heatmap(adata, var_names = names(topgenes_vals)[1:300], sortby = "latent_time", col_color = "Anno.Level.Fig.1")
