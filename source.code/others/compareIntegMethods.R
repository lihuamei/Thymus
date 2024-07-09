suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/trajectory_methods.R')

#-------------------------------------------------------------
# Load thymus single-cell RNA-seq data.

markers.save <- list()
merged.obj.1 <- readRDS(file.path('../3.results', '1.remove_doublet/data/merged.obj.clean.rds'))
DefaultAssay(merged.obj.1) <- "RNA"

SIG.REF <- read.table(SCIENCE_EXPR, header = T, sep = '\t') %>%
    .[!duplicated(.[, 1]), ] %>%
    `rownames<-`(.$index) %>%
    .[, -1] %>%
    .[read.table(SCIENCE_MARKERS, sep = ',', header = T) %>% as.matrix %>% as.vector %>% unique, ] %>%
    .[intersect(rownames(merged.obj.1), rownames(.)), ]

anno.title <- 'Level.0'
merged.obj.1 <- AddMetaData(merged.obj.1, metadata = Idents(merged.obj.1) , col.name = sprintf('Anno.%s', anno.title))

#-------------------------------------------------------------
# Update annotated cell types -> Level 1.

anno.title <- 'Level.1'
obj <- merged.obj.1
cell.names <- Idents(obj) %>% as.vector
cell.names[cell.names %in% levels(obj)[6 : 17]] <- 'T'
cell.names[cell.names == 'TEC'] <- 'Stromal'
cell.names <- factor(cell.names, levels = ANNO_SC_IDENT_LEVEL_1)
merged.obj.1 <- AddMetaData(merged.obj.1, metadata = cell.names %>% as.vector, col.name = sprintf('Anno.%s', anno.title))
Idents(obj) <- cell.names

obj$AgeStatusNew <- obj$AgeStatus %>% as.vector
obj$AgeStatusNew[obj$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

obj$AgeStatusNew <- factor(obj$AgeStatusNew, levels = c("Prental", "Children", "Adult", "Older"))

DimPlot(obj, label = TRUE) 
ggsave('revise.ms/reviewer.2/figs/UMAP.CCA.pdf', width = 6, height = 6)

DimPlot(obj, label = FALSE, raster = TRUE, group.by = 'AgeStatus', label.size = 6, pt.size=0.1)
ggsave('revise.ms/reviewer.2/figs/Anno.Split.by.AgeStatus.UMAP.pdf', width = 12, height = 6)

FeaturePlot(obj, features = c('HBB', 'HBG1', 'IGHG1', 'IGHG2', 'CD19', 'MS4A1', 'CST3', 'LYZ', 'IFITM3', 'ACTA2', 'CD3D', 'CD3E'), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.canonical.markers.expr.pdf', anno.title)), height = 8, width = 12)

gp <- DimPlot(obj, label = TRUE) + DimPlot(obj, label = TRUE, group.by = 'seurat_clusters') + DimPlot(obj, group.by = 'AgeStatusNew', cols = AGE_STATUS_EXTENT_COLOR[1 : 4])
ggsave('revise.ms/reviewer.2/figs/seurat.UMAP.pdf', width = 18, height = 5)

#-------------------------------------------------------------

obj$AgeStatusNew <- obj$AgeStatus %>% as.vector
obj$AgeStatusNew[obj$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

obj$AgeStatusNew <- factor(obj$AgeStatusNew, levels = c("Prental", "Children", "Adult", "Older"))

#-------------------------------------------------------------
# Run Honmony

require(harmony)
obj.harmon <- RunHarmony(obj, group.by.vars = 'orig.ident', plot_convergence = TRUE)
obj.harmon <- RunUMAP(obj.harmon, reduction = "harmony", dims = 1 : 30, seed.use = 8123456)
obj.harmon <- FindNeighbors(obj.harmon, reduction = "harmony", dims = 1 : 30) 

DefaultAssay(obj.harmon) <- 'integrated'
obj.harmon <- FindClusters(obj.harmon, resolution = 0.8)
obj.harmon$Anno.Level.1 <- cell.names

DefaultAssay(obj.harmon) <- 'RNA'
require(future)
options(future.globals.maxSize = 400000 * 1024^2)
plan("multiprocess", workers = 10)
de.markers <- FindAllMarkers(obj.harmon, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)
top20 <- de.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

anno <- c(
	'0' = 'T',
	'1' = 'T',
	'2' = 'T',
	'3' = 'T',
	'4' = 'T',
	'5' = 'T',
	'6' = 'T',
	'7' = 'T',
	'8' = 'T',
	'9' = 'T',
	'10' = 'T',
	'11' = 'B',
	'12' = 'T',
	'13' = 'T',
	'14' = 'T',
	'15' = 'B',
	'16' = 'T',
	'17' = 'T',
	'18' = 'Stromal',
	'19' = 'B',
	'20' = 'T',
	'21' = 'Myeloid',
	'22' = 'T',
	'23' = 'Stromal',
	'24' = 'Stromal',
	'25' = 'T',
	'26' = 'Plasma',
	'27' = 'B',
    '28' = 'T'
)	

obj.harmon <- RenameIdents(obj.harmon, anno)
saveRDS(obj.harmon, file = 'revise.ms/reviewer.2/data/obj.harmon.RDS')

obj.harmon <- readRDS('revise.ms/reviewer.2/data/obj.harmon.RDS')
cells <- subset(obj.harmon, seurat_clusters == '5')	%>% Cells
obj.xx <- subset(obj.harmon, cells = cells)
cells <- colnames(obj.xx)[FetchData(obj.xx, vars = c('HBG1')) > 0]	
obj.harmon$Harmony <- Idents(obj.harmon) %>% as.vector
obj.harmon$Harmony[cells] <- 'Ery'

obj.harmon$Harmony <- factor(obj.harmon$Harmony, levels = levels(obj))
obj.harmon$AgeStatusNew <- obj$AgeStatusNew
gp <- DimPlot(obj.harmon, label = TRUE, group.by = 'Harmony') + DimPlot(obj.harmon, label = TRUE, group.by = 'seurat_clusters') + DimPlot(obj.harmon, group.by = 'AgeStatusNew', cols = AGE_STATUS_EXTENT_COLOR[1 : 4]) 
ggsave('revise.ms/reviewer.2/figs/harmony.UMAP.pdf', width = 18, height = 5)

obj$harmony <- obj.harmon$Harmony 
obj$Seurat <- Idents(obj)
plot.df <- table(obj$harmony, obj$Seurat) %>% as.data.frame.matrix
plot.df <- plot.df[intersect(colnames(plot.df), rownames(plot.df)), ]
pheatmap::pheatmap(sweep(plot.df, 2, colSums(plot.df), '/'), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, fontsize = 14, filename = 'revise.ms/reviewer.2/figs/harmony.seurat.pdf', width = 5, height = 5)

obj.xx <- subset(obj, Seurat == 'Ery')
VlnPlot(obj.xx, feature = c('HBB'), group.by = 'harmony')
ggsave('revise.ms/reviewer.2/figs/harmony.ery.vln.pdf', width = 5, height = 4)

#-------------------------------------------------------------	
# Scanpy analysis
library(SeuratDisk)

obj$CellType <- Idents(obj) %>% as.vector
obj$orig.ident <- as.vector(obj$orig.ident)
SaveH5Seurat(obj, filename = "revise.ms/reviewer.2/data/thymus.h5Seurat")
Convert("revise.ms/reviewer.2/data/thymus.h5Seurat", dest = "h5ad")

import scanpy as sc
import pandas as pd

#ad = sc.read_h5ad('revise.ms/reviewer.2/data/thymus.h5ad')
ad = sc.read_h5ad('thymus.h5ad')
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
ad.raw = ad
ad = ad[:, ad.var.highly_variable]
sc.pp.regress_out(ad, ['S.Score', 'G2M.Score', 'percent.mt', 'percent.ribo'])	
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, svd_solver="arpack")
#sc.external.pp.bbknn(ad, batch_key="orig.ident")	
bbknn.bbknn(ad,batch_key='orig.ident')
sc.pp.pca(ad)
sc.pp.neighbors(ad, n_neighbors=10, n_pcs=30)
sc.tl.leiden(ad, resolution = 0.8, key_added = "leiden_0.8")
pd.DataFrame(ad.obsm['X_umap']).to_csv('umap.coord.csv', index = False)
pd.DataFrame(ad.obs).to_csv('meta.data.bbknn.csv')

sc.tl.rank_genes_groups(ad, "leiden", method="wilcoxon")
sc.tl.umap(ad)

#-------------------------------------------------------------
# Read BBKNN results	

meta.data <- read.csv('revise.ms/meta.data.bbknn.csv')
coord <- read.csv('revise.ms/umap.coord.csv')
rownames(coord) <- meta.data[, 1]
colnames(coord) <- c('BBKNN_1', 'BBKNN_2')
UMAP_coordinates_mat <- as(coord, "matrix")

# Create DimReducObject and add to object
obj[['BBKNN']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "BBKNN_", global = T, assay = "RNA")

obj$leiden_0.8 <- meta.data$leiden_0.8

require(future)
options(future.globals.maxSize = 400000 * 1024^2)
plan("multiprocess", workers = 20)
Idents(obj) <- obj$leiden_0.8
de.markers <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)
top20 <- de.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
anno <- c(
	'0' = 'T',
	'1' = 'T',
	'2' = 'T',
	'3' = 'T',
	'4' = 'B',
	'5' = 'T',
	'6' = 'T',
    '7' = 'T',
	'8' = 'T',
	'9' = 'Stromal',
    '10' = 'T',
	'11' = 'T',
	'12' = 'Myeloid',
    '13' = 'Stromal',
	'14' = 'Stromal',
	'15' = 'Plasma',
	'16' = 'Ery'
)

obj <- RenameIdents(obj, anno)
levels(obj) <- levels(obj.harmon$Harmony)
obj$bbknnRes <- Idents(obj)	
gp <- DimPlot(obj, label = TRUE, reduction = 'BBKNN') + DimPlot(obj, label = TRUE, group.by = 'leiden_0.8', reduction = 'BBKNN') + DimPlot(obj, group.by = 'AgeStatusNew', reduction = 'BBKNN', cols = AGE_STATUS_EXTENT_COLOR[1 : 4])
ggsave('revise.ms/reviewer.2/figs/BBKNN.UMAP.pdf', width = 18, height = 5)

plot.df <- table(obj$bbknnRes, obj$Seurat) %>% as.data.frame.matrix
plot.df <- plot.df[intersect(colnames(plot.df), rownames(plot.df)), ]
pheatmap::pheatmap(sweep(plot.df, 2, colSums(plot.df), '/'), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, fontsize = 14, filename = 'revise.ms/reviewer.2/figs/bbknn.seurat.pdf', width = 5, height = 5)

obj.xx <- subset(obj, Seurat == 'Ery')
VlnPlot(obj.xx, feature = c('HBB'), group.by = 'bbknnRes')
ggsave('revise.ms/reviewer.2/figs/bbknn.ery.vln.plot', width = 5, height = 4)


plot.df <- table(obj$bbknnRes, obj$harmony) %>% as.data.frame.matrix
plot.df <- plot.df[intersect(colnames(plot.df), rownames(plot.df)), ]
pheatmap::pheatmap(sweep(plot.df, 2, colSums(plot.df), '/'), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, fontsize = 14, filename = 'revise.ms/reviewer.2/figs/bbknn.Harmony.pdf', width = 5, height = 5)

