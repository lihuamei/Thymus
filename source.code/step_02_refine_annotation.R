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

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '2.fine_anno/data')
out.figs.dir <- file.path('../3.results', '2.fine_anno/figs')

dir.create(file.path('../3.results', '2.fine_anno'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

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

# FindAllMarkers  
Idents(obj) <- factor(merged.obj.1@meta.data$Anno.Level.1, levels = ANNO_SC_IDENT_LEVEL_1)
de.markers <- findAllMarkersPatched(obj, sprintf('sc.markers.%s.rds', anno.title), TRUE, anno.title, out.data.dir, min.pct = 0.25)
top10 <- de.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlotCustom(obj, top10$gene, flip = FALSE) + theme(axis.text.x = element_text(face = "italic"))
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.top10.Bubble.pdf', anno.title)), height = 4, width = 14)

names(ANNO_SC_COLOR_LEVEL_1) <- ANNO_SC_IDENT_LEVEL_1
gp1 <- DimPlot(merged.obj.1, label = T, raster = FALSE, cols = ANNO_SC_COLOR_LEVEL_1, group.by = sprintf('Anno.%s', anno.title), label.size = 6)
gp2 <- DimPlot(merged.obj.1, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Anno.UMAP.pdf', anno.title)), width = 12, height = 6)

gp3 <- DimPlot(merged.obj.1, label = F, raster = FALSE, cols = ANNO_SC_COLOR_LEVEL_1, split.by = 'AgeStatus', label.size = 6, group.by = sprintf('Anno.%s', anno.title))
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DotPlotCustom(obj, ANNO_SC_MARKERS %>% unlist %>% as.vector, flip = TRUE)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.canonical.markers.Bubble.pdf', anno.title)), height = 8, width = 5)

FeaturePlot(obj, features = c('HBB', 'HBG1', 'IGHG1', 'IGHG2', 'CD19', 'MS4A1', 'CST3', 'LYZ', 'IFITM3', 'ACTA2', 'CD3D', 'CD3E'), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.canonical.markers.expr.pdf', anno.title)), height = 8, width = 12)

markers.save[['Level.1']] <- de.markers
#-------------------------------------------------------------
# Refine annotation for B cells

anno.title <- 'Level.2'
obj.B <- subset(merged.obj.1, Anno.Level.1 == 'B')
obj <- runReClustering(obj.B, integed = F, res = 0.3, regress = T)

DefaultAssay(obj) <- "RNA"
DimPlot(obj, label = T, label.size = 6) + FeaturePlot(obj, feature = c('CD19', 'CD79A', 'CD79B', 'CD27', 'IGHD', 'IGHA1', 'IGHM', 'CD24'), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.B.markers.UMAP.pdf', anno.title)), width = 16, height = 8)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('%s.B.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

# Find marker genes
de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.B.rds', T, 'Cluster.B', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.3', 'Anno.Level.2.B', out.figs.dir, width = 16, height = 4.5)
obj@misc[['markers']] <- ovp.markers

root.cells <- subset(obj, idents = '2') %>% Cells %>% .[1 : 10]
psedo <- monocle3Analysis(obj, root.cells)
ggsave(file.path(out.figs.dir, 'Anno.Level.2.B.pseudotime.UMAP.pdf'), plot = psedo$gp,  height = 6, width = 6)

# Update annotation
new.cluster <- c('B_memory', 'B_memory', 'B_naive', 'B_trans', 'B_memory', 'B_memory')
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.1", "Anno.Level.2")
Idents(obj) <- factor(Idents(obj), levels = ANNO_SC_B_IDENT)

# UMAP plots	
sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_B_IDENT, sample.size = sample.size, prefix = sprintf('Anno.%s.B', anno.title), out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.B.model.pickle', anno.title)))

de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.B.rds', T, 'Anno.B', out.data.dir, min.pct = 0.25)
obj@meta.data$Anno.Level.2 <- Idents(obj) %>% as.vector
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.2', 'Anno.B', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['B.degs']] <- de.markers
merged.obj.1@misc[['B.markers']] <- ovp.markers

gp1 <- DimPlot(obj, label = T, raster = FALSE,, label.size = 5) + ggtitle('Fine annotation: B cells') 
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.B.UMAP.pdf', anno.title)), width = 12, height = 6)

DotPlotCustom(obj, ANNO_SC_B_MARKERS %>% unlist %>% as.vector, flip = TRUE)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.B.canonical.markers.Bubble.pdf', anno.title)), height = 6, width = 5)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.B.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.B.rds'))

markers.save[['Level.2']] <- de.markers
#--------------------------------------------------------------
# Refine annotation for Myeloid cells.

anno.title <- 'Level.3'
obj.Myeloid <- subset(merged.obj.1, Anno.Level.2 == 'Myeloid')
obj <- runReClustering(obj.Myeloid, integed = F, regress = T, res = 0.5)

DefaultAssay(obj) <- 'RNA'
de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.Myeloid.rds', T, 'Cluster.Myeloid', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.5', 'Anno.Myeloid', out.figs.dir, width = 18, height = 4)
obj@misc[['markers']] <- ovp.markers

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('%s.Myeloid.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DimPlot(obj, label = T, label.size = 6) + FeaturePlot(obj, feature = c('C1QA', 'C1QB', 'S100A8', "S100A9", 'HLA-DQB1', 'CLEC9A','LILRA4', "JCHAIN"), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Myeloid.markers.UMAP.pdf', anno.title)), width = 16, height = 8)

new.cluster <- c('Mac1', 'DC', 'Mac2', 'Mono1', 'Mono2', 'pDC')
names(new.cluster) <- levels(obj)
obj.bak <- obj <- RenameIdents(obj, new.cluster)
Idents(obj) <- factor(Idents(obj), levels = ANNO_SC_MYELOID_IDENT)

gp1 <- DimPlot(obj, label = T, raster = FALSE, label.size = 6)
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Myeloid.UMAP.pdf', anno.title)), width = 12, height = 6)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.2", "Anno.Level.3")

obj@meta.data$Anno.Level.3 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.Myeloid.rds', T, 'Anno.Myeloid', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.3', 'Anno.Myeloid', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['Myeloid.degs']] <- de.markers
merged.obj.1@misc[['Myeloid.markers']] <- ovp.markers

# Tranin model
sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_MYELOID_IDENT, sample.size = sample.size, prefix = 'Level.3.Myeloid', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.Myeloid.model.pickle', anno.title)))

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Myeloid.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.Myeloid.rds'))

markers.save[['Level.3']] <- de.markers
#-------------------------------------------------------------
# Refine annotation of TEC subset 

anno.title <- 'Level.4'
obj.Stromal <- subset(merged.obj.1, Anno.Level.3 == 'Stromal')
obj <- runReClustering(obj.Stromal, integed = F, res = 0.3, regress = F)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('%s.Stromal.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DefaultAssay(obj) <- 'RNA'
de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.Stromal.rds', T, 'Cluster.Stromal', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.3', 'Cluster.Stromal', out.figs.dir, top = 5, width = 18, height = 4)
obj@misc[['markers']] <- ovp.markers

DimPlot(obj, label = T, label.size = 6) + FeaturePlot(obj, feature = c("COL1A2", "COL1A1", "TOP2A", "UBE2C", "PECAM1", "PLVAP", "CD74", "HLA-DRB1", "ACTA2", "NOTCH3", "CCL25", "S100A14", "KRT19", "KRT17"), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Stromal.select.markers.expr.pdf', anno.title)), width = 16, height = 8)

new.cluster <- c('Fb1', 'TEC', 'Endo', 'Fb2', 'Fb3', 'VSMC', 'Lymph', 'TEC', 'Endo', 'Fb4', 'Fb_cycling', 'Endo')
names(new.cluster) <- levels(obj)
obj.bak <- obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.3", "Anno.Level.4")

gp1 <- DimPlot(obj.bak, label = T, raster = FALSE, label.size = 6) + ggtitle('Stromal annotation') + theme(plot.title = element_text(hjust = 0.5))
gp2 <- DimPlot(obj.bak, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus') 
gp1 + gp2
ggsave(file.path(out.figs.dir, 'Anno.Level.4.Stromal.anno.UMAP.pdf'), height = 6, width = 12)

obj@meta.data$Anno.Level.4 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.Stromal.rds', T, 'Anno.Stromal', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.4', 'Anno.Stromal', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['Stromal.degs']] <- de.markers
merged.obj.1@misc[['Stromal.markers']] <- ovp.markers

# Tranin model
sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_STROMAL_IDENT, sample.size = sample.size, prefix = 'Level.4.Stromal', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.Stromal.model.pickle', anno.title)))

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Stromal.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.Stromal.rds'))

markers.save[['Level.4']] <- de.markers
#-------------------------------------------------------------------
# Refine TEC subset

anno.title <- 'Level.5'
obj.sub <- runReClustering(subset(obj, idents = 'TEC'), integed = F, res = 0.5, regress = F)
gp3 <- DimPlot(obj.sub, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.TEC.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DefaultAssay(obj.sub) <- 'RNA'
de.markers <- findAllMarkersPatched(obj.sub, 'sc.markers.cluster.TEC.rds', T, 'Cluster.TEC', out.data.dir, min.pct = 0.25)
markerValidate(obj.sub, de.markers, SIG.REF, 'integrated_snn_res.0.5', 'Cluster.TEC', out.figs.dir, width = 18, height = 4, top2 = 50)
obj.sub@misc[['markers']] <- ovp.markers

DimPlot(obj.sub, label = T, raster = FALSE) + FeaturePlot(obj.sub, feature = c("CCL25", "PRSS16", "S100A14", "PAX1", "KRT14", "CCL19", "CLU", "KRT15", "MDK"), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.TEC.select.markers.expr.pdf', anno.title)), height = 8, width = 16)

new.cluster <- c('mTEC1', 'cTEC1', 'cTEC2', 'cTEC3', 'mTEC2', 'mTEC3', 'cTEC4')
names(new.cluster) <- levels(obj.sub)
obj.sub <- RenameIdents(obj.sub, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj.sub), merged.obj.1, "Anno.Level.4", "Anno.Level.5")
DimPlot(obj.sub, label = T)
ggsave(file.path(out.figs.dir, 'Anno.Level.5.TEC.anno.UMAP.pdf'), height = 6, width = 8)

gp1 <- DimPlot(obj.sub, label = T, raster = FALSE, label.size = 5)
gp2 <- DimPlot(obj.sub, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.TEC.anno.UMAP.pdf', anno.title)), width = 12, height = 6)

obj.sub@meta.data$Anno.Level.5 <- Idents(obj.sub)
de.markers <- findAllMarkersPatched(obj.sub, 'sc.markers.anno.TEC.rds', T, 'Anno.TEC', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj.sub, de.markers, SIG.REF, 'Anno.Level.5', 'Anno.TEC', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['TEC.degs']] <- de.markers
merged.obj.1@misc[['TEC.markers']] <- ovp.markers

sample.size <- table(Idents(obj.sub)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj.sub, ANNO_SC_TEC_IDENT, sample.size = sample.size, prefix = 'Level.5.TEC', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.TEC.model.pickle', anno.title)))

gp3 <- DimPlot(obj.sub, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.TEC.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj.sub, file.path(out.data.dir, 'obj.TEC.rds'))
markers.save[['Level.5']] <- de.markers

#------------------------------------------------------------
# Refine T cells

obj.T <- subset(merged.obj.1, Anno.Level.5 == 'T')
merged.obj.1 <- updateAnnoRes(Idents(obj.T), merged.obj.1, "Anno.Level.5", "Anno.Level.6")
saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))

#-------------------------------------------------------------
# Refine DN cells

anno.title <- "Level.7"
obj.DN <- subset(merged.obj.1, Anno.Level.6 %in% c('DN'))
obj <- runReClustering(obj.DN, integed = F, res = 0.5, regress = F)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DN.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DefaultAssay(obj) <- 'RNA'
de.markers <- findAllMarkersPatched(obj, 'sc.markers.T.cluster.DN.rds', T, 'Cluster.DN', out.data.dir, min.pct = 0.25)
markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.5', 'Cluster.T.DN', out.figs.dir, width = 18, height = 4, top2 = 50)
obj@misc[['markers']] <- ovp.markers

DimPlot(obj, label = T, raster = FALSE) + FeaturePlot(obj, feature = c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CDK1', 'RAG1', 'IGLL1'), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DN.markers.UMAP.pdf', anno.title)), height = 8, width = 16)

root.cells <- subset(obj, idents = '7') %>% Cells %>% .[1:2]
psedo <- monocle3Analysis(obj, root.cells)
ggsave(file.path(out.figs.dir, 'Anno.Level.7.DN.pseudotime.UMAP.pdf'), plot = psedo$gp,  height = 6, width = 6)

new.cluster <- c('DN_re1', 'DN_re2', 'DN_trans', 'DN_blast1', 'DN_ISP', 'DN_blast2', 'DN_mito', 'DN_early', 'DN_re3')
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.6", "Anno.Level.7")

gp1 <- DimPlot(obj, label = T, raster = FALSE, label.size = 5)
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DN.anno.UMAP.pdf', anno.title)), width = 12, height = 6)

obj@meta.data$Anno.Level.7 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.DN.rds', T, 'Anno.DN', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.7', 'Anno.DN', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['DN.degs']] <- de.markers
merged.obj.1@misc[['DN.markers']] <- ovp.markers

sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_DN_IDENT, sample.size = sample.size, prefix = 'Level.7.DN', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.DN.model.pickle', anno.title)))

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DN.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.DN.rds'))
markers.save[['Level.7']] <- de.markers

#------------------------------------------------------------
# Refine DP cells 

anno.title <- "Level.8"
obj.DP <- subset(merged.obj.1, Anno.Level.7 %in% c('DP1', 'DP2', 'DP3', 'DP4'))
DimPlot(obj.DP, label = T, raster = FALSE) + FeaturePlot(obj.DP, feature = c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CDK1', 'RAG1', 'IL2RA', 'ITM2A'), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.markers.UMAP.pdf', anno.title)), height = 8, width = 16)

obj <- subset(merged.obj.1, Anno.Level.7 %in% c('DP1', 'DP2'))
obj <- runReClustering(obj, integed = F, res = 0.3, regress = F)
DimPlot(obj, label = T, pt.size = 0.6, label.size = 6)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.P.Cluster.UMAP.pdf', anno.title)), width = 6, height = 6)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.P.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DefaultAssay(obj) <- 'RNA'
de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.T.DP.P.rds', T, 'Cluster.DP.P', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.3', 'Cluster.T.DP.P', out.figs.dir, width = 18, height = 4, top2 = 50)
obj@misc[['markers']] <- ovp.markers

root.cells <- subset(obj, idents = '6') %>% Cells %>% .[1 : 2]
psedo <- monocle3Analysis(obj, root.cells)
ggsave(file.path(out.figs.dir, 'Anno.Level.8.DP.P.pseudotime.UMAP.pdf'), plot = psedo$gp,  height = 6, width = 6)

new.cluster <- c('DP_blast1', 'DP_blast2', 'DP_blast3', 'DP_blast4', 'DP_blast5', 'DP_blast6', 'DP_blast7', 'DP_blast8')
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.7", "Anno.Level.8")

gp1 <- DimPlot(obj, label = T, raster = FALSE, label.size = 5)
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.P.anno.UMAP.pdf', anno.title)), width = 12, height = 6)

# Re-find all markers.
obj@meta.data$Anno.Level.8 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.DP.P.rds', T, 'Anno.DP.P', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.8', 'Anno.DP.P', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['DP.blast.degs']] <- de.markers
merged.obj.1@misc[['DP.blast.markers']] <- ovp.markers

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.P.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_DP_P_IDENT, sample.size = sample.size, prefix = 'Level.8.DP.P', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.DP.P.model.pickle', anno.title)))

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.DP.P.rds'))
markers.save[['Level.8']] <- de.markers

#----------------------------------------------------------------------
# DP re-arrangement

anno.title <- "Level.9"
obj.DP.R <- subset(merged.obj.1, Anno.Level.8 %in% c('DP3', 'DP4'))
obj <- runReClustering(obj.DP.R, integed = F, res = 0.3, regress = F)
DimPlot(obj, label = T, pt.size = 0.6, label.size = 6)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.R.Cluster.UMAP.pdf', anno.title)), width = 6, height = 6)
	
DefaultAssay(obj) <- 'RNA'
obj.sub <- subset(obj, idents = '6')
pred.res <- lr.model$predict(obj.sub %>% GetAssayData(.) %>% t) %>% `names<-`(colnames(obj.sub))
merged.obj.1 <- updateAnnoRes(pred.res, merged.obj.1, "Anno.Level.8", "Anno.Level.9")

obj <- subset(obj, idents = c('0', '1', '2', '3', '4', '5'))
gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.R.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.T.DP.R.rds', T, 'Cluster.DP.R', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.0.3', 'Cluster.T.DP.R', out.figs.dir, width = 18, height = 4, top2 = 50)
obj@misc[['markers']] <- ovp.markers

root.cells <- subset(obj, idents = '5') %>% Cells %>% .[1 : 2]
psedo <- monocle3Analysis(obj, root.cells)
ggsave(file.path(out.figs.dir, 'Anno.Level.9.DP.R.pseudotime.UMAP.pdf'), plot = psedo$gp,  height = 6, width = 6)

new.cluster <- c('DP_re1', 'DP_re2', 'DP_re3', 'DP_re4', 'DP_re5', 'DP_re6')
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.9", "Anno.Level.9")

gp1 <- DimPlot(obj, label = T, raster = FALSE, label.size = 5)
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.R.anno.UMAP.pdf', anno.title)), width = 12, height = 6)

# Re-find all markers.
obj@meta.data$Anno.Level.9 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.DP.R.rds', T, 'Anno.DP.R', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.9', 'Anno.DP.R', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['DP.R.degs']] <- de.markers
merged.obj.1@misc[['DP.R.markers']] <- ovp.markers

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.DP.R.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_DP_R_IDENT, sample.size = sample.size, prefix = 'Level.9.DP.R', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.DP.R.model.pickle', anno.title)))

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.DP.R.rds'))
markers.save[['Level.9']] <- de.markers

#------------------------------------------------------------
# Refine other T cells annotation

anno.title <- 'Level.10'
obj.T <- subset(merged.obj.1, Anno.Level.9 %in% c('CD4+ T', "NKT", "Treg", "CD8aa", "CD8+ T"))
obj <- runReClustering(obj.T, split.by = 'orig.ident', integed = T, res = 1.5, regress = T)

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Other.T.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

DefaultAssay(obj) <- 'RNA'
de.markers <- findAllMarkersPatched(obj, 'sc.markers.cluster.Other.T.rds', T, 'Cluster.Other.T', out.data.dir, min.pct = 0.5)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'integrated_snn_res.1.5', 'Cluster.Other.T', out.figs.dir, width = 18, height = 10, top = 5, top2 = 50)
obj@misc[['markers']] <- ovp.markers

DimPlot(obj, label = T, raster = FALSE) + 
	FeaturePlot(obj, feature = c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'S100A4', 'CCR7', 'SELL', 'FOXP3', "FXYD5", "NKG7"), cols = FEATURES_COLOR)
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Other.T.markers.UMAP.pdf', anno.title)), height = 8, width = 20)

new.cluster <- c(
	'CD4T', 'CD8T', 'CD8T', 'CD8aa.I', 'Treg.diff', 'CD8T_mem', 'T_agonist', 'CD4T_mem', 'T_apoptosis', 'CD8aa.II', 'CD4T', 'CD4T', 'CD4T_mem', 'Treg', 'CD8aa.I', 
	'DP.Q', 'DP.P', 'CD4T_mem', 'NKT', 'Treg.diff', 'CD8T', 'T_proliferating', 'ILC3'			
)
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
merged.obj.1 <- updateAnnoRes(Idents(obj), merged.obj.1, "Anno.Level.9", "Anno.Level.10")

# LR model
obj.sub <- subset(obj, idents = 'DP.Q')
lr.model <- py_load_object(file.path(out.data.dir, 'Level.9_lr.DP.R.model.pickle'))
pred.res <- lr.model$predict(obj.sub %>% GetAssayData(.) %>% t) %>% `names<-`(colnames(obj.sub))
merged.obj.1 <- updateAnnoRes(pred.res, merged.obj.1, "Anno.Level.10", "Anno.Level.10")

obj.sub <- subset(obj, idents = 'DP.P')
lr.model <- py_load_object(file.path(out.data.dir, sprintf('Level.8_lr.DP.P.model.pickle')))
pred.res <- lr.model$predict(obj.sub %>% GetAssayData(.) %>% t) %>% `names<-`(colnames(obj.sub))
merged.obj.1 <- updateAnnoRes(pred.res, merged.obj.1, "Anno.Level.10", "Anno.Level.10")

#----------------------------------------
obj <- subset(obj, idents = setdiff(levels(obj), c('DP.Q', 'DP.P')))
gp1 <- DimPlot(obj, label = T, raster = FALSE, label.size = 5)
gp2 <- DimPlot(obj, label = F, raster = FALSE, cols = SampleClassifyColors, group.by = 'AgeStatus')
gp1 + gp2
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.T.Others.anno.UMAP.pdf', anno.title)), width = 12, height = 6)

obj@meta.data$Anno.Level.10 <- Idents(obj)
de.markers <- findAllMarkersPatched(obj, 'sc.markers.anno.other.T.rds', T, 'Anno.Other.T', out.data.dir, min.pct = 0.25)
ovp.markers <- markerValidate(obj, de.markers, SIG.REF, 'Anno.Level.10', 'Anno.Other.T', out.figs.dir, width = 16, height = 4.5, top2 = 50)
merged.obj.1@misc[['Other.T.degs']] <- de.markers
merged.obj.1@misc[['Other.T.markers']] <- ovp.markers

gp3 <- DimPlot(obj, label = F, raster = FALSE, split.by = 'AgeStatus', label.size = 6, pt.size = 0.5)
gp3 + ggtitle('')
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.Other.T.Split.by.AgeStatus.UMAP.pdf', anno.title)), width = 12, height = 6)

sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_SC_OTHER_T_INDET, sample.size = sample.size, prefix = 'Level.10.Other.T', out.figs.dir = out.figs.dir)
py_save_object(lr.model, file = file.path(out.data.dir, sprintf('%s_lr.Other.T.model.pickle', anno.title)))

saveRDS(merged.obj.1@meta.data, file.path(out.data.dir, 'meta.data.rds'))
saveRDS(obj, file.path(out.data.dir, 'obj.Other.T.rds'))
markers.save[['Level.10']] <- de.markers

#-------------------------------------------
# save to output directory

obj <- merged.obj.1
Idents(obj) <- factor(obj@meta.data$Anno.Level.10, levels = ANNO_ENTIRE_IDNET)
new.cluster <- c(
		"Ery", "B_naive", "B_trans", "B_memory", "Plasma", "Mono", "Mono", "Mac", "Mac", "DC", "pDC", "Fb", "Fb", "Fb", "Fb", "Fb_cycling", "Endo", "Tuft", 
		"VSMC", "Lymph", "cTEC", "cTEC", "cTEC", "cTEC", "mTEC", "mTEC", "mTEC","DN_early", "DN_re", "DN_blast", "DN_blast", "DN_re", "DN_re", "DN_re", 
		"DN_re", "DN_re", "DP_blast", "DP_blast", "DP_blast", "DP_blast", "DP_blast", "DP_blast", "DP_blast", "DP_blast", "DP_re", "DP_re", 
		"DP_re", "DP_re", "DP_re","DP_re", "abT(entry)", "abT(entry)", "CD4T", "CD4T_mem", "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa","CD8aa", 
		"T_agonist", "T_apoptosis", "T_proliferating", "NKT", "ILC3")

names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
Idents(obj) <- factor(Idents(obj), levels = ANNO_ENTIRE_IDNET_FIG1)
obj@meta.data$Anno.Level.Fig.1 <- Idents(obj) %>% as.vector	

saveRDS(markers.save, file = file.path(out.data.dir, 'markers.save.rds'))
saveRDS(obj, file = file.path(out.data.dir, 'merged.obj.anno.rds'))
saveAsH5adFile(obj, 'merged.obj.anno', out.data.dir)

