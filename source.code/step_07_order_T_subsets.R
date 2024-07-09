suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(colorRamps))
suppressMessages(library(tidyverse))
suppressMessages(library(writexl))
suppressMessages(library(clusterProfiler))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(ggtree))
suppressMessages(library(ape))
suppressMessages(library(monocle))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/findMedullaClusters.R')
source('modules/trajectory_methods.R')
source('modules/tcrSubgroup.R')
source('modules/vis/step_07_order_T_subsets_vis.R')

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '6.clonetype/data')
out.figs.dir <- file.path('../3.results', '6.clonetype/figs')

dir.create(file.path('../3.results', '6.clonetype'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load spatial data

obj.st.lst <- readRDS('../3.results/4.proj_spatial/data/obj.sp.lst.rds')
obj.st.merged <- readRDS('../3.results/4.proj_spatial/data/obj.st.merged.rds')
anno <- read.table('configs/gencode.gene.info.v22.tsv', sep = '\t', header = T) %>% .[!duplicated(.[, 2]), ] %>% `rownames<-`(.[, 2])
anno <- read.table('configs/gencode.gene.info.v22.tsv', sep = '\t', header = T) %>% .[!duplicated(.[, 2]), ] %>% `rownames<-`(.[, 2])
coding.genes <- subset(anno, (gene_type %in% ANNO.AREAS) & (gene_status == 'KNOWN')) %>% .[, 'gene_name']
obj.st.lst <- lapply(obj.st.lst, function(xx) subset(xx, idents = setdiff(levels(xx), 'Contamination spots')))
obj.merged.sc <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')

anno.meta <- read.table('pbs/.tcell.pseudotime.xls', sep = '\t', header = TRUE, row.names = 1)
obj.merged.sc$dpt_pseudotime <- anno.meta[colnames(obj.merged.sc), 'dpt_pseudotime']
obj.merged.sc$AgeStatusNew <- as.vector(obj.merged.sc$AgeStatus)
obj.merged.sc$AgeStatusNew[obj.merged.sc$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'
obj.merged.sc$AgeStatusNew <- factor(obj.merged.sc$AgeStatusNew, levels = SampleClassifyNew %>% names %>% .[1 : 4])

#--------------------------------------------------------------
# Re-order cell types for DP-blast subsets.

obj.merged.sc <- mappCellNames(obj.merged.sc, c('DN_early'), c('DN_early'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DN_blast1', 'DN_blast2'), c('DN_blast', 'DN_blast'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DN_trans', 'DN_re3'), c('DN_re2', 'DN_re2'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DN_re1', 'DN_re2'), c('DN_re1', 'DN_re1'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DN_mito', 'DN_ISP'), c('ISP', 'ISP'), 'Anno.Level.10', 'Anno.Level.11')

xx.bak <- subset(obj.merged.sc, Anno.Level.11 %in% c('DN_early', 'DN_blast', 'DN_re1', 'DN_re2', 'ISP'))
xx <- readRDS('../3.results/2.fine_anno/data/obj.DN.rds')
xx <- xx[, colnames(xx.bak)]
xx$Anno.Level.11 <- factor(xx.bak$Anno.Level.11[colnames(xx)], levels = c('DN_early', 'DN_blast', 'DN_re1', 'DN_re2', 'ISP'))
DimPlot(xx, label = T, group.by = 'Anno.Level.11', cols = ANNO_TREND_IDENT_COLOR[levels(xx$Anno.Level.11)])
ggsave(file.path(out.figs.dir, 'DN.fine.UMAP.pdf'), width = 6, height = 6)

Idents(xx) <- xx@meta.data$Anno.Level.11
de.markers <- findAllMarkersPatched(xx, sprintf('sc.markers.%s.rds', 'DN'), TRUE, 'DN', out.data.dir, min.pct = 0.25)
top10 <- de.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
de.markers <- de.markers[de.markers$gene %in% coding.genes, ]

xx <- ScaleData(xx)
DotPlotCustom(xx, top10$gene %>% unique, flip = TRUE) + theme(axis.text.x = element_text(face = "italic"))
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.top10.Bubble.pdf', 'DN')), height = 10, width = 5)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
FeaturePlot(xx, feature = c('IGLL1', 'CDK1', 'RAG1', 'RAG2', 'CD8A'), cols = myPalette(100))
ggsave(file.path(out.figs.dir, 'dn.feature.plot.pdf'), width = 8, height = 8)

# Cell-cycle scoring
cell.cycle.genes <- readxl::read_excel('./configs/cellcyle.genes.xlsx')
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(G1.S = cell.cycle.genes %>% .[, 'G1/S'] %>% unlist %>% unique %>% .[!is.na(.)]), name = 'G1.S')
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(S = cell.cycle.genes %>% .[, 'S'] %>% unlist %>% unique %>% .[!is.na(.)]), name = 'S')
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(G2.M = cell.cycle.genes %>% .[, 'G2/M'] %>% unlist %>% unique %>% .[!is.na(.)]), name = 'G2.M')
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(M = cell.cycle.genes %>% .[, 'M'] %>% unlist %>% unique %>% .[!is.na(.)]), name = 'M')
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(M.G1 = cell.cycle.genes %>% .[, 'M/G1'] %>% unlist %>% unique %>% .[!is.na(.)]), name = 'M.G1')

obj <- subset(obj.merged.sc, Anno.Level.10 %in% c('DP_blast1', 'DP_blast2', 'DP_blast3', 'DP_blast4', 'DP_blast5', 'DP_blast6', 'DP_blast7', 'DP_blast8'))
meta.data <- obj@meta.data[, c('Anno.Level.10', 'G1.S1', 'S1', 'G2.M1', 'M1', 'M.G11')]
plot.data <- tidyr::gather(meta.data, 'Status', 'Score', -Anno.Level.10)
ggplot(plot.data, aes(x = Status, y = Score, fill = Anno.Level.10)) + geom_boxplot(outlier.shape = NA)

# Re-name DP-blast cells.
obj.merged.sc <- subset(obj.merged.sc, Anno.Level.10 != 'DP_blast8')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_blast3', 'DP_blast4'), rep('DP_blast1', 2), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_blast1'), c('DP_blast2'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_blast2'), c('DP_blast3'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_blast5'), c('DP_blast4'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_blast6', 'DP_blast7'), rep('DP_blast5', 2), 'Anno.Level.10', 'Anno.Level.11')

idents <- c('DP_blast1', 'DP_blast2', 'DP_blast3', 'DP_blast4', 'DP_blast5')
xx.bak <- subset(obj.merged.sc, Anno.Level.11 %in% idents)
xx <- readRDS('../3.results/2.fine_anno/data/obj.DP.P.rds')
xx <- xx[, colnames(xx.bak)]
xx$Anno.Level.11 <- factor(xx.bak$Anno.Level.11[colnames(xx)], levels = idents)
DimPlot(xx, label = T, group.by = 'Anno.Level.11', cols = ANNO_TREND_IDENT_COLOR[levels(xx$Anno.Level.11)])
ggsave(file.path(out.figs.dir, 'DP.fine.UMAP.pdf'), width = 6, height = 6)

Idents(xx) <- xx@meta.data$Anno.Level.11
de.markers <- findAllMarkersPatched(xx, sprintf('sc.markers.%s.rds', 'DP.P'), TRUE, 'DP.P', out.data.dir, min.pct = 0.25)
de.markers <- de.markers[de.markers$gene %in% coding.genes, ]
de.markers <- de.markers[!grepl('^MT-', rownames(de.markers)), ]
top10 <- de.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

xx <- ScaleData(xx)
DotPlotCustom(xx, top10$gene %>% unique, flip = TRUE) + theme(axis.text.x = element_text(face = "italic"))
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.top10.Bubble.pdf', 'DP.P')), height = 14, width = 5)

# Re-name DP-re cells.
obj.merged.sc <- subset(obj.merged.sc, Anno.Level.10 != 'DP_re5')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_re2'), rep('DP_re1', 3), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_re1'), c('DP_re2'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_re3'), c('DP_re4'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('DP_re6', 'DP_re4', 'DP_re5'), rep('DP_re3', 2), 'Anno.Level.10', 'Anno.Level.11')

idents <- c('DP_re1', 'DP_re2', 'DP_re3', 'DP_re4')
xx.bak <- subset(obj.merged.sc, Anno.Level.11 %in% idents)
xx <- readRDS('../3.results/2.fine_anno/data/obj.DP.R.rds')
xx <- xx[, colnames(xx.bak)]
xx$Anno.Level.11 <- factor(xx.bak$Anno.Level.11[colnames(xx)], levels = idents)
DimPlot(xx, label = T, group.by = 'Anno.Level.11', cols = ANNO_TREND_IDENT_COLOR[levels(xx$Anno.Level.11)])
ggsave(file.path(out.figs.dir, 'DP.Q.fine.UMAP.pdf'), width = 6, height = 6)

Idents(xx) <- xx@meta.data$Anno.Level.11
de.markers <- findAllMarkersPatched(xx, sprintf('sc.markers.%s.rds', 'DP.P.Q'), TRUE, 'DP.P.Q', out.data.dir, min.pct = 0.25)
de.markers <- de.markers[de.markers$gene %in% coding.genes, ]
de.markers <- de.markers[!grepl('^MT-|^RPS', rownames(de.markers)), ]
top10 <- de.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

xx <- ScaleData(xx)
DotPlotCustom(xx, top10$gene, flip = TRUE) + theme(axis.text.x = element_text(face = "italic"))
ggsave(file.path(out.figs.dir, sprintf('Anno.%s.top10.Bubble.pdf', 'DP.P.Q')), height = 12, width = 5)

plotRePseudo(xx.bak)
ggsave(file.path(out.figs.dir, 'dp.re.pseudo.pdf'), width = 5, height = 4)

# SP cells
obj.merged.sc <- mappCellNames(obj.merged.sc, c('abT(entry)1', 'abT(entry)2'), c('abT(entry)', 'abT(entry)'), 'Anno.Level.10', 'Anno.Level.11')
obj.merged.sc <- mappCellNames(obj.merged.sc, c('CD8aa.I', 'CD8aa.II'), c('CD8aa', 'CD8aa'), 'Anno.Level.10', 'Anno.Level.11')

saveRDS(obj.merged.sc, file = file.path(out.data.dir, 'obj.merged.sc.reformat.new.rds'))

#-------------------------------------------------------------
# Pseudotime

order.cells <- c(ANNO_TREND_IDENT[1 : 15], 'CD8T', 'CD4T')
saveRDS(data.frame(ANNO_TREND_IDENT_COLOR[order.cells]) %>%  `colnames<-`('Color'), file = '.color.T.rds')

obj.sub <- subset(obj.merged.sc, Anno.Level.11 %in% order.cells)
obj.sub$Anno.Level.11 <- as.vector(obj.sub$Anno.Level.11)
saveAsH5adFile(obj.sub, '.TCR.tcell.obj.diffmap', './')
system('python .plot_heatmap.py')

#------------------------------------------------------------
# Vilidation for DP-blast 1 to 5.

source('modules/.DoMultiBarHeatmap.R')
obj.blast <- subset(obj.merged.sc, Anno.Level.11 %in% c('DP_blast1', 'DP_blast2', 'DP_blast3', 'DP_blast4', 'DP_blast5'))
cycle.genes <- tidyr::gather(cell.cycle.genes)[, 2] %>% unlist(., use.names = F) %>% unique %>% .[!is.na(.)]
obj.blast <- ScaleData(obj.blast)
cellcyleHeatmap(obj.blast, cell.cycle.genes, out.figs.dir)

matureStatus(obj.merged.sc, out.figs.dir)	
#-------------------------------------------------------------
# T-subgroups differential trend

tcrSubgroupsTrend(obj.merged.sc, permu = 100)
ggsave(file.path(out.figs.dir, 'fit.percent.four.types_scatterplot.pdf'), width = 8, height = 10)

#-------------------------------------------------------------
# DEGs across T-subgroups

library(Scillus)
source('modules/DoMultiBarHeatmap.R')

order.cells <- c(ANNO_TREND_IDENT[1 : 15], 'CD8T', 'CD4T')
tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17 : 33]))
obj <- tcrSubgroup1(obj) %>% subset(., trb != '')
obj.1 <- subset(obj, t.cell_group %in% tcr.subtypes)
obj.1$Anno.Level.11 <- factor(obj.1$Anno.Level.11, levels = ANNO_TREND_IDENT)

Idents(obj.1) <- obj.1$t.cell_group
levels(obj.1) <- tcr.subtypes
obj.sub <- downSamplSeurat(obj.1, cnt = 1000, seed = 1234)
require(future)
options(future.globals.maxSize = 40000 * 1024^2)
plan("multiprocess", workers = 40)

obj.sub <- subset(obj.sub, (Anno.Level.11 %in% setdiff(levels(obj.sub$Anno.Level.11), c('T_proliferating', 'CD8aa'))))
de.markers <- FindAllMarkers(obj.sub, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)
de.markers <- de.markers[de.markers$gene %in% coding.genes, ]

obj.sub <- ScaleData(obj.sub)
pdf(file.path(out.figs.dir, 'deg.expr.four.subtypes_heatmap.pdf'), width = 10, height = 21)
plotHeamtmapForDEGs(obj.sub)
dev.off()

pdf(file.path(out.figs.dir, 'four.subgroups.de.markers_Venn.pdf'), width = 7, height = 7)
temp <- overlapDegFourSubgroups(de.markers)
grid.draw(temp)
dev.off()

#------------------------------------------------------------
# Score for top20 genes

top20 <- de.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
xx.genes <- split(top20$gene, top20$cluster)

obj.st.lst <- lapply(obj.st.lst, function(obj.st) {
    obj.st <- AddModuleScore(object = obj.st, features = xx.genes, name = 'T.Subgroups')
	obj.st
})

plot.tar <- c('T.Subgroups1', 'T.Subgroups2', 'T.Subgroups3', 'T.Subgroups4')
fourGroupsTrend(obj.st.lst, plot.tar)
ggsave(file.path(out.figs.dir, 't.subgourps.score.along.sp.distance.pdf'), width = 6, height = 4)

fourSubgroupScores(obj.merged.sc, out.figs.dir)

#------------------------------------------------------------
# Track clone

tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17 : 33]))
obj$Anno.Level.12 <- obj$Anno.Level.11
obj$Anno.Level.12[obj$Anno.Level.12 %in% ANNO_TREND_IDENT[16 : 26]] <- 'SP'
obj$Anno.Level.12 <- factor(obj$Anno.Level.12, levels = c(ANNO_TREND_IDENT[1 : 15], 'SP'))

obj <- tcrSubgroup(obj) %>% subset(., trb != '')
obj.1 <- subset(obj, t.cell_group %in% tcr.subtypes)
Idents(obj.1) <- obj.1$t.cell_group
obj.2 <- subset(obj.1, trb != '')
xx.1 <- subset(obj.2, idents = "TRA-TRB+TRD+")$trb %>% unique
xx.2 <- subset(obj.2, idents = "TRA+TRB+TRD+")$trb %>% unique
xx.3 <- subset(obj.2, idents = "TRA-TRB+TRD-")$trb %>% unique
xx.4 <- subset(obj.2, idents = "TRA+TRB+TRD-")$trb %>% unique
tcr.ovp <- Reduce(intersect, list(xx.1 = xx.1, xx.2 = xx.2, xx.3 = xx.3, xx.4 = xx.4))

obj.2.sub <- subset(obj.2, trb %in% tcr.ovp)
pp <- table(obj.2.sub$t.cell_group, obj.2.sub$Anno.Level.12) %>% t
trackClonalLine(pp)
ggsave(file.path(out.figs.dir, 'four.subgroups.tracking.clones_Line.pdf'), width = 7, height = 5)

travUsages(obj.1, out.figs.dir)
trbvUsages(obj.1, out.figs.dir)

#-------------------------------------------------------------
# Validation the ordered cell subsets.

merged.obj <- subset(obj.merged.sc, tcr == 1)

detach("package:monocle3", unload=TRUE)
obj <- subset(merged.obj, Anno.Level.11 %in% c("DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5")) 
Idents(obj) <- obj$Anno.Level.11
obj <- obj %>% downSamplSeurat(., cnt = 500)

sc.cds <- createMonocleObject(obj)
sc.cds <- detectGenes(sc.cds, min_expr = 3) 
sc.cds <- sc.cds[fData(sc.cds)$num_cells_expressed > 10, ]

cds <- sc.cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 

pData(cds)$Cluster <- pData(cds)$Anno.Level.11
diff.test.res <- differentialGeneTest(cds, fullModelFormulaStr = "~Cluster")

ordering.genes <- row.names(subset(diff.test.res, qval < 1e-5))
cds <- setOrderingFilter(cds, ordering.genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
plot_cell_trajectory(cds, color_by = "Cluster") 
ggsave(file.path(out.figs.dir, 'DP_blast.monocle2.celltype_scatter.pdf'), width = 6, height = 6)

cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(file.path(out.figs.dir, 'DP_blast.monocle2.celltype_Pseudotime.pdf'), width = 6, height = 6)

idents <- c("DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5")
xx.bak <- obj
xx <- readRDS('../3.results/2.fine_anno/data/obj.DP.P.rds')
xx <- xx[, colnames(obj)]
xx$Pseudotime <- pData(cds)[colnames(xx), 'Pseudotime']
xx$Anno.Level.11 <- factor(xx.bak$Anno.Level.11[colnames(xx)], levels = idents)
DimPlot(xx, label = T, group.by = 'Anno.Level.11', cols = ANNO_TREND_IDENT_COLOR[levels(xx$Anno.Level.11)]) +
FeaturePlot(xx, feature = 'Pseudotime', cols = myPalette(100))
ggsave(file.path(out.figs.dir, 'DP_re.monocle2.celltype_Pseudotime_new.pdf'), width = 10, height = 6)

#-----------------------------------------------------------------
# Pie-plot for DP-blasts

cols <- c('gray', '#4D71B2', '#4C0E81')
xx.counts <- table(obj$Anno.Level.11, obj$Phase) %>% as.data.frame.matrix %>% t

pdf(file.path(out.figs.dir, 'DP_blast1.G.M.S_pieplot.pdf'), width = 4, height = 4)
pie(xx.counts[, 1], labels = rownames(xx.counts), col = cols, main = 'DP_blast1')
dev.off()

pdf(file.path(out.figs.dir, 'DP_blast2.G.M.S_pieplot.pdf'), width = 4, height = 4)
pie(xx.counts[, 2], labels = rownames(xx.counts), col = cols, main = 'DP_blast2')
dev.off()

pdf(file.path(out.figs.dir, 'DP_blast3.G.M.S_pieplot.pdf'), width = 4, height = 4)
pie(xx.counts[, 3], labels = rownames(xx.counts), col = cols, main = 'DP_blast3')
dev.off()

#------------------------------------------------------------------
# For DP-re subsets

merged.obj <- subset(obj.merged.sc, tcr == 1)
obj <- subset(merged.obj, Anno.Level.11 %in% c("DP_re1", "DP_re2", "DP_re3", "DP_re4")) 
Idents(obj) <- obj$Anno.Level.11
obj <- obj %>% downSamplSeurat(., cnt = 1000)
sc.cds <- createMonocleObject(obj)
sc.cds <- detectGenes(sc.cds, min_expr = 3)
sc.cds <- sc.cds[fData(sc.cds)$num_cells_expressed > 10, ]

cds <- sc.cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

pData(cds)$Cluster <- pData(cds)$Anno.Level.11
diff.test.res <- differentialGeneTest(cds, fullModelFormulaStr = "~Cluster")

ordering.genes <- row.names(subset(diff.test.res, qval < 1e-2))
cds <- setOrderingFilter(cds, ordering.genes)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', residualModelFormulaStr = "~orig.ident")
plot_cell_trajectory(cds, color_by = "Cluster")
ggsave(file.path(out.figs.dir, 'DP_re.monocle2.celltype_scatter.pdf'), width = 6, height = 6)

cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(file.path(out.figs.dir, 'DP_re.monocle2.celltype_Pseudotime.pdf'), width = 6, height = 6)

idents <- c('DP_re1', 'DP_re2', 'DP_re3', 'DP_re4')
xx.bak <- subset(obj, Anno.Level.11 %in% idents)
xx <- readRDS('../3.results/2.fine_anno/data/obj.DP.R.rds')
xx <- xx[, colnames(xx.bak)]
xx$Pseudotime <- pData(cds)[colnames(xx), 'Pseudotime']
xx$Anno.Level.11 <- factor(xx.bak$Anno.Level.11[colnames(xx)], levels = idents)
DimPlot(xx, label = T, group.by = 'Anno.Level.11', cols = ANNO_TREND_IDENT_COLOR[levels(xx$Anno.Level.11)]) + 
FeaturePlot(xx, feature = 'Pseudotime', cols = myPalette(100))
ggsave(file.path(out.figs.dir, 'DP_blast.monocle2.celltype_Pseudotime_new.pdf'), width = 10, height = 6)

#------------------------------------------------------------------------------
# Relative distance

ratio <- 0.75
ratio.res <- lapply(obj.st.lst, function(obj.st) {
    obj.st.tmp <- subset(obj.st, HE.Labels %in% c("Cortex", c("Medulla_centric", "Medulla_hi", "Medulla_edge")))
    lapply(c('T.Subgroups1', 'T.Subgroups2', 'T.Subgroups3', 'T.Subgroups4'), function(gene) {
        kk.1 <- subset(obj.st.tmp, HE.Labels == "Cortex")@meta.data[, gene]
        kk.1 <- kk.1[kk.1 > quantile(kk.1, ratio)] %>% mean
        kk.2 = subset(obj.st.tmp, HE.Labels %in% c("Medulla_centric", "Medulla_hi"))@meta.data[, gene]
        kk.2 = kk.2[kk.2 > quantile(kk.2, ratio)] %>% mean
        max(kk.1, kk.2)
    }) %>% unlist(.) %>% `names<-`(c('T.Subgroups1', 'T.Subgroups2', 'T.Subgroups3', 'T.Subgroups4'))
}) %>% do.call(rbind, .)

ratio.plot <- cbind.data.frame(ratio.res, Sample = rownames(ratio.res)) %>% tidyr::gather(., 'Gene', 'Ratio', -Sample)
ratio.plot$Gene <- factor(ratio.plot$Gene, levels = c('T.Subgroups1', 'T.Subgroups2', 'T.Subgroups3', 'T.Subgroups4'))
ratio.plot$Sample <- factor(ratio.plot$Sample, levels = names(SAMPLE_COLORS))

ggplot(ratio.plot, aes(x = Gene, y = Ratio, color = Sample, group = Sample)) +
    geom_point(aes(fill=Sample), size=3) +
    geom_line(aes(color=Sample, linetype = Sample)) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab('Maturity score') + xlab('') +
    scale_color_manual(values = SAMPLE_COLORS) +
    scale_x_discrete(labels = c('TRA-TRB+TRD+(C)', 'TRA+TRB+TRD+(C)', 'TRA-TRB+TRD-(C)', 'TRA+TRB+TRD-(M)'))

ggsave(file.path(out.figs.dir, 'four.subgroups.Maturity.Cortex.medulla_Scatter.pdf'), width = 5, height = 6)

#-----------------------------------------------------------------
# CD8A vs. CD8B

obj.dn.dp <- subset(merged.obj.sc, idents = merged.obj.sc %>% levels %>% .[18 : 34])
obj.dn.dp.sub <- subset(obj.dn.dp, Anno.Level.11 %in% order.cells[1 : 12])
obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c('CD8A', 'CD8B')), CellType = obj.dn.dp$Anno.Level.10) %>% subset(., CellType %in% ANNO_ENTIRE_IDNET[51 : 65])
cd8a.avg <- aggregate(CD8A ~ CellType, data = obj.expr, FUN = 'mean')
cd8b.avg <- aggregate(CD8B ~ CellType, data = obj.expr, FUN = 'mean')
plot.expr <- cbind.data.frame(CellType = cd8a.avg[, 1], CD8A = cd8a.avg[, 2], CD8B = cd8b.avg[, 2])
ggplot(plot.expr, aes(x = CD8B, y = CD8A, color = CellType, label = CellType)) +
    geom_point(size = 4) +
    geom_text(hjust = 0, nudge_x = 0.05) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
    ggtitle('CD8A/CD8B Ratio')
ggsave(file.path(out.figs.dir, 'dn2sp_CD8A_CD8B_Ratio.scatter.pdf'), width = 6, height = 6)

#-----------------------------------------------------------------
# CD4 vs. CD8B

obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c('CD4', 'CD8B')), CellType = obj.dn.dp$Anno.Level.11) %>% subset(., CellType %in% c(rename.cells, "NKT", "ILC3"))
cd4.avg <- aggregate(CD4 ~ CellType, data = obj.expr, FUN = 'mean')
cd8b.avg <- aggregate(CD8B ~ CellType, data = obj.expr, FUN = 'mean')
plot.expr <- cbind.data.frame(CellType = cd4.avg[, 1], CD4 = cd4.avg[, 2], CD8B = cd8b.avg[, 2])
ggplot(plot.expr, aes(x = CD4, y = CD8B, color = CellType, label = CellType)) +
    geom_point(size = 4) +
    geom_text(hjust = 0, nudge_x = 0.05) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
    ggtitle('CD4/CD8B Ratio')
ggsave(file.path(out.figs.dir, 'sp_CD4_CD8B_Ratio.scatter.pdf'), width = 6, height = 6)

#---------------------------------------------------------------
# TRDV1 vs. TRDV2

obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c('TRDV1', 'TRDV2')), CellType = Idents(obj.dn.dp))
trdv1.avg <- aggregate(TRDV1 ~ CellType, data = obj.expr, FUN = 'mean')
trdv2.avg <- aggregate(TRDV2 ~ CellType, data = obj.expr, FUN = 'mean')
plot.expr <- cbind.data.frame(CellType = trdv1.avg[, 1], TRDV1 = trdv1.avg[, 2], TRDV2 = trdv2.avg[, 2])
ggplot(plot.expr, aes(x = TRDV1, y = TRDV2, color = CellType, label = CellType)) +
    geom_point(size = 4) +
    geom_text(hjust = 0, nudge_x = 0.002) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
    ggtitle('TRDV1/TRDV2 Ratio')
ggsave(file.path(out.figs.dir, 'dn2sp_TRDV1_TRDV2_Ratio.scatter.pdf'), width = 6, height = 6)

#---------------------------------------------------------------
# alpha chain vs. beta chain vs. delta chain

sc.dn.dp.sum <- obj.dn.dp$Anno.Level.11 %>% table %>% .[rename.cells] %>% as.vector
count.tra <- table(obj.dn.dp$Anno.Level.11, obj.dn.dp$tra != '') %>% as.data.frame.matrix %>% .[rename.cells, 2]
count.trb <- table(obj.dn.dp$Anno.Level.11, obj.dn.dp$trb != '') %>% as.data.frame.matrix %>% .[rename.cells, ]

plot.expr <- cbind.data.frame(CellType = rename.cells, TRA_chain = count.tra / sc.dn.dp.sum, TRB_chain = count.trb / sc.dn.dp.sum)
ggplot(plot.expr, aes(x = TRA_chain, y = TRB_chain, color = CellType, label = CellType)) +
    geom_point(size = 4) +
    geom_text(hjust = 0, nudge_x = 0.05) +
    theme_bw(base_size = 16) +
    theme(legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
    ggtitle('TRA_chain/TRB_chain Ratio')
ggsave(file.path(out.figs.dir, 'tra_vs._trb.chains_Ratio.scatter.pdf'), width = 6, height = 6)

#------------------------------------------------------------------------------
sessionInfo()
