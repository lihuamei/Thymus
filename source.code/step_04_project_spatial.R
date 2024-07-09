suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(ROGUE))
suppressMessages(library(colorRamps))
suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(ggtree))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/trajectory_methods.R')
source('modules/findMedullaClusters.R')
source('modules/vis/step_04_project_spatial.vis.R')
source('modules/removeOutlierSpots.R')

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '4.proj_spatial/data')
out.figs.dir <- file.path('../3.results', '4.proj_spatial/figs')

dir.create(file.path('../3.results', '4.proj_spatial'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load annotated single-cell RNA-seq data and spatial data.

obj.sc <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')
obj.st <- readSPMergedData()
obj.sp.lst <- readRDS('../3.results/2.fine_anno/data/obj.st.lst.rds')

#---------------------------------------------------------------
# Corr. SC vs. SP

umiCorrSpWithSt(obj.sc, obj.sp.lst)
ggsave(file.path(out.figs.dir, 'sp_vs_sc_UMI.scatter.pdf'), width = 12, height = 6)

countSpotsPerSample(obj.sp.lst)
ggsave(file.path(out.figs.dir, 'number.spots.per.sample_bar.pdf'), width = 6, height = 4)

showAllImages(obj.sp.lst)
ggsave(file.path(out.figs.dir, 'spatial.HE.images.pdf'), width = 8, height = 4)

QCPlots(obj.sp.lst)
ggsave(file.path(out.figs.dir, 'spatial.QC.pdf'), width = 8, height = 6)

#---------------------------------------------------------------
# ScalerSP algorithm

medulla.genes <- list(c('EBI3', 'CCL17', 'CCR7', 'CSF2RB', 'CCL21', 'CCL22', 'TNFRSF18', 'CCL27', 'CXCL10', 'CXCL9', 'MS4A1', 'LAMP3'))
p.cut <- c(T1 = 1e-4, T2 = 1e-2, T3 = 1e-5, T4 = 1e-5, T5 = 1e-2, T6 = 1e-4, T7 = 1e-4, T8 = 1e-2)	
module.size <- c(T1 = 11, T2 = 8, T3 = 10, T4 = 10, T5 = 10, T6 = 10, T7 = 10, T8 = 10)
remove.spots <- list(
	T1 = c('TTGTATCACACAGAAT-1', 'TTAATGCGAGGTAACT-1', "TTGTGTATGCCACCAA-1", "TGCTCGGCGAAACCCA-1", 'AACCTTTAAATACGGT-1'),
	T2 = c('ATACCACGGGCAACTT-1', 'TGTCCTAAGTCACCGC-1', 'AGTCAAGATGACACTT-1', 'CGGGAGCTTCAGTGTA-1'),
	#T4 = c('GTCCTTCTAGTGGGTT-1', 'CACCGATACACCGAGC-1', "CTGTTGGCTCTTCTGA-1", "GTTCAAATCAGATGTC-1"),
	T5 = c('TGATTCTGTCGCCGGT-1'),
	T6 = c('AAGAGCTCTTTATCGG-1'),
	T7 = c('ATGTTGATTAGAGACT-1', 'TACAAGGGCTTCTTTA-1', 'CCCAAGTCATTACACT-1', 'CGAGAGCTTTCACTAC-1', 'ATCTTATCGCACACCC-1'),
	T8 = c('CACCCTAACAAGATCT-1', 'CGATTAAATATCTCCT-1', 'CCACACTGAGATATTA-1', 'GTGGACGCATTTGTCC-1', 'TGCCGTGGATCGTCCT-1')
)
obj.sp.lst <- findMedullaClusters(obj.sp.lst, medulla.genes, out.figs.dir = out.figs.dir, module.size = module.size, p.cut = p.cut, remove.spots = remove.spots)
obj.sp.lst <- calcSpot2ModuleDist(obj.sp.lst)

obj.sp.lst <- searchOutliers(obj.sp.lst)
obj.sp.lst <- lapply(obj.sp.lst, function(obj) { subset(obj, Outlier.Spot == 'N')})

gp.lst <- assignSpots2Centre(obj.sp.lst)
ggsave(file = file.path(out.figs.dir, 'assign2centric_spots.pdf'),  width = 12, height = 8)

cortexMeduallaCorrHclust(obj.sp.lst, feature.genes = obj.st@assays$integrated@var.features)
ggsave(file.path(out.figs.dir, 'entire_modules_hclust.pdf'), width = 10, height = 8)

cortexMedullaSpotStat(obj.sp.lst)
ggsave(file.path(out.figs.dir, 'cortex_medulla_module_spot_counts.pdf'), width = 15, height = 12)

showExampleThyScalerSp(obj.sp.lst[['T2']])
ggsave(file.path(out.figs.dir, 'show.example.ScalerSP.Thy5.pdf'), width = 8, height = 4)

showExampleThyScalerSp(obj.sp.lst[['T8']], centric = "ATAGTTCCACCCACTC-1")
ggsave(file.path(out.figs.dir, 'show.example.ScalerSP.Thy7.pdf'), width = 8, height = 4)

gps <- lapply(obj.sp.lst, function(xx) {
	SpatialFeaturePlot(xx, feature = 'medulla.score1')	
}) 

ggarrange(plotlist = gps, ncol = 3, nrow = 3)
ggsave(file.path(out.figs.dir, 'medulla.score.HEs.pdf'), width = 8, height = 8)

saveRDS(obj.sp.lst, file = file.path(out.data.dir, 'obj.sp.lst.update.rds'))

#--------------------------------------------------------------
# Scoring cell types in ST spots using marker genes.

method <- 'ks' #ks, raw.score, wkde
smooth <- TRUE
lapply(obj.sp.lst, function(obj) {
	sn.name <- obj@images %>% names
	gp.lst <- lapply(ANNO_ENTIRE_IDNET_FIG1, function(cell) {
		scoreCellTypePlots(obj, cell, smooth = smooth, method = method, adjust = 1)
	})
	ggarrange(plotlist = gp.lst, ncol = 5, nrow = 7)
	ggsave(file.path(out.figs.dir, sprintf('%s.cell.type.score.%s_smooth.pdf', sn.name, method)), height = 40, width = 30)
})

cols <- ANNO_ENTIRE_COLOR_FIG1 %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)

order.cells <- Idents(obj.sc) %>% levels
order.cells <- c(order.cells[15 : 34], order.cells[1 : 14] %>% rev)
gp <- trendLinesWithDistPlot2(obj.sp.lst, plot.tar = order.cells, degree = 10)
gp <- gp + facet_wrap(~label, scale = 'free', strip.position="top", ncol = 6, nrow = 6) +
    theme(strip.text.x = element_text(size = 14), legend.position = 'none') +
    scale_linetype_manual(values=c(rep("solid", length(Idents(obj.sc) %>% levels)))) +
    scale_color_manual(values = cols) + ylab('Score')
ggsave(file.path(out.figs.dir, sprintf('%s.to.trend.C_J_M_trend2.pdf', 'celltypes')), height = 12, width = 12)

#--------------------------------------------------------------
# Co-exists of cell types.

pdf(file.path(out.figs.dir, 'celltype_coexist.pdf'), width = 8, height = 8)
xx <- coExistsCelltypes(obj.sp.lst)
dev.off()

pdf(file.path(out.figs.dir, 'celltype_coexist_with_number.pdf'), width = 28, height = 10)
yy <- xx %>% round(., 2)
for (idx.1 in 1 : dim(xx)[2]) {
    for (idx.2 in 1 : dim(xx)[2]) {
        if (idx.1 <= idx.2) yy[idx.1, idx.2] <- '.'
    }
}
grid.table(yy)
dev.off()

#---------------------------------------------------------------
# Project annotation to spatial data using integrated methods.

obj.merged <- readRDS(file.path(out.data.dir, 'obj.sc.st.merged.rds'))
obj.merged <- RunPCA(obj.merged, seed.use = 8123456, npcs = 50)
obj.merged <- RunUMAP(obj.merged, reduction = "pca", dims = 1:30, seed.use = 8123456)
obj.merged <- FindNeighbors(obj.merged, reduction = "pca", dims = 1:30)
obj.merged <- FindClusters(obj.merged, resolution = 1.5, random.seed = 8123456)
obj.merged$SeqType <- gsub('_.*', '', rownames(obj.merged@meta.data))
obj.merged$Anno.Level.Fig.1[is.na(obj.merged$Anno.Level.Fig.1)] <- 'Spatial_spot'
saveRDS(obj.merged, file = file.path(out.data.dir, 'obj.merged.new.rds'))

obj.merged <- readRDS(file = file.path(out.data.dir, 'obj.merged.new.rds'))
DimPlot(obj.merged, reduction = "umap", pt.size = 0.3, group.by = 'SeqType')
ggsave(file = file.path(out.figs.dir, "sc.st.merged.UMAP.pdf"), width = 8, height = 8)

DimPlot(obj.merged, reduction = "umap", pt.size = 0.3, group.by = 'SeqType', split.by = 'SeqType')
ggsave(file = file.path(out.figs.dir, "sc.st.merged.split.UMAP.pdf"), width = 12, height = 8)

#-------------------------------------------------------------
# Search patterns of cell types.

require(ape)
cell.patterns <- searchAllTrendPatterns(obj.sp.lst, query.genes = ANNO_ENTIRE_IDNET_FIG1, fine.zero = FALSE)
hc <- as.dist(1 - cor(cell.patterns)) %>% upgma(., method = "ward.D")
celltypePatterns(hc)
ggsave(file = file.path(out.figs.dir, "celltype.patterns.along.sp.dist_ggtree.pdf"), width = 6, height = 6)

#--------------------------------------------------------------
# Estimate cellular proportions based on integrated single-cell RNA-seq and spatial data.

pseud.coord.lst <- readRDS(file.path(out.data.dir, 'pseud.coord.lst.rds'))
lapply(names(pseud.coord.lst), function(sn.name) {
	print(sn.name)
	figs.dir <- file.path(out.figs.dir, 'sc.project')
	dir.create(figs.dir, showWarnings = FALSE)
	obj <- obj.sp.lst[[sn.name]]
	pseud.coord <- pseud.coord.lst[[sn.name]]
	projCellTypeToSpatial(pseud.coord, obj)
	ggsave(file = file.path(figs.dir, sprintf("%s.sc.project.sp.UMAP.pdf", sp2SampleName(sn.name))), width = 8, height = 6)
})

coExistAnalysis(pseud.coord.lst, focus.ct = c('mTEC', 'DC'))
ggsave(file.path(out.figs.dir, 'coexist.celltypes.pdf'), width = 10, height = 6)

# SC vs. SP
corrSpVsSCMappedCtypes(pseud.coord.lst, obj.sc)
ggsave(file = file.path(out.figs.dir, "sp.vs.sc.celltype.propotion.pdf"), width = 6, height = 6)

# Project cell type images
sn.name <- 'T8'
obj <- obj.sp.lst[[sn.name]]
centric.spots <- obj$Assign.Centric %>% unique
centric.coord <- Seurat::GetTissueCoordinates(object = obj@images[[1]])[centric.spots, ]
pseud.coord <- searchNearestNeigCentricSpots(centric.coord, pseud.coord.lst[[sn.name]])
show.centric <- "ATAGTTCCACCCACTC-1"
pseud.coord.sub <- subset(pseud.coord, Centric == show.centric)

projCellTypeToSpatial(pseud.coord.sub[sample(rownames(pseud.coord.sub), 3000), ], title = 'Thy7', size = 3, shape = 20) + theme(legend.position = 'none')
ggsave(file.path(out.figs.dir, 'ours.example.enlarge.pdf'), width = 6, height = 6)

saveRDS(pseud.coord.sub, file = file.path(out.data.dir, 'pseud.coord.sub.example.rds'))

# Distance distribution of cell-types 
res.dist <- lapply(names(obj.sp.lst), function(sn.name) {
	obj <- obj.sp.lst[[sn.name]]
	centric.spots <- obj$Assign.Centric %>% unique
	centric.coord <- Seurat::GetTissueCoordinates(object = obj@images[[1]])[centric.spots, ]
	pseud.coord <- searchNearestNeigCentricSpots(centric.coord, pseud.coord.lst[[sn.name]])	
}) %>% do.call(rbind, .)

distOfCellTypes(res.dist)
ggsave(file.path(out.figs.dir, 'cell.type.distance.spatial_boxplot.pdf'), width = 8, height = 4)

# 
res.bak <- showSpatialCellTypes(obj.sp.lst, pseud.coord.lst, c('mTEC', 'B_memory', 'Mono'), 'T8', sp.size = c(3, 0), depend = 'mTEC')
SpatialFeaturePlot(subset(obj.sp.lst[['T8']], cells = res.bak$spot)) + res.bak$gp
ggsave(file.path(out.figs.dir, 'Fb.spatial.T8.pdf'), width = 8, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, pseud.coord.lst, c('DC', 'Treg'), 'T8', sp.size = c(3, 3), depend = 'Treg')
SpatialFeaturePlot(subset(obj.sp.lst[['T8']], cells = res.bak$spot)) + res.bak$gp
ggsave(file.path(out.figs.dir, 'DC.spatial.T8.pdf'), width = 8, height = 5)

expr.lst <- readRDS('.expr.lst.rds')
cor.lst <-  calcCorsSpSC(obj.sp.lst, expr.lst)
corDistPlot(cor.lst)
ggsave(file.path(out.figs.dir, 'corr.spot.pesudo.sc_violin.pdf'), width = 6, height = 4)

#--------------------------------------------------------------
# Proj. using Seurat

obj.seu <- readRDS(file.path(out.data.dir, 'obj.SrtCT.rds'))

predictBySeurat(obj.seu)
ggsave(file.path(out.figs.dir, 'cell.type.proj.using.seurat.pdf'), width = 12, height = 8)

corrSpVsSCMappedCtypesSeurat(obj.seu, obj.sc)
ggsave(file = file.path(out.figs.dir, "sp.vs.sc.celltype.seurat.propotion.pdf"), width = 6, height = 6)

#--------------------------------------------------------------
# CellTrek analysis

obj.celltrek <- readRDS(file.path(out.data.dir, 'obj.celltrek.rds'))
celltrekProj(obj.celltrek)
ggsave(file.path(out.figs.dir, 'celltrek.proj.celltypes.pdf'), width = 8, height = 8)

corrSpVsSCMappedCtypesCellTrek(obj.celltrek, obj.sc)
ggsave(file = file.path(out.figs.dir, "sp.vs.sc.celltype.celltrek.propotion.pdf"), width = 8, height = 6)

#---------------------------------------------------------------
# CRAD analysis

obj.card <- readRDS('../3.results/4.proj_spatial/data/CARD.obj.lst.rds')
cardPiePlot(obj.card, out.figs.dir)

corrCardPiePlot(obj.card)	
ggsave(file.path(out.figs.dir, 'corr.CARD.pdf'), height = 12, width = 12)

sc.map.lst <- lapply(obj.card, function(xx) {
	sc.mapping <- CARD_SCMapping(xx, RILST, obj.sp.lst, shapeSpot = "Square", numCell = 20, ncore = 10)		
	sc.mapping
})
saveRDS(sc.map.lst, file = file.path(out.data.dir, 'sc.map.lst.rds'))

sc.map.card <- readRDS(file = file.path(out.data.dir, 'sc.map.lst.rds'))
cardMapSc(sc.map.lst, out.figs.dir)
corrSpVsSCMappedCtypesCARD(sc.map.lst, obj.sc)
ggsave(file = file.path(out.figs.dir, "sp.vs.sc.celltype.card.propotion.pdf"), width = 8, height = 6)

assign.spots <- searchCentricSpotForCARD(sc.map.lst, obj.sp.lst)	
saveRDS(assign.spots, file = file.path(out.data.dir, 'assign.spots.for.card.rds'))

sc.map.T8 <- sc.map.lst[['T8']]
show.centric <- "ATAGTTCCACCCACTC-1"
card.coord <- as.data.frame(colData(sc.map.T8))
card.coord$Assign <- assign.spots[['T8']]
card.coord.sub <- subset(card.coord, Assign == show.centric) %>% rename(CT = 'cellType')
projCellTypeToSpatial(card.coord.sub[sample(rownames(card.coord.sub), 3000), ], title = 'Thy7', size = 3, shape = 20) + theme(legend.position = 'none')
ggsave(file.path(out.figs.dir, 'card.example.enlarge.pdf'), width = 6, height = 6)

#--------------------------------------------------------------
# 
sc.map.lst <- readRDS(file = file.path(out.data.dir, 'pseud.coord.lst.rds'))
cardMapScSplit(sc.map.lst, obj.sp.lst, out.figs.dir, label = 'Ours')

assign.spots <- searchCentricSpotForCARD(sc.map.lst, obj.sp.lst)
sc.map.T8 <- sc.map.lst[['T8']]
show.centric <- "ATAGTTCCACCCACTC-1"
card.coord <- as.data.frame(colData(sc.map.T8))	
card.coord$Assign <- assign.spots[['T8']]
card.coord.sub <- subset(card.coord, Assign == show.centric) %>% rename(CT = 'cellType')
projCellTypeToSpatial(card.coord.sub[sample(rownames(card.coord.sub), 3000), ], title = 'Thy7', size = 3, shape = 20) + theme(legend.position = 'none')
ggsave(file.path(out.figs.dir, 'proj2sp.example.enlarge.pdf'), width = 6, height = 6)

require(alluvial)

sc.map.lst.card <- readRDS(file = file.path(out.data.dir, 'sc.map.lst.rds'))
sankyPlot(sc.map.lst, sc.map.lst.card, obj.sp.lst) 
distCellTypeDot(sc.map.lst, obj.sp.lst)

cardMapScCompare(sc.map.lst, sc.map.lst.card, obj.sp.lst, out.figs.dir, ctypes = 'cTEC', label = 'Compare', xx.filter = 'CARD')	
cardMapScCompare(sc.map.lst, sc.map.lst.card, obj.sp.lst, out.figs.dir, ctypes = 'cTEC', label = 'Compare', xx.filter = 'Ours')

coExistAnalysis(pseud.coord.lst, sc.map.lst, focus.ct = NULL, out.data.dir = out.data.dir)
ggsave(file.path(out.figs.dir, 'coexist.new.celltypes.pdf'), width = 10, height = 6)

markerExpresion(sc.map.lst, obj.sp.lst, out.figs.dir, ctypes = ANNO_ENTIRE_IDNET_FIG1)
ggsave(file.path(out.figs.dir, 'ours.marker.score.celltypes.pdf'), width = 10, height = 4)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('CD8aa', 'B_memory', 'DC'), 'T8', sp.size = c(3, 0), depend = 'DC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'CD8aa.B_memory.DC.spatial.T8.pdf'), width = 6, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('CD8aa', 'B_memory', 'DC'), 'T2', sp.size = c(3, 0), depend = 'DC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'CD8aa.B_memory.DC.spatial.T2.pdf'), width = 6, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('cTEC', 'Fb'), 'T8', sp.size = c(2, 0), depend = 'cTEC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'cTEC.Fb.spatial.T8.pdf'), width = 6, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('cTEC', 'Fb'), 'T2', sp.size = c(2, 0), depend = 'cTEC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'cTEC.Fb.spatial.T2.pdf'), width = 6, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('cTEC', 'Fb', 'DP_re'), 'T8', sp.size = c(1.5, 0), depend = 'cTEC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'cTEC.DP_re.spatial.T8.pdf'), width = 6, height = 5)

res.bak <- showSpatialCellTypes(obj.sp.lst, sc.map.lst, c('cTEC', 'Fb', 'DP_re'), 'T2', sp.size = c(1.5, 0), depend = 'cTEC', bool = F)
res.bak$gp
ggsave(file.path(out.figs.dir, 'cTEC.DP_re.spatial.T2.pdf'), width = 6, height = 5)

#---------------------------------------------------------------
# Projection of T cells.

gps.lst <- tcrGroupsProj(pseud.coord.lst, obj.sc, obj.sp.lst)
lapply(names(gps.lst), function(sn.name) {
	gp <- gps.lst[[sn.name]]
	ggsave(file.path(out.figs.dir, sprintf('%s_tcell.subgroups_plots.pdf', sn.name)), plot = gp, width = 10, height = 8)		
})

#--------------------------------------------------------------
# Project using marker genes

de.markers <- read.table('../3.results/3.anno_assess/data/Anno.entire.T.de.markers.xls', sep = '\t', header = TRUE)
top20 <- de.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
markers <- split(top50$gene, top50$cluster)

obj.sp.lst <- lapply(obj.sp.lst, function(obj.st) {
	obj.st <- AddModuleScore(object = obj.st, features = markers)
	colnames(obj.st@meta.data)[(length(colnames(obj.st@meta.data)) - length(markers) + 1) : length(colnames(obj.st@meta.data))] <- names(markers)
	obj.st
}) 

obj.sp.lst <- lapply(obj.sp.lst, function(xx) {
    figs.dir <- file.path(out.figs.dir, xx@images %>% names)
    dir.create(figs.dir, showWarnings = FALSE)

    for (cell in names(markers)) {
        g1 <- SpatialPlot(xx, features = cell, image.alpha = 1, alpha = 0) + NoLegend()
        g2 <- SpatialPlot(xx, features = cell, image.alpha = 0, alpha = 1)
        g1 + g2
        ggsave(file.path(figs.dir, sprintf('%s_%s_markers.pdf', xx@images %>% names, cell)), width = 8, height = 5)
    }
    return(xx)
})
saveRDS(obj.sp.lst, file = file.path(out.data.dir, 'obj.sp.lst.proj.markers.RDS'))

#---------------------------------------------------------------
# Extract high confidence modules from spatial compartment.

spots.module.cnt <- lapply(obj.st.lst, function(xx) {
    table(xx@meta.data$Assign.Centric, xx@meta.data$HE.Labels) %>% as.data.frame.matrix %>% `rownames<-`(paste0(xx@images %>% names, '_', rownames(.)))
}) %>% do.call(rbind, .) %>% .[, 1 : 4]
count.df <- cbind.data.frame(Cortex = spots.module.cnt[, 'Cortex'], Medulla = rowSums(spots.module.cnt[, c('Medulla_edge', 'Medulla_hi')]))
rownames(count.df) <- gsub('T.*\\.', '', rownames(count.df))
count.df <- count.df[order(rownames(count.df)), ]

plot.df <- sweep(count.df, 1, rowSums(count.df), '/') %>% cbind.data.frame(., Label = rownames(.))
ggplot(tidyr::gather(plot.df, 'key', 'value', -Label), aes(x = Label, y = value, fill = key)) + geom_bar(position="stack", stat="identity") +
    theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab('Fraction') +
    scale_fill_manual(values = CM_COLORS %>% unlist %>% .[c('Cortex', 'Medulla')]) + xlab('') +
    geom_hline(yintercept=plot.df$Medulla %>% mean, linetype="dashed", color = "yellow") +
    scale_y_continuous(breaks = c(0.00, 0.18, 0.25, 0.50, 0.75, 1.00), expand = c(0, 0))
ggsave(file.path(out.figs.dir, 'cortex_medulla_module_relative_fraction.pdf'), width = 15, height = 7)

#-----------------------------------------------------------
# Number of spots in each compartment

plot.df.2 <- cbind.data.frame(count.df, Summary = rowSums(count.df), Label = rownames(count.df))
ggplot(plot.df.2, aes(x = Label, y = Summary, fill = 'blue', alpha = 0.5)) + geom_bar(stat = "identity") +
    geom_text(aes(label = Summary, y = Summary), vjust = -.5) +
    theme_classic(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'none'
    ) + ylab('Number of spots in each compartment') + xlab('') +
    scale_y_continuous(expand = c(0, 0), limit = c(0, 750), breaks = c(0, floor(median(plot.df.2$Summary)), 250, 500, 750)) +
    geom_hline(yintercept=floor(median(plot.df.2$Summary)), linetype="dashed", color = "blue")
ggsave(file.path(out.figs.dir, 'cortex_medulla_module_spot_counts.pdf'), width = 15, height = 8)

#-----------------------------------------------------------
# Correlation across compartments

obj.st.lst.sub <- lapply(obj.st.lst, function(xx) {
    xx <- subset(xx, HE.Labels != "Medulla_lo")
    xx@meta.data$Class <- xx@meta.data$HE.Labels
    xx@meta.data$Class[xx@meta.data$Class %in% c('Medulla_edge', 'Medulla_hi', "Medulla_centric")] <- 'Medulla'
    xx@meta.data$Class[xx@meta.data$Class %in% c("Cortex")] <- 'Cortex'
    xx@meta.data$Class.new <- paste0(xx@meta.data$Assign.Centric, '_', xx@meta.data$Class)
    xx
})

expr.corr <- lapply(obj.st.lst.sub, function(xx) {
    AverageExpression(xx, group.by = 'Class.new')$Spatial
}) %>% do.call(cbind, .) %>% cor(.)

anno.df <- lapply(obj.st.lst.sub, function(xx) {
    cbind.data.frame(Class.new = xx@meta.data$Class.new, Class = xx@meta.data$Class, Centric = xx@meta.data$Assign.Centric, Sample = xx@images %>% names) %>% .[!duplicated(.[, 1]), ]
}) %>% do.call(rbind, .) %>% `rownames<-`(.[, 'Class.new']) %>% .[colnames(expr.corr), ]

hc <- as.dist(1 - expr.corr) %>% upgma(., method="ward.D")
gt <- ggtree(hc, layout="circular", size = 0.5)
gt <- gt %<+% anno.df + geom_tippoint(aes(color = Class), alpha = 0.5, size = 3) + scale_color_manual(values = CM_COLORS[c('Cortex', 'Medulla')]) +
    guides(colour = guide_legend(override.aes = list(size = 8, title = 'HE.Position')))
gt1 <- gheatmap(gt, anno.df[, "Sample", drop = FALSE], width = 0.1, color = "black") + scale_fill_manual(values = SAMPLE_COLORS) +
    guides(colour = guide_legend(override.aes = list(title = 'SampleID')))
ggsave(file.path(out.figs.dir, 'entire_modules_hclust.pdf'), width = 10, height = 8)

