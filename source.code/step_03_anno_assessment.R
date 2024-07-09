suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(ROGUE))
suppressMessages(library(colorRamps))
suppressMessages(library(tidyverse))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/trajectory_methods.R')
source('modules/vis/step_03_anno_assessment.vis.R')

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '3.anno_assess/data')
out.figs.dir <- file.path('../3.results', '3.anno_assess/figs')

dir.create(file.path('../3.results', '3.anno_assess'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load data and Gene stats

cols <- ANNO_ENTIRE_COLOR_FIG1 %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
merged.obj <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')
merged.obj <- mappCellNames(merged.obj, 'Tuft', 'Fb', 'Anno.Level.Fig.1', 'Anno.Level.Fig.1')
merged.obj@meta.data$Anno.Level.Fig.1 <- factor(merged.obj@meta.data$Anno.Level.Fig.1, levels = ANNO_ENTIRE_IDNET_FIG1)
Idents(merged.obj) <- merged.obj@meta.data$Anno.Level.Fig.1
saveRDS(merged.obj, file = '../3.results/2.fine_anno/data/merged.obj.anno.new.rds')

numberOfGeneDetecedPerCellType(merged.obj)
ggsave(file.path(out.figs.dir, 'number.of.detected.genes.per.celltype.pdf'), width = 12, height = 5)

numberOfDetectedUMIs(merged.obj)
ggsave(file.path(out.figs.dir, 'number.of.detected.UMIs.per.celltype.pdf'), width = 12, height = 5)

#--------------------------------------------------------------
# rogue analysis

if (!file.exists(file.path(out.data.dir, 'rogue.res.rds'))) {
	rogue.res <- rogue(GetAssayData(merged.obj), labels = merged.obj@meta.data$Anno.Level.Fig.1, samples = merged.obj@meta.data$orig.ident, platform = "UMI", span = 0.6)
	saveRDS(rogue.res, file = file.path(out.data.dir, 'rogue.res.rds'))
} else {
	rogue.res <- readRDS(file.path(out.data.dir, 'rogue.res.rds'))
}
rogueSocrePerCellType(rogue.res)	
ggsave(file.path(out.figs.dir, 'rogue.boxplot.pdf'), width = 10, height = 5)

#--------------------------------------------------------------
# Logistic regression

obj <- merged.obj
sample.size <- table(Idents(obj)) %>% min %>% {. * 0.5} %>% as.integer
lr.model <- trainLRmodel(obj, ANNO_ENTIRE_IDNET_FIG1, sample.size = sample.size, prefix = 'Anno.Fig1.LR', out.figs.dir = out.figs.dir, wt = 10, ht = 7, color = cols)
pred.res <- lr.model$predict(obj %>% GetAssayData(.) %>% t) %>% `names<-`(colnames(obj))
probs <- table(pred.res, Idents(obj)) %>% as.matrix %>% sweep(., 1, rowSums(.), '/') %>% .[colnames(.), ]

pdf(file.path(out.figs.dir, 'sc.anno_logistic_predicted_similarity.pdf'), width = 8, height = 8)
pheatmap(probs, cluster_rows=F, cluster_cols=F, color = colorRampPalette (c ('grey80', "white", "red")) (100), border_color=NA)
dev.off()

#-------------------------------------------------------------
# Hcluster analysis

avg.expr <- AverageExpression(merged.obj)$RNA
cv <- apply(avg.expr, 1, function(x) sd(x) / mean(x))
avg.expr.sub <- cv[order(cv) %>% rev][1 : 15000] %>% names %>% avg.expr[., ]

pdf(file.path(out.figs.dir, 'sc.anno_hclust.pdf'), width = 8, height = 8)
pheatmap(cor(avg.expr.sub), clustering_method = 'ward.D', border_color=NA, color = colorRampPalette (c('grey80', "white", "red")) (100))
dev.off()

#--------------------------------------------------------------
# Correlation between ours vs. science

science.markers <- read.table('../4.extdata/science_thymus_markers.csv', sep = ',', header = T)
expr.refs <- read.table('../4.extdata/science_thymus_gene_expression.xls', sep = '\t', header = TRUE) %>% .[!duplicated(.[, 1]), ] %>% `rownames<-`(.[, 1]) %>% .[, -1]
expr.ours <- AverageExpression(merged.obj)$RNA %>% { log2(1 + .) }

ovp.genes <- intersect(rownames(expr.refs), rownames(expr.ours))
marker.genes <- unlist(science.markers, use.names = FALSE) %>% intersect(ovp.genes, .)
corr.ref <- apply(expr.ours[marker.genes, ], 2, function(x) {
	apply(expr.refs[marker.genes, ], 2, function(y) {
    	cor(x, y)
    })
})

pdf(file.path(out.figs.dir, 'sc.anno_ours.vs.science.pdf'), width = 8, height = 8)
pheatmap(corr.ref, clustering_method = 'ward.D', border_color=NA, color = colorRampPalette (c('grey80', "white", "red")) (100))
dev.off()

#-----------------------------------------------------------
# R(O/E)

obj <- merged.obj
obj$AgeStatusNew <- obj$AgeStatus %>% as.vector
obj$AgeStatusNew[obj$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

N.o <- table(obj$AgeStatusNew, obj$Anno.Level.Fig.1) %>% t
res.chisq <- chisq.test(N.o)
R.oe <- (res.chisq$observed)/(res.chisq$expected)

ROEPlot(R.oe)
ggsave(file.path(out.figs.dir, 'roe.indedx.across.celltypes_Dotplot.pdf'), width = 10, height = 3.5)

#-----------------------------------------------------------
# R(O/E) for science data

obj.science <- readLoadScience()
obj <- obj.science
obj$Anno_level_fig1 <- factor(obj$Anno_level_fig1, levels = order.cell.science)
N.o <- table(obj$AgeStatus, obj$Anno_level_fig1) %>% t
res.chisq <- chisq.test(N.o)
R.oe <- (res.chisq$observed)/(res.chisq$expected)

ROEPlot(R.oe)
ggsave(file.path(out.figs.dir, 'roe.indedx.science.across.celltypes_Dotplot.pdf'), width = 12, height = 4)

#------------------------------------------------------------
# Project to UMAP

obj <- merged.obj	
obj$AgeStatusNew <- obj$AgeStatus %>% as.vector
obj$AgeStatusNew[obj$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

markers <- readRDS(file.path('../3.results/2.fine_anno/data', 'markers.save.rds'))
de.markers <- readRDS(file.path('../3.results/3.anno_assess/data', 'sc.markers.anno.entire.rds'))
var.genes <- de.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) %>% .[, 'gene'] %>% unlist %>% as.vector
xx <- read.table(SCIENCE_MARKERS, sep = ',', header = T) %>% as.list
obj.sub <- obj[c(var.genes, xx$CD8.T, xx$CD8.Tmem, xx$CD4.T, xx$CD4.Tmem, unlist(xx, use.names = FALSE)) %>% unique, ]
obj.sub@meta.data$Anno.Level.Fig.1 <- as.vector(obj.sub@meta.data$Anno.Level.Fig.1)
saveAsH5adFile(obj.sub, 'merged.obj.UMAP', out.data.dir)

system('python anno_UMAP.py')
file.remove(file.path(out.data.dir, 'merged.obj.UMAP.h5ad'))
file.remove(file.path(out.data.dir, 'merged.obj.UMAP.h5Seurat'))

# Show UMAP with Seurat.
Convert('../3.results/3.anno_assess/data/merged.obj.UMAP.Scanpy.h5ad', ".merged.obj.UMAP.Scanpy.h5seurat", overwrite = TRUE)
obj <- LoadH5Seurat(".merged.obj.UMAP.Scanpy.h5seurat")
obj@meta.data$Anno.Level.Fig.1 <- factor(obj@meta.data$Anno.Level.Fig.1, levels = ANNO_ENTIRE_IDNET_FIG1)
obj <- mappCellNames(obj, obj$Anno.Level.Fig.1 %>% levels, as.character(1 : length(ANNO_ENTIRE_IDNET_FIG1)), 'Anno.Level.Fig.1', 'Anno.Level.Cluster')
annoUMAPRes(obj)
ggsave(file.path(out.figs.dir, 'thymus.anno_UMAP.pdf'), width = 12, height = 6)

# Split by age status.
obj.1 <- merged.obj
obj.1$AgeStatusNew <- obj.1$AgeStatus %>% as.vector
obj.1$AgeStatusNew[obj.1$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'
obj$AgeStatusNew <- obj.1$AgeStatusNew
obj$AgeStatusNew <- factor(obj$AgeStatusNew, levels = c("Prental", "Children", "Adult", "Older"))
DimPlot(obj, group.by = 'AgeStatusNew', cols = AGE_STATUS_EXTENT_COLOR[1 : 4], pt.size = 0.5)
ggsave(file.path(out.figs.dir, 'thymus.anno.splitby.agestatus_UMAP.pdf'), height = 6, width = 8)

#------------------------------------------------
# Supplementary figures

ovp.markers <- validateMarkers(de.markers, label = 'major', out.figs.dir = out.figs.dir)
marker.genes <- list(
	Ery = c("HBG1", "HBG2"),
	B_naive = c("CD19", "CD79A", "IGHD"),
	B_trans = c("MS4A1", "CD24"),
	B_memory = c("IGHA1", 'CD27'),
	Plasma = c("IGHG1", "IGHG2"),
	Mono = c("S100A8", "S100A9"),
	Mac = c("C1QA", "C1QC"),
	DC = c("HLA-DPA1", "HLA-DQB1"),
	pDC = c("IL3RA", "LILRA4"),
	Fb = c("DCN", "COL1A1"),
	Fb_cycling = c("UBE2C", "TOP2A"),
	Endo = c("ACKR1", "PLVAP"),
	VSMC = c("TAGLN", "ACTA2"),
	Lymph = c("TFF3", "NTS"),
	cTEC = c("CCL25", "GNG11"),
	mTEC = c("KRT19", "KRT17"),
	DN_early = c("PTPRC", "CD3D", "CD3G", "CD4", "CD8A", "CD8B", "TRDC", 'TRAC', 'TRBC1', "IGLL1", "SMIM24"),
	DN_blast = c("TOP2A", "TYMS"),
	DN_re = c("PTCRA", 'RAG1', 'RAG2'),
	DP_blast = c("TOP2A", "UBE2C"),
	DP_re = c("RAG1", "RAG2"),
	`abT(entry)` = c("TOX2", "CCR9"),
	CD4T = c("CD4", "SELL"),
	CD4T_mem = c("IL7R", "CCR7"),
	Treg.diff = c("FOXP3", "IFITM1"),
	Treg = c("FOXP3", "IL32"),
	CD8T = c("GZMM", "ABLIM1"),
	CD8T_mem = c("CCL4", "GZMK"),
	CD8aa = c("CD27", "ID3"),
	T_agonist = c("TNFRSF1B", "LCP2"),
	T_apoptosis = c("BCL2L11"),
	T_proliferating = c("TYMS", "UBE2T"),
	NKT = c("NKG7", "KLRB1"),
	ILC3 = c("IL4I1", "SLC16A3", "TNFRSF25")
)
bubbleMarkers(merged.obj, marker.genes)
ggsave(file.path(out.figs.dir, 'marker.genes.bubble.pdf'), width = 10, height = 15)

#-------------------------------------
# UMI, Genes, and Mito percentage plots

gex_ids <- list.files('../1.data/cellranger.results/RNAseq')
source("../transfer/XQ/data_sc_st/add.tcr.funs.new.R")
datapath <- '../1.data/cellranger.results/RNAseq/'
seurat.obj.list <- list()

for(i in 1 : length(gex_ids)){
    data.dir = file.path(datapath, gex_ids[i], "filtered_feature_bc_matrix")
    counts <- Read10X(data.dir = data.dir)
    project <- paste("Thy", i, sep = "")
    seurat.obj <- CreateSeuratObject(counts = counts, project = project, min.cells = 3, min.features = 200)
    seurat.obj.list[[i]] <- seurat.obj
}
merged.obj.tumor <- merge(seurat.obj.list[[1]], y = seurat.obj.list[-1], add.cell.ids = gex_ids, project = "Thymus")

source("../transfer/XQ/data_sc_st/seurat_QC.R")
counts <- merged.obj.tumor@assays$RNA@counts
UMIPlot(counts, out.name = "total_umi_distr.pdf", out.dir = out.figs.dir)
GenePlot(counts, min.cut = 200, max.cut = 4500, out.name = "total_genes_distr.pdf", out.dir = out.figs.dir)
features <- rownames(counts)
pattern = "^MT-"
mito.features <- grep(pattern, features, value = T)
MitoGenePlot(counts, mito.features, max.cut = 0.2, out.name = "mito_genes_distr.pdf", out.dir = out.figs.dir)
pattern = "^RP[SL]"
ribo.features <- grep(pattern, features, value = T)
MitoGenePlot(counts, ribo.features, max.cut = 0.5, out.name = "ribo_genes_distr.pdf", out.dir = out.figs.dir)

mgenePlotCellsPerGene(counts, out.figs.dir, color = 'darkgreen', main.title = 'Gene expressed per cell')
mcellPlotUmisPerCell(counts, out.figs.dir, main.title = 'Total UMIs per cell')
mcellPlotMitoProp(counts, out.figs.dir, patten = '^MT-', max.cut = 0.2, color = 'darkorchid1')
mcellPlotMitoProp(counts, out.figs.dir, patten = '^RP[SL]', color = 'cyan4', show.cutoff = T, main.title = 'Fraction of ribosomal gene expression per cell', file.name = 'total_ribosome_distr.pdf', max.cut = 0.5, scale2 = -0.2)

#--------------------------------------------------------
# Major Cell types.

merged.obj$Anno.Level.15 <- merged.obj$Anno.Level.1
obj.sub <- subset(merged.obj, Anno.Level.1 == 'T')
merged.obj.1 <- updateAnnoRes(obj.sub$Anno.Level.0, merged.obj, "Anno.Level.15", "Anno.Level.15")
merged.obj <- mappCellNames(
	merged.obj.1,
	c("DP4", "CD8aa", "abT(entry)2", 'DP2', 'NKT', 'abT(entry)1', 'CD8+ T', 'DN', "DP3", "DP1", "CD4+ T", "Treg"),
	c('DP_re', 'SP', 'abT(entry)', 'DP_blast', 'SP', 'abT(entry)', 'SP', 'DN', 'DP_re', 'DP_blast', 'SP', 'SP'),
	'Anno.Level.15',
	'Anno.Level.16'
)

DimPlot(merged.obj, group.by = 'Anno.Level.16', label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'major.anno.res_UMAP.pdf'), width = 8, height = 8)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.B.rds'))
DimPlot(obj, label = T) + + NoLegend()
ggsave(file.path(out.figs.dir, 'B.anno.res_UMAP.pdf'), width = 6, height = 6)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.Myeloid.rds'))
DimPlot(obj, label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'Myeloid.anno.res_UMAP.pdf'), width = 6, height = 6)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.TEC.rds'))
DimPlot(obj, label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'TEC.anno.res_UMAP.pdf'), width = 6, height = 6)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.Stromal.rds'))
DimPlot(obj, label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'Stromal.anno.res_UMAP.pdf'), width = 6, height = 6)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.Other.T.rds'))
DimPlot(obj, label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'Other.T.anno.res_UMAP.pdf'), width = 6, height = 6)

obj <- readRDS(file.path('../3.results/2.fine_anno/data', 'obj.DN.rds'))
new.cluster <- c('DN_re', 'DN_re', 'DN_re', 'DN_blast', 'DN_re', 'DN_blast', 'DN_re', 'DN_early', 'DN_re')
names(new.cluster) <- levels(obj)
obj <- RenameIdents(obj, new.cluster)
DimPlot(obj, label = T) + NoLegend()
ggsave(file.path(out.figs.dir, 'DN.T.anno.res_UMAP.pdf'), width = 6, height = 6)
