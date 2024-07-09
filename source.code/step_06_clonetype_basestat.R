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

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/findMedullaClusters.R')
source('modules/tcrSubgroup.R')
source('modules/vis/step_08_clonetype_basestat_vis.R')

#--------------------------------------------------------------
# Create output directories

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
coding.genes <- subset(anno, gene_type == 'protein_coding') %>% .[, 'gene_name'] %>% intersect(., rownames(obj.st.lst[[1]]))
obj.merged.sc <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')

merged.obj <- obj.merged.sc %>% tcrSubgroup(trab.cut = 1, trd.cut = 1) 
tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")

#---------------------------------------------------------------
# Statistics of T cell subtypes 

trg.stat <- trgStat(merged.obj)
baseStatOfEightSubgroups(merged.obj, trg.stat)
ggsave(file.path(out.figs.dir, 't.cell.subgroup.stat.Bar_with_Pie.pdf'), width = 6, height = 6)

# Validate using science thymus dataset
obj.science.t <- LoadH5Seurat("../4.extdata/HTA08.v01.A06.Science_human_tcells_new.h5seurat")
obj.science.t$tcr <- 1
obj.science.t <- obj.science.t %>% tcrSubgroup(trab.cut = 1, trd.cut = 1)
trg.stat.ref <- trgStat(obj.science.t)

#---------------------------------------------------------------
# four subtypes percentage.

pdf(file.path(out.figs.dir, 'four.subgroups.percentage_Pie.pdf'), width = 6, height = 6)
vennDiagForTSubgroups(merged.obj)
dev.off()

pdf(file.path(out.figs.dir, 'four.subgroups.percentage.science_Pie.pdf'), width = 6, height = 6)
vennDiagForTSubgroups(obj.science.t)
dev.off()

dn.sort <- readRDS('../4.extdata/GSE195812_Sort1_DN123.rds')
isp.dp2sp.sort <- readRDS('../4.extdata/GSE195812_Sort2_ISP_DP_SP.rds')
dn.sort <- dn.sort[intersect(rownames(dn.sort), rownames(isp.dp2sp.sort)), ]
isp.dp2sp.sort <- isp.dp2sp.sort[intersect(rownames(dn.sort), rownames(isp.dp2sp.sort)), ]

merged.extra.objs <- merge(dn.sort, y = isp.dp2sp.sort, project = "T") 
merged.extra.objs <- merged.extra.objs[, !is.na(merged.extra.objs$CDR3)]

merged.extra.objs$tcr <- 1
merged.extra.objs <- merged.extra.objs %>% tcrSubgroup
pdf(file.path(out.figs.dir, 'four.subgroups.percentage.bioxv_Pie.pdf'), width = 6, height = 6)
vennDiagForTSubgroups(merged.extra.objs)
dev.off()

#---------------------------------------------------------------- 
# four subtypes percentage of each sample.

fourSubgroupsPerSample(merged.obj, tcr.subtypes)
ggsave(file.path(out.figs.dir, 'four.subgroups.percentage.samples_Barplot.pdf'), width = 10, height = 6)

obj <- subset(merged.obj, t.cell_group %in% tcr.subtypes)
obj$AgeStatusNew <- obj$AgeStatus %>% as.vector
obj$AgeStatusNew[obj$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

N.o <- table(obj$AgeStatusNew, obj$t.cell_group) %>% t
res.chisq <- chisq.test(N.o)
R.oe <- (res.chisq$observed)/(res.chisq$expected)
R.oe <- R.oe[tcr.subtypes, AGE_STATUS_EXTENT[1 : 4]]

ROEPlot(R.oe) + scale_size_continuous(limits = c(0,2))
ggsave(file.path(out.figs.dir, 'roe.indedx.across.subgroups_Dotplot.pdf'), width = 8, height = 3.5)

#----------------------------------------------------------------
# four subtypes compare with age status

corrBetweenAgeAndPropSubgroup(merged.obj)	
ggsave(file.path(out.figs.dir, 'four.subgroups.percentage.vs.age_Scatter.pdf'), width = 7, height = 6)

#----------------------------------------------------------------
# Number of clonotype and number of unique clonotypes.

numberOfClonetypes(merged.obj)
ggsave(file.path(out.figs.dir, 'Number.of.clonotypes.in.each.sample_barplot.pdf'), width = 14, height = 6)

# Clone size distribution
obj <- merged.obj
Idents(obj) <- obj$orig.ident
divs.df <- mclapply(1 : 100, function(idx) {
    print(idx)
    obj.sub = downSamplSeurat(obj, cnt = 1000, seed = idx)
    divers.tcr <- lapply(merged.obj$orig.ident %>% levels, function(sn) {
        obj.sub <- subset(obj.sub, orig.ident == sn)
        pp <- obj.sub$tcr_frequency / sum(obj.sub$tcr_frequency)
        sum(pp * log2(pp)) * (-1)
    }) %>% unlist %>% as.vector %>% `names<-`(names(tcr.count))
}, mc.cores = 40) %>% do.call(rbind, .)

tcrDivShannonEntro(divs.df)
ggsave(file.path(out.figs.dir, 'tcr.diversity.along.with.age_boxplot.pdf'), width = 5, height = 6)

#------------------------------------------------------------------
# Clonotypes vs. age status.

obj <- merged.obj
Idents(obj) <- obj$orig.ident

plot.df <- parallel::mclapply(1 : 100, function(idx) {
	print(idx)
	obj.sub <- downSamplSeurat(obj, cnt = 3000, seed = idx)	
	obj.sub <- subset(obj.sub, tcr == 1)
	cnt.tcr.unq <- lapply(obj.sub$orig.ident %>% levels, function(xx) subset(obj.sub, orig.ident == xx)$tcr_cdr3s_aa %>% unique %>% length) %>% unlist
	plot.df <- cbind.data.frame(Count = cnt.tcr.unq, Type = c(rep('Prental', 4), rep('Children', 8), rep('Adult', 2), rep('Older', 2)), X = 1 : 16)
}, mc.cores = 10) %>% do.call(rbind, .)

tcrDiversityPermutation(plot.df)
ggsave(file.path(out.figs.dir, 'uniq.tcr.from.prental2adult_scatterplot.pdf'), width = 7, height = 6)

#------------------------------------------------------------------
# Distribution of clonotype abundance.

cloneSizePlot(merged.obj)
ggsave(file.path(out.figs.dir, 'clone.size.dist_pieplot.pdf'), width = 8, height = 8)

merged.obj$CloneSize <- merged.obj$tcr_frequency
merged.obj$CloneSize[merged.obj$CloneSize >=3] <- '>=3'

obj <- merges.obj
Idents(obj) <- obj$CloneSize
de.markers <- FindAllMarkers(obj, only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
top20 <- de.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% data.frame
cloneSizeProjectSp(obj.st.lst, top20)
ggsave(file.path(out.figs.dir, 'clone.size.project.spatial.pdf'), width = 15, height = 8)

obj.st.lst <- lapply(obj.st.lst, function(xx) {
    xx <- AddModuleScore(xx, features = split(top20$gene, top20$cluster))
})

trendLinesWithDistPlot1(obj.st.lst, plot.tar = c('Cluster1', 'Cluster2', 'Cluster3')) + scale_color_manual(values = c('cornsilk2', 'darkorange', 'darkred'))
ggsave(file.path(out.figs.dir, 'clone.size.fit.spatial.pdf'), width = 6, height = 6)

#-------------------------------------------------------------------------------
# TRA, TRB, and TRAB chain.

obj.bak <- LoadH5Seurat(".merged.obj.UMAP.Scanpy.h5seurat")
uamp.new <- obj.bak[["umap"]]@cell.embeddings %>% `colnames<-`(c('umapn_1', 'umapn_2'))
obj.new <- obj.merged.sc

uamp.new.sub <- uamp.new[intersect(colnames(obj.new), rownames(uamp.new)), ]
obj.new <- obj.new[, rownames(uamp.new.sub)]
obj.new[["umapn"]] <- CreateDimReducObject(embeddings = uamp.new.sub, key = "UMAPN_", assay = DefaultAssay(obj.new))
obj <- subset(obj.new, Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17 : 33])

obj <- subset(obj, tcr == 1)	
obj <- subset(obj, Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17 : 33])
trb.seq <- obj$tcr_cdr3s_aa %>% strsplit(., ';') %>% lapply(., function(xx) { xx %>% grep('TRB:', .)  %>% xx[.] %>% paste0(., collapse = ';') })
obj$trb <- trb.seq %>% unlist(.)
tra.seq <- obj$tcr_cdr3s_aa %>% strsplit(., ';') %>% lapply(., function(xx) { xx %>% grep('TRA:', .)  %>% xx[.] %>% paste0(., collapse = ';') })
obj$tra <- tra.seq %>% unlist(.)

alphaBetaChainUMAP(obj)
ggsave(file.path(out.figs.dir, 'alpha.beta.rearrangement.T.subsets.UMAP.pdf'), width = 8, height = 6)

tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
obj <- obj %>% tcrSubgroup(trd.cut = 1) %>% subset(., t.cell_group %in% tcr.subtypes) 
cells.1 <- subset(obj, Anno.Level.Fig.1 == 'DN_re' & tcr.subtypes == "TRA-TRB+TRD-") %>% Cells

obj.1 <- subset(obj, Anno.Level.Fig.1 %in% c(ANNO_ENTIRE_IDNET_FIG1[17:23], 'CD8T'))
Idents(obj.1) <- obj.1$t.cell_group
levels(obj.1) <- tcr.subtypes
obj.sub <- downSamplSeurat(obj.1, cnt = 2000, seed = 1234)
DimPlot(obj.sub, group.by = 't.cell_group', cols = T.CELL.GROUPS.COLORS, pt.size = 0.3) + 
	DimPlot(obj.sub, group.by = 'Anno.Level.Fig.1', cols = ANNO_ENTIRE_COLOR_FIG1[17 : 33], pt.size = 0.3)

ggsave(file.path(out.figs.dir, 'four.T.subsets.UMAP.pdf'), width = 12, height = 6)

#-------------------------------------------------------------------------------
# Mitochondrial score of T cell subsets.

mt.genes <- rownames(obj)[grep('^MT-', rownames(obj))]
obj <- AddModuleScore(obj, features = list(xx = mt.genes), name = 'mito.score')
FeaturePlot(obj, feature = 'mito.score1', col = c("grey50", "lightgrey", "red"), pt.size = 0.5)
ggsave(file.path(out.figs.dir, 'mito.score.T.subsets.UMAP.pdf'), width = 8, height = 6)

#-------------------------------------------------------------------------------
# TRD gene expression

trd.genes <- rownames(sc.dn.dp)[grep('^TRD[VJC]', rownames(sc.dn.dp))]
sc.dn.dp <- AddModuleScore(sc.dn.dp, features = list(xx = trd.genes), name = 'trd.score')
FeaturePlot(sc.dn.dp, feature = 'trd.score1', col = c("grey50", "lightgrey", "red"), pt.size = 0.5)

#-------------------------------------------------------------------------------
# TOP2A and RAG1

expr <- FetchData(obj, vars = c('TOP2A', 'RAG1'))
obj$P.Q <- log2(exp(expr[, 2]) / exp(expr[, 1]))
FeaturePlot(obj, feature = 'P.Q', col = CustomPalette(low = "#209CA0", high = "red", mid = "white", k = 50), pt.size = 0.3) + ggtitle('Proliferation / Recombination')
ggsave(file.path(out.figs.dir, 'prolifer.recombin.T.subsets.UMAP.pdf'), width = 8, height = 6)

#-----------------------------------------------------------------
# TRB share between samples

obj <- merged.obj
Idents(obj) <- obj$orig.ident
obj.sub <- downSamplSeurat(obj, cnt = 2000)

obj.list <- SplitObject(obj.sub, split.by = 'orig.ident') %>% .[names(SampleInfos)]
clone.mat <- lapply(obj.list, function(obj.1) {
    trb.1 <- obj.1$trb
    jacard.idx <- lapply(obj.list, function(obj.2) {
        trb.2 <- obj.2$trb
        intersect(trb.1, trb.2) %>% length
    }) %>% as.data.frame
}) %>% do.call(rbind, .)

for(i in 1 : ncol(clone.mat)){
	for (j in 1 : ncol(clone.mat))
		if (i == j) {
			clone.mat[i, j] <- NA
	}
}
cols <- colorRampPalette(c("#67001f", "#d6604d", "#f7f7f7", "#4393c3", "#053061") %>% rev )(1024)
pheatmap(
	log2(1 + clone.mat), 
	color = cols, 
	cluster_rows = T, 
	cluster_cols = T, 
	display_numbers = clone.mat, 
	main = 'Repertoire overlap', 
	legend_labels = c('4', '16', '64', '256'),
	filename = file.path(out.figs.dir, 'trb.share.in.samples_heatmap.pdf'),
	height = 8,
	width = 8
)

#------ TRA
clone.mat <- lapply(obj.list, function(obj.1) {
    trb.1 <- obj.1$tra
    jacard.idx <- lapply(obj.list, function(obj.2) {
        trb.2 <- obj.2$tra
        intersect(trb.1, trb.2) %>% length
    }) %>% as.data.frame
}) %>% do.call(rbind, .)

for(i in 1 : ncol(clone.mat)){
    for (j in 1 : ncol(clone.mat))
        if (i == j) {
            clone.mat[i, j] <- NA
    }
}

pheatmap(
    log2(1 + clone.mat),
    color = cols,
    cluster_rows = T,
    cluster_cols = T,
    display_numbers = clone.mat,
    main = 'Repertoire overlap',
    legend_labels = c('4', '16', '64', '256'),
    filename = file.path(out.figs.dir, 'tra.share.in.samples_heatmap.pdf'),
    height = 8,
    width = 8
)

#----------CDR3
clone.mat <- lapply(obj.list, function(obj.1) {
    trb.1 <- obj.1$tcr_cdr3s_aa
    jacard.idx <- lapply(obj.list, function(obj.2) {
        trb.2 <- obj.2$tcr_cdr3s_aa
        intersect(trb.1, trb.2) %>% length
    }) %>% as.data.frame
}) %>% do.call(rbind, .)

for(i in 1 : ncol(clone.mat)){
    for (j in 1 : ncol(clone.mat))
        if (i == j) {
            clone.mat[i, j] <- NA
    }
}

pheatmap(
    log2(1 + clone.mat),
    color = cols,
    cluster_rows = T,
    cluster_cols = T,
    display_numbers = clone.mat,
    main = 'Repertoire overlap',
    legend_labels = c('4', '16', '64', '256'),
    filename = file.path(out.figs.dir, 'cdr3.share.in.samples_heatmap.pdf'),
    height = 8,
    width = 8
)

#----------------------------------------------------------------------------------
# TCR statistics

obj.dn.dp <- subset(merged.obj.sc, idents = merged.obj.sc %>% levels %>% .[18 : 34])
count.mat <- table(obj.dn.dp$Anno.Level.10, obj.dn.dp$t.cell_group) %>% as.data.frame.matrix %>% .[ANNO_ENTIRE_IDNET[28 : 65], ]
t.sc.count <- table(sc.dn.dp$Anno.Level.10)[ANNO_ENTIRE_IDNET[28 : 65]]

pdf(file.path(out.figs.dir, 'number_of_cloned_cells.barplot.pdf'), width = 14, height = 6)
rotate_x(count.mat %>% rowSums %>% { . / t.sc.count } , count.mat %>% rowSums %>% names, 45, ylab = 'Number of cloned cells')
dev.off()

pdf(file.path(out.figs.dir, 'fraction_of_cloned_cells.barplot.pdf'), width = 14, height = 6)
rotate_x(count.mat %>% rowSums %>% { . / t.sc.count } , count.mat %>% rowSums %>% names, 45, ylab = 'Fraction of cloned cells')
dev.off()

annot <- as.data.frame(SampleClassify %>% unlist)
colnames(annot) <- 'Sample'
annot[, 'Var1'] <- gsub('\\d', '', rownames(annot))
rownames(annot) <- annot$Sample
annot <- annot[, -1, drop = FALSE]
anno.colors <- list(Var1 = SampleClassifyColors)

yy <- table(obj.dn.dp$orig.ident, obj.dn.dp$Anno.Level.10) / (1 + table(sc.dn.dp$orig.ident, sc.dn.dp$Anno.Level.10))
pheatmap(yy[, ANNO_ENTIRE_IDNET[28 : 65]], cluster_rows = F, cluster_cols = F, border_color=NA, annotation_row = annot, annotation_colors = anno.colors,
        main = 'Fraction of cloned cells in each sample',
        filename = file.path(out.figs.dir, 'fraction_of_cloned_cells.heatmap.pdf'), height = 6, width = 14)

xx <- table(obj.dn.dp$orig.ident, obj.dn.dp$Anno.Level.10)
xx <- sweep(xx, 1, rowSums(xx), '/')
pheatmap(yy[, ANNO_ENTIRE_IDNET[28 : 65]], cluster_rows = F, cluster_cols = F, border_color=NA, annotation_row = annot, annotation_colors = anno.colors,
        main = 'Fraction of cloned cells across cell types'
        filename = file.path(out.figs.dir, 'fraction_of_cloned_cells_across_celltypes.heatmap.pdf'), height = 6,  width = 14)

#----------------------------------------------------------------------------------

sessionInfo()
