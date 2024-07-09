suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(colorRamps))
suppressMessages(library(tidyverse))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')
source('modules/findMedullaClusters.R')
source('modules/vis/step_05_spatial_distscore_vis.R')

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '5.spdist_score/data')
out.figs.dir <- file.path('../3.results', '5.spdist_score/figs')

dir.create(file.path('../3.results', '5.spdist_score'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#---------------------------------------------------------------
# Load spatial data

obj.st.lst <- readRDS('../3.results/4.proj_spatial/data/obj.sp.lst.rds')
anno <- read.table('configs/gencode.gene.info.v22.tsv', sep = '\t', header = T) %>% .[!duplicated(.[, 2]), ] %>% `rownames<-`(.[, 2])
coding.genes <- subset(anno, (gene_type %in% ANNO.AREAS) & (gene_status == 'KNOWN')) %>% .[, 'gene_name']
obj.st.lst <- lapply(obj.st.lst, function(obj) obj[intersect(rownames(obj), coding.genes), ] )

# GLM test
obj.st.lst <- glmFitTest(obj.st.lst, tar.genes = rownames(obj.st.lst[[1]]))
saveRDS(obj.st.lst, file = file.path(out.data.dir, 'obj.st.lst.with.fit.rds'))

obj.st.lst <- readRDS(file.path(out.data.dir, 'obj.st.lst.with.fit.rds'))
#-------------------------------------------------------------
# Visulation marker genes

markers.to.plot <- c("PTPRC", "CD3E", "CD3G", "TRDC", "TRAC", "CD8A", "CD8B", "GNG4", "CD4", "TOX2", "FOXP3", "CTLA4", "GNG8", "KLRB1", "KLRD1", "MS4A1", "CD19", "JCHAIN", "IGHG1", "LAMP3", "LYZ", "C1QC", "HLA-DPA1", "HLA-DPB1", "EPCAM", "CCL19", "CCL25", "PECAM1", "ACTA2")

gp1 <- plotBoxplots(obj.st.lst, markers.to.plot) + ylab('Row-normalized expression') + stat_summary(fun.y=mean, geom="point", size=0.5, color="red", fill = 'red', aes(group = factor(Label)), position=position_dodge(width=0.75)) + coord_flip()
ggsave(file.path(out.figs.dir, 'markers.to.boxplot.C_J_M.pdf'), height = 14, width = 6)

gp.lst <- lapply(markers.to.plot, function(gene){
	plot.df <- trendLinesWithDistPlot2(obj.st.lst, plot.tar = gene) + xlab('') + ylab('') + theme(legend.position = 'none') + ggtitle(gene) + scale_fill_manual(values = 'grey')
})

ggarrange(plotlist=gp.lst, widths = c(6, 5))
ggsave(file.path(out.figs.dir, 'markers.to.trend.C_J_M.pdf'), height = 8, width = 14)

#------------------------------------------------------------
# Dynamic-related genes

markers.to.plot <- c('SELPLG', 'CXCR4', 'CXCL12', 'CCR7',  'CCL19', 'CCL21', 'CCR9', 'CCL25', 'SELL')
trendLinesWithDistPlot2(obj.st.lst, plot.tar = markers.to.plot)
ggsave(file.path(out.figs.dir, 'development_sel.genes.to.trend.C_J_M.pdf'), height = 6, width = 8)

gp1 <- plotBoxplots(obj.st.lst, markers.to.plot) + ylab('Row-normalized expression') + 
	stat_summary(fun.y=mean, geom="point", size=0.5, color="red", fill = 'red', aes(group = factor(Label)), position=position_dodge(width=0.75)) + coord_flip()
ggsave(file.path(out.figs.dir, 'development_sel.genes.to.boxplot.C_J_M.pdf'), height = 10, width = 6)

#-------------------------------------------------------------
# Select significant genes

merged.pvals <- combinePvalsByStouffer(obj.st.lst)
p.sig.gene <- subset(merged.pvals, FDR.1.merged < 0.01 & FDR.2.merged < 0.01 & FDR.3.merged < 0.01 & FDR.4.merged < 0.01)
saveRDS(merged.pvals, file = file.path(out.data.dir, 'merged.pvals.rds'))

coefs.df <- lapply(obj.st.lst, function(obj) {
	xx@misc$GLM_FIT_TEST[rownames(p.sig.gene), ]		
}) %>% `names<-`(sp2SampleName(names(obj.st.lst))) %>% .[paste0('Thy', 5 : 12)] %>% do.call(cbind, .) %>% .[, !grepl('FDR', colnames(.))]

p.sig.gene.combine <- cbind.data.frame(p.sig.gene, coefs.df)
write.table(p.sig.gene.combine, file.path(out.data.dir, 'expr_vs_distance.sig.gene.xls'), sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

#---------------------------------------------------
# DEGs analysis between C-M-J

obj.merged <- Reduce(function(x,y) merge(x, y, add.cell.ids = c(x@project.name,y@project.name)) , obj.st.lst)
obj.merged <- subset(obj.merged, HE.Labels != 'Medulla_lo')
obj.merged@meta.data$Class <- obj.merged@meta.data$HE.Labels
obj.merged@meta.data$Class[obj.merged@meta.data$Class %in% c("Medulla_hi", "Medulla_centric", "Medulla_edge")] <- 'Medulla'

Idents(obj.merged) <- obj.merged@meta.data$Class
as.loom(obj.merged, filename = file.path(out.data.dir, "obj.sp.merged.loom"), verbose = FALSE)

require(future)
options(future.globals.maxSize = 500000 * 1024^2)
plan("multiprocess", workers = 20)

degs.roc <- FindAllMarkers(obj.merged, test.use = 'roc', only.pos = TRUE, group.by = 'Class', return.thresh = 0.60)
degs.wic <- FindAllMarkers(obj.merged, only.pos = TRUE, group.by = 'Class', return.thresh = 0.01)
openxlsx::write.xlsx(list(Markers = degs.wic), file = file.path(out.data.dir, 'cortext.medulla.degs.xlsx'))

degs.roc <- subset(degs.roc, cluster %in% c('Medulla', 'Cortex'))
degs.wic <- subset(degs.wic, cluster %in% c('Medulla', 'Cortex'))
degs.wic.sub <- subset(degs.wic, avg_log2FC >= 0.25 & p_val_adj < 0.01) 

pdf(file = file.path(out.figs.dir, "dist.vary.genes.vs.degs_venn.pdf"), width = 5, height = 3)
temp <- ovpGenesSig(rownames(p.sig.gene), degs.wic.sub$gene)
grid.draw(temp)
dev.off()

ovp.sig.genes <- intersect(rownames(p.sig.gene), degs.wic.sub$gene)
top.genes <- intersect(degs.wic.sub %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% .[, 'gene'] %>% unlist, ovp.sig.genes)

top5GenesTrend(obj.st.lst, top.genes)
ggsave(file.path(out.figs.dir, 'top5.genes.along.spot.distance.pdf'), width = 7, height = 4)	

obj.merged <- ScaleData(obj.merged)
pdf(file.path(out.figs.dir, 'medulla.cortex.degs.pdf'), width = 6, height = 10)
doheatmapTop30Deg(obj.merged, topn = 30)
dev.off()

ovp.sig.genes <- intersect(rownames(p.sig.gene), degs.wic.sub$gene)
obj.sub <- obj.merged[ovp.sig.genes, ]
saveRDS(obj.sub, file = file.path(out.data.dir, 'obj.sp.for.auc.rds'))

#------------------------------------------------------
# Identify distance vary genes.

p.sig.gene.sub <- p.sig.gene
cls.genes <- clusterHyoerVariabGenesSp(obj.st.lst, p.sig.gene.sub)
gene.cls <- cls.genes$genes
expr.df <- cls.genes$expr
idx <- cls.genes$genes[names(cls.genes$genes) == 'TRDC']
up.genes <- gene.cls[gene.cls != idx] %>% names
down.genes <- gene.cls[gene.cls == idx] %>% names

anno.res <- gannoRes(up.genes, down.genes)
pdf(file.path(out.figs.dir, 'Sig.genes.along_distance.heatmap.pdf'), width = 10, height = 12)
hp <- plotComplexHeatmapWithAnno(expr.df, anno.res$ego.up, anno.res$ego.down, idx.2 = c(1, 16, 52), idx.1 = c(1, 4, 6))
draw(hp, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()
	
lapply(names(obj.st.lst), function(sn) {
	obj <- obj.st.lst[[sn]]
	SpatialPlot(obj, alpha = 0) + NoLegend()
	ggsave(file.path(out.figs.dir, sprintf('%s.HE.image.pdf', sn)), width = 5, height = 5)
})

#------------------------------------------------------
# Show marker genes

top50 <- degs.wic.sub %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
obj.merged <- ScaleData(obj.merged)

library("viridis")
library(writexl)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(
	subset(obj.merged, Class != 'Junction'), 
	group.by = 'Class', 
	features = top50[, 'gene'] %>% unlist,
	group.color = CM_COLORS[c('Medulla', 'Cortex')],
#angle = 45,
	hjust = 0,
	size = 4
) + scale_fill_gradientn(colours = rev(mapal)) 

ggsave(file.path(out.figs.dir, 'medulla_cortex_roc_genes.heatmap.pdf'), width = 8, height = 12)
write_xlsx(list(ROC = degs.roc), file.path(out.data.dir, "degs_cortex_medulla.xlsx"))

#------------------------------------------------------
# Functional analysis

library(clusterProfiler)
gene <- bitr(degs.roc$gene, fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
colnames(degs.roc)[8] <- 'SYMBOL'
data_all <- degs.roc %>% inner_join(gene,by="SYMBOL")

data_all_sort <- data_all %>% arrange(desc(avg_log2FC))
geneList <- data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$SYMBOL

go_gmt <- read.gmt("configs/c5.all.v7.4.symbols.gmt")
gsea <- GSEA(geneList, TERM2GENE = go_gmt)

#-----------

ego.medulla <- clusterProfiler::enrichGO(
	gene          = subset(data_all, cluster == 'Medulla')$ENTREZID %>% unique,
	OrgDb         = org.Hs.eg.db,
	ont           = "BP",
	pAdjustMethod = "BH",
	pvalueCutoff  = 0.01,
	qvalueCutoff  = 0.05,
	readable      = TRUE
)
ego.cortex <- clusterProfiler::enrichGO(
    gene          = subset(data_all, cluster == 'Cortex')$ENTREZID %>% unique,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
)

anno <- cbind.data.frame(Category = 'BP', rbind.data.frame(head(ego.medulla %>% data.frame, n = 4), head(ego.cortex %>% data.frame, n = 4)))
anno.sub <- anno[, c(1, 2, 3, 7, 9)] %>% `colnames<-`(c('Category', 'ID', 'Term', 'adj_pval', 'Genes'))
anno.sub$Genes <- gsub('/', ', ', anno.sub$Genes)

library(GOplot)
genes.tmp <- degs.wic
colnames(genes.tmp) <- c('p_val', 'logFC', 'pct.1', 'pct.2', 'adj.P.Val', 'cluster', 'ID')
EC <- list(david = anno.sub, genelist = genes.tmp)
circ <- circle_dat(EC$david, EC$genelist)

pdf(file.path(out.figs.dir, 'degs_go_anno.pdf'), width = 13, height = 8)
GOplot::GOCircle(circ, nsub = 7, rad2 = 3.5, rad1 = 2.5)
dev.off()

trendLinesWithDistPlot2(obj.st.lst, c('CD99', 'CD1E', 'TRBC1', 'STMN1', 'MZB1', 'CD1B', 'RAG1', 'HMGB2', 'PTP4A2', 'H2AFZ', 'CD8A', 'CD8B'))
ggsave(file.path(out.figs.dir, 'top_roc_genes_expr_dist_trend.pdf'), width = 7, heigth = 6)

#----------------------------------------------------
# Co-expressed analyses

degs.wic.cortex <- subset(degs.wic, cluster == 'Cortex') 
avg.entire.paterns <- searchAllTrendPatterns(obj.st.lst, obj.st.lst[[1]] %>% rownames, win = 30, GLM = T)
cd99.paterns.corrs <- identifySigCorrPatterns(avg.entire.paterns, 'CD99')

cd99.paterns.corrs$Degs <- 'N'
cd99.paterns.corrs[rownames(cd99.paterns.corrs) %in% degs.wic.cortex$gene, 'Degs'] <- 'Y'
cd99.paterns.corrs <- cbind.data.frame(Symbol = rownames(cd99.paterns.corrs), cd99.paterns.corrs)
writexl::write_xlsx(list(ROC = cd99.paterns.corrs), file.path(out.data.dir, "CD99_Coexpress_genes.xlsx"))
saveRDS(cd99.paterns.corrs, file.path(out.data.dir, "CD99_Coexpress_genes.rds"))
