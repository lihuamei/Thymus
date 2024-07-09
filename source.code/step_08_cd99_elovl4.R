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

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path('../3.results', '7.cd99_elovl4/data')
out.figs.dir <- file.path('../3.results', '7.cd99_elovl4/figs')

dir.create(file.path('../3.results', '7.cd99_elovl4'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load spatial data

obj.st.lst <- readRDS('../3.results/4.proj_spatial/data/obj.sp.lst.rds')
obj.st.merged <- readRDS('../3.results/4.proj_spatial/data/obj.st.merged.rds')
anno <- read.table('configs/gencode.gene.info.v22.tsv', sep = '\t', header = T) %>% .[!duplicated(.[, 2]), ] %>% `rownames<-`(.[, 2])
coding.genes <- subset(anno, gene_type == 'protein_coding') %>% .[, 'gene_name'] %>% intersect(., rownames(obj.st.lst[[1]]))
obj.st.lst <- lapply(obj.st.lst, function(xx) subset(xx, idents = setdiff(levels(xx), 'Contamination spots')))
obj.sc <- obj.merged.sc <- readRDS('../3.results/2.fine_anno/data/merged.obj.anno.rds')

#-----------------------------------------------------------
# CD99  vs. TRCA and TRCB

TRAC.expr <- TRAC.expr.bak <- lapply(obj.st.lst, function(obj) {
    expr.sub <- AverageExpression(obj)$Spatial
	expr.sub <- cbind.data.frame(
        CD99 = expr.sub['CD99', ],
        TRAC = expr.sub[grep(("^TRA[VJC]"), rownames(obj), value = T), ] %>% colMeans(.),
        TRBC1 = expr.sub[grep(("^TRB[VJC]"), rownames(obj), value = T), ] %>% colMeans(.)
    ) %>% t
	expr.sub <- rbind.data.frame(expr.sub, Cluster = colnames(expr.sub), Sample = obj@images %>% names)
    colnames(expr.sub) <- paste0(obj@images %>% names, '.', colnames(expr.sub))
    expr.sub
}) %>% do.call(cbind, .) %>% t %>% as.data.frame

TRAC.expr[, 'CD99'] <- as.numeric(TRAC.expr[, 'CD99']) %>% scale
TRAC.expr[, 'TRAC'] <- as.numeric(TRAC.expr[, 'TRAC']) %>% scale
TRAC.expr[, 'TRBC1'] <- as.numeric(TRAC.expr[, 'TRBC1']) %>% scale

plot.df <- tidyr::gather(TRAC.expr, 'Gene', 'Expr', -Cluster, -Sample)
plot.df.new <- aggregate(Expr ~ Cluster + Gene, data = plot.df, FUN = 'mean')
plot.df.new$Cluster <- factor(plot.df.new$Cluster, levels = setdiff(SP_ANNO_ORDER, 'Contamination spots')) 
ggplot(plot.df.new, aes(x = Cluster, y = Expr, color = Gene, group = Gene)) + 
	geom_point(aes(fill=Gene), size=3) +
    geom_line(aes(color=Gene, linetype = Gene)) +
	theme_classic(base_size = 16) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab('Scaled expression') + xlab('')
ggsave(file.path(out.figs.dir, 'TRAC_TRBC_vs._CD99_expr.pdf'), width = 7, height = 5)

TRAC.expr.bak[, 'CD99'] <- as.numeric(TRAC.expr.bak[, 'CD99']) %>% { log2(1 + .) }
TRAC.expr.bak[, 'TRAC'] <- as.numeric(TRAC.expr.bak[, 'TRAC']) %>% { log2(1 + .) }
TRAC.expr.bak[, 'TRBC1'] <- as.numeric(TRAC.expr.bak[, 'TRBC1']) %>% { log2(1 + .) }
TRAC.expr.bak$Cluster <- factor(TRAC.expr.bak$Cluster, levels = setdiff(SP_ANNO_ORDER, 'Contamination spots'))
ggplot(data = TRAC.expr.bak, aes(x = CD99, y = TRAC, color = Sample, shape = Cluster)) + 
	geom_point() + 
	scale_color_manual(values = SAMPLE_COLORS) + 
	theme_classic(base_size = 16) + 
	geom_smooth(method = "lm", se = TRUE, aes(group = 1)) + 
	stat_cor(method = "pearson", aes(group = 1)) + xlab('CD99 expression') + ylab('TRAC expression') + 
	scale_shape_manual(values = 1 : (TRAC.expr.bak$Sample %>% unique %>% length))
ggsave(file.path(out.figs.dir, 'TRAC_vs._CD99_expr_scatter.pdf'), width = 7, height = 6)

ggplot(data = TRAC.expr.bak, aes(x = CD99, y = TRBC1, color = Sample, shape = Cluster)) +
    geom_point() +
    scale_color_manual(values = SAMPLE_COLORS) +
    theme_classic(base_size = 16) +
    geom_smooth(method = "lm", se = TRUE, aes(group = 1)) +
    stat_cor(method = "pearson", aes(group = 1)) + xlab('CD99 expression') + ylab('TRBC1 expression') + 
	scale_shape_manual(values = 1 : (TRAC.expr.bak$Sample %>% unique %>% length))
ggsave(file.path(out.figs.dir, 'TRBC1_vs._CD99_expr_scatter.pdf'), width = 7, height = 6)

#------------------------------------------------------------------------------
# TRBV vs. CD99 on the basis of single-cell data

obj.sc.sub <- subset(obj.merged.sc, idents = c("DN_early", "DN_blast", "DN_re", "DP_blast", "DP_re", "abT(entry)"))
obj.sc.sub <- subset(obj.sc.sub, orig.ident %in% (obj.sc.sub@meta.data$orig.ident %>% levels %>% .[5 : 12]))

TRAC.expr <- TRAC.expr.bak <- lapply(levels(obj.sc.sub@meta.data$orig.ident) %>% .[5 : 12], function(sn) {
    obj <- subset(obj.sc.sub, orig.ident == sn)
	expr.sub <- AverageExpression(obj)$RNA
    expr.sub <- cbind.data.frame(
        CD99 = expr.sub['CD99', ],
        TRAC = expr.sub[grep(("^TRA[VJC]"), rownames(obj), value = T), ] %>% colMeans(.),
        TRBC1 = expr.sub[grep(("^TRB[VJC]"), rownames(obj), value = T), ] %>% colMeans(.)
    ) %>% t
    expr.sub <- rbind.data.frame(expr.sub, Cluster = colnames(expr.sub), Sample = SC_MAP_SP[[sn]] %>% SP_NAMES[.])
    colnames(expr.sub) <- paste0(obj@images %>% names, '.', colnames(expr.sub))
    expr.sub
}) %>% do.call(cbind, .) %>% t %>% as.data.frame

TRAC.expr[, 'CD99'] <- as.numeric(TRAC.expr[, 'CD99']) %>% scale
TRAC.expr[, 'TRAC'] <- as.numeric(TRAC.expr[, 'TRAC']) %>% scale
TRAC.expr[, 'TRBC1'] <- as.numeric(TRAC.expr[, 'TRBC1']) %>% scale

plot.df <- tidyr::gather(TRAC.expr, 'Gene', 'Expr', -Cluster, -Sample)
plot.df.new <- aggregate(Expr ~ Cluster + Gene, data = plot.df, FUN = 'mean')
plot.df.new$Cluster <- factor(plot.df.new$Cluster, levels = c("DN_early", "DN_blast", "DN_re", "DP_blast", "DP_re", "abT(entry)"))

ggplot(plot.df.new, aes(x = Cluster, y = Expr, color = Gene, group = Gene)) +
    geom_point(aes(fill=Gene), size=3) +
    geom_line(aes(color=Gene, linetype = Gene)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab('Scaled expression') + xlab('')
ggsave(file.path(out.figs.dir, 'TRAC_TRBC_vs._CD99_expr_sc.pdf'), width = 7, height = 5)

TRAC.expr.bak[, 'CD99'] <- as.numeric(TRAC.expr.bak[, 'CD99']) %>% { log2(1 + .) }
TRAC.expr.bak[, 'TRAC'] <- as.numeric(TRAC.expr.bak[, 'TRAC']) %>% { log2(1 + .) }
TRAC.expr.bak[, 'TRBC1'] <- as.numeric(TRAC.expr.bak[, 'TRBC1']) %>% { log2(1 + .) }

TRAC.expr.bak$Cluster <- factor(TRAC.expr.bak$Cluster, levels = c("DN_early", "DN_blast", "DN_re", "DP_blast", "DP_re", "abT(entry)"))
ggplot(data = TRAC.expr.bak, aes(x = CD99, y = TRAC, color = Sample, shape = Cluster)) +
    geom_point() +
    scale_color_manual(values = SAMPLE_COLORS) +
    theme_classic(base_size = 16) +
    geom_smooth(method = "lm", se = TRUE, aes(group = 1)) +
    stat_cor(method = "pearson", aes(group = 1)) + xlab('CD99 expression') + ylab('TRAC expression') +
    scale_shape_manual(values = 1 : (TRAC.expr.bak$Sample %>% unique %>% length))
ggsave(file.path(out.figs.dir, 'TRAC_vs._CD99_sc_expr_scatter.pdf'), width = 7, height = 6)

ggplot(data = TRAC.expr.bak, aes(x = CD99, y = TRBC1, color = Sample, shape = Cluster)) +
    geom_point() +
    scale_color_manual(values = SAMPLE_COLORS) +
    theme_classic(base_size = 16) +
    geom_smooth(method = "lm", se = TRUE, aes(group = 1)) +
    stat_cor(method = "pearson", aes(group = 1)) + xlab('CD99 expression') + ylab('TRBC1 expression') +
    scale_shape_manual(values = 1 : (TRAC.expr.bak$Sample %>% unique %>% length))
ggsave(file.path(out.figs.dir, 'TRBC1_vs._CD99_sc_expr_scatter.pdf'), width = 7, height = 6)
	
#--------------------------------------------------------

genes <- rownames(obj.sc)

tra <- grep(("^TRA[VJC]"), genes, value = T)
trb <- grep(("^TRB[VJC]"), genes, value = T)
trg <- grep(("^TRG[VJC]"), genes, value = T)
trd <- grep(("^TRD[VJC]"), genes, value = T)

obj.sc <- AddModuleScore(obj.sc, features = list(xx = tra, yy = trb, zz = trd))
xx = subset(obj.sc, idents = c("DN_early","DN_blast","DN_re","DP_blast","DP_re","abT(entry)","CD8T", 'CD4T'))
xx <- subset(xx, tcr == 1)

xx$AgeStatus.New <- as.vector(xx$AgeStatus)
xx$AgeStatus.New[xx$orig.ident %in% c('Thy15', 'Thy16')] <- 'Older'

yy = FetchData(xx, vars = c('CD99', 'ELOVL4', 'Age', 'Cluster1', 'Anno.Level.Fig.1', 'orig.ident', 'AgeStatus.New')) %>% subset(., Anno.Level.Fig.1 %in% c('DP_blast', 'DP_re'))
vv = subset(yy, Anno.Level.Fig.1 == 'DP_re')
plot.vv <- aggregate(ELOVL4 ~ orig.ident, vv, 'mean')
plot.vv$Xaes <- 1 : nrow(plot.vv)
plot.vv$Age <- unlist(SampleInfos)
ggplot(plot.vv, aes(x = Xaes, y = ELOVL4)) +
    geom_point(size = 5, aes(color = orig.ident)) +
    geom_smooth(method = loess, color = 'red', linetype = 'dashed') +
    scale_x_continuous(breaks = 1 : nrow(plot.vv), labels = as.vector(plot.vv$orig.ident)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') +
    ylab('Average expression of ELOVL4 in DP_re cells') + scale_color_manual(values = SAMPLE.COLORS_SC[1 : 16]) +
    geom_text(aes(label = plot.vv$Age), check_overlap = TRUE) #+
#ggpubr::stat_cor()
ggsave(file.path(out.figs.dir, 'ELOVL4.exp.vs.dp.re.age.pdf'), width = 5, height = 4)

ggplot(yy, aes(x = Anno.Level.Fig.1, y = CD99)) +
    geom_violin(aes(color = Anno.Level.Fig.1, fill = Anno.Level.Fig.1)) +
    geom_boxplot(width = 0.8, color = 'black', fill = 'white', alpha = 0.5, outlier.shape = NA) +
    stat_compare_means(method = 't.test') +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')]) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')])

ggsave(file.path(out.figs.dir, 'cd99.score.vs.dp.pdf'), width = 5, height = 5)

ggplot(yy, aes(x = Anno.Level.Fig.1, y = Cluster1)) +
    geom_violin(aes(color = Anno.Level.Fig.1, fill = Anno.Level.Fig.1)) +
    geom_boxplot(width = 0.8, color = 'black', fill = 'white', alpha = 0.5, outlier.shape = NA) +
    stat_compare_means(method = 't.test') +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')]) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')])

ggsave(file.path(out.figs.dir, 'tra.score.vs.dp.pdf'), width = 5, height = 5)

# CD99 vs. tra/b/d  

yy = FetchData(xx, vars = c('CD99', 'Anno.Level.Fig.1'))
kk = aggregate(CD99 ~ Anno.Level.Fig.1, yy, 'mean')

yy = FetchData(xx, vars = c('Cluster1', 'Anno.Level.Fig.1'))
vv = aggregate(Cluster1 ~ Anno.Level.Fig.1, yy, 'mean')
plot.df <- cbind.data.frame(kk, TRA = vv[, 'Cluster1'])
ggplot(plot.df, aes(x = scale(CD99), y = scale(TRA), label = Anno.Level.Fig.1)) +
    geom_point(aes(color = Anno.Level.Fig.1), size = 3) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[plot.df[, 1]]) +
    theme_bw(base_size = 12) +
    geom_smooth(method = lm) +
    theme(panel.grid = element_blank(), legend.position = 'none') +
    ggpubr::stat_cor(method = "spearman") +
    geom_text(hjust=0, vjust=0)
ggsave(file.path(out.figs.dir, 'tra.score.vs.cd99.pdf'), width = 5, height = 5)

yy = FetchData(xx, vars = c('Cluster2', 'Anno.Level.Fig.1'))
vv = aggregate(Cluster2 ~ Anno.Level.Fig.1, yy, 'mean')
plot.df <- cbind.data.frame(kk, TRA = vv[, 'Cluster2'])
ggplot(plot.df, aes(x = scale(CD99), y = scale(TRA), label = Anno.Level.Fig.1)) +
    geom_point(aes(color = Anno.Level.Fig.1), size = 3) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[plot.df[, 1]]) +
    theme_bw(base_size = 12) +
    geom_smooth(method = lm) +
    theme(panel.grid = element_blank(), legend.position = 'none') +
    ggpubr::stat_cor(method = "spearman") +
    geom_text(hjust=0, vjust=0)
ggsave(file.path(out.figs.dir, 'trb.score.vs.cd99.pdf'), width = 5, height = 5)

yy = FetchData(xx, vars = c('Cluster3', 'Anno.Level.Fig.1'))
vv = aggregate(Cluster3 ~ Anno.Level.Fig.1, yy, 'mean')
plot.df <- cbind.data.frame(kk, TRA = vv[, 'Cluster3'])
ggplot(plot.df, aes(x = scale(CD99), y = scale(TRA), label = Anno.Level.Fig.1)) +
    geom_point(aes(color = Anno.Level.Fig.1), size = 3) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[plot.df[, 1]]) +
    theme_bw(base_size = 12) +
    geom_smooth(method = lm, linetype = 'dashed', se = F) +
    theme(panel.grid = element_blank(), legend.position = 'none') +
    ggpubr::stat_cor(method = "spearman") +
    geom_text(hjust=0, vjust=0) + ggtitle('Ours')
ggsave(file.path(out.figs.dir, 'trd.score.vs.cd99.pdf'), width = 5, height = 5)

#------------------------------------------------------
# science

obj.science.t <- LoadH5Seurat("../4.extdata/HTA08.v01.A06.Science_human_tcells_new.h5seurat")
obj.science.t <- AddModuleScore(obj.science.t, features = list(xx = trd, yy = tra))
xx = subset(obj.science.t, cell.types %in% c('DN(early)', "DN(P)", 'DN(Q)', 'DP(P)', 'DP(Q)', 'αβT(entry)', 'CD4+T', 'CD8+T'))
yy = FetchData(xx, vars = c('Cluster1', 'cell.types'))
vv = aggregate(Cluster1 ~ cell.types, yy, 'mean')

yy = FetchData(xx, vars = c('CD99', 'cell.types'))
kk = aggregate(CD99 ~ cell.types, yy, 'mean')

plot.df.1 <- cbind.data.frame(kk, TRA = vv[, 'Cluster1'])
cols <-  ANNO_ENTIRE_COLOR_FIG1[plot.df[, 1]] %>% `names<-`(c('DN(early)', "DN(P)", 'DN(Q)', 'DP(P)', 'DP(Q)', 'αβT(entry)', 'CD8+T', 'CD4+T'))
ggplot(plot.df.1, aes(x = scale(CD99), y = scale(TRA), label = cell.types)) +
    geom_point(aes(color = cell.types), size = 3) +
    scale_color_manual(values = cols) +
    theme_bw(base_size = 12) +
    geom_smooth(method = lm, linetype = 'dashed', se = F) +
    theme(panel.grid = element_blank(), legend.position = 'none') +
    ggpubr::stat_cor(method = "spearman") +
    geom_text(hjust=0, vjust=0) + ggtitle('Park et al. (2020, Science)')
ggsave(file.path(out.figs.dir, 'trd.score.science.vs.cd99.pdf'), width = 5, height = 5)

yy = FetchData(xx, vars = c('CD99', 'Cluster2', 'cell.types')) %>% subset(., cell.types %in% c('DP(P)', 'DP(Q)'))
ggplot(yy, aes(x = cell.types, y = CD99)) +
    geom_violin(aes(color = cell.types, fill = cell.types)) +
    geom_boxplot(width = 0.8, color = 'black', fill = 'white', alpha = 0.5, outlier.shape = NA) +
    stat_compare_means(method = 't.test') +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')] %>% `names<-`(c('DP(P)', 'DP(Q)'))) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')] %>% `names<-`(c('DP(P)', 'DP(Q)')))

ggsave(file.path(out.figs.dir, 'science.cd99.score.vs.dp.pdf'), width = 5, height = 5)

ggplot(yy, aes(x = cell.types, y = Cluster2)) +
    geom_violin(aes(color = cell.types, fill = cell.types)) +
    geom_boxplot(width = 0.8, color = 'black', fill = 'white', alpha = 0.5, outlier.shape = NA) +
    stat_compare_means(method = 't.test') +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')] %>% `names<-`(c('DP(P)', 'DP(Q)'))) +
    scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1[c('DP_blast', 'DP_re')] %>% `names<-`(c('DP(P)', 'DP(Q)')))

ggsave(file.path(out.figs.dir, 'science.tra.score.vs.dp.pdf'), width = 5, height = 5)

#------------------------------------------------------

yy = FetchData(xx, vars = c('CD99', 'ELOVL4', 'Age', 'Cluster1', 'Anno.Level.Fig.1', 'orig.ident', 'AgeStatus.New')) %>% subset(., Anno.Level.Fig.1 %in% c('DP_re'))
plot.vv <- aggregate(ELOVL4 ~ orig.ident, yy, 'mean')
plot.vv$Xaes <- 1 : nrow(plot.vv)
plot.vv$Age <- unlist(SampleInfos)
plot.vv <- cbind.data.frame(plot.vv, Group = 'ELOVL4')
colnames(plot.vv)[2] <- 'Score'
plot.vv$Score <- scale(plot.vv$Score)

plot.vv.1 <- aggregate(Cluster1 ~ orig.ident, yy, 'mean')
plot.vv.1$Xaes <- 1 : nrow(plot.vv)
plot.vv.1$Age <- unlist(SampleInfos)
plot.vv.1 <- cbind.data.frame(plot.vv.1, Group = 'Cluster1')
colnames(plot.vv.1)[2] <- 'Score'
plot.vv.1$Score <- scale(plot.vv.1$Score)

plot.df <- rbind.data.frame(plot.vv, plot.vv.1)
plot.df.com <- cbind.data.frame(plot.vv, TRA = plot.vv.1[, 2])

ggplot(plot.df, aes(x = Xaes, y = Score, color = Group)) +
    geom_point(size = 5, alpha = 0.5) +
    geom_smooth(method = loess, aes(color = Group), se = FALSE, linetype = 'dashed') +
    scale_x_continuous(breaks = 1 : nrow(plot.vv), labels = as.vector(plot.vv$Age)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab('') +
    ylab('Average expression of ELOVL4 in DP_re cells') +
    scale_color_manual(values = c('ELOVL4' = '#DA650D', Cluster1 = '#8480BB'))
ggsave(file.path(out.figs.dir, 'ELOVL4.tra.pdf'), width = 7, height = 3)

# CD99
plot.vv <- aggregate(CD99 ~ orig.ident, yy, 'mean')
plot.vv$Xaes <- 1 : nrow(plot.vv)
plot.vv$Age <- unlist(SampleInfos)
plot.vv <- cbind.data.frame(plot.vv, Group = 'CD99')
colnames(plot.vv)[2] <- 'Score'
plot.vv$Score <- scale(plot.vv$Score)

plot.vv.1 <- aggregate(Cluster1 ~ orig.ident, yy, 'mean')
plot.vv.1$Xaes <- 1 : nrow(plot.vv)
plot.vv.1$Age <- unlist(SampleInfos)
plot.vv.1 <- cbind.data.frame(plot.vv.1, Group = 'Cluster1')
colnames(plot.vv.1)[2] <- 'Score'
plot.vv.1$Score <- scale(plot.vv.1$Score)

plot.df <- rbind.data.frame(plot.vv, plot.vv.1)
plot.df.com <- cbind.data.frame(plot.vv, TRA = plot.vv.1[, 2])


ggplot(plot.df, aes(x = Xaes, y = Score, color = Group)) +
    geom_point(size = 5, alpha = 0.5) +
    geom_smooth(method = loess, aes(color = Group), se = FALSE, linetype = 'dashed') +
    scale_x_continuous(breaks = 1 : nrow(plot.vv), labels = as.vector(plot.vv$Age)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab('') +
    ylab('Average expression of CD99 in DP_re cells') +
    scale_color_manual(values = c('CD99' = '#22A079', Cluster1 = '#8480BB'))
ggsave(file.path(out.figs.dir, 'CD99.tra.pdf'), width = 7, height = 3)

# science

obj.science.t <- LoadH5Seurat("../4.extdata/HTA08.v01.A05.Science_human_figX.h5Seurat")
#obj.science.t <- NormalizeData(obj.science.t)

obj.science.t <- AddModuleScore(obj.science.t, features = list(xx = tra))

yy = FetchData(obj.science.t, vars = c('CD99', 'ELOVL4', 'Anno_level_5', 'Age', 'Cluster1')) %>% subset(.,Anno_level_5 %in% c('DP(Q)'))
order.sample <- c("7w", "8w", "9w", "10w", "11w", "12w", "13w", "14w", "16w", "17w", "3m", "6m", "10m", "15m", "30m", "13y", "24y", "35y")
yy$Cluster1 <- yy$Cluster1
yy$Age <- factor(yy$Age, levels = order.sample)
plot.vv <- aggregate(CD99 ~ Age, vv, 'mean')
plot.vv$Xaes <- 1 : nrow(plot.vv)
colnames(plot.vv)[2] <- 'Score'
plot.vv$Score <- scale(plot.vv$Score)

yy$Age <- factor(yy$Age, levels = order.sample)
plot.vv.1 <- aggregate(Cluster1 ~ Age, yy, 'mean')
plot.vv.1$Xaes <- 1 : nrow(plot.vv.1)
colnames(plot.vv.1)[2] <- 'Score'
plot.vv.1$Score <- scale(plot.vv.1$Score)

plot.df.sci <- plot.df <- cbind.data.frame(plot.vv, TRA = plot.vv.1[, 2])
plot.df <- rbind.data.frame(cbind.data.frame(plot.df.sci, SRC = 'Sci'), cbind.data.frame(plot.df.com[, colnames(plot.df.sci)], SRC = 'Ours'))

ggplot(plot.df, aes(x = TRA, y = Score)) +
    geom_point(size = 3, aes(color = SRC)) +
    geom_smooth(method = lm, color = 'red', linetype = 'dashed') +
    #scale_x_continuous(breaks = 1 : nrow(plot.vv), labels = as.vector(plot.vv$Age)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = 'none') +
    xlab('') + ggpubr::stat_cor() +
    ylab('Average expression of CD99 in DP_re cells') + ggpubr::stat_cor()
    #geom_text(label = plot.vv$Age, color = 'black')

ggsave(file.path(out.figs.dir, 'CD99.exp.vs.dp.re.age.scatter.pdf'), width = 4, height = 4) #+

plot.vv <- aggregate(ELOVL4 ~ Age, yy, 'mean')
plot.vv$Xaes <- 1 : nrow(plot.vv)
colnames(plot.vv)[2] <- 'Score'
plot.vv$Score <- scale(plot.vv$Score)

yy$Age <- factor(yy$Age, levels = order.sample)
plot.vv.1 <- aggregate(Cluster1 ~ Age, yy, 'mean')
plot.vv.1$Xaes <- 1 : nrow(plot.vv.1)
colnames(plot.vv.1)[2] <- 'Score'
plot.vv.1$Score <- scale(plot.vv.1$Score)

plot.df.sci <- plot.df <- cbind.data.frame(plot.vv, TRA = plot.vv.1[, 2])

plot.df <- rbind.data.frame(cbind.data.frame(plot.df.sci, SRC = 'Sci'), cbind.data.frame(plot.df.com[, colnames(plot.df.sci)], SRC = 'Ours'))

ggplot(plot.df, aes(x = TRA, y = Score)) +
    geom_point(size = 3, aes(color = SRC)) +
    geom_smooth(method = lm, color = 'red', linetype = 'dashed') +
    #scale_x_continuous(breaks = 1 : nrow(plot.vv), labels = as.vector(plot.vv$Age)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab('') +
    ylab('Average expression of ELOVL4 in DP_re cells') +
    #geom_text(label = plot.vv$Age, color = 'black') +
    ggpubr::stat_cor()

ggsave(file.path(out.figs.dir, 'ELOVL4.exp.vs.dp.re.age.scatter.pdf'), width = 4, height = 4) #+

#------------------------------------------------------
sessionInfo()
