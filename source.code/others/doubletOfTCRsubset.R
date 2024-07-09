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

source("./modules/utils.R")
source("./modules/global_params.R")
source("./modules/seurat_methods.R")
source("./modules/findMedullaClusters.R")
source("./modules/trajectory_methods.R")
source("./modules/tcrSubgroup.R")
source("./modules/vis/step_07_order_T_subsets_vis.R")

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path("../3.results", "8.clonetype/data.update")
out.figs.dir <- file.path("../3.results", "8.clonetype/figs.update")

dir.create(file.path("../3.results", "8.clonetype"))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#-------------------------------------------------------------
# Pseudotime

obj.merged.sc <- readRDS(file.path("../../3.results/8.clonetype/data.update", "obj.merged.sc.reformat.new.rds"))
# obj <- subset(obj.merged.sc, tcr == 1)
obj.xx <- subset(obj.merged.sc, Anno.Level.11 %in% c("DN_re1", "DN_blast", "DN_re2", "DN_early", "ISP"))
cc <- FetchData(obj.xx, slot = "counts", vars = c("CD34", "CD38", "CD1A"))
cut.off <- 0
dn1 <- subset(cc, CD34 > cut.off & CD38 <= cut.off & CD1A <= cut.off) %>% rownames()
dn2 <- subset(cc, CD34 > cut.off & CD38 > cut.off & CD1A <= cut.off) %>% rownames()
dn3 <- subset(cc, CD34 <= cut.off & CD38 > cut.off & CD1A > cut.off) %>% rownames()
obj.xx$DN <- "-"
obj.xx$DN[rownames(obj.xx@meta.data) %in% dn1] <- "DN1"
obj.xx$DN[rownames(obj.xx@meta.data) %in% dn2] <- "DN2"
obj.xx$DN[rownames(obj.xx@meta.data) %in% dn3] <- "DN3"

vv <- table(obj.xx$DN, obj.xx$Anno.Level.11) %>%
    as.data.frame.matrix() %>%
    .[c(2, 3, 4), ]
vv <- sweep(vv, 1, rowSums(vv), "/")
vv <- sweep(vv, 2, colSums(vv), "/")
vv <- vv[, c("DN_early", "DN_re2", "DN_blast", "DN_re1", "ISP")]
pheatmap::pheatmap(vv, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, filename = file.path("reviewer.1/figs/dn1.to.3.pdf"), width = 5, height = 2.6)

#------------------------------------------------------------
# Vilidation for DP-blast 1 to 5.
cell.cycle.genes <- readxl::read_excel("../configs/cellcyle.genes.xlsx")
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(G1.S = cell.cycle.genes %>% .[, "G1/S"] %>% unlist() %>% unique() %>% .[!is.na(.)]), name = "G1.S")
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(S = cell.cycle.genes %>% .[, "S"] %>% unlist() %>% unique() %>% .[!is.na(.)]), name = "S")
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(G2.M = cell.cycle.genes %>% .[, "G2/M"] %>% unlist() %>% unique() %>% .[!is.na(.)]), name = "G2.M")
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(M = cell.cycle.genes %>% .[, "M"] %>% unlist() %>% unique() %>% .[!is.na(.)]), name = "M")
obj.merged.sc <- AddModuleScore(obj.merged.sc, feature = list(M.G1 = cell.cycle.genes %>% .[, "M/G1"] %>% unlist() %>% unique() %>% .[!is.na(.)]), name = "M.G1")
VlnPlot(obj.blast, features = c("G1.S1", "S1", "G2.M1", "M1", "M.G11"), group.by = "Anno.Level.11", pt.size = 0)
ggsave(file.path("reviewer.1/figs/dp.blast.pdf"), width = 10, height = 6)

source("modules/.DoMultiBarHeatmap.R")
obj.blast <- subset(obj.merged.sc, Anno.Level.11 %in% c("DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5"))
cycle.genes <- tidyr::gather(cell.cycle.genes)[, 2] %>%
    unlist(., use.names = F) %>%
    unique() %>%
    .[!is.na(.)]
obj.blast <- ScaleData(obj.blast)
cellcyleHeatmap(obj.blast, cell.cycle.genes, out.figs.dir)

matureStatus(obj.merged.sc, out.figs.dir)
#-------------------------------------------------------------
# T-subgroups differential trend

tcrSubgroupsTrend(obj.merged.sc, permu = 100)
ggsave(file.path(out.figs.dir, "fit.percent.four.types_scatterplot.pdf"), width = 8, height = 10)

#-------------------------------------------------------------
# Doublet analysis of the four T cell subgroups

order.cells <- c(ANNO_TREND_IDENT[1:15], "CD8T", "CD4T")
tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17:33]))
obj <- tcrSubgroup2(obj) %>% subset(., trb != "")
obj.1 <- subset(obj, t.cell_group %in% tcr.subtypes)
obj.1$Anno.Level.11 <- factor(obj.1$Anno.Level.11, levels = ANNO_TREND_IDENT)

Idents(obj.1) <- obj.1$t.cell_group
levels(obj.1) <- tcr.subtypes

# Run scrublet tool for detecting doublets.

runScrublet <- function(obj.merged) {
    require(reticulate)
    scr <- import("scrublet")
    scr.scores <- parallel::mclapply(obj.merged$orig.ident %>% unique(), function(sn) {
        obj <- subset(obj.merged, orig.ident == sn)
        count.mat <- obj %>%
            GetAssayData(., slot = "count") %>%
            as.matrix() %>%
            t()
        scrub <- scr$Scrublet(count.mat)
        scr.score <- scrub$scrub_doublets()[[1]]
        scr.score <- scr.score %>% `names<-`(colnames(obj))
        scr.doublet <- scrub$scrub_doublets()[[2]]
        scr.doublet <- scr.doublet %>% `names<-`(colnames(obj))
        cbind.data.frame(DoubletScore = scr.score, Doublet = scr.doublet) %>% `rownames<-`(colnames(obj))
        # scr.score
    }, mc.cores = 1) %>% do.call(rbind, .)
    head(scr.scores)
    obj.merged <- AddMetaData(obj.merged, scr.scores, col.name = c("DoubletScore", "Doublet"))
    obj.merged
}

obj.1 <- runScrublet(obj.1)
VlnPlot(obj.1, features = c("DoubletScore"), pt.size = 0, cols = T.CELL.GROUPS.COLORS) + geom_boxplot(width = 0.1, outlier.shape = NA) + NoLegend()
ggsave("./revise.ms/reviewer.1/figs/doublet.score.tcells.pdf", width = 5, height = 4)

obj.vv <- LoadH5Seurat(".merged.obj.UMAP.Scanpy.h5seurat")
obj.vv <- obj.vv[, colnames(obj.1)]
obj.vv@meta.data <- obj.1@meta.data
DimPlot(obj.vv, cells.highlight = subset(obj.vv, Doublet == TRUE) %>% Cells())
ggsave("./revise.ms/reviewer.1/figs/doublet.UMAP.tcells.pdf", width = 6, height = 5)

library(DoubletFinder)
library(future)

options(future.globals.maxSize = 1000000 * 1024^2)
obj.lsts <- SplitObject(obj.1, split.by = "orig.ident")
cells.singlet <- lapply(obj.lsts, function(obj) {
    plan("multiprocess", workers = 10)
    obj <- NormalizeData(object = obj, normalization.method = "LogNormalize")
    obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(object = obj, features = rownames(x = obj), vars.to.regress = c("nCount_RNA", "percent.mt", "G2M.Score", "S.Score"))
    obj <- RunPCA(object = obj, features = VariableFeatures(object = obj), verbose = FALSE)
    obj <- FindNeighbors(object = obj, dims = 1:20)
    removeDoublets(obj)
})

saveRDS(cells.singlet, file = "./revise.ms/reviewer.1/data/cells.singlet.RDS")

obj.1$DoubletFinder <- "doublet"
obj.1$DoubletFinder[rownames(obj.1@meta.data) %in% unlist(cells.singlet)] <- "singlet"

# Doublet dection

db <- import("doubletdetection")
scr.scores <- parallel::mclapply(obj.1$orig.ident %>% unique(), function(sn) {
    obj <- subset(obj.1, orig.ident == sn)
    count.mat <- obj %>%
        GetAssayData(., slot = "count") %>%
        as.matrix() %>%
        t()
    clf <- db$BoostClassifier()
    scr.doublet <- clf$fit(count.mat)$predict()
    scr.score <- clf$doublet_score()
    scr.score <- scr.score %>% `names<-`(colnames(obj))
    scr.doublet <- scr.doublet %>% `names<-`(colnames(obj))
    cbind.data.frame(DoubletScore = scr.score, Doublet = scr.doublet) %>% `rownames<-`(colnames(obj))
}, mc.cores = 1) %>% do.call(rbind, .)

obj.1 <- AddMetaData(obj.1, scr.scores, col.name = c("DoubletScore", "Doublet"))
obj.vv@meta.data <- obj.1@meta.data
DimPlot(obj.vv, cells.highlight = subset(obj.vv, Doublet == 1) %>% Cells())
ggsave("./revise.ms/reviewer.1/figs/doublet.UMAP.tcells.doubletdetection.pdf", width = 6, height = 5)

plot.df <- table(obj.vv$Doublet, obj.vv$t.cell_group)
plot.df <- sweep(plot.df, 2, colSums(plot.df), "/") %>% as.data.frame()
plot.df$Var2 <- factor(plot.df$Var2, levels = tcr.subtypes)
ggplot(plot.df, aes(x = Var2, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity")
ggsave("./revise.ms/reviewer.1/figs/doublet.barplot.doubletdetection.pdf", width = 5, height = 4)

#-----------------------------------------------------------------
sessionInfo()
