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
suppressMessages(library(monocle))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(enrichplot))

#--------------------------------------------------------------
# Load own modules

source("modules/utils.R")
source("modules/global_params.R")
source("modules/seurat_methods.R")
source("modules/visualization.R")
source("modules/findMedullaClusters.R")
source("modules/trajectory_methods.R")

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path("../3.results", "4.proj_spatial/data")
out.figs.dir <- file.path("../3.results", "4.proj_spatial/figs")

dir.create(file.path("../3.results", "4.proj_spatial"))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load spatial data and single-cell RNA-seq data for normal samples.

anno <- read.table("configs/gencode.gene.info.v22.tsv", sep = "\t", header = T) %>%
    .[!duplicated(.[, 2]), ] %>%
    `rownames<-`(.[, 2])
coding.genes <- subset(anno, (gene_type %in% ANNO.AREAS) & (gene_status == "KNOWN")) %>% .[, "gene_name"]
obj.merged.sc <- readRDS(file.path(out.data.dir, "obj.sc.Fb.RDS"))
cells.1 <- subset(obj.merged.sc, idents = levels(obj.merged.sc)[c(1:15, 33)])
cells.1 <- subset(cells.1, tcr == 0) %>% Cells()
cells.2 <- subset(obj.merged.sc, idents = levels(obj.merged.sc)[16:32])
cells.2 <- subset(cells.2, tcr == 1) %>% Cells()
obj.merged.sc <- subset(obj.merged.sc, cells = c(cells.1, cells.2))

obj.merged.sc <- downSamplSeurat(obj.merged.sc, cnt = 1000)

#----------------------------------------------------------------
# Merged data and set levels.
cellChatAnlysisNew <- function(obj, group.by = "Anno.Level.Fig.1", db.chat = CellChatDB.human) {
    data.input <- GetAssayData(obj, slot = "data")
    meta.data <- obj@meta.data
    cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = group.by)
    cellchat <- addMeta(cellchat, meta = meta.data)
    cellchat <- setIdent(cellchat, ident.use = group.by)
    cellchat@DB <- db.chat

    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 10)
    cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.1, thresh.fc = 0.1)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 2)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    return(cellchat)
}

require(CellChat)
require(patchwork)
options(stringsAsFactors = FALSE)

options(future.globals.maxSize = 40000 * 1024^2)
sub.cellx <- c("DN_early", "DN_blast", "DN_re", "DP_blast", "DP_re", "abT(entry)", "Fb", "VSMC", "cTEC", "mTEC", "CD4T_mem", "CD8T_mem")
obj.sub <- subset(obj.merged.sc, Anno.Level.Fig.1 %in% sub.cellx)
Idents(obj.sub) <- obj.sub$Anno.Level.Fig.1

out.data.dir <- "pbs/fb.revise"
obj.sub$CellChat <- factor(obj.sub$Anno.Level.Fig.1)
cellchat.fb <- cellChatAnlysisNew(obj.sub, group.by = "CellChat", db.chat = CellChatDB.human)

#----------------------------------------------------------------
# Variation analysis

sm.pct <- seq(50, 100, 5) / 100
parallel::mclapply(sm.pct, function(pct) {
    vv <- downSamplSeurat(obj.sub, percent = pct)
    vv$CellChat <- factor(vv$Anno.Level.Fig.1)
    cellchat.fb <- cellChatAnlysisNew(vv, group.by = "CellChat", db.chat = CellChatDB.human)
    saveRDS(cellchat.fb, file = file.path(out.data.dir, sprintf("cellchat.fb.sub.%g.rds", pct)))
}, mc.cores = 5)
