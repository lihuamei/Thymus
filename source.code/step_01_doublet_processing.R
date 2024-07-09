suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(future))
suppressMessages(library(DoubletFinder))

#--------------------------------------------------------------
# Load own modules

source('modules/utils.R')
source('modules/global_params.R')
source('modules/seurat_methods.R')
source('modules/visualization.R')

#--------------------------------------------------------------
# Create save directories

options(future.globals.maxSize = 10000 * 1024^2)
out.data.dir <- file.path('../3.results', '1.remove_doublet/data')
out.figs.dir <- file.path('../3.results', '1.remove_doublet/figs')

dir.create(file.path('../3.results', '1.remove_doublet'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#-------------------------------------------------------------
# Load thymus single-cell RNA-seq data.

merged.obj.1 <- readRDS('../transfer/merged.obj.1.rds')
new.cluster.ids <- c(
	"DP4", "DP4", "DP4", "DP4", "abT(entry)2", "DP4", "DP2", "DP4", "CD8+ T", "DP2", "DN", 
	"B", "CD8+ T", "abT(entry)1", "CD8aa", "B", "NKT", "DP1", "DP3", "Treg", "TEC", "DP2", "B", "Myeloid", 
	"TEC", "CD4+ T", "CD4+ T", "DP2", "DP4", "TEC", "TEC", "Plasma", "Ery"
)
names(new.cluster.ids) <- levels(merged.obj.1)
merged.obj.1 <- RenameIdents(merged.obj.1, new.cluster.ids)
Idents(merged.obj.1) <- factor(Idents(merged.obj.1), levels = ANNO_SC_RAW_INDENT)
merged.obj.1 <- mapAgesStatus(merged.obj.1, SampleClassify)
merged.obj.1@meta.data <- merged.obj.1@meta.data %>% mutate(Age = recode(orig.ident, !!!(SampleInfos)))
DefaultAssay(merged.obj.1) <- "RNA"

obj.lsts <- SplitObject(merged.obj.1, split.by = "orig.ident") 
cells.singlet <- lapply(obj.lsts, function(obj) {
	plan("multiprocess", workers = 10)
	obj <- NormalizeData(object = obj, normalization.method = "LogNormalize")
	obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
	obj <- ScaleData(object = obj, features = rownames(x = obj), vars.to.regress = c("nCount_RNA", "percent.mt", "G2M.Score", "S.Score"))
	obj <- RunPCA(object = obj, features = VariableFeatures(object = obj), verbose = FALSE)
	obj <- FindNeighbors(object = obj, dims = 1 : 20)
	removeDoublets(obj)	
}) 

# Remove doublet cells and save to the directory.
saveRDS(cells.singlet, file.path(out.data.dir, 'Singlet.Cells.rds'))
cells.singlet <- readRDS(file.path(out.data.dir, 'Singlet.Cells.rds'))
dublet.cells <- setdiff(Cells(merged.obj.1), unlist(cells.singlet))
message(sprintf('Number of doublet cells: %g', length(dublet.cells)))

merged.obj.2 <- merged.obj.1[, unlist(cells.singlet)]
saveRDS(merged.obj.2, file.path(out.data.dir, 'merged.obj.clean.rds'))

# Highlight show doublet cells.

raw.counts  <- table(merged.obj.1@meta.data$orig.ident)
clean.counts <- table(merged.obj.2@meta.data$orig.ident)
tab.df <- cbind.data.frame(Singlet = clean.counts, Doublet = raw.counts - clean.counts)[, c(1, 2, 4)] %>% `rownames<-`(.[, 1]) %>% .[, -1] %>% `colnames<-`(c('Singlet', 'Doublet'))

my.cols <- c('gray', 'red')
DefaultAssay(merged.obj.1) <- 'integrated'
merged.obj.1 <- merged.obj.1 %>% RunTSNE(dims = 1 : 30)
gp1 <- DimPlot(merged.obj.1, cells.highlight = dublet.cells, reduction = 'tsne') + scale_color_manual(values = my.cols, labels = c("Singlet", "Doublet"))
gp2 <- ggtexttable(tab.df, rows = rownames(tab.df), theme = ttheme("mOrange"))
gp1 + gp2
ggsave(file.path(out.figs.dir, 'Doublets.Finder.UMAP.pdf'), width = 9, height = 6)

DimPlot(merged.obj.1, cells.highlight = dublet.cells, split.by = 'AgeStatus', pt.size = 0.5, reduction = 'tsne') + 
	scale_color_manual(values = my.cols, labels = c("Singlet", "Doublet"))
ggsave(file.path(out.figs.dir, 'Doublets.Finder.Split.by.AgeStatus.UMAP.pdf'), width = 15, height = 6)

#------------------------------------------------------

sessionInfos()
