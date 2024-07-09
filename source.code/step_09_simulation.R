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

out.data.dir <- file.path('../3.results', 'simu/data')
out.figs.dir <- file.path('../3.results', 'simu/figs')

dir.create(file.path('../3.results', 'simu'))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load annotated single-cell RNA-seq data and spatial data.

obj.sc <- readRDS('../3.results/2.fine_anno/data.bak/merged.obj.anno.new.rds')
obj.sp.lst <- readRDS('../3.results/6.function_zone/data/obj.st.lst.rds')
anno <- read.table('configs/gencode.gene.info.v22.tsv', sep = '\t', header = T) %>% .[!duplicated(.[, 2]), ] %>% `rownames<-`(.[, 2])
coding.genes <- subset(anno, gene_type == 'protein_coding') %>% .[, 'gene_name'] %>% intersect(., rownames(obj.sp.lst[[1]]))

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

#----------------------------------------------------------
# Deconvolution for spot

require(CARD)
obj.sc.sub <- downSamplSeurat(obj.sc, cnt = 1000)
CARD.obj.t8 <- lapply(obj.sp.lst['T8'], function(xx) {
    sp.counts <- GetAssayData(xx, slot = 'count')
    sc.counts <- GetAssayData(obj.sc.sub, slot = 'count')
    sc.meta <- obj.sc.sub@meta.data
    sc.meta$Anno.Level.Fig.1 <- as.vector(sc.meta$Anno.Level.Fig.1)
    sc.meta$Anno.Level.Fig.1[sc.meta$Anno.Level.Fig.1 == 'abT(entry)'] <- 'abT_entry'
    sp.loc <- GetTissueCoordinates(object = xx@images[[1]]) %>% `colnames<-`(c('x', 'y'))
    sc.meta$sampleInfo <- 'ALL'
    CARD.obj <- createCARDObject(
        sc_count = sc.counts,
        sc_meta = sc.meta,
        spatial_count = sp.counts,
        spatial_location = sp.loc,
        ct.varname = "Anno.Level.Fig.1",
        ct.select = unique(sc.meta$Anno.Level.Fig.1),
        sample.varname = 'sampleInfo',
        minCountGene = 100,
        minCountSpot = 5
    )
    CARD.obj <- CARD_deconvolution(CARD_object = CARD.obj)
    CARD.obj
}) %>% `names<-`('T8')

saveRDS(CARD.obj.t8, file = file.path(out.data.dir, 'card.decon.t8.rds'))

getWeightForCell <- function(sc_eset,B){
   count = GetAssayData(sc_eset, slot = 'count')
   count = count[,colSums(count) > 0]
   count = sweep(count,2,colSums(count),"/")
   count = count[rownames(count) %in% rownames(B),]
   count = count[match(rownames(B),rownames(count)),]
   Mean_Cell <- sapply(1:ncol(count),function(icell){
        mod1 = nnls::nnls(as.matrix(B),as.matrix(count[,icell]))
        return(mod1$x)}
        )
   Mean_Cell = t(Mean_Cell)
   rownames(Mean_Cell) = colnames(count)
   colnames(Mean_Cell) = colnames(B)
   return(Mean_Cell)
}

simSpots <- function(sp.ref, sc.prop, top.n = 10) {
    cor.mat <- cor(t(sp.ref), t(sc.prop))
	sp.sim <- lapply(rownames(cor.mat), function(sp) {
        cor.vec <- cor.mat[sp, ]
        colnames(cor.mat)[order(cor.vec) %>% rev]
    }) %>% `names<-`(rownames(cor.mat))
    return(sp.sim)
}

#-----------------------------------------------------------
# Simulation based on Thy7 samples

obj.sp <- obj.sp.lst[['T8']]
obj.sc.train <- downSamplSeurat(obj.sc, cnt = 500)
ovp.genes <- intersect(rownames(obj.sc.train), rownames(obj.sp))
obj.sp <- obj.sp[ovp.genes, ]
obj.sc.train <- obj.sc.train[ovp.genes, ]
sp.ref <- CARD.obj.t8$T8@Proportion_CARD
saveRDS(obj.sc.train, file = file.path(out.data.dir, 'obj.sc.train.rds'))

sc.prop <- getWeightForCell(obj.sc.train, CARD.obj.t8$T8@algorithm_matrix$B)
sc.prop <- sweep(sc.prop, 1, rowSums(sc.prop), '/')
sim.res <- simSpots(sp.ref, sc.prop)	
for (ncell in c(5, 10, 15, 20, 25, 30, 35, 40)) {
	ref.ctpes <- Idents(obj.sc.train)
	prop.vec <- rep(0, levels(obj.sc) %>% length) %>% `names<-`(levels(obj.sc))	
	true.props <- lapply(sim.res, function(sim.vec) {
		ctype <- ref.ctpes[sim.vec]
		sim.prop <- table(ctype) / length(ctype)	
		prop.vec[names(sim.prop)] <- sim.prop
		prop.vec
	}) %>% do.call(rbind, .) %>% as.data.frame
	saveRDS(true.props, file = file.path(out.data.dir, sprintf('true.props.%g.rds', ncell)))

	sc.expr <- GetAssayData(obj.sc.train, slot = 'count')
	sim.expr <- lapply(sim.res, function(sim.vec) {
		sc.expr[, sim.vec] %>% rowSums(.)		
	}) %>% do.call(cbind, .) %>% as.data.frame
	idxes <- sample(nrow(sim.expr), floor(nrow(sim.expr) * ncell * 0.01))
	rand.shuffle <- sim.expr[idxes, ]
	rand.shuffle <- rand.shuffle[, sample(ncol(rand.shuffle))]
	sim.expr[idxes, ] <- rand.shuffle
	obj.sim <- CreateSeuratObject(counts = sim.expr, project = 'Sim', assay = 'Spatial') %>% NormalizeData
	obj.sim@images <- obj.sp@images
	obj.sim@meta.data <- obj.sp@meta.data
	saveRDS(obj.sim, file = file.path(out.data.dir, sprintf('obj.T8.sim.%g.rds', ncell)))
}
#-----------------------------------------------------------
# SrtCT

true.prop.5 <- readRDS(file.path(out.data.dir, 'true.props.5.rds'))
true.prop.10 <- readRDS(file.path(out.data.dir, 'true.props.10.rds'))
true.prop.15 <- readRDS(file.path(out.data.dir, 'true.props.15.rds'))
true.prop.20 <- readRDS(file.path(out.data.dir, 'true.props.20.rds'))
true.prop.25 <- readRDS(file.path(out.data.dir, 'true.props.25.rds'))
true.prop.30 <- readRDS(file.path(out.data.dir, 'true.props.30.rds'))
true.prop.35 <- readRDS(file.path(out.data.dir, 'true.props.35.rds'))
true.prop.40 <- readRDS(file.path(out.data.dir, 'true.props.40.rds'))
true.props <- list(
	x5 = true.prop.5, 
	x10 = true.prop.10, 
	x15 = true.prop.15, 
	x20 = true.prop.20, 
	x25 = true.prop.25, 
	x30 = true.prop.30, 
	x35 = true.prop.35, 
	x40 = true.prop.40
)
obj.srtct <- readRDS(file.path(out.data.dir, 'obj.SrtCT.rds'))
pred.srtct <- lapply(obj.srtct, function(xx) xx@meta.data[, 'predicted.id']) %>% do.call(cbind, .) %>% as.data.frame %>% `rownames<-`(rownames(true.prop.5))

srtct.acc <- lapply(1 : 8, function(pidx) {
    tt <- true.props[[pidx]] * 4 * 5
    nums <- colSums(tt)
    cut.off <- pidx * 5 * 0.1
    true.lst <- lapply(rownames(tt), function(sp) {
        colnames(tt)[tt[sp, ] >= 2]
    }) %>% `names<-`(rownames(tt))
    pred.lst <- split(obj.srtct[[pidx]]$predicted.id, obj.srtct[[pidx]] %>% colnames)
    names(pred.lst)[names(pred.lst) == "abT_entry"] <- "abT(entry)"
	pred.lst <- lapply(pred.lst, function(ff) table(ff)[table(ff) >= 1] %>% names)
    accurate <- lapply(names(true.lst), function(ct) {
        length(intersect(true.lst[[ct]], pred.lst[[ct]])) / length(true.lst[[ct]])
    }) %>% unlist
}) %>% do.call(cbind, .) %>% as.data.frame %>% `colnames<-`(c('x5', 'x10', 'x15', 'x20', 'x25', 'x30', 'x35', 'x40'))

#-------------------------------------------------------------------
# CellTrek
obj.celltrek <- readRDS(file.path(out.data.dir, 'obj.celltrek.rds'))
sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c('x', 'y'))
celltrek.maplist <- lapply(1 : 8, function(kk) {
	celltrek.loc <- GetTissueCoordinates(obj.celltrek[[kk]])
	map.sp <- parallel::mclapply(1 : nrow(celltrek.loc), function(idx) {
		coord.1 <- celltrek.loc[idx, ] %>% unlist 		
		dist.vec <- apply(sp.loc, 1, function(coord.2) {
			df <- cbind.data.frame(coord.1, coord.2) %>% t
			dist(df) 
		})	%>% { rownames(sp.loc)[order(.)[1]] }
	}, mc.cores = 20) %>% unlist %>% `names<-`(rownames(celltrek.loc))
	
	sc.ctypes <- Idents(obj.sc)
	ctypes <- sc.ctypes[obj.celltrek[[kk]]@meta.data[, 'id_raw']]
	cbind.data.frame(Spot = map.sp, CellTypes = ctypes, RawID = obj.celltrek[[kk]]@meta.data[, 'id_raw'])
})
saveRDS(celltrek.maplist, file = file.path(out.data.dir, 'celltrek.maplist.rds'))

celltrek.maplist <- readRDS(file.path(out.data.dir, 'celltrek.maplist.rds'))
celltrek.acc <- lapply(1 : 8, function(pidx) {
    tt <- true.props[[pidx]] * 4 * 5
    nums <- colSums(tt)
    cut.off <- pidx * 5 * 0.1
    true.lst <- lapply(rownames(tt), function(sp) {
        colnames(tt)[tt[sp, ] >= 2]
    }) %>% `names<-`(rownames(tt))
    pred.lst <- split(celltrek.maplist[[pidx]]$CellTypes, celltrek.maplist[[pidx]]$Spot)
    names(pred.lst)[names(pred.lst) == "abT_entry"] <- "abT(entry)"
	pred.lst <- lapply(pred.lst, function(ff) table(ff)[table(ff) >= 1] %>% names)
    accurate <- lapply(names(true.lst), function(ct) {
        length(intersect(true.lst[[ct]], pred.lst[[ct]])) / length(true.lst[[ct]])
    }) %>% unlist
}) %>% do.call(cbind, .) %>% as.data.frame %>% `colnames<-`(c('x5', 'x10', 'x15', 'x20', 'x25', 'x30', 'x35', 'x40'))

#-------------------------------------------------------------------
# TSO-hismap

obj.tsohis <- readRDS(file.path(out.data.dir, 'sc.map.lst.simu.tso.rds'))
obj.sp <- obj.sp.lst[['T8']]
coord <- GetTissueCoordinates(obj.sp)
obj.tsohis <- lapply(obj.tsohis, function(xx) {
	qry.label <- paste0(colData(xx)$centerx, '_', colData(xx)$centery)	
	ref.label <- paste0(coord[, 1], '_', coord[, 2])
	ref.ser <- rownames(coord) %>% `names<-`(ref.label)
	spots <- ref.ser[ref.label]
	xx$SpotName <- spots %>% as.vector
	xx
})

card.acc <- lapply(1 : 8, function(pidx) {
    tt <- true.props[[pidx]] * 4 * 5
    nums <- colSums(tt)
    cut.off <- pidx * 5 * 0.1
	true.lst <- lapply(rownames(tt), function(sp) {
        colnames(tt)[tt[sp, ] >= 2]
    }) %>% `names<-`(rownames(tt))
    pred.lst <- split(obj.tsohis[[pidx]]$CT, obj.tsohis[[pidx]]$SpotName)
	pred.lst <- lapply(pred.lst, function(ff) table(ff)[table(ff) >= 1] %>% names)
    accurate <- lapply(names(true.lst), function(ct) {
		length(intersect(true.lst[[ct]], pred.lst[[ct]])) / length(true.lst[[ct]])
    }) %>% unlist
}) %>% do.call(cbind, .) %>% as.data.frame %>% `colnames<-`(c('x5', 'x10', 'x15', 'x20', 'x25', 'x30', 'x35', 'x40'))

#------------------------------------------------------------------
# TSO-hismap

obj.tsohis <- readRDS(file.path(out.data.dir, 'sc.map.lst.simu.tso.rds'))
tsohis.acc <- lapply(1 : 8, function(pidx) {
    tt <- true.props[[pidx]] * 4 * 5
    nums <- colSums(tt)
    cut.off <- pidx * 5 * 0.1
    true.lst <- lapply(rownames(tt), function(sp) {
        colnames(tt)[tt[sp, ] >= 2]
    }) %>% `names<-`(rownames(tt))
    pred.lst <- split(obj.tsohis[[pidx]]$CT, obj.tsohis[[pidx]]$SpotName)
    pred.lst <- lapply(pred.lst, function(ff) table(ff)[table(ff) >= 1] %>% names)
    accurate <- lapply(names(true.lst), function(ct) {
        length(intersect(true.lst[[ct]], pred.lst[[ct]])) / length(true.lst[[ct]])
    }) %>% unlist
}) %>% do.call(cbind, .) %>% as.data.frame %>% `colnames<-`(c('x5', 'x10', 'x15', 'x20', 'x25', 'x30', 'x35', 'x40'))

#-------------------------------------------------------------------
# Plot

srtct.acc.p <- cbind.data.frame(srtct.acc, Method = 'SrtCT', CT = rownames(srtct.acc)) %>% tidyr::gather(., 'Group', 'Acc', -Method, -CT)
celltrek.acc.p <- cbind.data.frame(celltrek.acc, Method = 'CellTrek', CT = rownames(srtct.acc)) %>% tidyr::gather(., 'Group', 'Acc', -Method, -CT)
tsohis.acc.p <- cbind.data.frame(tsohis.acc, Method = 'TSO-hismap', CT = rownames(srtct.acc)) %>% tidyr::gather(., 'Group', 'Acc', -Method, -CT)
card.acc.p <- cbind.data.frame(card.acc, Method = 'CARD', CT = rownames(srtct.acc)) %>% tidyr::gather(., 'Group', 'Acc', -Method, -CT)

plot.data <- rbind.data.frame(srtct.acc.p, celltrek.acc.p, tsohis.acc.p, card.acc.p)
plot.data$Method <- factor(plot.data$Method, levels = c("SrtCT", "CellTrek", "CARD", "TSO-hismap"))

plt <- ggstatsplot::ggbetweenstats(
	data = subset(plot.data, Group == 'x20'),
	x = Method,
	p.adjust.method = 'none',
	y = Acc
)
ggsave(file.path(out.figs.dir, 'sim.performace.20.pdf'), width = 5, height = 3)

plt <- ggstatsplot::ggbetweenstats(
	data = subset(plot.data, Group == 'x15'),
	x = Method,
	p.adjust.method = 'none',
	y = Acc
)
ggsave(file.path(out.figs.dir, 'sim.performace.15.pdf'), width = 5, height = 3)

plt <- ggstatsplot::ggbetweenstats(
    data = subset(plot.data, Group == 'x10'),
    x = Method,
	p.adjust.method = 'none',
    y = Acc
)
ggsave(file.path(out.figs.dir, 'sim.performace.10.pdf'), width = 5, height = 3)

plt <- ggstatsplot::ggbetweenstats(
    data = subset(plot.data, Group == 'x5'),
    x = Method,
	p.adjust.method = 'none',
    y = Acc
)
ggsave(file.path(out.figs.dir, 'sim.performace.5.pdf'), width = 5, height = 3)


plt <- ggstatsplot::ggbetweenstats(
    data = subset(plot.data, Group == 'x35'),
    x = Method,
	p.adjust.method = 'none',
    y = Acc
)
ggsave(file.path(out.figs.dir, 'sim.performace.35.pdf'), width = 5, height = 3)

