#-----------------------------------------------------
# UMI correlation between SP and ST

umiCorrSpWithSt <- function(obj.sc, obj.sp.lst) {
    figs.dir <- file.path(out.figs.dir, "sc_cor_sp")
    dir.create(figs.dir, showWarnings = FALSE)

    gp.lst <- lapply(names(SC_MAP_SP), function(sn) {
        obj.sc.sub <- subset(obj.sc, orig.ident == sn)
        obj.sc.umi <- GetAssayData(obj.sc.sub, slot = "count") %>% rowSums(.)
        obj.st.sub <- obj.sp.lst[[SP_NAMES[[SC_MAP_SP[[sn]]]]]]
        obj.st.umi <- GetAssayData(obj.st.sub, slot = "count") %>% rowSums(.)

        ovp.genes <- intersect(names(obj.sc.umi), names(obj.st.umi))
        plot.dat <- cbind.data.frame(SC = obj.sc.umi[ovp.genes], SP = obj.st.umi[ovp.genes])
        gp <- ggplot(plot.dat, aes(x = log2(SC), y = log2(SP))) +
            geom_point(color = "blue", alpha = 0.3, size = 1.5, shape = 16) +
            geom_smooth(method = lm, color = "red", size = 0.8, alpha = 0.6) +
            theme_classic(base_size = 12) +
            stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method = "spearman") +
            ggtitle(sprintf("%s", sn))
        ggsave(file.path(figs.dir, sprintf("%s_vs_%s.scatter.pdf", sn, SP_NAMES[[SC_MAP_SP[[sn]]]])), width = 6, height = 6)
        gp
    })
    ggarrange(plotlist = gp.lst, ncol = 4, nrow = 2)
}

#---------------------------------------------------
# Count spots of each sample.

countSpotsPerSample <- function(obj.sp.lst) {
    plot.df <- lapply(obj.sp.lst, dim) %>%
        as.data.frame() %>%
        .[2, ] %>%
        as.data.frame() %>%
        t() %>%
        `colnames<-`("Spots")
    plot.df <- as.data.frame(plot.df)
    plot.df$Sample <- sp2SampleName(rownames(plot.df))
    plot.df$Sample <- factor(plot.df$Sample, levels = names(SC_MAP_SP))
    cols <- SAMPLE_COLORS %>% `names<-`(sp2SampleName(names(SAMPLE_COLORS)))
    ggplot(plot.df, aes(x = Sample, y = Spots, fill = Sample), alpha = 0.5) +
        geom_bar(stat = "identity") +
        theme_bw(base_size = 14) +
        theme(
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 3100)) +
        xlab("Spatial transcriptomics samples") +
        ylab("Number of spots") +
        scale_fill_manual(values = cols) +
        geom_hline(yintercept = 2000, linetype = "dotted", color = "black", size = 1.5) +
        geom_hline(yintercept = 2500, linetype = "dotted", color = "red", size = 1.5) +
        geom_text(aes(label = Spots), vjust = -0.5) +
        scale_x_discrete(labels = paste0(plot.df$Sample %>% levels(), " (", SampleInfos[plot.df$Sample %>% levels()] %>% unlist(), ")"))
}

#--------------------------------------------------
# Images plot

showAllImages <- function(obj.sp.lst) {
    names(obj.sp.lst) <- sp2SampleName(names(obj.sp.lst))
    obj.sp.lst <- obj.sp.lst[names(SC_MAP_SP)]
    gps <- lapply(names(obj.sp.lst), function(sn) {
        obj <- obj.sp.lst[[sn]]
        SpatialFeaturePlot(obj, feature = "CCR9", alpha = 0) + NoLegend() +
            ggtitle(sn) +
            theme(plot.title = element_text(hjust = 0.5))
    }) %>% `names<-`(names(obj.sp.lst))
    ggarrange(plotlist = gps, nrow = 2, ncol = 4)
}

#--------------------------------------------------
# QC plot

QCPlots <- function(objs) {
    names(objs) <- sp2SampleName(names(objs))
    objs <- objs[names(SC_MAP_SP)]
    meta.infos <- lapply(names(objs), function(sn) {
        obj <- objs[[sn]]
        meta <- cbind.data.frame(obj@meta.data, Label = sn)
    }) %>% do.call(rbind, .)

    cols <- SAMPLE_COLORS %>% `names<-`(sp2SampleName(names(SAMPLE_COLORS)))
    meta.infos$Label <- factor(meta.infos$Label, levels = names(SC_MAP_SP))
    gp1 <- ggplot(meta.infos, aes(x = Label, y = log2(1 + nCount_Spatial), fill = Label)) +
        geom_violin() +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        xlab("") +
        ylab("Number of UMIs") +
        scale_fill_manual(values = cols) +
        scale_y_continuous(labels = c("0", "16", "256", "4096", "16384", ""), breaks = c(0, 4, 8, 12, 14, 16)) +
        ggtitle("UMI")

    gp2 <- ggplot(meta.infos, aes(x = Label, y = log2(1 + nFeature_Spatial), fill = Label)) +
        geom_violin() +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        xlab("") +
        ylab("Number of genes") +
        scale_fill_manual(values = cols) +
        scale_y_continuous(labels = c("", "16", "256", "4096"), breaks = c(0, 4, 8, 12)) +
        ggtitle("Genes")

    gp3 <- ggplot(meta.infos, aes(x = nCount_Spatial, color = Label)) +
        geom_density(size = 1.0) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(color = "black"),
            legend.position = "none"
        ) +
        ggtitle("UMI") +
        ylab("Density") +
        scale_color_manual(values = cols) +
        scale_x_continuous(labels = c("0", "20k", "40k", "60k"), breaks = c(0, 20000, 40000, 60000))

    gp4 <- ggplot(meta.infos, aes(x = nFeature_Spatial, color = Label)) +
        geom_density(size = 1.0) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(color = "black")
        ) +
        ggtitle("Gene") +
        ylab("Density") +
        scale_color_manual(values = cols) +
        scale_x_continuous(labels = c("0", "2k", "4k", "6k", "8k"), breaks = c(0, 2000, 4000, 6000, 8000))
    ggarrange(gp1, gp2, gp3, gp4, nrow = 2, ncol = 2)
}

#-------------------------------------------------------------
# Plot patterns of cell types.

celltypePatterns <- function(hc) {
    require(ggtree)
    require(ggnewscale)
    gt <- ggtree(hc, layout = "fan", size = 0.5) + geom_tiplab(size = 4, aes(angle = 90), offset = 0.4)
    gt <- gt %<+% anno.df + geom_tippoint(aes(color = cellType), alpha = 1, size = 5) +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        theme(legend.position = "none")
}

#------------------------------------------------------------
# SC vs. SP correlation

corrSpVsSCMappedCtypes <- function(pseud.coord, obj.sc) {
    pseud.coord.df <- do.call(rbind, pseud.coord)
    count.sp.ctypes <- pseud.coord.df$cellType %>% table()
    count.sc.ctypes <- Idents(obj.sc) %>% table()
    count.sp.ctypes <- count.sp.ctypes[names(count.sc.ctypes)]
    count.sp.ctypes <- count.sp.ctypes / sum(count.sp.ctypes)
    count.sc.ctypes <- count.sc.ctypes / sum(count.sc.ctypes)
    plot.df <- cbind.data.frame(SP = count.sp.ctypes %>% as.vector(), SC = count.sc.ctypes %>% as.vector(), CellType = names(count.sp.ctypes))
    plot.df$CellType <- factor(plot.df$CellType, levels = ANNO_ENTIRE_IDNET_FIG1)
    ggplot(plot.df, aes(x = SC, y = SP)) +
        geom_point(size = 4, aes(color = CellType)) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        geom_smooth(method = "lm", se = TRUE) +
        stat_cor(method = "pearson") +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        xlab("Single-cell data") +
        ylab("Spatial-spot data (with cell segmentation)")
}

#----------------------------------------------------------
# Project cell-type to spatial data.

projCellTypeToSpatial <- function(pseud.coord, obj = NULL, title = NULL, size = 0.3, shape = 21) {
    idxes <- which(!duplicated(pseud.coord$cellType))
    pseud.coord$Label <- ""
    pseud.coord$Label[idxes] <- pseud.coord$cellType[idxes]
    g1 <- ggplot(pseud.coord, aes(x = y, y = x, color = cellType), alpha = 0.3) +
        geom_point(size = size, shape = shape) +
        ggrepel::geom_text_repel(data = subset(pseud.coord, Label != ""), aes(label = Label), color = "white", box.padding = 1.2) +
        ggplot2::scale_y_reverse() +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = "black"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        ylab("") +
        xlab("") +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        guides(colour = guide_legend(override.aes = list(size = 4)))

    if (is.null(title)) {
        g1 <- g1 + ggtitle(sp2SampleName(obj@images %>% names()))
    } else {
        g1 <- g1 + ggtitle(title)
    }
    return(g1)
}

#------------------------------------------------------
# Medulla vs. Cortex

cortexMeduallaCorrHclust <- function(obj.sp.lst, feature.genes) {
    obj <- lapply(obj.sp.lst, function(xx) xx[feature.genes, ])
    obj <- lapply(obj, function(xx) {
        xx <- subset(xx, HE.Labels != "Medulla_lo")
        xx@meta.data$Class <- xx@meta.data$HE.Labels
        xx@meta.data$Class[xx@meta.data$Class %in% c("Medulla_edge", "Medulla_hi", "Medulla_centric")] <- "Medulla"
        xx@meta.data$Class[xx@meta.data$Class %in% c("Cortex")] <- "Cortex"
        xx@meta.data$Class.new <- paste0(xx@images %>% names(), "_", xx@meta.data$Assign.Centric, "_", xx@meta.data$Class)
        xx
    })
    expr.corr <- lapply(obj, function(xx) {
        AverageExpression(xx, group.by = "Class.new")$Spatial
    }) %>%
        do.call(cbind, .) %>%
        cor(.)

    anno.df <- lapply(obj, function(xx) {
        cbind.data.frame(
            Class.new = xx@meta.data$Class.new,
            Class = xx@meta.data$Class, Centric = xx@meta.data$Assign.Centric, Sample = sp2SampleName(xx@images %>% names())
        ) %>% .[!duplicated(.[, 1]), ]
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(.[, "Class.new"]) %>%
        .[colnames(expr.corr), ]
    anno.df$Sample <- factor(anno.df$Sample, levels = names(SC_MAP_SP))
    cols <- SAMPLE_COLORS
    names(cols) <- sp2SampleName(names(SAMPLE_COLORS))
    hc <- as.dist(1 - expr.corr) %>% upgma(., method = "ward.D")
    gt <- ggtree(hc, layout = "circular", size = 0.5)
    gt <- gt %<+% anno.df + geom_tippoint(aes(color = Class), alpha = 0.5, size = 4) + scale_color_manual(values = CM_COLORS[c("Cortex", "Medulla")]) +
        guides(colour = guide_legend(override.aes = list(size = 12, title = "HE.Position")))
    gt1 <- gheatmap(gt, anno.df[, "Sample", drop = FALSE], width = 0.1, color = "black") + scale_fill_manual(values = cols[names(SC_MAP_SP)]) +
        guides(colour = guide_legend(override.aes = list(title = "SampleID")))
}

#----------------------------------------------------
# Cortex and Medulla spots

cortexMedullaSpotStat <- function(obj.sp.lst) {
    spots.module.cnt <- lapply(obj.sp.lst, function(xx) {
        table(xx@meta.data$Assign.Centric, xx@meta.data$HE.Labels) %>%
            as.data.frame.matrix() %>%
            `rownames<-`(paste0(xx@images %>% names(), "_", rownames(.)))
    }) %>%
        do.call(rbind, .) %>%
        .[, 1:4]
    count.df <- cbind.data.frame(Cortex = spots.module.cnt[, "Cortex"], Medulla = rowSums(spots.module.cnt[, c("Medulla_edge", "Medulla_hi")]))
    rownames(count.df) <- gsub("T.*\\.", "", rownames(count.df))
    count.df <- count.df[order(rownames(count.df)), ]

    plot.df <- sweep(count.df, 1, rowSums(count.df), "/") %>% cbind.data.frame(., Label = rownames(.))
    sn.names <- sp2SampleName(gsub("_.*", "", plot.df$Label))
    plot.df$SN <- sp2SampleName(gsub("_.*", "", plot.df$Label))
    plot.df <- plot.df[str_order(plot.df$SN, numeric = T), ]
    plot.df$Label <- paste0(plot.df$SN, "_", gsub("T.*_", "", plot.df$Label))
    plot.df$Label <- factor(plot.df$Label, levels = plot.df$Label)
    gp1 <- ggplot(tidyr::gather(plot.df, "key", "value", -Label, -SN), aes(x = Label, y = value, fill = key)) +
        geom_bar(position = "stack", stat = "identity") +
        theme_classic(base_size = 16) +
        theme(
            axis.text.x = element_blank(),
            legend.position = "top"
        ) +
        ylab("Fraction") +
        scale_fill_manual(values = CM_COLORS %>% unlist() %>% .[c("Cortex", "Medulla")]) +
        xlab("") +
        geom_hline(yintercept = plot.df$Medulla %>% mean(), linetype = "dashed", color = "yellow") +
        scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), expand = c(0, 0))

    # Number of spots in each compartment
    plot.df.2 <- cbind.data.frame(count.df, Summary = rowSums(count.df), Label = rownames(count.df))
    sn.names <- sp2SampleName(gsub("_.*", "", plot.df.2$Label))
    plot.df.2$SN <- sp2SampleName(gsub("_.*", "", plot.df.2$Label))
    plot.df.2 <- plot.df.2[str_order(plot.df.2$SN, numeric = T), ]
    plot.df.2$Label <- paste0(plot.df.2$SN, "_", gsub("T.*_", "", plot.df.2$Label))
    plot.df.2$Label <- factor(plot.df.2$Label, levels = plot.df.2$Label)

    gp2 <- ggplot(plot.df.2, aes(x = Label, y = Summary, fill = "blue", alpha = 0.5)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Summary, y = Summary), vjust = -.5) +
        theme_classic(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "none"
        ) +
        ylab("Number of spots in each compartment") +
        xlab("") +
        scale_y_continuous(expand = c(0, 0), limit = c(0, 750), breaks = c(0, floor(median(plot.df.2$Summary)), 250, 500, 750)) +
        geom_hline(yintercept = floor(median(plot.df.2$Summary)), linetype = "dashed", color = "blue")
    ggarrange(gp1, gp2, nrow = 2, heights = c(2, 3.5))
}

#----------------------------------------------
# Assign spots to differnet centric

assignSpots2Centre <- function(obj.st.lst) {
    names(obj.st.lst) <- sp2SampleName(names(obj.st.lst))
    obj.st.lst <- obj.st.lst[names(SC_MAP_SP)]
    gp.lst <- lapply(obj.st.lst, function(xx) {
        img.name <- xx@images %>% names()
        gp1 <- SpatialPlot(xx, group.by = "Assign.Centric", image.alpha = 0) + NoLegend() + ggtitle(sp2SampleName(img.name))

        cols <- c(Cortex = "#00B6EB", Medulla_centric = "green", Medulla_edge = "white", Medulla_hi = "#FB61D7")
        centric.spots <- subset(xx, Centric == "Y") %>% Cells()
        xx@meta.data[centric.spots, "HE.Labels"] <- "Medulla_centric"
        low.cells <- subset(xx, HE.Labels == "Medulla_lo") %>% Cells()
        xx@meta.data[low.cells, "HE.Labels"] <- "Cortex"
        xx.sub <- subset(xx, HE.Labels != "Cortex")
        xx.sub$HE.Labels <- factor(xx.sub$HE.Labels, levels = c("Medulla_edge", "Medulla_hi", "Medulla_centric"))
        cols <- cols[c("Medulla_edge", "Medulla_hi", "Medulla_centric")]
        gp2 <- SpatialPlot(xx.sub, group.by = "HE.Labels", image.alpha = 1, cols = cols) + guides(fill = guide_legend(override.aes = list(size = 3)))
        ggarrange(gp2, gp1, ncol = 2)
    })
    ggarrange(plotlist = gp.lst, ncol = 2, nrow = 4)
}

#---------------------------------------------
# Show an example for Thy5

showExampleThyScalerSp <- function(obj, centric = "GAGTACGGGTATACAA-1") {
    xx <- obj
    img.name <- xx@images %>% names()
    gp1 <- SpatialPlot(xx, group.by = "Assign.Centric", image.alpha = 0) + NoLegend() + ggtitle(sp2SampleName(img.name))

    cols <- c(Cortex = "#00B6EB", Medulla_centric = "green", Medulla_edge = "white", Medulla_hi = "#FB61D7")
    centric.spots <- subset(xx, Centric == "Y") %>% Cells()
    xx@meta.data[centric.spots, "HE.Labels"] <- "Medulla_centric"
    low.cells <- subset(xx, HE.Labels == "Medulla_lo") %>% Cells()
    xx@meta.data[low.cells, "HE.Labels"] <- "Cortex"
    xx.sub.cells.1 <- subset(xx, HE.Labels != "Cortex") %>% Cells()
    xx.sub.cells.2 <- subset(xx, HE.Labels == "Cortex" & Assign.Centric %in% centric) %>% Cells()
    xx.sub <- xx[, c(xx.sub.cells.1, xx.sub.cells.2)]
    xx.sub$HE.Labels <- factor(xx.sub$HE.Labels, levels = c("Cortex", "Medulla_edge", "Medulla_hi", "Medulla_centric"))
    cols <- cols[c("Cortex", "Medulla_edge", "Medulla_hi", "Medulla_centric")]
    gp2 <- SpatialPlot(xx.sub, group.by = "HE.Labels", image.alpha = 1, cols = cols) + guides(fill = guide_legend(override.aes = list(size = 3)))
    gp3 <- SpatialFeaturePlot(obj, feature = "medulla.score1", image.alpha = 0)
    gp3 + gp2 + gp1
}

#---------------------------------------------
# Scoring cell types using marker genes.

scoreCellTypePlots <- function(obj, cell.type, smooth = FALSE, method = "ks", adjust = 1) {
    require(Nebulosa)
    w <- obj@meta.data[, cell.type]
    x <- Seurat::GetTissueCoordinates(object = obj@images[[1]])
    if (smooth == TRUE) {
        w <- Nebulosa:::calculate_density(w, x, method, adjust)
    }
    edge.spots <- ifelse(obj$HE.Labels == "Medulla_edge", "Edge", "Other")
    plot.dat <- cbind.data.frame(x, expr = w, HE = edge.spots)
    plot.dat$expr <- scale(plot.dat$expr)
    cut.off <- quantile(plot.dat$expr, c(0.01, 0.99))
    plot.dat$expr[plot.dat$expr >= cut.off[2]] <- cut.off[2]
    plot.dat$expr[plot.dat$expr <= cut.off[1]] <- cut.off[1]

    ggplot(plot.dat, aes(x = imagecol, y = imagerow)) +
        geom_point(size = 3, aes(fill = expr), shape = 21, stroke = 0) +
        ggplot2::scale_y_reverse() +
        theme_bw() +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        ylab("") +
        xlab("") +
        scale_fill_viridis_c(option = "viridis") +
        # scale_fill_gradient2('expr', low = "#762A83", mid = "grey", high = "#1B7837") +
        ggtitle(cell.type) +
        geom_point(data = subset(plot.dat, HE == "Edge"), aes(color = HE, shape = HE)) +
        scale_color_manual(values = "white")
}

#-----------------------------------------------
# Distance distribution of cell-types.

distOfCellTypes <- function(res.dist) {
    order.celltypes <- aggregate(Centrici.Dist ~ cellType, res.dist, "median") %>% {
        .[order(.$Centrici.Dist), ]
    }
    res.dist$cellType <- factor(res.dist$cellType, levels = order.celltypes$cellType)
    ggplot(res.dist, aes(x = cellType, y = Centrici.Dist, color = cellType)) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        ) +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        ylab("Distance to medullary centroid") +
        xlab("")
}

#------------------------------------------------
# Co-localization of cell-types.

coExistsCelltypes <- function(obj.sp.lst) {
    corr.merged.1 <- lapply(obj.sp.lst, function(xx) {
        ff <- smoothSpatialData(xx, ANNO_ENTIRE_IDNET_FIG1, method = "wkde")
        ff
    }) %>% do.call(rbind, .)

    corr.merged <- lapply(obj.sp.lst, function(xx) {
        cell.fracs <- xx@meta.data[, ANNO_ENTIRE_IDNET_FIG1] %>% scale()
    }) %>% do.call(rbind, .) # %>% cor(.)

    corr.merged[, "T_proliferating"] <- corr.merged.1[, "T_proliferating"]
    corr.merged <- cor(corr.merged)
    corr.merged <- corr.merged[setdiff(rownames(corr.merged), c("Ery")), setdiff(rownames(corr.merged), c("Ery"))]
    corrplot(corr.merged, method = "color", order = "hclust", tl.col = "black", hclust.method = "ward.D2", insig = "pch", addrect = 3, rect.col = "blue", col = colorRampPalette(c("grey80", "white", "red"))(100))
}

#-------------------------------------------------
# Projection using Seurat

predictBySeurat <- function(obj.lst) {
    gps <- lapply(obj.lst, function(obj) {
        ovp.cells <- intersect(ANNO_ENTIRE_IDNET_FIG1, obj$predicted.id %>% unique())
        obj$predicted.id <- factor(obj$predicted.id, levels = ovp.cells)
        SpatialPlot(obj, group.by = "predicted.id", cols = ANNO_ENTIRE_COLOR_FIG1[ovp.cells], image.alpha = 0) +
            guides(fill = guide_legend(override.aes = list(size = 3))) +
            ggtitle(sp2SampleName(obj@images %>% names()))
    })
    ggarrange(plotlist = gps, ncol = 3, nrow = 3)
}

#-----------------------------------------------
# Co-localization analysis

showSpatialCellTypes <- function(obj.sp.lst, pseud.coord.lst, cell.types, sn.name, sp.size = c(3, 4.0), depend = NULL, down.sample = NULL, percent = 0.5, bool = T) {
    obj.sp <- obj.sp.lst[[sn.name]]
    pseud.coord <- pseud.coord.lst[[sn.name]]
    if (!bool) {
        pseud.coord <- pseud.coord %>%
            colData() %>%
            rename("HE.Labels" = "HE", CT = "cellType", SpotName = "spot") %>%
            as.data.frame()
        colnames(pseud.coord)[1:2] <- c("x.spot", "y.spot")
    }
    pseud.edges <- subset(pseud.coord, HE == "Medulla_edge") %>% .[!duplicated(.$spot), ]
    pseud.edges$x <- pseud.edges$x.spot
    pseud.edges$y <- pseud.edges$y.spot
    pseud.edges$cellType <- "Edge"

    img <- obj.sp@images[[sn.name]]@image
    pseud.coord.sub <- subset(pseud.coord, cellType %in% cell.types)
    pseud.coord.sub$x <- pseud.coord.sub$x.spot
    pseud.coord.sub$y <- pseud.coord.sub$y.spot
    if (!is.null(depend)) {
        dep.spots <- subset(pseud.coord.sub, cellType == depend)$spot %>% unique()
        pseud.coord.sub <- subset(pseud.coord.sub, spot %in% dep.spots)
    }
    df.dups <- pseud.coord.sub[c("x.spot", "y.spot")]
    pseud.parent <- pseud.coord.sub[!duplicated(df.dups), ]
    pseud.parent$x <- pseud.parent$x.spot
    pseud.parent$y <- pseud.parent$y.spot
    pseud.parent$cellType <- "Parent"
    plot.df <- rbind.data.frame(pseud.coord.sub, pseud.edges)
    plot.df <- rbind.data.frame(plot.df, pseud.parent)

    plot.df.sub <- subset(plot.df, !cellType %in% c("Parent", "Edge"))
    # plot.df.sub <- plot.df
    spot.size <- c(rep(sp.size[1], length(cell.types)), sp.size[2]) %>% `names<-`(c(cell.types, "Parent"))
    spot.shape <- c(rep(21, length(cell.types)), 1) %>% `names<-`(c(cell.types, "Parent"))

    if (!is.null(down.sample)) {
        idx.rand <- which(plot.df.sub$cellType == down.sample)
        plot.df.sub.tmp <- subset(plot.df.sub, cellType != down.sample)
        int.index <- sample(idx.rand, floor(length(idx.rand) * percent))
        plot.df.sub <- rbind.data.frame(plot.df.sub.tmp, plot.df.sub[int.index, ])
    }

    gp1 <- ggplot(plot.df.sub, aes(x = y, y = x, fill = cellType), alpha = 0.3) +
        annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        geom_point(aes(shape = cellType, size = cellType)) +
        ggplot2::scale_y_reverse() +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = "black"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_size_manual(values = c(spot.size, Edge = 0.5)) +
        scale_shape_manual(values = c(spot.shape, Edge = 0)) +
        scale_fill_manual(values = c(ANNO_ENTIRE_COLOR_FIG1[cell.types], Parent = "grey", Edge = "white")) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        ggtitle(sp2SampleName(obj.sp@images %>% names())) +
        xlab("") +
        ylab("")
    return(list(gp = gp1, spot = pseud.coord.sub$spot %>% unique()))
}

#------------------------------------------
# Compare spatial expression vs. Single-cell expressin

celltrekProj <- function(obj.celltrek) {
    gps <- lapply(obj.celltrek, function(xx) {
        SpatialPlot(xx, group.by = "Anno.Level.Fig.1", image.alpha = 0, cols = ANNO_ENTIRE_COLOR_FIG1) + NoLegend() + ggtitle(sp2SampleName(xx@images %>% names()))
    })
    ggarrange(plotlist = gps, ncol = 3, nrow = 3)
}

#------------------------------------------------
#

calcCorsSpSC <- function(obj.sp.lst, expr.lst) {
    cor.lst <- lapply(obj.sp.lst, function(xx) {
        ref.mat <- GetAssayData(xx)
        qry.mat <- expr.lst[[xx@images %>% names()]]
        ovp.genes <- intersect(rownames(ref.mat), rownames(qry.mat))
        qry.mat <- qry.mat[ovp.genes, colnames(ref.mat)]
        ref.mat <- ref.mat[ovp.genes, ]
        cors <- lapply(colnames(ref.mat), function(cell) {
            print(cell)
            cor(qry.mat[, cell], ref.mat[, cell])
        }) %>%
            unlist() %>%
            as.vector()
        cors
    })
    cor.lst
}

corDistPlot <- function(cor.lst) {
    cor.df <- lapply(names(cor.lst), function(xx) {
        cbind.data.frame(Cor = cor.lst[[xx]], SN = sp2SampleName(xx))
    }) %>% do.call(rbind, .)
    cor.df$SN <- factor(cor.df$SN, levels = names(SC_MAP_SP))
    ggplot(cor.df, aes(x = SN, y = Cor, color = SN)) +
        geom_violin() +
        geom_boxplot(width = 0.1, aes(color = SN), outlier.shape = NA) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        scale_color_manual(values = SAMPLE_COLORS %>% `names<-`(sp2SampleName(names(SAMPLE_COLORS)))) +
        ylab("Corr.\n(Expression of spot vs. segmentian spot)")
}

#----------------------------------------------
# Corr with CellTrek

corrSpVsSCMappedCtypesCellTrek <- function(obj.celltrek, obj.sc) {
    qry.vec <- rep(0, length(ANNO_ENTIRE_IDNET_FIG1)) %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
    tmp.bak <- lapply(obj.celltrek, function(xx) xx@meta.data[, c("Anno.Level.Fig.1")]) %>%
        unlist() %>%
        table()

    qry.vec[names(tmp.bak)] <- tmp.bak
    qry.vec <- qry.vec / sum(qry.vec)

    count.sc.ctypes <- Idents(obj.sc) %>% table()
    count.sc.ctypes <- count.sc.ctypes / sum(count.sc.ctypes)
    plot.df <- cbind.data.frame(SP = qry.vec %>% as.vector(), SC = count.sc.ctypes[names(qry.vec)] %>% as.vector(), CellType = names(qry.vec))
    plot.df$CellType <- factor(plot.df$CellType, levels = ANNO_ENTIRE_IDNET_FIG1)
    ggplot(plot.df, aes(x = SC, y = SP)) +
        geom_point(size = 4, aes(color = CellType)) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        geom_smooth(method = "lm", se = TRUE) +
        stat_cor(method = "pearson") +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        xlab("Single-cell data") +
        ylab("Projected using CellTrek") +
        guides(color = guide_legend(nrow = 6, byrow = TRUE))
}

corrSpVsSCMappedCtypesSeurat <- function(obj.seu, obj.sc) {
    qry.vec <- rep(0, length(ANNO_ENTIRE_IDNET_FIG1)) %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
    tmp.bak <- lapply(obj.seu, function(xx) xx@meta.data[, c("predicted.id")]) %>%
        unlist() %>%
        table()

    qry.vec[names(tmp.bak)] <- tmp.bak
    qry.vec <- qry.vec / sum(qry.vec)

    count.sc.ctypes <- Idents(obj.sc) %>% table()
    count.sc.ctypes <- count.sc.ctypes / sum(count.sc.ctypes)
    plot.df <- cbind.data.frame(SP = qry.vec %>% as.vector(), SC = count.sc.ctypes[names(qry.vec)] %>% as.vector(), CellType = names(qry.vec))
    plot.df$CellType <- factor(plot.df$CellType, levels = ANNO_ENTIRE_IDNET_FIG1)
    ggplot(plot.df, aes(x = SC, y = SP)) +
        geom_point(size = 4, aes(color = CellType)) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        geom_smooth(method = "lm", se = TRUE) +
        stat_cor(method = "pearson") +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        xlab("Single-cell data") +
        ylab("Projected using CellTrek")
}

#----------------------------------------------
# Boxplot of socre

validScoreOfCellTypes <- function(obj.sub, ovp.spots, cell.types, genes = NULL) {
    if (!is.null(genes)) {
        bak <- FetchData(obj.sub, vars = genes)
        cc <- cbind.data.frame(bak, HE.Labels = obj.sub$HE.Labels)
        cell.types <- genes
    } else {
        cc <- obj.sub@meta.data[ovp.spots, c(cell.types, "HE.Labels")]
    }
    plot.df <- tidyr::gather(cc, "CellType", "Score", -HE.Labels)
    plot.df$HE.Labels[plot.df$HE.Labels %in% c("Medulla_hi", "Medulla_edge", "Medulla_centric")] <- "Medulla"
    plot.df$HE.Labels[plot.df$HE.Labels != "Medulla"] <- "Cortex"
    plot.df$CellType <- factor(plot.df$CellType, levels = cell.types)
    dodge <- position_dodge(width = 0.1)
    ggplot(plot.df, aes(x = HE.Labels, y = Score, fill = HE.Labels, color = HE.Labels)) +
        geom_violin() +
        geom_boxplot(outlier.shape = NA, width = 0.1, color = "black") +
        theme_bw(base_size = 12) +
        facet_grid(. ~ CellType) +
        stat_compare_means(aes(group = HE.Labels), label = "p.signif") +
        scale_fill_manual(values = CM_COLORS[c("Cortex", "Medulla")]) +
        scale_color_manual(values = CM_COLORS[c("Cortex", "Medulla")]) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#------------------------------------------------------
# CARD pie-plot

corrCardPiePlot <- function(obj.card, hc.order = T) {
    devtools::load_all("../../../softwares/CARD-master")
    ref.mat <- lapply(obj.card, function(xx) {
        xx@Proportion_CARD
    }) %>% do.call(rbind, .)
    proportion <- ref.mat[, intersect(names(ANNO_ENTIRE_COLOR_FIG1), colnames(ref.mat))]
    proportion <- proportion[, colnames(proportion)]
    cor_CARD <- cor(as.matrix(proportion))
    colors <- c("#91a28c", "white", "#8f2c37")

    p <- suppressMessages(ggcorrplot(cor_CARD,
        hc.order = hc.order,
        outline.color = "white",
        tl.srt = 60,
        tl.cex = 18,
        lab_size = 7,
        colors = colors
    ) +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
            axis.text = element_text(size = 10),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(0.45, "cm")
        ) +
        coord_fixed() +
        ggtitle("Correlation") + theme(plot.title = element_text(size = 14, face = "bold")))
}

cardPiePlot <- function(obj.card, out.figs.dir) {
    devtools::load_all("../../../softwares/CARD-master")
    parallel::mclapply(names(obj.card), function(sn.name) {
        print(sn.name)
        xx <- obj.card[[sn.name]]
        pdf(file.path(out.figs.dir, sprintf("%s.CARD.Pieplot.pdf", sn.name)), width = 8, height = 8)
        fp <- CARD.visualize.pie(xx@Proportion_CARD, xx@spatial_location, colors = ANNO_ENTIRE_COLOR_FIG1, rc = 3.0)
        print(fp)
        dev.off()
    }, mc.cores = 8)
}

cardMapSc <- function(sc.map.lst, out.figs.dir) {
    require(SingleCellExperiment)
    lapply(names(sc.map.lst), function(sn.name) {
        scMapping <- sc.map.lst[[sn.name]]
        MapCellCords <- as.data.frame(colData(scMapping))
        count_SC <- assays(scMapping)$counts

        df <- MapCellCords
        p9 <- ggplot(df, aes(x = y, y = x, colour = CT)) +
            geom_point(size = 0.3, shape = 21) +
            scale_colour_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
            ggplot2::scale_y_reverse() +
            theme( # plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                panel.background = element_rect(colour = "black", fill = "black"),
                plot.background = element_rect(colour = "black", fill = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                legend.title = element_text(size = 13, face = "bold"),
                legend.text = element_text(size = 12),
                legend.key = element_rect(colour = "transparent", fill = "white"),
                legend.key.size = unit(0.45, "cm"),
                strip.text = element_text(size = 15, face = "bold")
            )
        ggsave(file = file.path(out.figs.dir, sprintf("%s.sc.project.sp.CRAD.UMAP.pdf", sp2SampleName(sn.name))), width = 8, height = 6)
    })
}

corrSpVsSCMappedCtypesCARD <- function(sc.map.lst, obj.sc) {
    qry.vec <- rep(0, length(ANNO_ENTIRE_IDNET_FIG1)) %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
    require(SingleCellExperiment)
    tmp.bak <- lapply(names(sc.map.lst), function(sn.name) {
        scMapping <- sc.map.lst[[sn.name]]
        MapCellCords <- as.data.frame(colData(scMapping))
        MapCellCords$CT
    }) %>%
        unlist() %>%
        table()
    names(tmp.bak)[which(names(tmp.bak) == "abT_entry")] <- "abT(entry)"

    qry.vec[names(tmp.bak)] <- tmp.bak
    qry.vec <- qry.vec / sum(qry.vec)

    count.sc.ctypes <- Idents(obj.sc) %>% table()
    count.sc.ctypes <- count.sc.ctypes / sum(count.sc.ctypes)
    plot.df <- cbind.data.frame(SP = qry.vec %>% as.vector(), SC = count.sc.ctypes[names(qry.vec)] %>% as.vector(), CellType = names(qry.vec))
    plot.df$CellType <- factor(plot.df$CellType, levels = ANNO_ENTIRE_IDNET_FIG1)
    ggplot(plot.df, aes(x = SC, y = SP)) +
        geom_point(size = 4, aes(color = CellType)) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        geom_smooth(method = "lm", se = TRUE) +
        stat_cor(method = "pearson") +
        scale_color_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
        xlab("Single-cell data") +
        ylab("Projected using CARD") +
        guides(color = guide_legend(nrow = 6, byrow = TRUE))
}

#------------------------------------------------------
# Compare proportions between Proj2Spatial and CARD.

compareProps <- function(obj.sc, sc.map.lst, pseud.coord.lst) {
    ref.prop <- table(obj.sc$Anno.Level.Fig.1, obj.sc$orig.ident)[, paste0("Thy", 5:12)] %>% as.data.frame.matrix()
    ref.prop <- sweep(ref.prop, 2, colSums(ref.prop), "/")
    prop.vec <- rep(0, length(ANNO_ENTIRE_COLOR_FIG1)) %>% `names<-`(names(ANNO_ENTIRE_COLOR_FIG1))
    card.df <- lapply(names(sc.map.lst), function(sn.name) {
        xx <- sc.map.lst[[sn.name]]
        xx <- as.data.frame(colData(xx))
        xx.bak <- (xx$CT %>% table()) / sum(xx$CT %>% table())
        idx <- which(names(xx.bak) == "abT_entry")
        names(xx.bak)[idx] <- "abT(entry)"
        prop.bak <- prop.vec
        prop.bak[names(xx.bak)] <- xx.bak
        yy <- cbind.data.frame(Prop = prop.bak, Class = "CARD", Sample = sp2SampleName(sn.name), CellType = names(prop.bak))
    }) %>% do.call(rbind, .)

    ours.df <- lapply(names(pseud.coord.lst), function(sn.name) {
        xx <- pseud.coord.lst[[sn.name]]
        xx.bak <- (xx$cellType %>% table()) / sum(xx$cellType %>% table())
        prop.bak <- prop.vec
        prop.bak[names(xx.bak)] <- xx.bak
        yy <- cbind.data.frame(Prop = prop.bak, Class = "Proj2Spatial", Sample = sp2SampleName(sn.name), CellType = names(prop.bak))
    }) %>% do.call(rbind, .)
    ref.df <- tidyr::gather(cbind.data.frame(ref.prop, CellType = rownames(ref.prop)), "Sample", "Prop", -CellType) %>% cbind.data.frame(., Class = "Ref")
    ref.df <- ref.df[, colnames(ours.df)]

    plot.df <- rbind.data.frame(card.df, ours.df, ref.df)
    plot.df.new <- lapply(plot.df$CellType %>% unique(), function(cell) {
        df.sub <- subset(plot.df, CellType == cell)
        df.sub$Prop <- scale(df.sub$Prop)
        df.sub
    }) %>% do.call(rbind, .)
    ggplot(plot.df.new, aes(x = CellType, y = Prop, color = Class)) +
        geom_boxplot(outlier.shape = NA) +
        theme_classic(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
        )
}

#----------------------------------------------------------------

cardMapScSplit <- function(sc.map.lst, obj.sp.lst, out.figs.dir, ctypes = ANNO_ENTIRE_IDNET_FIG1, label = "Ours") {
    require(SingleCellExperiment)
    lapply(names(sc.map.lst), function(sn.name) {
        scMapping <- sc.map.lst[[sn.name]]
        MapCellCords <- as.data.frame(colData(scMapping))
        img <- obj.sp.lst[[sn.name]]@images[[sn.name]]@image

        df <- MapCellCords
        df$CT[df$CT == "abT_entry"] <- "abT(entry)"
        lapply(ctypes, function(ct) {
            p9 <- ggplot(subset(df, CT == ct), aes(x = y, y = x, colour = CT)) +
                annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                geom_point(size = 0.6, shape = 21) +
                scale_colour_manual(values = ANNO_ENTIRE_COLOR_FIG1) +
                ggplot2::scale_y_reverse() +
                theme(
                    panel.background = element_rect(colour = "black", fill = "black"),
                    plot.background = element_rect(colour = "black", fill = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    legend.title = element_text(size = 13, face = "bold"),
                    legend.text = element_text(size = 12),
                    legend.key = element_rect(colour = "transparent", fill = "white"),
                    legend.key.size = unit(0.45, "cm"),
                    strip.text = element_text(size = 15, face = "bold")
                )
            ggsave(file = file.path(out.figs.dir, sprintf("%s.%s.%s.sc.project.sp.CRAD.UMAP.pdf", label, sp2SampleName(sn.name), ct)), width = 8, height = 6)
        })
    })
}


sankyPlot <- function(sc.map.lst, sc.map.lst.card, obj.sp.lst) {
    ours <- statCMCTypes(sc.map.lst, obj.sp.lst)
    ours <- ours %>%
        sweep(., 2, colSums(.), "/") %>%
        cbind.data.frame(., Region = rownames(.))
    card <- statCMCTypes(sc.map.lst.card, obj.sp.lst, group = "card")
    card <- card %>%
        sweep(., 2, colSums(.), "/") %>%
        cbind.data.frame(., Region = rownames(.))

    ours <- ours %>%
        tidyr::gather(., "CellType", "Proportion", -Region) %>%
        cbind.data.frame(., Source = "Ours") %>%
        .[rowSums(is.na(.)) == 0, ]
    card <- card %>%
        tidyr::gather(., "CellType", "Proportion", -Region) %>%
        cbind.data.frame(., Source = "CARD") %>%
        .[rowSums(is.na(.)) == 0, ]
    plot.df <- rbind.data.frame(ours, card)
    plot.df <- subset(plot.df, CellType %in% intersect(ours$CellType %>% unique(), card$CellType %>% unique()))

    library(ggalluvial)
    plot.df$CellType <- factor(plot.df$CellType, levels = names(ANNO_ENTIRE_COLOR_FIG1))

    plot.df$CellType <- factor(plot.df$CellType, levels = names(ANNO_ENTIRE_COLOR_FIG1) %>% rev())
    plot.df$Region <- factor(plot.df$Region, levels = c("Cortex", "Medulla_edge", "Medulla") %>% rev())
    plot.sub <- subset(plot.df, Source == "Ours")
    cols <- lapply(plot.sub$CellType, function(xx) ANNO_ENTIRE_COLOR_FIG1[xx]) %>%
        unlist() %>%
        as.vector()
    plot.sub$Col <- cols


    pdf(file.path(out.figs.dir, "ours.sanky.plot.pdf"), width = 6.5, height = 12)
    alluvial(
        plot.sub[, c("CellType", "Region")],
        freq = plot.sub[, "Proportion"],
        col = plot.sub$Col
    )
    dev.off()

    plot.sub <- subset(plot.df, Source == "CARD")
    cols <- lapply(plot.sub$CellType, function(xx) ANNO_ENTIRE_COLOR_FIG1[xx]) %>%
        unlist() %>%
        as.vector()
    plot.sub$Col <- cols

    pdf(file.path(out.figs.dir, "CARD.sanky.plot.pdf"), width = 6.5, height = 12)
    alluvial(
        plot.sub[, c("CellType", "Region")],
        freq = plot.sub[, "Proportion"],
        col = plot.sub$Col
    )
    dev.off()
}

cardMapScCompare <- function(sc.map.lst, sc.map.lst.card, obj.sp.lst, out.figs.dir, ctypes = ANNO_ENTIRE_IDNET_FIG1, label = "Compare", xx.filter = "CARD") {
    require(SingleCellExperiment)
    lapply(names(sc.map.lst), function(sn.name) {
        scMapping <- sc.map.lst[[sn.name]]
        MapCellCords <- as.data.frame(colData(scMapping))
        img <- obj.sp.lst[[sn.name]]@images[[sn.name]]@image
        MapCellCords.ref <- sc.map.lst.card[[sn.name]] %>%
            colData(.) %>%
            as.data.frame()

        df <- MapCellCords %>% reformatCARD(., obj.sp.lst[[sn.name]])
        df$CT[df$CT == "abT_entry"] <- "abT(entry)"

        MapCellCords.ref <- MapCellCords.ref %>% reformatCARD(., obj.sp.lst[[sn.name]])
        MapCellCords.ref$CT[MapCellCords.ref$CT == "abT_entry"] <- "abT(entry)"

        df.edges <- subset(df, HE.Labels == "Medulla_edge")
        df.edges$Source <- "Ours"
        df.edges <- df.edges[!duplicated(df.edges$centerSPOT), ]
        df.edges$x <- df.edges$centerx
        df.edges$y <- df.edges$centery
        df.edges$Source <- "Edge"

        lapply(ctypes, function(ct) {
            print(ct)
            plot.df.1 <- cbind.data.frame(subset(df, CT %in% ct), Source = "Ours")
            plot.df.2 <- cbind.data.frame(subset(MapCellCords.ref, CT %in% ct), Source = "CARD")
            plot.df <- rbind.data.frame(plot.df.1, plot.df.2)

            plot.df.new <- rbind.data.frame(cbind.data.frame(plot.df, Shape = "C"), cbind.data.frame(df.edges, Shape = "R"))

            plot.df.new$CT.NEW <- plot.df.new$CT
            plot.df.new$CT.NEW[plot.df.new$Source == "Edge"] <- "Edge"
            p9 <- ggplot(subset(plot.df.new, Source != xx.filter), aes(x = y, y = x, fill = CT.NEW, color = Source)) +
                annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                geom_point(aes(shape = Shape, size = Shape)) +
                scale_color_manual(values = c(Ours = "black", CARD = "black", Edge = "grey")) +
                scale_fill_manual(values = c(ANNO_ENTIRE_COLOR_FIG1[ct], Edge = "white")) +
                ggplot2::scale_y_reverse() +
                theme(
                    panel.background = element_rect(colour = "black", fill = "black"),
                    plot.background = element_rect(colour = "black", fill = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    legend.title = element_text(size = 13, face = "bold"),
                    legend.text = element_text(size = 12),
                    legend.key = element_rect(colour = "transparent", fill = "white"),
                    legend.key.size = unit(0.45, "cm"),
                    strip.text = element_text(size = 15, face = "bold")
                ) +
                scale_shape_manual(values = c(C = 21, R = 24)) +
                scale_size_manual(values = c(2, 1))
            ggsave(file = file.path(out.figs.dir, sprintf("%s.%s.%s.%s.UMAP.pdf", label, sp2SampleName(sn.name), paste0(ct, collapse = "_"), xx.filter)), width = 7, height = 6)
        })
    })
}

markerExpresion <- function(sc.map.lst, obj.sp.lst, out.figs.dir, ctypes = ANNO_ENTIRE_IDNET_FIG1) {
    lapply(names(sc.map.lst), function(sn) {
        sc.map <- sc.map.lst[[sn]] %>%
            colData() %>%
            as.data.frame()
        meta.info <- obj.sp.lst[[sn]]@meta.data
        sc.map <- reformatCARD(sc.map, obj.sp.lst[[sn]])
        plot.df <- lapply(sc.map$CT %>% unique(), function(ct) {
            print(ct)
            if (ct == "abT_entry") ct <- "abT(entry)"
            sc.map.bak <- subset(sc.map, CT == ct)
            df.tmp <- cbind.data.frame(meta.info[, ct, drop = F], CT = ct, Group = "N")
            df.tmp[sc.map.bak$SpotName %>% unique(), "Group"] <- "Y"
            df.tmp <- df.tmp %>% `colnames<-`(c("Score", "CT", "Group"))
            df.tmp
        }) %>% do.call(rbind, .)
        plot.df <- cbind.data.frame(plot.df, SN = sn)
        plot.df
    }) %>% do.call(rbind, .) -> plot.res
    plot.res$CT <- factor(plot.res$CT, levels = intersect(ctypes, plot.res$CT %>% unique()))
    gp <- ggplot(plot.res, aes(x = CT, y = Score, fill = Group, color = Group)) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(hjust = 1, angle = 45)
        ) +
        stat_compare_means(label = "p.signif") +
        xlab("") +
        ylab("Signature score") +
        scale_fill_manual(values = c("grey", "grey")) +
        scale_color_manual(values = c(Y = "red", N = "blue"))
    gp
}
