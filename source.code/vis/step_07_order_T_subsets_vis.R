tcrSubgroupsTrend <- function(obj.merged.sc, permu = 100) {
    obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17:33]))
    obj.1 <- tcrSubgroup1(obj, trd.cut = 1) %>% subset(., trb != "")
    tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
    obj.2 <- subset(obj.1, trb != "" & (t.cell_group %in% tcr.subtypes))
    obj.2 <- subset(obj.2, Anno.Level.11 %in% c(ANNO_TREND_IDENT[1:15], "CD8T", "CD4T"))
    obj.2$Anno.Level.11 <- factor(obj.2$Anno.Level.11, levels = c(ANNO_TREND_IDENT[1:15], "CD8T", "CD4T"))

    order.cells <- c(ANNO_TREND_IDENT[1:15], "CD8T", "CD4T")
    Idents(obj.2) <- obj.2$t.cell_group
    percent.df <- lapply(1:permu, function(idx) {
        print(idx)
        obj.2.tmp <- downSamplSeurat(obj.2, cnt = 2000, seed = idx)
        xx <- table(obj.2.tmp$Anno.Level.11, obj.2.tmp$t.cell_group)
        sweep(xx, 1, rowSums(xx), "/")[order.cells, tcr.subtypes]
    }) %>% do.call(cbind, .)

    percent.df.new <- lapply(tcr.subtypes, function(xx) {
        percent.df[, which(colnames(percent.df) == xx)] %>% rowMeans()
    }) %>%
        do.call(cbind, .) %>%
        `colnames<-`(tcr.subtypes)

    res <- lapply(1:nrow(percent.df.new), function(idx) {
        if (idx == nrow(percent.df.new)) {
            tmp <- rbind.data.frame(percent.df.new[idx, ]) %>%
                `colnames<-`(tcr.subtypes) %>%
                `rownames<-`(c(rownames(percent.df.new)[idx]))
            tmp.bak <- cbind.data.frame(tmp, CellType = rownames(tmp), Idx = idx)
            return(tmp.bak)
        } else {
            tmp <- rbind.data.frame(percent.df.new[idx, ], colMeans(percent.df.new[idx:(idx + 1), ])) %>%
                `colnames<-`(tcr.subtypes) %>%
                `rownames<-`(c(rownames(percent.df.new)[idx], "xx"))
        }
        tmp.bak <- cbind.data.frame(tmp, CellType = rownames(tmp), Idx = c(idx, (idx * 2 + 1) / 2))
    }) %>% do.call(rbind, .)

    plot.df <- tidyr::gather(res, "Class", "Proportion", -CellType, -Idx)
    plot.df.sub <- subset(plot.df, CellType != "xx")
    plot.df.sub$Class <- factor(plot.df.sub$Class, levels = tcr.subtypes)

    plot.df.ht <- sweep(percent.df.new, 2, apply(percent.df.new, 2, max), "/")
    plot.df.ht <- plot.df.ht %>%
        as.data.frame() %>%
        cbind.data.frame(., CellType = rownames(.))
    plot.df.ht <- tidyr::gather(plot.df.ht, "Class", "Proportion", -CellType)

    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space = "Lab")
    plot.df.ht$CellType <- factor(plot.df.ht$CellType, levels = order.cells)
    plot.df.ht$Class <- factor(plot.df.ht$Class, levels = tcr.subtypes %>% rev())
    gp1 <- ggplot(plot.df.ht, aes(CellType, Class, fill = Proportion)) +
        geom_tile() +
        theme_bw() +
        theme_classic() +
        scale_fill_gradientn(colours = myPalette(100)) +
        easy_remove_axes(
            which = c("both"),
            what = c("ticks", "title", "line"),
            teach = FALSE
        ) +
        theme(legend.position = "top", axis.text.x = element_blank())


    gp2 <- ggplot(plot.df.sub, aes(x = Idx, y = Proportion)) +
        geom_point(aes(color = CellType), size = 3) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 9), color = "grey40", size = 1.2, linetype = "solid") +
        facet_wrap(~Class, nrow = 4, scale = "free_y") +
        scale_x_continuous(limits = c(1, nrow(percent.df.new)), breaks = seq(1, nrow(percent.df.new), 1), labels = rownames(percent.df.new)) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_line(color = "#90E669", linetype = "dashed", size = 0.25),
            panel.grid.major = element_blank(),
            legend.position = "none"
        ) +
        xlab("") +
        scale_color_manual(values = ANNO_TREND_IDENT_COLOR %>% `names<-`(order.cells))
    ggarrange(gp1, gp2, nrow = 2, heights = c(1, 2.5))
}

#-----------------------------------------------
# GEGs heatmap

plotHeamtmapForDEGs <- function(obj.sub) {
    obj.sub$Cluster <- Idents(obj.sub)
    obj.sub$AgeStatus <- factor(obj.sub$AgeStatus, levels = names(SampleClassifyColors))
    plot_heatmap(
        obj.sub,
        de.markers,
        n = 25,
        sort_var = c("Cluster", "Anno.Level.Fig.1"),
        anno_var = c("Cluster", "Anno.Level.Fig.1", "S.Score", "G2M.Score", "AgeStatusNew", "tra_sum", "trb_sum", "trd_sum", "trg_sum", "orig.ident") %>% rev(.),
        anno_colors = list(
            T.CELL.GROUPS.COLORS,
            ANNO_ENTIRE_COLOR_FIG1[setdiff(levels(obj.sub$Anno.Level.Fig.1)[17:33], c("T_proliferating", "CD8aa"))],
            "Reds",
            "Greens",
            AGE_STATUS_EXTENT_COLOR[1:4],
            "BuPu",
            "Purples",
            "Greys",
            "Oranges",
            SAMPLE.COLORS_SC[1:16]
        ) %>% rev(.)
    )
}

#-----------------------------------------------
# Trend of score

fourGroupsTrend <- function(obj.st.lst, plot.tar) {
    cols <- T.CELL.GROUPS.COLORS %>% `names<-`(plot.tar)
    trendLinesWithDistPlot1(obj.st.lst, plot.tar = plot.tar) +
        scale_color_manual(values = cols, labels = names(T.CELL.GROUPS.COLORS)) +
        ylab("Signature score")
}

#-----------------------------------------------
# Track clonal

trackClonalLine <- function(pp) {
    kk <- apply(pp, 2, function(xx) cumsum(rev(xx)))
    kk.frac <- sweep(kk, 2, apply(kk, 2, max), "/") %>%
        as.data.frame() %>%
        add_column(CellType = rownames(.), Idex = 1:nrow(kk))
    plot.df <- reshape2::melt(kk.frac, id.vars = c("CellType", "Idex"))
    plot.df$CellType <- factor(plot.df$CellType, levels = rownames(kk.frac))
    ggplot(plot.df, aes(x = Idex, y = value, color = variable)) +
        geom_smooth(method = "gam", se = FALSE) +
        theme_bw(base_size = 14) +
        xlab("CellType") +
        ylab("Cumulative frequency distribution") +
        scale_color_manual(values = T.CELL.GROUPS.COLORS) +
        ylim(0, 1) +
        xlab("") +
        scale_x_continuous(limits = c(1, nrow(kk)), breaks = seq(1, nrow(kk), 1), labels = rownames(kk)) +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            panel.grid.major = element_line(color = "grey", size = 0.5, linetype = "dotted"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        guides(color = guide_legend(title = "Subgroup")) +
        ggtitle("Tracking clones with amplification\nin four subgroups")
}

#------------------------------------------------
# DP-blast plot

cellcyleHeatmap <- function(obj.blast, cell.cycle.genes, out.figs.dir) {
    Idents(obj.blast) <- obj.blast$Anno.Level.11
    levels(obj.blast) <- levels(obj.blast) %>% .[order(.)]

    genes <- cell.cycle.genes %>%
        unlist(use.names = T) %>%
        as.data.frame() %>%
        .[!duplicated(.[, 1]), , drop = F] %>%
        .[!is.na(.[, 1]), , drop = F]
    genes <- cbind.data.frame(genes, Status = rownames(genes))
    genes[, 2] <- gsub("^G1/S.*", "G1/S", genes[, 2])
    genes[, 2] <- gsub("^S.*", "S", genes[, 2])
    genes[, 2] <- gsub("^G2/M.*", "G2/M", genes[, 2])
    genes[, 2] <- gsub("^M\\d+", "M", genes[, 2])
    genes[, 2] <- gsub("^M/G1.*", "M/G1", genes[, 2])
    DEG.color <- c(`G1/S` = "#277221", `S` = "#4D71B2", `G2/M` = "#4C0E81", `M` = "#FF0800", `M/G1` = "#CD853F")
    pdf(file.path(out.figs.dir, "cellcycle.heatmap.dp.blast.pdf"), width = 8, height = 10)
    plot_heatmapNew(
        obj.blast,
        genes,
        sort_var = c("Anno.Level.11", "dpt_pseudotime"),
        anno_var = c("Anno.Level.11"),
        anno_colors = list(
            ANNO_TREND_IDENT_COLOR[levels(obj.blast)]
        ),
        DEG.color = DEG.color,
        hm_limit = c(0, 0.5, 1.0),
        hm_colors = c("grey", "white", "red")
    )
    dev.off()

    obj.blast <- AddModuleScore(obj.blast, features = list(xx = genes[, 1], yy = sample(rownames(obj.blast), 500)))
    plot.data <- obj.blast@meta.data[, c("dpt_pseudotime", "Cluster1", "Cluster2", "Anno.Level.11")]
    plot.data <- plot.data[with(plot.data, order(Anno.Level.11, dpt_pseudotime)), ]
    plot.df <- tidyr::gather(plot.data, "Class", "Score", -dpt_pseudotime, -Anno.Level.11)
    ggplot(data = plot.df, aes(x = rep(1:nrow(obj.blast@meta.data), 2), y = Score, color = Class)) +
        geom_smooth(se = TRUE, formula = y ~ splines::ns(x, 10), method = lm) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
        ) +
        xlab("") +
        ylab("Signature score") +
        scale_color_manual(values = c("red", "grey"))
    ggsave(file.path(out.figs.dir, "cellcycle.score.trend.line.pdf"), width = 8, height = 3)
}

#--------------------------------------------------
# Signature score of four T subgroups.

fourSubgroupScores <- function(obj.merged.sc, out.figs.dir) {
    require(future)
    options(future.globals.maxSize = 500000 * 1024^2)
    plan("multiprocess", workers = 20)

    tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")
    obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17:33]))
    obj <- tcrSubgroup1(obj, trd.cut = 0) %>% subset(., trb != "")
    obj <- subset(obj, t.cell_group %in% tcr.subtypes)

    Idents(obj) <- obj$t.cell_group
    levels(obj) <- tcr.subtypes

    de.markers <- FindAllMarkers(obj, only.pos = TRUE)
    de.markers <- de.markers[de.markers$gene %in% coding.genes, ]

    top20 <- de.markers %>%
        group_by(cluster) %>%
        top_n(n = 100, wt = avg_log2FC)
    xx.genes <- split(top20$gene, top20$cluster)

    lapply(names(obj.st.lst), function(sn.name) {
        obj.st <- obj.st.lst[[sn.name]]
        obj.st <- AddModuleScore(object = obj.st, features = xx.genes, name = "T.Subgroups")
        colnames(obj.st@meta.data) <- gsub("T.Subgroups1", "TRA-TRB+TRD+", colnames(obj.st@meta.data))
        colnames(obj.st@meta.data) <- gsub("T.Subgroups2", "TRA+TRB+TRD+", colnames(obj.st@meta.data))
        colnames(obj.st@meta.data) <- gsub("T.Subgroups3", "TRA-TRB+TRD-", colnames(obj.st@meta.data))
        colnames(obj.st@meta.data) <- gsub("T.Subgroups4", "TRA+TRB+TRD-", colnames(obj.st@meta.data))
        SpatialFeaturePlot(obj.st, features = tcr.subtypes, image.alpha = 0)
        ggsave(file.path(out.figs.dir, sprintf("%s.signature.scores.pdf", sp2SampleName(sn.name))), width = 8, height = 8)
    })
}

#-------------------------------------------------
# Venn diagram of differential expressed genes.

overlapDegFourSubgroups <- function(de.markers) {
    de.markers.sub <- subset(de.markers, p_val_adj < 0.01)
    de.markers.lst <- split(de.markers.sub$gene, de.markers.sub$cluster) %>% .[names(T.CELL.GROUPS.COLORS)]
    temp <- venn.diagram(
        de.markers.lst,
        filename = NULL,
        fill = T.CELL.GROUPS.COLORS,
        col = "white",
        category.names = c(paste0(names(de.markers.lst), "\n(", lapply(de.markers.lst, length) %>% as.vector(), ")")),
        cat.pos = c(-8, 5, 350, 350)
    )
    return(temp)
}

#--------------------------------------------------
# tcrExamples

pairwiseArrow <- function(plot.df.sub) {
    xx.start <- subset(plot.df.sub, Cluster == "C1")
    xx.end <- subset(plot.df.sub, Cluster == "C2")
    df.pos <- lapply(xx.start$Group %>% as.vector() %>% unique(), function(xx.seq) {
        tmp.pos.end <- xx.end[which(xx.end$Group == xx.seq), ]
        tmp.pos.start <- xx.start[which(xx.start$Group == xx.seq), ]
        lapply(1:dim(tmp.pos.start)[1], function(pp.1.idx) {
            pp.1 <- tmp.pos.start[pp.1.idx, ]
            cc <- lapply(1:dim(tmp.pos.end)[1], function(pp.2.idx) {
                pp.2 <- tmp.pos.end[pp.2.idx, ]
                c(pp.1[1], pp.1[2], pp.2[1], pp.2[2])
            }) %>%
                do.call(rbind, .) %>%
                data.frame() %>%
                `colnames<-`(c("x1", "y1", "x2", "y2"))
            cc <- cbind.data.frame(cc, TCR = tmp.pos.end$Group)
        }) %>%
            do.call(rbind, .) %>%
            data.frame()
    }) %>%
        do.call(rbind, .) %>%
        data.frame()
    df.pos <- apply(df.pos, 1, unlist) %>%
        t() %>%
        as.data.frame()
    df.pos[, 1] <- as.numeric(df.pos[, 1])
    df.pos[, 2] <- as.numeric(df.pos[, 2])
    df.pos[, 3] <- as.numeric(df.pos[, 3])
    df.pos[, 4] <- as.numeric(df.pos[, 4])
    df.pos
}

#--------------------------------------------------
# Mature status

matureStatus <- function(obj.merged.sc, out.figs.dir) {
    obj <- subset(obj.merged.sc, (tcr == 1) & (Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17:33]))
    obj <- tcrSubgroup1(obj) %>% subset(., trb != "")

    library(ggrepel)
    obj.dn.dp <- obj
    obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c("CD4", "CD8B")), CellType = obj.dn.dp$Anno.Level.11)
    cd4.avg <- aggregate(CD4 ~ CellType, data = obj.expr, FUN = "mean")
    cd8b.avg <- aggregate(CD8B ~ CellType, data = obj.expr, FUN = "mean")
    plot.expr <- cbind.data.frame(CellType = cd4.avg[, 1], CD4 = cd4.avg[, 2], CD8B = cd8b.avg[, 2])
    ggplot(plot.expr, aes(x = CD4, y = CD8B, label = CellType)) +
        geom_point(size = 4, aes(color = CellType)) +
        geom_text_repel(hjust = 0, nudge_x = 0.02) +
        theme_bw(base_size = 16) +
        theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
        ggtitle("CD4/CD8B Ratio") +
        scale_color_manual(values = ANNO_TREND_IDENT_COLOR) +
        geom_smooth(method = "lm", linetype = "dashed", alpha = 0.5)
    ggsave(file.path(out.figs.dir, "sp_CD4_CD8B_Ratio.scatter.pdf"), width = 8, height = 8)

    obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c("TRDV1", "TRDV2"), slot = "count"), CellType = obj.dn.dp$Anno.Level.11)
    obj.expr[, 1:2] <- obj.expr[, 1:2] > 0
    trdv1.avg <- aggregate(TRDV1 ~ CellType, data = obj.expr, FUN = "sum") %>% as.data.frame.matrix()
    trdv1.avg[, 2] <- trdv1.avg[, 2] / as.vector(table(obj.expr$CellType))
    trdv2.avg <- aggregate(TRDV2 ~ CellType, data = obj.expr, FUN = "sum") %>% as.data.frame.matrix()
    trdv2.avg[, 2] <- trdv2.avg[, 2] / as.vector(table(obj.expr$CellType))

    plot.expr <- cbind.data.frame(CellType = trdv1.avg[, 1], TRDV1 = trdv1.avg[, 2], TRDV2 = trdv2.avg[, 2])
    ggplot(plot.expr, aes(x = TRDV1, y = TRDV2, label = CellType)) +
        geom_point(size = 4, aes(color = CellType)) +
        geom_text_repel(hjust = 0, nudge_x = 0.002) +
        theme_bw(base_size = 16) +
        theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
        ggtitle("TRDV1/TRDV2 Ratio") +
        scale_color_manual(values = ANNO_TREND_IDENT_COLOR) +
        geom_smooth(method = "lm", linetype = "dashed", alpha = 0.5)
    ggsave(file.path(out.figs.dir, "dn2sp_TRDV1_TRDV2_Ratio.scatter.pdf"), width = 8, height = 8)

    genes <- rownames(obj.dn.dp)
    tra <- grep(("^TRA[VJC]"), genes, value = T)
    trb <- grep(("^TRB[VJC]"), genes, value = T)

    obj.dn.dp <- AddModuleScore(obj.dn.dp, features = list(xx = tra, yy = trb))
    obj.expr <- cbind.data.frame(FetchData(obj.dn.dp, vars = c("Cluster1", "Cluster2")), CellType = obj.dn.dp$Anno.Level.11)
    cd4.avg <- aggregate(Cluster1 ~ CellType, data = obj.expr, FUN = "mean")
    cd8b.avg <- aggregate(Cluster2 ~ CellType, data = obj.expr, FUN = "mean")
    plot.expr <- cbind.data.frame(CellType = cd4.avg[, 1], TCRa = cd4.avg[, 2], TCRb = cd8b.avg[, 2])
    ggplot(plot.expr, aes(x = TCRa, y = TCRb, label = CellType)) +
        geom_point(size = 4, aes(color = CellType)) +
        geom_text_repel(hjust = 0, nudge_x = 0.0) +
        theme_bw(base_size = 16) +
        theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) +
        ggtitle("TRA_chain/TRB_chain Ratio") +
        scale_color_manual(values = ANNO_TREND_IDENT_COLOR) +
        geom_smooth(method = "lm", linetype = "dashed", alpha = 0.5)
    ggsave(file.path(out.figs.dir, "tra_vs._trb.chains_Ratio.scatter.pdf"), width = 8, height = 8)
}

#------------------------------------

plotRePseudo <- function(xx.bak) {
    plot.df <- cbind.data.frame(dpt_pseudotime = xx.bak$dpt_pseudotime, CellType = xx.bak$Anno.Level.11)
    ggplot(plot.df, aes(x = CellType, y = dpt_pseudotime, color = CellType)) +
        geom_boxplot(outlier.shape = NA) +
        ylim(0.25, 0.26) +
        stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed") +
        scale_color_manual(values = ANNO_TREND_IDENT_COLOR[xx.bak$Anno.Level.11 %>%
            unique() %>%
            .[order(.)]]) +
        theme_bw(base_size = 12) +
        theme_classic()
}

travUsages <- function(obj.1, out.figs.dir) {
    genes <- rownames(obj.1)
    trav <- grep(("^TRA[V]"), genes, value = T)
    trav <- tra[str_order(trav, numeric = T)]
    traj <- grep(("^TRA[J]"), genes, value = T)
    traj <- traj[str_order(traj, numeric = T)] %>% rev()
    tra <- c(trav, traj)

    obj.x <- subset(obj.1, idents = c("TRA+TRB+TRD+"))
    obj.y <- subset(obj.1, idents = c("TRA+TRB+TRD-"))

    xx.1 <- FetchData(obj.x, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)
    xx.2 <- FetchData(obj.y, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)

    xx.df <- cbind.data.frame(xx.1, xx.2) %>% `colnames<-`(c("TRA+TRB+TRD+", "TRA+TRB+TRD-"))
    res.chisq <- chisq.test(xx.df)
    R.oe <- (res.chisq$observed) / (res.chisq$expected)
    R.oe <- R.oe[rowSums(is.na(R.oe)) == 0, ]

    plot.df <- R.oe %>%
        as.data.frame.matrix() %>%
        rownames_to_column(var = "CellType") %>%
        tidyr::gather(., "Status", "Roe", -CellType)
    plot.df$CellType <- factor(plot.df$CellType, levels = rownames(R.oe))
    plot.df$Color <- plot.df$Roe >= 1
    ggplot(plot.df) +
        geom_point(aes(x = CellType, y = Status, size = Roe, color = Color)) +
        scale_size(range = c(min(plot.df$Roe), max(plot.df$Roe))) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linetype = "dotted"),
            legend.position = "top"
        ) +
        scale_color_manual(values = TUMOR_NORMAL %>% as.vector() %>% rev(), label = c("Depletion", "Enrichment")) +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("")
    ggsave(file.path(out.figs.dir, "trav.usages.pdf"), width = 17, height = 3)
}

trbvUasges <- function(obj.1, out.figs.dir) {
    genes <- rownames(obj.1)
    trav <- grep(("^TRB[V]"), genes, value = T)
    trav <- trav[str_order(trav, numeric = T)]
    traj <- grep(("^TRB[J]"), genes, value = T)
    traj <- traj[str_order(traj, numeric = T)] %>% rev()
    tra <- c(trav, traj)

    obj.x <- subset(obj.1, idents = c("TRA-TRB+TRD+"))
    obj.y <- subset(obj.1, idents = c("TRA+TRB+TRD+"))
    obj.w <- subset(obj.1, idents = c("TRA-TRB+TRD-"))
    obj.z <- subset(obj.1, idents = c("TRA+TRB+TRD-"))

    xx.1 <- FetchData(obj.x, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)
    xx.2 <- FetchData(obj.y, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)
    xx.3 <- FetchData(obj.w, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)
    xx.4 <- FetchData(obj.z, vars = tra) %>%
        {
            . > 0
        } %>%
        colSums(.)

    xx.df <- cbind.data.frame(xx.1, xx.2, xx.3, xx.4) %>% `colnames<-`(tcr.subtypes)
    res.chisq <- chisq.test(xx.df)
    R.oe <- (res.chisq$observed) / (res.chisq$expected)
    R.oe <- R.oe[rowSums(is.na(R.oe)) == 0, ]

    plot.df <- R.oe %>%
        as.data.frame.matrix() %>%
        rownames_to_column(var = "CellType") %>%
        tidyr::gather(., "Status", "Roe", -CellType)
    plot.df$CellType <- factor(plot.df$CellType, levels = rownames(R.oe))
    plot.df$Color <- plot.df$Roe >= 1
    plot.df$Status <- factor(plot.df$Status, levels = tcr.subtypes %>% rev())

    ggplot(plot.df) +
        geom_point(aes(x = CellType, y = Status, size = Roe, color = Color)) +
        scale_size(range = c(min(plot.df$Roe), max(plot.df$Roe))) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linetype = "dotted"),
            legend.position = "top"
        ) +
        scale_color_manual(values = TUMOR_NORMAL %>% as.vector() %>% rev(), label = c("Depletion", "Enrichment")) +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        xlab("") +
        ylab("")

    ggsave(file.path(out.figs.dir, "trbv.usages.pdf"), width = 17, height = 4)
}
