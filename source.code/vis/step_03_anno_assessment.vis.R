#-------------------------------------------------------------
# Number of gene detected per cell type.

numberOfGeneDetecedPerCellType <- function(merged.obj) {
    ggplot(merged.obj@meta.data, aes(x = Anno.Level.Fig.1, y = nFeature_RNA, fill = "#94B3DE")) +
        geom_boxplot(outlier.colour = "red") +
        geom_boxplot(aes(color = "#94B3DE"), fill = NA, fatten = NULL, outlier.alpha = 0) +
        theme_classic(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "none"
        ) +
        xlab("") +
        ylab("Number of Genes detected") +
        ylim(0, 5000) +
        scale_fill_manual(values = "#94B3DE") +
        scale_color_manual(values = "#94B3DE")
}

#-------------------------------------------------------------
# Number of detected UMIs

numberOfDetectedUMIs <- function(merged.obj) {
    ggplot(merged.obj@meta.data, aes(x = Anno.Level.Fig.1, y = nCount_RNA, fill = "#94B3DE")) +
        geom_boxplot(outlier.colour = "red") +
        geom_boxplot(aes(color = "#94B3DE"), fill = NA, fatten = NULL, outlier.alpha = 0) +
        theme_classic(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "none"
        ) +
        xlab("") +
        ylab("Number of UMIs detected") +
        ylim(0, 20000) +
        scale_fill_manual(values = "#94B3DE") +
        scale_color_manual(values = "#94B3DE")
}

#-------------------------------------------------------------
# Rouge score of each cell type.

rogueSocrePerCellType <- function(rogue.res) {
    rogue.pf <- rogue.res %>% tidyr::gather(key = clusters, value = ROGUE)
    rogue.pf$clusters <- factor(rogue.pf$clusters, levels = ANNO_ENTIRE_IDNET_FIG1)

    ggplot(data = rogue.pf, aes(clusters, ROGUE)) +
        geom_boxplot(aes(fill = clusters), outlier.shape = NA) +
        geom_point(color = "#FF3E96", size = 1.5) +
        theme_bw() +
        theme(
            axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 13, colour = "black"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        labs(x = "", y = "ROGUE score") +
        scale_fill_manual(values = cols)
}

#----------------------------------------------------------
# ROE plot

ROEPlot <- function(R.oe) {
    plot.df <- R.oe %>%
        as.data.frame.matrix() %>%
        rownames_to_column(var = "CellType") %>%
        tidyr::gather(., "Status", "Roe", -CellType)
    plot.df$Status <- factor(plot.df$Status, levels = c("Prental", "Children", "Adult", "Older") %>% rev())
    plot.df$CellType <- factor(plot.df$CellType, levels = rownames(R.oe))
    plot.df$Color <- plot.df$Roe >= 1
    ggplot(plot.df) +
        geom_point(aes(x = CellType, y = Status, size = Roe, color = Color)) +
        scale_size(range = c(min(plot.df$Roe), max(plot.df$Roe))) +
        theme_bw(base_size = 14) +
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
}

#------------------------------------------------------------
# UMAP Seurat

annoUMAPRes <- function(merged.obj) {
    cols <- ANNO_ENTIRE_COLOR_FIG1 %>% `names<-`(1:34)
    p1 <- DimPlot(merged.obj, pt.size = 0.5, group.by = "Anno.Level.Cluster", cols = cols, label = T, label.size = 4.5) +
        NoLegend() + ggtitle(sprintf("Thymus annotation (n = %g)", dim(obj)[2]))
    leg.pos <- rbind.data.frame(
        data.frame(x = 1, y = (1:17) / 2, Col = cols[1:17], Cell = merged.obj$Anno.Level.Fig.1 %>% levels() %>% .[1:17] %>% rev(), Label = rev(1:17) %>% as.character()),
        data.frame(x = 1.2, y = (1:17) / 2, Col = cols[18:34], Cell = merged.obj$Anno.Level.Fig.1 %>% levels() %>% .[18:34] %>% rev(), Label = rev(18:34) %>% as.character())
    )
    leg.pos$Cell <- factor(leg.pos$Cell, levels = merged.obj$Anno.Level.Fig.1 %>% levels())
    leg.pos$Label <- factor(leg.pos$Label, levels = as.character(1:34))
    p2 <- ggplot(leg.pos, aes(x = x, y = y)) +
        geom_point(size = 7, shape = 21, aes(fill = Label, color = Label)) +
        theme_classic(base_size = 12) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none"
        ) +
        geom_text(aes(x = x, y = y, label = Label), size = 5) +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        geom_text(aes(x = x + 0.02, y = y, label = Cell), size = 5, hjust = 0) +
        xlim(1, 1.3)
    ggarrange(p1, p2, ncol = 2, widths = c(3, 1.5))
}

#---------------------------------------------------
# Marker bubble plot

bubbleMarkers <- function(merged.obj, marker.genes) {
    DotPlot(merged.obj, feature = unlist(marker.genes) %>% unique(), cols = c("lightgrey", "red")) + coord_flip() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom"
        ) + xlab("") + ylab("") +
        theme(panel.border = element_rect(colour = "black", size = 0.8)) +
        annotate(
            ymin = -1.05, ymax = 0,
            xmin = c(-Inf, 2, 11, 19, 33, 67),
            xmax = c(2, 11, 19, 33, 67, Inf),
            geom = "rect",
            fill = ANNO_SC_COLOR_LEVEL_1.New[1:6],
            color = ANNO_SC_COLOR_LEVEL_1.New[1:6],
        ) +
        annotate(
            xmin = 71, xmax = 72.5,
            ymin = c(-Inf, 1, 5, 9, 16, 33),
            ymax = c(1, 5, 9, 16, 33, Inf),
            geom = "rect",
            fill = ANNO_SC_COLOR_LEVEL_1.New[1:6],
            color = ANNO_SC_COLOR_LEVEL_1.New[1:6],
        )
}

#----------------------------------------------------
mcellPlotUmisPerCell <- function(mat, out.figs.dir, min.umis.cutoff = 500, max.umis.cutoff = 20000, bin.for.cutoff = 50, main.title = NULL, ...) {
    uc <- Matrix::colSums(mat)
    pdf(file.path(out.figs.dir, "total_umi_distr.pdf"), width = 6, height = 6)
    h <- hist(
        log2(uc),
        200,
        xlab = "Total cell UMIs (log2)",
        ylab = "Number of cells",
        main = main.title,
        col = "blue",
        border = "blue",
        cex.lab = 1.2,
        cex.axis = 1.2,
        cex.main = 1.1,
        cex.sub = 1.2,
        ...
    )
    abline(v = log2(min.umis.cutoff), col = "darkgreen", lty = 4, lwd = 2.0)
    abline(v = log2(max.umis.cutoff), col = "red", lty = 4, lwd = 2.0)
    text(log2(min.umis.cutoff) + 1, max(h$count) / 2, min.umis.cutoff, pos = 2, cex = 1.2)
    text(log2(max.umis.cutoff) + 1, max(h$count) / 2, max.umis.cutoff, pos = 2, cex = 1.2)

    text(log2(1 + min.umis.cutoff), max(h$count) * 0.75, sprintf("Max = %g\nMin = %g\nMean = %g", max(uc), min(uc), round(mean(uc))), pos = 4)

    dev.off()
    return(min.umis.cutoff)
}

mgenePlotCellsPerGene <- function(mat, out.figs.dir, main.title = NULL, color = "green", max.cut = 5000, min.cut = 200, ...) {
    uc <- Matrix::colSums(mat > 0)
    pdf(file.path(out.figs.dir, "total_gene_distr.pdf"), width = 6, height = 6)
    h <- hist(
        log2(uc + 1),
        200,
        xlab = "Number of gene expressed (log2)",
        ylab = "Number of cells",
        main = main.title,
        col = color,
        border = color,
        cex.lab = 1.2,
        cex.axis = 1.2,
        cex.main = 1.1,
        cex.sub = 1.2,
        ...
    )
    abline(v = log2(1 + max.cut), col = "red", lty = 4, lwd = 2.0)
    text(log2(max.cut + 1), max(h$count) / 2, max.cut, pos = 2, cex = 1.2)
    abline(v = log2(1 + min.cut) + 0.5, col = "darkgreen", lty = 4, lwd = 2.0)
    text(log2(min.cut + 1), max(h$count) / 2, min.cut, pos = 2, cex = 1.2)

    text(8, max(h$count) * 0.75, sprintf("Max = %g\nMin = %g\nMean = %g", max(uc), min(uc), round(mean(uc))), pos = 4)
    dev.off()
}

mcellPlotMitoProp <- function(mat, out.figs.dir, patten = "^MT-", color = "red",
                              main.title = "Fraction of mitochondrial gene expression per cell", file.name = "total_mito_distr.pdf", max.cut = 0.6, scale1 = 1, scale2 = 1, show.cutoff = TRUE) {
    mito.genes <- grep(patten, rownames(mat), v = TRUE, perl = TRUE)
    uc <- Matrix::colSums(mat)
    uc <- Matrix::colSums(mat[mito.genes, ]) / uc
    pdf(file.path(out.figs.dir, file.name), width = 6, height = 6)
    h <- hist(
        uc,
        200,
        xlab = "Fraction",
        ylab = "Number of cells",
        main = main.title,
        col = color,
        border = color,
        cex.lab = 1.2,
        cex.axis = 1.2,
        cex.main = 1.1,
        cex.sub = 1.2,
    )
    if (show.cutoff) {
        abline(v = max.cut, col = "red", lty = 4, lwd = 2.0)
        text(max.cut + 0.1 * scale1, max(h$count) / 2, max.cut, pos = 2, cex = 1.2)
    }
    text(max.cut + 0.05 * scale2, max(h$count) * 0.75, sprintf("Max = %g\nMin = %g\nMean = %g", round(max(uc), 3), round(min(uc), 3), round(mean(uc), 3)), pos = 4)
    dev.off()
}
