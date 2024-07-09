tcrGroupsProj <- function(pseud.coord.lst, obj.sc, obj.sp.lst) {
    obj <- obj.sc
    obj.tcr <- tcrSubgroup(obj)
    tcr.subtypes <- c("TRA-TRB+TRD+", "TRA+TRB+TRD+", "TRA-TRB+TRD-", "TRA+TRB+TRD-")

    res.lst <- lapply(pseud.coord.lst, function(xx) {
        yy.names <- gsub("SC_", "", xx$sc)
        lapply(tcr.subtypes, function(tc) {
            cells.1 <- subset(obj.tcr, t.cell_group == tc) %>% Cells()
            ovp.names <- intersect(yy.names, cells.1)
            xx.sub <- xx[yy.names %in% ovp.names, ]
            yy <- cbind.data.frame(xx.sub, obj.tcr@meta.data[gsub("SC_", "", xx.sub$sc), ])
            yy.tcr <- subset(yy, Anno.Level.Fig.1 %in% ANNO_ENTIRE_IDNET_FIG1[17:33])
            yy.tcr$Anno.Level.Fig.1 <- factor(yy.tcr$Anno.Level.Fig.1, levels = ANNO_ENTIRE_IDNET_FIG1[17:33])
            yy.tcr$Group <- tc
            return(yy.tcr)
        }) %>%
            do.call(rbind, .) %>%
            as.data.frame()
    })
    figs.dir <- file.path(out.figs.dir, "sc.project")
    gps <- lapply(names(res.lst), function(sn.name) {
        img <- obj.sp.lst[[sn.name]]@images[[sn.name]]@image
        pseud.coord <- res.lst[[sn.name]]
        pseud.coord.bak <- pseud.coord.lst[[sn.name]]
        edge <- subset(pseud.coord.bak, HE == "Medulla_edge")
        edge.sub <- edge[!duplicated(edge[, c("x.spot", "y.spot")]), ]
        edge.sub$x <- edge.sub[, "x.spot"]
        edge.sub$y <- edge.sub[, "y.spot"]
        edge.sub$Group <- "Edge"
        edge.sub[, setdiff(colnames(pseud.coord), colnames(edge.sub))] <- "NA"
        edge.sub <- edge.sub[, colnames(pseud.coord)]

        plot.df <- rbind.data.frame(pseud.coord, edge.sub)
        spot.shape <- c(rep(21, length(T.CELL.GROUPS.COLORS) + 1), 8) %>% `names<-`(c(names(T.CELL.GROUPS.COLORS), "Edge"))

        gps.sub <- lapply(names(T.CELL.GROUPS.COLORS), function(xx) {
            plot.df.sub <- subset(pseud.coord, Group == xx & (cellType %in% ANNO_ENTIRE_IDNET_FIG1[c(17:24, 27, 28)]))
            plot.df.sub <- rbind.data.frame(plot.df.sub, edge.sub)
            gp <- ggplot(plot.df.sub, aes(x = y, y = x, fill = Group, shape = Group), alpha = 0.3) +
                annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                geom_point(size = 1, shape = 21) +
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
                xlab("") + # facet_wrap(~Group, nrow = 2) +
                scale_fill_manual(values = c(T.CELL.GROUPS.COLORS, Edge = "white")) +
                scale_shape_manual(values = spot.shape) +
                guides(fill = guide_legend(override.aes = list(size = 4))) +
                ggtitle(xx)
            return(gp)
        })
        gp <- ggarrange(plotlist = gps.sub, nrow = 2, ncol = 2)
        ggsave(file = file.path(figs.dir, sprintf("%s.sc.project.sp.UMAP.pdf", sp2SampleName(sn.name))), width = 10, height = 8)
    }) %>% `names<-`(names(res.lst))
}
