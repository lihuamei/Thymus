trgStat <- function(merged.obj) {
    trg.stat <- lapply(T.CELL.GROUPS, function(tg) {
        if (!(tg %in% (merged.obj$t.cell_group %>% unique()))) {
            xx <- c(0, 0) %>% `names<-`(c(FALSE, TRUE))
            return(xx)
        }
        obj.tmp <- subset(merged.obj, t.cell_group == tg)
        xx <-
            {
                obj.tmp$trg_sum > 2
            } %>% table()
        if (length(xx) > 1) {
            return(xx)
        }
        if (F %in% names(xx)) {
            xx <- c(xx, 0) %>% `names<-`(c(FALSE, TRUE))
        } else {
            xx <- c(0, xx) %>% `names<-`(c(FALSE, TRUE))
        }
        return(xx)
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(T.CELL.GROUPS) %>%
        `colnames<-`(c("TRG-", "TRG+"))
    trg.stat
}

#----------------------------------------------------
# Statistics the count of eight T subgroups.

baseStatOfEightSubgroups <- function(merged.obj, trg.stat) {
    cnt.all <- rep(0, length(T.CELL.GROUPS)) %>% `names<-`(T.CELL.GROUPS)
    cnt.all.bak <- table(merged.obj$t.cell_group)
    cnt.all[names(cnt.all.bak)] <- cnt.all.bak
    tgroup.stat <- cbind.data.frame(Group = cnt.all[rownames(trg.stat)], trg.stat)
    tgroup.stat <- cbind.data.frame("Group.Var1" = rownames(tgroup.stat), tgroup.stat)
    colnames(tgroup.stat)[2] <- "Group.Freq"
    tgroup.stat <- tgroup.stat[order(tgroup.stat$Group.Freq) %>% rev(), ]
    tgroup.stat$Group.Var1 <- factor(tgroup.stat$Group.Var1, levels = unique(tgroup.stat$Group.Var1) %>% rev())
    g1 <- ggplot(
        data = tgroup.stat,
        aes(x = Group.Var1, y = Group.Freq, fill = Group.Var1)
    ) +
        geom_bar(stat = "identity") +
        coord_flip() +
        xlab("Subgroup of T cells") +
        ylab("Count") +
        theme_bw(base_size = 12) +
        theme(legend.position = "none", panel.grid.major = element_line(color = "gray", size = 0.5, linetype = 2), panel.grid.minor = element_blank()) +
        scale_x_discrete(labels = paste0(tgroup.stat[, 1], "(", tgroup.stat[, 2], ")") %>% rev()) +
        scale_fill_manual(values = T.CELL.GROUPS.EIGHT.COLORS)

    percent.df <- sweep(tgroup.stat[, c("TRG-", "TRG+")], 1, rowSums(tgroup.stat[, c("TRG-", "TRG+")]), "/") %>% cbind.data.frame(., Group = rownames(.))
    percent.df <- tidyr::gather(percent.df, "Type", "Percentage", -Group)
    percent.df$Group <- factor(percent.df$Group, levels = tgroup.stat[, 1] %>% unique())
    g2 <- ggplot(percent.df, aes(x = "", y = Percentage, fill = Type)) +
        geom_bar(width = 1, stat = "identity", color = "white") +
        coord_polar("y", start = 0) +
        facet_wrap(~Group, nrow = 8) +
        theme_void() +
        theme(strip.text.x = element_blank()) +
        scale_fill_manual(values = c("gray", "red")) +
        geom_text(aes(label = paste0(round(Percentage * 100, 2), "%")), position = position_stack(vjust = 0.5), color = "black", size = 3)
    g1 + g2
}

#--------------------------------------------------
# Venn Diagram plot for T subgroups.

vennDiagForTSubgroups <- function(obj) {
    percent.df <- obj$t.cell_group %>%
        table() %>%
        as.data.frame() %>%
        `colnames<-`(c("Subgroup", "Freq"))
    percent.df[, 3] <- percent.df[, 2]
    percent.df[, 2] <- percent.df[, 2] / sum(percent.df[, 2])
    idxes <- which(percent.df$Freq <= 0.01)
    percent.df.sub <- percent.df[which(percent.df$Freq > 0.01), ]
    percent.df.sub[, 1] <- as.vector(percent.df.sub[, 1])
    pie.df <- rbind.data.frame(percent.df.sub, Other = c("Others", 1 - sum(percent.df.sub[, 2]), percent.df[idxes, 3] %>% sum()))
    pie.df[, 2] <- as.numeric(pie.df[, 2])
    pie.df[, 3] <- as.integer(pie.df[, 3])
    pie.df$Label <- sprintf("%s\n(%3.1f%s, n = %d)", pie.df[, 1], 100 * pie.df[, 2], "%", pie.df[, 3])
    pie.df <- pie.df[order(pie.df[, 2]), ]
    pie.df[, 1] <- factor(pie.df[, 1], levels = pie.df[, 1])
    pie(
        pie.df[, 3],
        clockwise = TRUE,
        label = pie.df[, "Label"],
        border = "white",
        cex = 1.0,
        radius = 1.0,
        col = c(T.CELL.GROUPS.COLORS, "Others" = "grey")[pie.df[, 1] %>% levels()],
        cex.main = 1.0,
        main = sprintf("Percentage and number of cells for T subgroups\n(N = %d)", sum(pie.df[, 3]))
    )
}

#----------------------------------------------------
# Subgroups of T cells for each sample

fourSubgroupsPerSample <- function(merged.obj, tcr.subtypes) {
    merged.obj.sub <- subset(merged.obj, t.cell_group %in% tcr.subtypes)
    count.df <- table(merged.obj.sub$orig.ident, merged.obj.sub$t.cell_group) %>% as.data.frame.matrix(.)
    percent.df <- sweep(count.df, 1, rowSums(count.df), "/") %>%
        cbind.data.frame(., Sample = rownames(.)) %>%
        reshape2::melt(.) %>%
        `colnames<-`(c("Sample", "Subgroup", "Percentage"))
    percent.df$Subgroup <- factor(percent.df$Subgroup, levels = tcr.subtypes)
    percent.df$Sample <- factor(percent.df$Sample, levels = names(SampleInfos))

    gp1 <- ggplot(percent.df, aes(x = Sample, y = Percentage, fill = Subgroup)) +
        geom_bar(position = "stack", stat = "identity") +
        theme_bw(base_size = 16) +
        scale_fill_manual(values = T.CELL.GROUPS.COLORS) +
        scale_x_discrete(expand = c(0, 0), labels = paste0(percent.df$Sample %>% levels(), ":(", SampleInfos %>% unlist() %>% as.vector(), ")")) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
            legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ylab("Proportion of T subgroups") +
        xlab("")
    gp1
}

#---------------------------------------------------
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


#------------------------------------------------
# Correlation between age and proportions of T subgroups

corrBetweenAgeAndPropSubgroup <- function(merged.obj) {
    merged.obj.sub <- subset(merged.obj, t.cell_group %in% tcr.subtypes)
    count.df <- table(merged.obj.sub$orig.ident, merged.obj.sub$t.cell_group) %>% as.data.frame.matrix(.)
    percent.df <- sweep(count.df, 1, rowSums(count.df), "/") %>%
        cbind.data.frame(., Age = 1:nrow(.), Sample = rownames(.)) %>%
        reshape2::melt(., id.vars = c("Age", "Sample"))
    percent.df$Sample <- factor(percent.df$Sample, levels = names(SampleInfos))
    percent.df$variable <- factor(percent.df$variable, levels = tcr.subtypes)

    ggplot(percent.df, aes(x = Age, y = value)) +
        geom_point(aes(fill = Sample), shape = 21, size = 3) +
        facet_wrap(~variable, scale = "free") +
        geom_smooth(method = "lm") +
        stat_cor(method = "pearson") +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        ) +
        xlab("-> From Prental to Adult") +
        ylab("Proportion") +
        scale_fill_manual(values = SAMPLE.COLORS[1:16])
}

#-------------------------------------------------
# Number of TCR clonetypes

numberOfClonetypes <- function(merged.obj) {
    cnt.tcr <- table(merged.obj$tcr, merged.obj$orig.ident) %>%
        as.data.frame.matrix() %>%
        t() %>%
        .[, 1]
    cnt.tcr.unq <- lapply(merged.obj$orig.ident %>% levels(), function(xx) {
        subset(merged.obj, orig.ident == xx)$tcr_cdr3s_aa %>%
            unique() %>%
            length()
    }) %>% unlist()
    plot.df <- rbind.data.frame(cbind.data.frame(Count = cnt.tcr, Class = "TCR-entire"), cbind.data.frame(Count = cnt.tcr.unq, Class = "TCR-unique"))
    plot.df <- cbind.data.frame(plot.df, Sample = rep(names(cnt.tcr), 2))
    plot.df$Sample <- factor(plot.df$Sample, levels = names(cnt.tcr))

    ggplot(plot.df, aes(x = Sample, y = Count, fill = Class, alpha = 0.8)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_bw(base_size = 16) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_x_discrete(expand = c(0, 0), labels = paste0(plot.df$Sample %>% levels(), ":(", SampleInfos %>% unlist() %>% as.vector(), ")")) +
        scale_y_continuous(expand = c(0, 0), limit = c(0, 9500)) +
        xlab("") +
        ylab("Clonotypes") +
        geom_text(aes(label = Count, group = Class), size = 4, position = position_dodge(width = 1), vjust = -0.5) +
        geom_hline(yintercept = 5000, linetype = "dashed", color = "gray13", size = 0.8) +
        ggtitle("Number of clonotypes\nNumber of unique clonotypes in each sample") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
}

#--------------------------------------------------
# TCR diversity analysis

tcrDiversityPermutation <- function(plot.df) {
    plot.df$Type <- factor(plot.df$Type, levels = names(AGE_STATUS_EXTENT_COLOR[1:4]))
    data_summary <- function(data, varname, groupnames) {
        require(plyr)
        summary_func <- function(x, col) {
            c(
                mean = mean(x[[col]], na.rm = TRUE),
                sd = sd(x[[col]], na.rm = TRUE)
            )
        }
        data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
        data_sum <- rename(data_sum, c("mean" = varname))
        return(data_sum)
    }
    plot.df.new <- data_summary(plot.df, "Count", c("Type", "X"))
    ggplot(plot.df.new, aes(x = X, y = Count)) +
        geom_point(aes(color = Type)) +
        geom_errorbar(aes(ymin = Count - sd, ymax = Count + sd, color = Type), width = .2, position = position_dodge(.9)) +
        geom_smooth(method = "lm", se = TRUE, color = "darkblue", size = 0.8, fullrange = F) +
        theme_classic(base_size = 16) +
        scale_color_manual(values = AGE_STATUS_EXTENT_COLOR[1:4]) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        stat_cor(method = "pearson") +
        xlab("") +
        ylab("Number of unique TCRs (100 permutations)") +
        xlab("Prental -> Adult")
}

#----------------------------------------------
# TCR diversity using Shannon-entropy

tcrDivShannonEntro <- function(divs.df) {
    plot.df <- tidyr::gather(as.data.frame(divs.df), "Sample", "Diversity")
    plot.df$Status <- "Children"
    plot.df$Status[plot.df[, "Sample"] %in% names(SampleInfos)[1:4]] <- "Prental"
    plot.df$Status[plot.df[, "Sample"] %in% names(SampleInfos)[13:14]] <- "Adult"
    plot.df$Status[plot.df[, "Sample"] %in% names(SampleInfos)[15:16]] <- "Older"
    plot.df$Status <- factor(plot.df$Status, levels = names(AGE_STATUS_EXTENT_COLOR[1:4]))

    my.comparison <- list(c("Prental", "Children"), c("Children", "Adult"), c("Adult", "Older"))
    ggplot(plot.df, aes(x = Status, y = Diversity, fill = Status)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = AGE_STATUS_EXTENT_COLOR[1:4]) +
        geom_jitter(aes(color = Sample), size = 1.2, alpha = 0.5) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ylab("Diversity of TCR repetoire (100 permutations)") +
        xlab("") +
        scale_color_manual(values = SAMPLE.COLORS[1:16]) +
        stat_compare_means(comparisons = my.comparison) +
        guides(colour = guide_legend(override.aes = list(size = 4)))
}

#------------------------------------------
# Alpha/Beta analysis

alphaBetaChainUMAP <- function(obj) {
    obj$Class <- ""
    obj$Class[obj$trb != ""] <- "TRB"
    obj$Class[obj$tra != ""] <- "TRA"
    obj$Class[(obj$trb != "") & (obj$tra != "")] <- "TRAB"

    meta.data <- obj@meta.data
    meta.data <- cbind.data.frame(meta.data, obj[["umapn"]]@cell.embeddings)
    meta.data$Class <- as.vector(meta.data$Class)
    meta.data$Class[meta.data$Class == ""] <- "Other"
    meta.data$Class <- factor(meta.data$Class, levels = c("TRA", "TRB", "TRAB", "Other"))

    ggplot(meta.data, aes(x = UMAPN_1, y = UMAPN_2, color = Class)) +
        geom_point(aes(size = Class)) +
        scale_size_manual(values = c(TRA = 1, TRB = 0.5, TRAB = 0.5, Other = 0.1)) +
        scale_color_manual(values = c(TRA = "#8032BC", TRB = "#005297", TRAB = "#00B17F", Other = "grey")) +
        theme_classic() +
        ggtitle("Alpha/Beta rearrangement")
}

#-------------------------------------------------
# Clone size distribution

cloneSizePlot <- function(merged.obj) {
    freq.df <- cbind.data.frame(Count = merged.obj$tcr_frequency, Sample = merged.obj$orig.ident)
    group.vec <- merged.obj$tcr_frequency
    group.vec[group.vec >= 3] <- ">=3"
    freq.df <- cbind.data.frame(freq.df, Group = group.vec)

    freq.df$Group <- factor(freq.df$Group, levels = c("1", "2", ">=3"))
    ggplot(freq.df, aes(x = 1, fill = Group)) +
        geom_bar(position = "fill") +
        coord_polar(theta = "y") +
        facet_wrap(~Sample) +
        theme_bw() +
        ylab("") +
        xlab("") +
        scale_fill_manual(values = c("cornsilk2", "darkorange", "darkred")) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        guides(fill = guide_legend(title = "Clone size"))
}

#-------------------------------------------------
# Clone Size project to spatial data.

cloneSizeProjectSp <- function(obj.st.lst, top20) {
    gps <- lapply(obj.st.lst, function(xx) {
        xx <- AddModuleScore(xx, features = split(top20$gene, top20$cluster))
        gp.1 <- SpatialFeaturePlot(xx, feature = "Cluster1", image.alpha = 0)
        gp.2 <- SpatialFeaturePlot(xx, feature = "Cluster2", image.alpha = 0)
        gp.3 <- SpatialFeaturePlot(xx, feature = "Cluster3", image.alph = 0)
        ggarrange(gp.1, gp.2, gp.3, nrow = 3, ncol = 1)
    })
    names(gps) <- sp2SampleName(names(gps))
    gps <- gps[names(SC_MAP_SP)]
    ggarrange(plotlist = gps, nrow = 1, ncol = 8)
}
