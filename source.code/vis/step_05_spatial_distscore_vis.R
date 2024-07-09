#----------------------------------------------
# Clustering high variable genes based on the spatial-seq data

clusterHyoerVariabGenesSp <- function(obj.st.lst, p.sig.gene, sample.size = 10000, up.cut = 2) {
    dist.lst <- distScaleToOne(obj.st.lst)
    p.sig.gene <- gsub("_", "-", rownames(p.sig.gene))

    expr.sub.dis <- lapply(names(obj.st.lst), function(lab) {
        expr.sub <- FetchData(obj.st.lst[[lab]], vars = p.sig.gene)
        expr.sub.dis <- cbind.data.frame(expr.sub, distance = dist.lst[[lab]])
        colnames(expr.sub.dis) <- gsub("-|/", "_", colnames(expr.sub.dis))
        expr.sub.dis
    }) %>%
        do.call(rbind, .) %>%
        as.data.frame()

    plot.df <- lapply(colnames(expr.sub.dis)[1:(dim(expr.sub.dis)[2] - 1)], function(gene) {
        form <- formula(sprintf("%s ~ ns(Distance, 4)", gene))
        mod <- glm(form, data = expr.sub.dis, family = gaussian)
        preds <- predict(mod, newdata = expr.sub.dis, type = "link", se.fit = TRUE)
        upr <- preds$fit + (1.96 * preds$se.fit)
        lwr <- preds$fit - (1.96 * preds$se.fit)
        fit <- preds$fit
        dat <- cbind.data.frame(label = rep(gene, dim(expr.sub.dis)[1]), Distance = expr.sub.dis[, "Distance"], fit = fit, lwr = lwr, upr = upr)
        dat <- dat[order(dat$Distance), ]
    })
    expr.df <- lapply(plot.df, function(xx) {
        xx$fit
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(colnames(expr.sub.dis)[1:(dim(expr.sub.dis)[2] - 1)]) %>%
        `colnames<-`(rownames(expr.sub.dis))

    expr.df <- expr.df[, sample(1:dim(expr.df)[2], 5000) %>% .[order(.)]]
    expr.df[is.na(expr.df)] <- 0
    expr.df <- expr.df %>%
        t() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    expr.df[is.na(expr.df)] <- 0

    expr.df[expr.df > up.cut] <- up.cut
    expr.df[expr.df < -up.cut] <- -up.cut

    pg.cls <- pheatmap::pheatmap(expr.df, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = F, show_colnames = F, silent = T)
    gene.cls <- cutree(pg.cls$tree_row, k = 2)
    return(list(genes = gene.cls, expr = expr.df))
}

#------------------------------------------------------
# Overlap genes

ovpGenesSig <- function(g1, g2) {
    require(VennDiagram)
    temp <- venn.diagram(
        list(A = g1, B = g2),
        fill = c("red", "green"),
        alpha = c(0.5, 0.5),
        cex = 1,
        cat.fontface = 2,
        lty = 2,
        filename = NULL
    )
    return(temp)
}

#------------------------------------------------------
# GO annotation

gannoRes <- function(up.genes, down.genes) {
    require(org.Hs.eg.db)
    up.genes.map <- clusterProfiler::bitr(up.genes, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
    down.genes.map <- clusterProfiler::bitr(down.genes, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

    ego.up <- clusterProfiler::enrichGO(
        gene = up.genes.map$ENTREZID %>% unique(),
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ego.down <- clusterProfiler::enrichGO(
        gene = down.genes.map$ENTREZID %>% unique(),
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = TRUE
    )
    return(list(ego.up = ego.up, ego.down = ego.down))
}

#-------------------------------------------------------
# Plot Heatmap for annotation results.

plotComplexHeatmapWithAnno <- function(expr.df, ego.up, ego.down, idx.1 = c(1, 2, 4), idx.2 = c(1, 2, 5)) {
    library(ComplexHeatmap)
    library(circlize)
    myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
    myBreaks <- seq(-4, 4, length.out = 100)
    ncut.medulla <- ceiling(ncol(expr.df) * 0.25)
    label.sample <- c(rep("Medulla", ncut.medulla), rep("Cortex", ncol(expr.df) - ncut.medulla))
    dis.vec <- 1:ncol(expr.df)
    anno.cols <- cbind.data.frame(Region = label.sample, Distance = (dis.vec - 1) / max(dis.vec - 1)) %>% `rownames<-`(colnames(expr.df))
    anno.cols$Distance[anno.cols$Region == "Cortex"] <- anno.cols$Distance[anno.cols$Region == "Cortex"] %>% rev()
    anno.cols$Distance[anno.cols$Region == "Medulla"] <- anno.cols$Distance[anno.cols$Region == "Medulla"] %>% rev()
    fo <- HeatmapAnnotation(Region = anno.cols$Region, Distance = anno.cols$Distance, col = list(Region = c(Cortex = "#00B6EB", Medulla = "#FB61D7")))
    stage.genes <- c("ST18", "HIVEP3", "RGPD3", "SMPD3", "AQP3", "RORC", "SATB1", "TOX2", "HOPX", "MEF2C", "TCF7")
    ovp.genes.1 <- intersect(rownames(expr.df), stage.genes)
    ovp.genes.2 <- intersect(rownames(expr.df), c(medulla.genes[[1]], "FOXP3"))
    genes.markes <- c(ovp.genes.1, ovp.genes.2)
    at.pos <- lapply(genes.markes, function(xx) which(rownames(expr.df) == xx)) %>% unlist()
    hp <- Heatmap(
        expr.df,
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        col = colorRamp2(myBreaks, myCol),
        left_annotation = rowAnnotation(
            foo = anno_block(
                gp = gpar(fill = 2:3),
                labels = c("Cluster1", "Cluster2"),
                labels_gp = gpar(col = "white", fontsize = 14)
            )
        ),
        row_km = 2,
        row_gap = unit(0, "mm"),
        right_annotation = rowAnnotation(
            foo3 = anno_mark(at = at.pos, labels = genes.markes),
            foo = anno_block(
                gp = gpar(fill = 2:3),
                labels = c(data.frame(ego.up)[idx.1[1], 2], data.frame(ego.down)[idx.2[1], 2]),
                labels_gp = gpar(col = "white", fontsize = 14)
            ),
            foo1 = anno_block(
                gp = gpar(fill = 2:3),
                labels = c(data.frame(ego.up)[idx.1[2], 2], data.frame(ego.down)[idx.2[2], 2]),
                labels_gp = gpar(col = "white", fontsize = 14)
            ),
            foo2 = anno_block(
                gp = gpar(fill = 2:3),
                labels = c(data.frame(ego.up)[idx.1[3], 2], data.frame(ego.down)[idx.2[3], 2]),
                labels_gp = gpar(col = "white", fontsize = 14)
            )
        ),
        top_annotation = fo,
        heatmap_legend_param = list(direction = "horizontal"),
        column_split = anno.cols$Region
    )
    return(hp)
}
