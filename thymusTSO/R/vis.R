#' @title fitDistLinesByWindows
#'
#' @description Curves of gene expression with spatial distance.
#' @param obj.st.lst A list of spatial seurat objects.
#' @param plot.tar Vector of gene names or metadata column names to plot.
#' @param win Integer specifying the window size for averaging data points.
#' @return ggplot object.
#' @export fitDistLinesByWindows

fitDistLinesByWindows <- function(obj.st.lst, plot.tar, win = 10) {
    dist.lst <- distScaleToOne(obj.st.lst)
    plot.df <- lapply(names(obj.st.lst), function(sn) {
        xx.obj <- obj.st.lst[[sn]]
        meta.cols <- intersect(colnames(xx.obj@meta.data), plot.tar)
        features <- intersect(rownames(xx.obj), plot.tar)

        plot.df <- data.frame(Distance = dist.lst[[sn]][rownames(xx.obj@meta.data), ])
        if (length(meta.cols) > 0) {
            plot.df <- cbind.data.frame(xx.obj@meta.data[, meta.cols, drop = FALSE] %>% as.data.frame(), plot.df)
        }
        if (length(features) > 0) {
            plot.df <- cbind.data.frame(FetchData(xx.obj, vars = features, slot = "data"), plot.df)
        }
        return(plot.df)
    }) %>% do.call(rbind, .)

    lab.cnts <- lapply(obj.st.lst, function(xx.obj) {
        cnts <- xx.obj@meta.data$HE.Labels %>% table()
        medulla <- cnts[names(cnts) %in% c("Medulla_centric", "Medulla_edge", "Medulla_hi")] %>% sum()
        return(c(all = sum(cnts), medulla = medulla))
    }) %>% do.call(rbind, .)

    percent <- sum(lab.cnts[, 2]) / sum(lab.cnts[, 1])
    xx <- plot.df[order(plot.df[, "Distance"]), ]
    avg.dat.2.1 <- apply(xx[, 1:(dim(xx)[2])] %>% as.data.frame(), 2, function(x) {
        zoo::rollapply(x, width = win, by = win - 1, FUN = mean)
    }) %>% `colnames<-`(colnames(xx)[1:(dim(xx)[2])])

    if (length(colnames(avg.dat.2.1)) > 1) avg.dat.2.1 <- avg.dat.2.1[, c(plot.tar[plot.tar %in% colnames(avg.dat.2.1)], colnames(xx)[ncol(xx)])]
    avg.dat.2.1[, 1:(ncol(avg.dat.2.1) - 1)] <- avg.dat.2.1[, 1:(ncol(avg.dat.2.1) - 1)] %>% scale()
    ctypes <- plot.tar[plot.tar %in% colnames(avg.dat.2.1)]
    avg.dat.2.plot <- tidyr::gather(avg.dat.2.1 %>% as.data.frame() %>% .[, c(ctypes, "Distance")], "CellType", "Score", -Distance)
    colnames(avg.dat.2.plot) <- c("cell_id", "label", "value")
    avg.dat.2.plot$label <- factor(avg.dat.2.plot$label, levels = plot.tar[plot.tar %in% unique(avg.dat.2.plot$label)])

    gp <- ggplot(avg.dat.2.plot, aes(x = cell_id, y = value, color = label)) +
        theme_classic(base_size = 14) +
        geom_smooth(se = TRUE, method = lm, formula = y ~ splines::ns(x, 10)) +
        theme(strip.background = element_blank(), strip.text.x = element_text(size = 0)) +
        labs(y = "Average expression", x = "Spot order along distance") +
        geom_vline(xintercept = 1 * percent, color = "blue", lty = "dashed", lwd = 0.5) +
        xlim(0, 1)
    return(gp)
}

#' @title fitDistLinesByGlm
#'
#' @description Fits generalized linear models (GLMs) to gene expression data against spatial distance and visualizes the fitted curves with confidence intervals.
#' @param obj.st.lst A list of spatial seurat objects.
#' @param plot.tar Vector of gene names or metadata column names to plot.
#' @param degree Integer specifying the degree of splines for smoothing.
#' @return ggplot object.
#' @export fitDistLinesByGlm

fitDistLinesByGlm <- function(obj.st.lst, plot.tar, degree = 4) {
    dist.lst <- distScaleToOne(obj.st.lst)
    plot.df <- lapply(names(obj.st.lst), function(sn) {
        xx.obj <- obj.st.lst[[sn]]
        meta.cols <- intersect(colnames(xx.obj@meta.data), plot.tar)
        features <- intersect(rownames(xx.obj), plot.tar)

        plot.df <- data.frame(Distance = dist.lst[[sn]][rownames(xx.obj@meta.data), ])
        if (length(meta.cols) > 0) {
            plot.df <- cbind.data.frame(xx.obj@meta.data[, meta.cols, drop = FALSE] %>% as.data.frame() %>% scale(.), plot.df)
        }
        if (length(features) > 0) {
            plot.df <- cbind.data.frame(FetchData(xx.obj, vars = features, slot = "data"), plot.df)
        }
        return(plot.df)
    }) %>% do.call(rbind, .)

    plot.df <- subset(plot.df, Distance <= 1)

    lab.cnts <- lapply(obj.st.lst, function(xx.obj) {
        cnts <- xx.obj@meta.data$HE.Labels %>% table()
        medulla <- cnts[names(cnts) %in% c("Medulla_centric", "Medulla_edge", "Medulla_hi")] %>% sum()
        return(c(all = sum(cnts), medulla = medulla))
    }) %>% do.call(rbind, .)

    percent <- sum(lab.cnts[, 2]) / sum(lab.cnts[, 1])
    plot.df <- plot.df[order(plot.df[, "Distance"]), ]
    colnames(plot.df) <- gsub("-|/|\\(|\\)", "_", colnames(plot.df))
    plot.df.new <- lapply(colnames(plot.df)[1:(dim(plot.df)[2] - 1)], function(gene) {
        form <- formula(sprintf("%s ~ ns(Distance, %g)", gene, degree))
        mod <- glm(form, data = plot.df, family = gaussian)
        preds <- predict(mod, newdata = plot.df, type = "link", se.fit = TRUE)
        upr <- preds$fit + (1.96 * preds$se.fit)
        lwr <- preds$fit - (1.96 * preds$se.fit)
        fit <- preds$fit
        dat <- cbind.data.frame(label = rep(gene, dim(plot.df)[1]), Distance = plot.df[, "Distance"], fit = fit, lwr = lwr, upr = upr)
    }) %>% do.call(rbind, .)

    plot.tar <- gsub("-|/|\\(|\\)", "_", plot.tar)
    plot.df.new$label <- factor(plot.df.new$label, levels = plot.tar[plot.tar %in% unique(plot.df.new$label)])
    gp <- ggplot(plot.df.new, aes(x = Distance, y = fit, color = label)) +
        geom_ribbon(aes(ymin = lwr, ymax = upr, color = NULL, fill = label), alpha = 0.1) +
        geom_line(aes(y = fit, linetype = label), size = 1) +
        theme_classic() +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_text(size = 0)
        ) +
        labs(y = "Average expression", x = "Spot order along distance") +
        geom_vline(xintercept = quantile(plot.df.new$Distance, percent), color = "blue", lty = "dashed", lwd = 0.8)
    return(gp)
}
