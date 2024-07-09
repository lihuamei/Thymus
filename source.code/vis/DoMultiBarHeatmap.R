# This function copy from https://github.com/satijalab/seurat/issues/2201.

suppressPackageStartupMessages({
    library(rlang)
    library(grid)
    library(ComplexHeatmap)
    library(circlize)
})

DoMultiBarHeatmap <- function(object,
                              features = NULL,
                              cells = NULL,
                              group.by = "ident",
                              additional.group.by = NULL,
                              additional.group.sort.by = NULL,
                              cols.use = NULL,
                              group.bar = TRUE,
                              disp.min = -2.5,
                              disp.max = NULL,
                              slot = "scale.data",
                              assay = NULL,
                              label = TRUE,
                              size = 5.5,
                              hjust = 0,
                              angle = 45,
                              raster = TRUE,
                              draw.lines = TRUE,
                              lines.width = NULL,
                              group.bar.height = 0.02,
                              combine = TRUE) {
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = object)[cells]
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    features <- features %||% VariableFeatures(object = object)
    ## Why reverse???
    features <- rev(x = unique(x = features))
    disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
        yes = 2.5, no = 6
    )
    possible.features <- rownames(x = GetAssayData(
        object = object,
        slot = slot
    ))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if (length(x = features) == 0) {
            stop(
                "No requested features found in the ", slot,
                " slot for the ", assay, " assay."
            )
        }
        warning(
            "The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                collapse = ", "
            )
        )
    }

    if (!is.null(additional.group.sort.by)) {
        if (any(!additional.group.sort.by %in% additional.group.by)) {
            bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
            additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
            if (length(x = bad.sorts) > 0) {
                warning(
                    "The following additional sorts were omitted as they were not a subset of additional.group.by : ",
                    paste(bad.sorts, collapse = ", ")
                )
            }
        }
    }

    data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
        object = object,
        slot = slot
    )[features, cells, drop = FALSE])))

    object <- suppressMessages(expr = StashIdent(
        object = object,
        save.name = "ident"
    ))
    group.by <- group.by %||% "ident"
    groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
    plots <- list()
    for (i in group.by) {
        data.group <- data
        if (!is_null(additional.group.by)) {
            additional.group.use <- additional.group.by[additional.group.by != i]
            if (!is_null(additional.group.sort.by)) {
                additional.sort.use <- additional.group.sort.by[additional.group.sort.by != i]
            } else {
                additional.sort.use <- NULL
            }
        } else {
            additional.group.use <- NULL
            additional.sort.use <- NULL
        }

        group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]

        for (colname in colnames(group.use)) {
            if (!is.factor(x = group.use[[colname]])) {
                group.use[[colname]] <- factor(x = group.use[[colname]])
            }
        }

        if (draw.lines) {
            lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                0.0025)
            placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) *
                lines.width), FUN = function(x) {
                return(Seurat:::RandomName(length = 20))
            })
            placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
            group.levels <- list()
            group.levels[[i]] <- levels(x = group.use[[i]])
            for (j in additional.group.use) {
                group.levels[[j]] <- levels(x = group.use[[j]])
                placeholder.groups[[j]] <- NA
            }

            colnames(placeholder.groups) <- colnames(group.use)
            rownames(placeholder.groups) <- placeholder.cells

            group.use <- sapply(group.use, as.vector)
            rownames(x = group.use) <- cells

            group.use <- rbind(group.use, placeholder.groups)

            for (j in names(group.levels)) {
                group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
            }

            na.data.group <- matrix(
                data = NA, nrow = length(x = placeholder.cells),
                ncol = ncol(x = data.group), dimnames = list(
                    placeholder.cells,
                    colnames(x = data.group)
                )
            )
            data.group <- rbind(data.group, na.data.group)
        }

        order_expr <- paste0("order(", paste(c(i, additional.sort.use), collapse = ","), ")")
        group.use <- with(group.use, group.use[eval(parse(text = order_expr)), , drop = F])

        plot <- Seurat:::SingleRasterMap(
            data = data.group, raster = raster,
            disp.min = disp.min, disp.max = disp.max, feature.order = features,
            cell.order = rownames(x = group.use), group.by = group.use[[i]]
        )

        if (group.bar) {
            pbuild <- ggplot_build(plot = plot)
            group.use2 <- group.use
            cols <- list()
            na.group <- Seurat:::RandomName(length = 20)
            for (colname in rev(x = colnames(group.use2))) {
                if (colname == i) {
                    colid <- paste0("Identity (", colname, ")")
                } else {
                    colid <- colname
                }

                # Default
                cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))

                # Overwrite if better value is provided
                if (!is_null(cols.use[[colname]])) {
                    req_length <- length(x = levels(group.use))
                    if (length(cols.use[[colname]]) < req_length) {
                        warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
                    } else {
                        if (!is_null(names(cols.use[[colname]]))) {
                            if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
                            } else {
                                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse = ","), ") are not represented.")
                            }
                        } else {
                            cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
                        }
                    }
                }

                # Add white if there's lines
                if (draw.lines) {
                    levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
                    group.use2[placeholder.cells, colname] <- na.group
                    cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
                }
                names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])

                y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
                y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
                y.max <- y.pos + group.bar.height * y.range
                pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

                plot <- suppressMessages(plot +
                    annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]), xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                    annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                    coord_cartesian(ylim = c(0, y.max), clip = "off"))

                if ((colname == i) && label) {
                    x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                    x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
                    group.use$x <- x.divs
                    label.x.pos <- tapply(
                        X = group.use$x, INDEX = group.use[[colname]],
                        FUN = median
                    ) * x.max
                    label.x.pos <- data.frame(
                        group = names(x = label.x.pos),
                        label.x.pos
                    )
                    plot <- plot + geom_text(
                        stat = "identity",
                        data = label.x.pos, aes_string(
                            label = "group",
                            x = "label.x.pos"
                        ), y = y.max + y.max *
                            0.03 * 0.5, angle = angle, hjust = hjust,
                        size = size
                    )
                    plot <- suppressMessages(plot + coord_cartesian(ylim = c(
                        0,
                        y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) *
                            size
                    ), clip = "off"))
                }
            }
        }
        plot <- plot + theme(line = element_blank())
        plots[[i]] <- plot
    }
    if (combine) {
        plots <- CombinePlots(plots = plots)
    }
    return(plots)
}

# This function copy from Scillus package (see details from https://scillus.netlify.app/)

plot_heatmap <- function(dataset,
                         markers,
                         sort_var = c("seurat_clusters"),
                         n = 8,
                         anno_var,
                         anno_colors,
                         hm_limit = c(-2, 0, 2),
                         hm_colors = c("#4575b4", "white", "#d73027"),
                         row_font_size = 12,
                         fc = T,
                         DEG.color = T.CELL.GROUPS.COLORS,
                         cluster_rows = F,
                         left_annotation = T,
                         row_names_side = "left",
                         add.features = NULL,
                         shuffle = F) {
    mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
    mat <- MinMax(mat, min = -2.5, max = 2.5)
    if (is.data.frame(markers)) {
        # genes <- get_top_genes(dataset, markers, n)
        if (fc) {
            topn <- markers %>%
                group_by(cluster) %>%
                top_n(n = n, wt = avg_log2FC)
        } else {
            markers$p_val <- 1 - markers$p_val
            topn <- markers %>%
                group_by(cluster) %>%
                top_n(n = n, wt = p_val)
        }
        topn <- topn[!duplicated(topn$gene), c("cluster", "gene")]
        genes <- topn$gene
        ctype <- topn$cluster
    } else if (is.character(markers)) {
        genes <- markers
    } else {
        stop("Incorrect input of markers")
    }
    print(genes)
    genes <- c(genes, add.features)
    print(genes)
    mat <- mat[match(genes, rownames(mat)), ]

    if (shuffle) {
        xx.rand <- sample(colnames(mat))
        saveRDS(xx.rand, file = ".xx.rand.order.rds")
        mat <- mat[, xx.rand]
    }
    if (shuffle) {
        anno <- dataset@meta.data[xx.rand, ] %>%
            rownames_to_column(var = "barcode") %>%
            arrange(!!!syms(sort_var))
    } else {
        anno <- dataset@meta.data %>%
            rownames_to_column(var = "barcode") %>%
            arrange(!!!syms(sort_var))
    }
    mat <- t(mat)
    mat <- mat[match(anno$barcode, rownames(mat)), ]
    mat <- t(mat)
    annos <- list()
    for (i in seq_along(1:length(anno_var))) {
        err_msg <- paste("Incorrect specification for annotation colors for", anno_var[i])
        value <- anno[[anno_var[i]]]
        if (is.numeric(value)) {
            if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != "qual"])) {
                n <- brewer.pal.info[anno_colors[[i]], ]["maxcolors"][[1]]
                pal <- brewer.pal(n = n, name = anno_colors[[i]])

                col_fun <- colorRamp2(
                    c(min(value), stats::median(value), max(value)),
                    c(pal[2], pal[(n + 1) / 2], pal[n - 1])
                )
            } else if (length(anno_colors[[i]]) == 3 & all(are_colors(anno_colors[[i]]))) {
                col_fun <- colorRamp2(
                    c(min(value), stats::median(value), max(value)),
                    anno_colors[[i]]
                )
            } else {
                stop(err_msg)
            }

            ha <- HeatmapAnnotation(
                a = anno[[anno_var[i]]],
                col = list(a = col_fun),
                border = TRUE,
                annotation_label = anno_var[i]
            )
        } else {
            l <- levels(factor(anno[[anno_var[i]]]))

            if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
                col <- set_colors(anno_colors[[i]], length(l))
            } else if (length(anno_colors[[i]]) >= length(l) & all(are_colors(anno_colors[[i]]))) {
                col <- anno_colors[[i]]
            } else {
                stop(err_msg)
            }

            names(col) <- l
            col <- col[!is.na(names(col))]
            col <- list(a = col)

            ha <- HeatmapAnnotation(
                a = anno[[anno_var[i]]],
                col = col,
                border = TRUE,
                annotation_label = anno_var[i]
            )
        }
        names(ha) <- anno_var[i]

        annos[[i]] <- ha
    }

    annos <- do.call(c, annos)

    annos@gap <- rep(unit(1, "mm"), length(annos))
    if (left_annotation) {
        ht <- Heatmap(mat,
            cluster_rows = cluster_rows,
            cluster_columns = FALSE,
            heatmap_legend_param = list(
                direction = "horizontal",
                legend_width = unit(6, "cm"),
                title = "Expression"
            ),
            col = colorRamp2(hm_limit, hm_colors),
            show_column_names = FALSE,
            row_names_side = row_names_side,
            row_names_gp = gpar(fontsize = row_font_size),
            left_annotation = rowAnnotation(DEGs = ctype, col = list(DEGs = DEG.color), show_legend = FALSE),
            top_annotation = annos
        )
    } else {
        ht <- Heatmap(mat,
            cluster_rows = cluster_rows,
            cluster_columns = FALSE,
            heatmap_legend_param = list(
                direction = "horizontal",
                legend_width = unit(6, "cm"),
                title = "Expression"
            ),
            col = colorRamp2(hm_limit, hm_colors),
            show_column_names = FALSE,
            row_names_side = row_names_side,
            row_names_gp = gpar(fontsize = row_font_size),
            top_annotation = annos
        )
    }

    draw(ht,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "right"
    )
}

get_top_genes <- function(dataset, markers, n) {
    int_features <- rownames(dataset@assays$RNA@scale.data)
    df <- markers %>%
        filter(.data$gene %in% int_features) %>%
        arrange(desc(.data$avg_log2FC)) %>%
        group_by(.data$cluster) %>%
        filter(row_number() <= n) %>%
        arrange(.data$cluster)
    return(df$gene)
}

are_colors <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
}
