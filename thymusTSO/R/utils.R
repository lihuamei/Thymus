#' @title medullaScore
#'
#' @description Calculates medulla scores for each spot in spatial seurat data objects using Seurat's method.
#' @param obj.st.lst A list of spatial seurat data objects.
#' @param medulla.genes A vector of genes associated with the medulla.
#' @return Updated obj.st.lst with added medulla scores.

medullaScore <- function(obj.st.lst, medulla.genes) {
    obj.st.lst <- lapply(obj.st.lst, function(obj.st) {
        obj.st <- AddModuleScore(object = obj.st, features = medulla.genes, name = "medulla.score")
    })
    return(obj.st.lst)
}

#' @title fitNorm
#'
#' @description Fit medulla score distribution and corresponding parameters.
#' @param obj.st.lst A list of spatial seurat data objects.
#' @return est.params Vector of estimated parameters from the normal distribution fitting.

fitNorm <- function(obj.st.lst) {
    x <- lapply(obj.st.lst, function(obj) {
        obj@meta.data$medulla.score1 %>% scale()
    }) %>% unlist()
    dor <- density(x, kernel = "gaussian")
    dist.max <- dor$x[which.max(dor$y)]
    dist.offset <- x - dist.max
    tmp.dist <- c(dist.offset[dist.offset <= 0], abs(dist.offset[dist.offset < 0])) + dist.max
    dist.norm.fit <- fitdistr(tmp.dist, "normal")
    est.params <- as.vector(dist.norm.fit$estimate)
    return(est.params)
}

#' @title signifTestByZtest
#'
#' @description Z-test to identify spots with significant Medulla scores.
#' @param obj Spatial seurat object data.
#' @param est.params A vector of estimated parameters of normal distribution (mean and standard deviation).
#' @return Updated Seurat object with adjusted p-values.

signifTestByZtest <- function(obj, est.params) {
    scale.score <- obj@meta.data$medulla.score1 %>% scale()
    p.raw <- pnorm(scale.score, mean = est.params[1], sd = est.params[2], lower.tail = FALSE)
    p.adj <- p.adjust(p.raw, method = "BH")
    obj@meta.data$medulla.padj <- p.adj
    obj
}

#' @title calcSpotsDist
#'
#' @description Identified significant spots located in medulla areas and calculate the distances between spots.
#' @param obj Spatial Seurat object data containing spatial and medulla score information.
#' @param p.cut Cutoff value to determine significance of spots based on adjusted p-values, default: 1e-5.
#' @param call.xgb Call trained XGBoost model or not. Default: FALSE.
#' @return A list containing Eulerian distance matrix and significant spots.

calcSpotsDist <- function(obj, p.cut = 1e-5, call.xgb = FALSE) {
    if (call.xgb) {
        xcv <- system.file("data/xcv.RDS", package = "thymusTSO") %>% readRDS(.)
        test.data <- GetAssayData(test.obj) %>%
            as.data.frame() %>%
            t() %>%
            .[, feature.genes]
        sig.spots <- predict(bst, test.data) %>%
            {
                which(. > 0.5)
            } %>%
            colnames(obj)[.]
    } else {
        sig.spots <- subset(obj, medulla.padj < p.cut) %>% Cells()
    }
    image.coord <- GetTissueCoordinates(object = obj@images[[names(obj@images)]])
    elu.dist <- dist(image.coord) %>% as.matrix()
    return(list(dist = elu.dist, sig.spots = sig.spots))
}

#' @title candEdgeSpots
#'
#' @description Determine the edge spots of each module.
#' @param elu.dist.sig Distance matrix between significant spots and entire spots.
#' @param dist.mat Distance matrix of all spots.
#' @return A character vector of candidate edge spots.

candEdgeSpots <- function(elu.dist.sig, dist.mat) {
    nonsig.cells <- setdiff(colnames(elu.dist.sig), rownames(elu.dist.sig))
    adj.sig.spots <- lapply(nonsig.cells, function(idx) {
        ff <- dist.mat[, idx]
        knn.names <- ff[order(ff)][2:7] %>% names()
        empt.spots <- sum(knn.names %in% rownames(elu.dist.sig))
        flag <- ifelse(empt.spots == 6, 1, 0)
        return(flag)
    }) %>% unlist(.)

    sig.spots <- c(rownames(elu.dist.sig), nonsig.cells[adj.sig.spots == 1])
    nonsig.cells <- setdiff(colnames(elu.dist.sig), sig.spots)
    cand.flags <- apply(elu.dist.sig, 1, function(ff) {
        knn.names <- ff[order(ff)][2:7] %>% names()
        empt.spots <- sum(knn.names %in% nonsig.cells)
        flag <- ifelse(empt.spots > 0, 1, 0)
        if (empt.spots == 6) flag <- -1
        return(flag)
    })
    cand.edges <- rownames(elu.dist.sig)[cand.flags == 1]
    return(cand.edges)
}

#' @title removeLowConfModules
#'
#' @description Identify low confidence modules of spots using traversal strategy.
#' @param elu.dist.sig Distance matrix between significant spots and other spots.
#' @param module.size Minimum module size for a confident module, default: 10.
#' @return A list containing removed and remaining spots.

removeLowConfModules <- function(elu.dist.sig, module.size = 10) {
    sig.spots <- rownames(elu.dist.sig)
    apply(elu.dist.sig, 1, function(ff) {
        flag <- rep(0, dim(elu.dist.sig)[2]) %>% `names<-`(colnames(elu.dist.sig))
        knn.names <- ff[order(ff)][2:7] %>% names()
        idx.names <- knn.names[which(knn.names %in% sig.spots)]
        if (length(idx.names) > 0) flag[idx.names] <- 1
        return(flag)
    }) %>% t() -> adj.mat

    nodes.lst <- list()
    g1 <- from_adj_matrix(adj.mat[sig.spots, sig.spots])
    for (node in sig.spots) {
        idx <- which(g1$nodes_df[, "label"] == node) %>% g1$nodes_df[., "id"]
        con.nodes <- g1 %>% get_all_connected_nodes(node = idx)
        if (!(node %in% unlist(nodes.lst))) {
            if (any(is.na(con.nodes))) {
                nodes.lst[[node]] <- node
            } else {
                nodes.lst[[node]] <- c(node, g1$nodes_df[con.nodes, "label"])
            }
        } else {
            next
        }
    }
    remove.spots <- names(nodes.lst)[lapply(nodes.lst, length) < module.size] %>%
        nodes.lst[.] %>%
        unlist(., use.names = FALSE)
    remain.spots <- names(nodes.lst)[lapply(nodes.lst, length) >= module.size] %>% nodes.lst[.]
    return(list(remove = remove.spots, remain = remain.spots))
}

#' @title updateModuleCenter
#'
#' @description Update center spot of each remained module.
#' @param obj Spatial Seurat object data containing spatial and image coordinates.
#' @param modules List of remained significant modules identified.
#' @param edge.spots List of edge spots to exclude from center calculation.
#' @return Updated modules with updated center spots.

updateModuleCenter <- function(obj, modules, edge.spots) {
    image.coord <- GetTissueCoordinates(object = obj@images[[names(obj@images)]])
    center.spots <- lapply(modules, function(sub.mod) {
        tmp.spots <- setdiff(sub.mod, edge.spots)
        if (length(tmp.spots) > 0) sub.mod <- tmp.spots
        sub.mod.coord <- image.coord[sub.mod, ]
        if (length(sub.mod) > 1) {
            sub.mod.update <- rbind.data.frame(sub.mod.coord, Center = colMeans(sub.mod.coord))
            dist(sub.mod.update) %>%
                as.matrix() %>%
                .["Center", ] %>%
                .[order(.)] %>%
                .[-which(names(.) == "Center")] %>%
                .[1] %>%
                names()
        } else {
            center.spot <- sub.mod
        }
    })
    names(modules) <- unlist(center.spots, use.names = FALSE)
    return(modules)
}

#' @title unpdateSigSpots
#'
#' @description Update significance spots based on Eulerian distance matrix.
#' @param elu.dist Eulerian distance matrix between spots.
#' @param sig.spots Current significance spots identified.
#' @return Updated vector of significance spots.

unpdateSigSpots <- function(elu.dist, sig.spots) {
    flags <- lapply(1:nrow(elu.dist), function(idx) {
        if (rownames(elu.dist)[idx] %in% sig.spots) {
            return(1)
        }
        ff <- elu.dist[idx, ]
        knn.names <- ff[order(ff)][2:7] %>% names()
        idx.names <- knn.names[which(knn.names %in% sig.spots)]
        flag <- ifelse(length(idx.names) >= 4, 1, 0)
        return(flag)
    }) %>% `names<-`(rownames(elu.dist))
    sig.spots <- flags[flags == 1] %>% names()
}

#' @title removeEdgeSpotsWithOutInternalSpots
#'
#' @description Remove low confidence edge spots based on edge distance matrix.
#' @param edge.dist Edge distance matrix.
#' @param labels A vector of labels of spots.
#' @return A vector of removed spots.

removeEdgeSpotsWithOutInternalSpots <- function(edge.dist, labels) {
    edge.dist.sub <- edge.dist[, names(labels)]
    ordered.spots <- apply(edge.dist.sub, 1, function(kk) kk[order(kk)] %>% names()) %>% t()
    bool.vals <- lapply(rownames(ordered.spots), function(sp) {
        tmp.bak <- ordered.spots[sp, ]
        spot.names <- tmp.bak[2:7] %>% labels[.]
        spot.names %in% c("Medulla_hi", "Medulla_centric") %>% sum()
    }) %>%
        unlist(.) %>%
        `names<-`(rownames(ordered.spots))
    bool.vals[bool.vals == 0] %>% names()
}

#' @title calcSpot2ModuleDist
#'
#' @description Calculate the distance between non-significant spots and modules.
#' @param obj.st.lst A list of spatial Seurat objects containing spatial data.
#' @return Updated list of spatial Seurat objects with calculated distances.

calcSpot2ModuleDist <- function(obj.st.lst) {
    lapply(obj.st.lst, function(obj) {
        image.coord <- GetTissueCoordinates(object = obj@images[[names(obj@images)]])
        elu.dist <- dist(image.coord) %>% as.matrix()
        centric.spots <- subset(obj, Centric == "Y") %>% Cells(.)
        elu.sub <- elu.dist[, centric.spots]
        assig.modules <- apply(elu.sub, 1, which.min) %>% centric.spots[.]
        assig.mindist <- apply(elu.sub, 1, min)
        meta.data <- cbind.data.frame(Centric = assig.modules, Mindist = assig.mindist)
        obj <- AddMetaData(obj, metadata = meta.data, col.name = c("Assign.Centric", "Distance"))
        return(obj)
    }) -> obj.st.lst
    return(obj.st.lst)
}

#' @title glmFitTest
#'
#' @description Inferring the gene or cell type associated with distance significance using generalized linear models.
#' @param obj.st.lst A list of spatial Seurat objects containing spatial data.
#' @param tar.genes Target genes to be tested for association with distance.
#' @param other.meta Other metadata column names to be tested for association with distance.
#' @param query.name Column name which contains the distance metric, default: "Distance".
#' @param ncores Number of cores to use for parallel processing, default: 8.
#' @return A list of results (sig.df) containing coefficients from the generalized linear models.
#' @export glmFitTest

glmFitTest <- function(obj.st.lst, tar.genes = NULL, other.meta = NULL, query.name = "Distance", ncores = 8) {
    if (is.null(tar.genes) & is.null(other.meta)) {
        stop("At least one of the tar.genes and meta.data parameters must be assigned values")
    }
    obj.st.lst <- parallel::mclapply(obj.st.lst, function(obj) {
        dat <- data.frame(distance = obj@meta.data[, query.name])
        if (!is.null(tar.genes)) {
            tar.df.1 <- FetchData(obj, vars = tar.genes, slot = "data") %>% as.data.frame()
            colnames(tar.df.1) <- tar.genes
            dat <- cbind.data.frame(tar.df.1, dat)
        }
        if (!is.null(other.meta)) {
            tar.df.2 <- obj@meta.data[, other.meta] %>% as.data.frame()
            colnames(tar.df.2) <- other.meta
            dat <- cbind.data.frame(tar.df.2, dat)
        }
        colnames(dat) <- gsub("-|/", "_", colnames(dat))
        coef.df <- sapply(colnames(dat)[1:(dim(dat)[2] - 1)], function(idx) {
            message(sprintf("%s -> %s", obj@images %>% names(), idx))
            formu <- formula(sprintf("formula = %s ~ ns(distance, 4)", idx))
            coef.val <- tryCatch(
                {
                    coef.val <- glm(formu, data = dat, family = gaussian) %>%
                        summary(.) %>%
                        .$coefficients %>%
                        .[2:5, c(1, 4)] %>%
                        t() %>%
                        as.vector()
                },
                error = function(e) {
                    coef.val <- c(NA, NA, NA, NA, NA, NA, NA, NA)
                }
            )
            if (length(coef.val) == 0) coef.val <- c(NA, NA, NA, NA, NA, NA, NA, NA)
            return(coef.val)
        }) %>%
            `names<-`(colnames(dat)[1:(dim(dat)[2] - 1)]) %>%
            t() %>%
            as.data.frame()
        colnames(coef.df) <- c("Spline.1.coef", "Spline.1.pvals", "Spline.2.coef", "Spline.2.pvals", "Spline.3.coef", "Spline.3.pvals", "Spline.4.coef", "Spline.4.pvals")
        coef.df <- cbind.data.frame(
            coef.df,
            FDR.1 = p.adjust(coef.df[, 2], method = "BH"),
            FDR.2 = p.adjust(coef.df[, 4], method = "BH"),
            FDR.3 = p.adjust(coef.df[, 6], method = "BH"),
            FDR.4 = p.adjust(coef.df[, 8], method = "BH")
        )
        obj@misc[["GLM_FIT_TEST"]] <- coef.df
        obj
    }, mc.cores = getOption("mc.cores", ncores))
    return(obj.st.lst)
}

#' @title combinePvalsByStouffer
#'
#' @description Combine P values using the Stouffer's method.
#' @param obj.st.lst A list of spatial Seurat objects containing spatial data.
#' @return Combined P values with adjusted FDR values.

combinePvalsByStouffer <- function(obj.st.lst) {
    p.lst.1 <- lapply(obj.st.lst, function(obj) {
        Spline.1.merged <- obj@misc[["GLM_FIT_TEST"]][, "Spline.1.pvals"]
    })
    p.lst.2 <- lapply(obj.st.lst, function(obj) {
        Spline.2.merged <- obj@misc[["GLM_FIT_TEST"]][, "Spline.2.pvals"]
    })
    p.lst.3 <- lapply(obj.st.lst, function(obj) {
        Spline.3.merged <- obj@misc[["GLM_FIT_TEST"]][, "Spline.3.pvals"]
    })
    p.lst.4 <- lapply(obj.st.lst, function(obj) {
        Spline.3.merged <- obj@misc[["GLM_FIT_TEST"]][, "Spline.4.pvals"]
    })
    combine.pvals.1 <- do.call(scran::combinePValues, c(p.lst.1, method = "z"))
    combine.pvals.2 <- do.call(scran::combinePValues, c(p.lst.2, method = "z"))
    combine.pvals.3 <- do.call(scran::combinePValues, c(p.lst.3, method = "z"))
    combine.pvals.4 <- do.call(scran::combinePValues, c(p.lst.4, method = "z"))
    combine.pvals <- cbind.data.frame(
        FDR.1.merged = p.adjust(combine.pvals.1, method = "BH"),
        FDR.2.merged = p.adjust(combine.pvals.2, method = "BH"),
        FDR.3.merged = p.adjust(combine.pvals.3, method = "BH"),
        FDR.4.merged = p.adjust(combine.pvals.4, method = "BH")
    )
    rownames(combine.pvals) <- rownames(obj.st.lst[[1]]@misc[["GLM_FIT_TEST"]])
    return(combine.pvals)
}

#' @title distScaleToOne
#'
#' @description Scale distances within each module to a range of [0, 1].
#' @param obj.st.lst A list of spatial Seurat objects.
#' @return A list of scaled distance matrices for each object in obj.st.lst.

distScaleToOne <- function(obj.st.lst) {
    dist.lst <- lapply(obj.st.lst, function(obj) {
        modules.centric <- obj@meta.data$Assign.Centric %>% unique()
        distance <- obj@meta.data[, "Distance", drop = FALSE]
        for (kk in modules.centric) {
            obj.sub <- subset(obj, Assign.Centric == kk)
            max.dist <- quantile(obj.sub@meta.data$Distance, 0.99)
            distance[Cells(obj.sub), "Distance"] <- obj.sub@meta.data$Distance / max.dist
        }
        return(distance)
    })
    return(dist.lst)
}
