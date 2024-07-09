tcrSubgroup <- function(obj, trab.cut = 0, trd.cut = 1, seq = T) {
    merged.obj <- subset(obj, tcr == 1)
    if (seq) {
        trb.seq <- merged.obj$tcr_cdr3s_aa %>%
            strsplit(., ";") %>%
            lapply(., function(xx) {
                xx %>%
                    grep("TRB:", .) %>%
                    xx[.] %>%
                    paste0(., collapse = ";")
            })
        merged.obj$trb <- trb.seq %>% unlist(.)
        tra.seq <- merged.obj$tcr_cdr3s_aa %>%
            strsplit(., ";") %>%
            lapply(., function(xx) {
                xx %>%
                    grep("TRA:", .) %>%
                    xx[.] %>%
                    paste0(., collapse = ";")
            })
        merged.obj$tra <- tra.seq %>% unlist(.)
    }
    merged.obj@meta.data$t.cell_group <- "others"

    exp <- merged.obj@assays$RNA@counts
    genes <- rownames(exp)
    tra <- grep(("^TRA[VJC]"), genes, value = T)
    trb <- grep(("^TRB[VJC]"), genes, value = T)
    trg <- grep(("^TRG[VJC]"), genes, value = T)
    trd <- grep(("^TRD[VJC]"), genes, value = T)

    tra_exp <- FetchData(merged.obj, vars = tra, slot = "counts")
    trb_exp <- FetchData(merged.obj, vars = trb, slot = "counts")
    trg_exp <- FetchData(merged.obj, vars = trg, slot = "counts")
    trd_exp <- FetchData(merged.obj, vars = trd, slot = "counts")

    tra_exp$sum <- apply(tra_exp, 1, sum)
    trb_exp$sum <- apply(trb_exp, 1, sum)
    trg_exp$sum <- apply(trg_exp, 1, sum)
    trd_exp$sum <- apply(trd_exp, 1, sum)

    merged.obj@meta.data$tra_sum <- tra_exp$sum
    merged.obj@meta.data$trb_sum <- trb_exp$sum
    merged.obj@meta.data$trg_sum <- trg_exp$sum
    merged.obj@meta.data$trd_sum <- trd_exp$sum


    merged.obj_1 <- subset(merged.obj, tra_sum > trab.cut & trb_sum == 0 & trd_sum == 0, slot = "counts")
    merged.obj_1_cells <- Cells(merged.obj_1)
    merged.obj@meta.data[as.vector(merged.obj_1_cells), "t.cell_group"] <- "TRA+TRB-TRD-"

    merged.obj_2 <- subset(merged.obj, tra_sum > trab.cut & trb_sum > trab.cut & trd_sum == 0, slot = "counts")
    merged.obj_2_cells <- Cells(merged.obj_2)
    merged.obj@meta.data[as.vector(merged.obj_2_cells), "t.cell_group"] <- "TRA+TRB+TRD-"

    merged.obj_3 <- subset(merged.obj, tra_sum > trab.cut & trb_sum > trab.cut & trd_sum > trd.cut, slot = "counts")
    merged.obj_3_cells <- Cells(merged.obj_3)
    merged.obj@meta.data[as.vector(merged.obj_3_cells), "t.cell_group"] <- "TRA+TRB+TRD+"

    merged.obj_4 <- subset(merged.obj, tra_sum == 0 & trb_sum > trab.cut & trd_sum > trd.cut, slot = "counts")
    merged.obj_4_cells <- Cells(merged.obj_4)
    merged.obj@meta.data[as.vector(merged.obj_4_cells), "t.cell_group"] <- "TRA-TRB+TRD+"

    merged.obj_5 <- tryCatch(
        {
            merged.obj_5 <- subset(merged.obj, tra_sum == 0 & trb_sum == 0 & trd_sum > trd.cut, slot = "counts")
        },
        error = function(e) {
            return(NULL)
        }
    )
    if (!is.null(merged.obj_5)) {
        merged.obj_5_cells <- Cells(merged.obj_5)
        merged.obj@meta.data[as.vector(merged.obj_5_cells), "t.cell_group"] <- "TRA-TRB-TRD+"
    }
    merged.obj_6 <- tryCatch(
        {
            merged.obj_6 <- subset(merged.obj, tra_sum == 0 & trb_sum == 0 & trd_sum == 0, slot = "counts")
        },
        error = function(e) {
            return(NULL)
        }
    )
    if (!is.null(merged.obj_6)) {
        merged.obj_6_cells <- Cells(merged.obj_6)
        merged.obj@meta.data[as.vector(merged.obj_6_cells), "t.cell_group"] <- "TRA-TRB-TRD-"
    }

    merged.obj_7 <- tryCatch(
        {
            merged.obj_7 <- subset(merged.obj, tra_sum == 0 & trb_sum > trab.cut & trd_sum == 0, slot = "counts")
        },
        error = function(e) {
            return(NULL)
        }
    )
    if (!is.null(merged.obj_7)) {
        merged.obj_7_cells <- Cells(merged.obj_7)
        merged.obj@meta.data[as.vector(merged.obj_7_cells), "t.cell_group"] <- "TRA-TRB+TRD-"
    }

    merged.obj_8 <- tryCatch(
        {
            merged.obj_8 <- subset(merged.obj, tra_sum > trab.cut & trb_sum == 0 & trd_sum > trd.cut, slot = "counts")
        },
        error = function(e) {
            return(NULL)
        }
    )
    if (!is.null(merged.obj_8)) {
        merged.obj_8_cells <- Cells(merged.obj_8)
        merged.obj@meta.data[as.vector(merged.obj_8_cells), "t.cell_group"] <- "TRA+TRB-TRD+"
    }
    return(merged.obj)
}
