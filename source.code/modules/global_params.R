#----------------------------------------------------------------------
# Reference expression and markers

SCIENCE_EXPR <- "../4.extdata/science_thymus_gene_expression.xls"
SCIENCE_MARKERS <- "../4.extdata/science_thymus_markers.csv"

ANNO.AREAS <- c("protein_coding", "TEC", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene")

SIG.REF <- read.table(SCIENCE_EXPR, header = T, sep = "\t") %>%
    .[!duplicated(.[, 1]), ] %>%
    `rownames<-`(.$index) %>%
    .[, -1] %>%
    .[read.table(SCIENCE_MARKERS, sep = ",", header = T) %>%
        as.matrix() %>%
        as.vector() %>%
        unique(), ]

REF.MARKERS <- read.table(SCIENCE_MARKERS, sep = ",", header = T)

library(RColorBrewer)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
COLOR.VECS <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) %>% unique(.)

order.cell.science <- c(
    "DN", "DP", "αβT(entry)", "CD8αα", "γδT", "Treg",
    "T(agonist)", "CD4+T", "CD8+T", "CD4+Tmem", "CD8+Tmem", "NK",
    "ILC3", "NKT", "ETP", "NMP", "Mono", "Mac", "DC1", "DC2",
    "aDC", "pDC", "B_pro/pre", "B_naive", "B_memory",
    "B_plasma", "Mast", "Ery", "Mgk", "cTEC", "mcTEC",
    "mTEC(I)", "mTEC(II)", "mTEC(III)", "mTEC(IV)", "TEC(myo)",
    "TEC(neuro)", "Fb_1", "Fb_2", "Fb_cycling", "VSMC", "Endo", "Lymph"
)

#---------------------------------------------------------------------
# Annotation

ANNO_SC_RAW_INDENT <- rev(c("DN", "DP1", "DP2", "DP3", "DP4", "abT(entry)1", "abT(entry)2", "CD8+ T", "CD8aa", "CD4+ T", "Treg", "NKT", "Stromal", "Myeloid", "B", "Plasma", "Ery"))
ANNO_SC_RAW_COLOR <- rev(c("#ffccff", "grey50", "grey65", "grey75", "grey80", "#800000", "#cc0000", "dark orange", "red", "#00ffff", "purple", "hotpink", "yellow", "#40ff00", "blue", "#00008B", "#A0522D"))

ANNO_SC_IDENT_LEVEL_1 <- rev(c("T", "Stromal", "Myeloid", "B", "Plasma", "Ery"))
ANNO_SC_COLOR_LEVEL_1 <- rev(c("grey", "yellow", "#40ff00", "blue", "cornflowerblue", "#A0522D"))

ANNO_SC_MYELOID_IDENT <- c("Mac1", "Mac2", "Mono1", "Mono2", "DC", "pDC")
ANNO_SC_STROMAL_IDENT <- c("Fb1", "Fb2", "Fb3", "Fb4", "Fb_cycling", "Endo", "Tuft", "VSMC", "Lymph", "TEC")
ANNO_SC_TEC_IDENT <- c("cTEC1", "cTEC2", "cTEC3", "cTEC4", "mTEC1", "mTEC2", "mTEC3")

ANNO_SC_DN_IDENT <- c("DN_early", "DN_trans", "DN_blast1", "DN_blast2", "DN_re1", "DN_re2", "DN_re3", "DN_ISP", "DN_mito")
ANNO_SC_DP_P_IDENT <- c("DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5", "DP_blast6", "DP_blast7", "DP_blast8")

ANNO_SC_DP_R_IDENT <- c("DP_re1", "DP_re2", "DP_re3", "DP_re4", "DP_re5", "DP_re6")
ANNO_SC_OTHER_T_INDET <- c("CD4T", "CD4T_mem", "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa.I", "CD8aa.II", "T_agonist", "T_apoptosis", "T_proliferating", "NKT", "ILC3")

ANNO_ENTIRE_IDNET <- c("Ery", "B_naive", "B_trans", "B_memory", "Plasma", "Mono1", "Mono2", "Mac1", "Mac2", "DC", "pDC", "Fb1", "Fb2", "Fb3", "Fb4", "Fb_cycling", "Endo", "Tuft", "VSMC", "Lymph", "cTEC1", "cTEC2", "cTEC3", "cTEC4", "mTEC1", "mTEC2", "mTEC3", "DN_early", "DN_trans", "DN_blast1", "DN_blast2", "DN_re1", "DN_re2", "DN_re3", "DN_ISP", "DN_mito", "DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5", "DP_blast6", "DP_blast7", "DP_blast8", "DP_re1", "DP_re2", "DP_re3", "DP_re4", "DP_re5", "DP_re6", "abT(entry)1", "abT(entry)2", "CD4T", "CD4T_mem", "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa.I", "CD8aa.II", "T_agonist", "T_apoptosis", "T_proliferating", "NKT", "ILC3")

ANNO_ENTIRE_IDNET_FIG1 <- c("Ery", "B_naive", "B_trans", "B_memory", "Plasma", "Mono", "Mac", "DC", "pDC", "Fb", "Fb_cycling", "Endo", "VSMC", "Lymph", "cTEC", "mTEC", "DN_early", "DN_blast", "DN_re", "DP_blast", "DP_re", "abT(entry)", "CD4T", "CD4T_mem", "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa", "T_agonist", "T_apoptosis", "T_proliferating", "NKT", "ILC3")

SAMPLE_COLOR_FIG1 <- c(
    "#CBE6D4", "#D8A1BF", "#9BC6DA", "#FEECF0", "#C5B4A1", "#EBEED1", "#D2ECF1", "#CBAAAA", "#9FD1B6",
    "#B2BFDD", "#F1B3CC", "#C81864", "#B4B0B9", "#BBBDBA", "#C1C9AC", "#99B8B3", "#CEDAEC",
    "#CDC1DD", "#A5CADF", "#D2C7CC", "#FACCAC", "#B8999F", "#ACD9B6", "#E0E3AE", "#D8E5F2",
    "#E3DAEA", "#CCB9B7", "#A6DDE4", "#E2C9B5", "#EEABAC", "#EAC9DE", "#E4D3D1", "#A4B2A5", "#D2E9D1", "#FAD3D3"
)

ANNO_ENTIRE_COLOR_FIG1 <- COLOR.VECS[1:length(ANNO_ENTIRE_IDNET_FIG1)]
names(ANNO_ENTIRE_COLOR_FIG1) <- ANNO_ENTIRE_IDNET_FIG1

T.CELL.GROUPS <- c("TRA+TRB+TRD-", "TRA-TRB+TRD-", "TRA+TRB+TRD+", "TRA-TRB+TRD+", "TRA+TRB-TRD-", "TRA+TRB-TRD+", "TRA-TRB-TRD-", "TRA-TRB-TRD+")
T.CELL.GROUPS.COLORS <- c(`TRA-TRB+TRD+` = "#66C2A5", `TRA+TRB+TRD+` = "#FC8D62", `TRA-TRB+TRD-` = "#8DA0CB", `TRA+TRB+TRD-` = "#E78AC3")
T.CELL.GROUPS.EIGHT.COLORS <- c(T.CELL.GROUPS.COLORS[T.CELL.GROUPS[1:4]], COLOR.VECS[60:63]) %>% `names<-`(T.CELL.GROUPS)

notch.genes <- cogena::gmt2list("./configs/NORTCH.gmt")
DIFFERENTIAL_STATE_MARKERS <- list(
    NOTCH = c("TCF4", "NOTCH1", "HES4", "PTCRA", "DTX1"),
    Proliferation = c("TYMS", "TCF19", "ESCO2", "TOP2A", "HMMR"),
    POS = c("CD1D", "FOXN1", "H2-DMA", "ITPKB", "PTPN2", "SKINT1L", "SRF", "STK11", "TOX", "RAG1", "RAG2", "CD8A", "CD8B", "CD4"),
    NEG = c("AIRE", "ATG5", "CCR7", "CD28", "FAS", "GLI3", "MINK1", "SPN", "ZAP70", "CD3E", "PTPRC"),
    Surface = c("CD4", "CD8A", "CD8B", "CD3D", "CD3E", "CD3G")
)

#----------------------------------------------------------------------
# SampleInfos

SampleInfos <- list(
    Thy1 = "13w",
    Thy2 = "17w",
    Thy3 = "21w",
    Thy4 = "23w",
    Thy5 = "3m",
    Thy6 = "4m",
    Thy7 = "5m",
    Thy8 = "9m",
    Thy9 = "10m",
    Thy10 = "1y",
    Thy11 = "2y",
    Thy12 = "2y",
    Thy13 = "29y",
    Thy14 = "39y",
    Thy15 = "60y",
    Thy16 = "67y"
)

SampleClassify <- list(
    Prental = c("Thy1", "Thy2", "Thy3", "Thy4"),
    Children = c("Thy5", "Thy6", "Thy7", "Thy8", "Thy9", "Thy10", "Thy11", "Thy12"),
    Adult = c("Thy13", "Thy14", "Thy15", "Thy16")
)

SC_MAP_SP <- list(
    Thy5 = "A6",
    Thy6 = "A2",
    Thy7 = "A7",
    Thy8 = "A1",
    Thy9 = "A8",
    Thy10 = "A3",
    Thy11 = "A4",
    Thy12 = "A5"
)

SampleClassifyColors <- c(Prental = "#83534C", Children = "#CC78B1", Adult = "#217AB6")

HEATMAP_COLORS <- colorRampPalette(c("grey80", "white", "red"))(100)

# RE-Cluster T cells

ANNO_TCR_IDENT <- c(
    "DN_early", "DN_blast", "DN_re", "ISP", "DP_blast1", "DP_blast2", "DP_re1", "DP_re2", "DP_re3", "abT(entry)", "CD4T", "CD4T_mem",
    "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa", "T_agonist", "T_apoptosis"
)

ANNO_TCR_IDENT_COLOR <- c(COLOR.VECS[64:69] %>% `names<-`(setdiff(ANNO_TCR_IDENT, names(ANNO_ENTIRE_COLOR_FIG1))), ANNO_ENTIRE_COLOR_FIG1[intersect(ANNO_TCR_IDENT, names(ANNO_ENTIRE_COLOR_FIG1))])[ANNO_TCR_IDENT]

ANNO_TREND_IDENT <- c("DN_early", "DN_trans", "DN_blast", "DN_re", "ISP", "DP_blast1", "DP_blast2", "DP_blast3", "DP_blast4", "DP_blast5", "DP_re1", "DP_re2", "DP_re3", "DP_re4", "abT(entry)", "CD4T", "CD4T_mem", "Treg.diff", "Treg", "CD8T", "CD8T_mem", "CD8aa", "T_agonist", "T_apoptosis", "T_proliferating", "NKT")
ANNO_TREND_IDENT_COLOR <- c(COLOR.VECS[(70 - length(ANNO_TREND_IDENT) + 1):70] %>% `names<-`(setdiff(ANNO_TREND_IDENT, names(ANNO_ENTIRE_COLOR_FIG1))), ANNO_ENTIRE_COLOR_FIG1[intersect(ANNO_TREND_IDENT, names(ANNO_ENTIRE_COLOR_FIG1))])[ANNO_TREND_IDENT]

ANNO_TREND_IDENT_COLOR["DN_re1"] <- "#EED2EE"
ANNO_TREND_IDENT_COLOR["DP_blast3"] <- "#D15FEE"
ANNO_TREND_IDENT_COLOR["DP_blast2"] <- "#CDB7B5"

FEATURES_COLOR <- c("gray90", "#4077B3", "#419EAB", "#6DC6A4", "#E5ECAD", "#FAF9C3", "#EE793E", "#F1783D", "#96234B")

#-----------------------------------------------------------------------
# Spatial

SP_NAMES <- c(A1 = "T1", A2 = "T7", A3 = "T5", A4 = "T6", A5 = "T3", A6 = "T2", A7 = "T8", A8 = "T4")
sp2SampleName <- function(T.names) {
    tmp.rev <- names(SP_NAMES) %>% `names<-`(SP_NAMES)
    tmp.rev.sub <- tmp.rev[T.names]
    cc <- names(SC_MAP_SP) %>% `names<-`(SC_MAP_SP)
    cc <- cc[tmp.rev.sub] %>% `names<-`(names(tmp.rev.sub))
    cc
}

medulla.genes <- list(c("EBI3", "CCL17", "CCR7", "CSF2RB", "CCL21", "CCL22", "TNFRSF18", "CCL27", "CXCL10", "CXCL9", "MS4A1", "LAMP3"))
CM_COLORS <- c(Cortex = "#00B6EB", Junction = "#A68CFF", Medulla = "#FB61D7")
SAMPLE_COLORS <- c(T1 = "blue1", T2 = "orange", T3 = "darkgreen", T4 = "navy", T5 = "brown1", T6 = "darkmagenta", T7 = "grey", T8 = "darkkhaki")
