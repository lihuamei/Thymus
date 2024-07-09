suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(caTools))
suppressMessages(library(ROGUE))
suppressMessages(library(colorRamps))
suppressMessages(library(tidyverse))

#--------------------------------------------------------------
# Load own modules

source("modules/utils.R")
source("modules/global_params.R")
source("modules/seurat_methods.R")
source("modules/visualization.R")
source("modules/trajectory_methods.R")
source("modules/vis/step_03_anno_assessment.vis.R")

#--------------------------------------------------------------
# Create save directories

out.data.dir <- file.path("../3.results", "3.anno_assess/data")
out.figs.dir <- file.path("../3.results", "3.anno_assess/figs")

dir.create(file.path("../3.results", "3.anno_assess"))
dir.create(out.data.dir, showWarnings = FALSE)
dir.create(out.figs.dir, showWarnings = FALSE)

#--------------------------------------------------------------
# Load data and Gene stats

cols <- ANNO_ENTIRE_COLOR_FIG1 %>% `names<-`(ANNO_ENTIRE_IDNET_FIG1)
merged.obj <- readRDS("../3.results/2.fine_anno/data/merged.obj.anno.rds")
merged.obj@meta.data$Anno.Level.Fig.1 <- factor(merged.obj@meta.data$Anno.Level.Fig.1, levels = ANNO_ENTIRE_IDNET_FIG1)
Idents(merged.obj) <- merged.obj@meta.data$Anno.Level.Fig.1

openxlsx::write.xlsx(table(merged.obj$Anno.Level.Fig.1, merged.obj$orig.ident), file = "./revise.ms/reviewer.3/data/stat.cells.per.sample.xlsx", rowNames = TRUE)

xx <- table(merged.obj$Anno.Level.Fig.1, merged.obj$orig.ident)
xx[xx < 10] <- 0
xx[xx > 0] <- 1

my_palette <- alpha(colorRampPalette(c("skyblue", "red"))(25), 0.5)
pheatmap::pheatmap(xx, cluster_cols = FALSE, cluster_rows = TRUE, color = my_palette, border_color = "white", filename = "./revise.ms/reviewer.3/figs/cell.cnt.per.sample.pdf", width = 6, height = 6)

plot.df.1 <- cbind(G = rowSums(xx) / 16, L = 1 - rowSums(xx) / 16, CellTyep = rownames(xx), Group = "<=10") %>% as.data.frame()

xx <- table(merged.obj$Anno.Level.Fig.1, merged.obj$orig.ident)
xx[xx < 20] <- 0
xx[xx > 0] <- 1
plot.df.2 <- cbind(G = rowSums(xx) / 16, L = 1 - rowSums(xx) / 16, CellTyep = rownames(xx), Group = "<=20") %>% as.data.frame()


xx <- table(merged.obj$Anno.Level.Fig.1, merged.obj$orig.ident)
xx[xx < 30] <- 0
xx[xx > 0] <- 1
plot.df.3 <- cbind(G = rowSums(xx) / 16, L = 1 - rowSums(xx) / 16, CellTyep = rownames(xx), Group = "<=30") %>% as.data.frame()

# plot.df <- rbind.data.frame(plot.df.1, plot.df.2, plot.df.3)

order.cells <- rownames(plot.df.1)[order(as.numeric(plot.df.1[, 1]))]
plot.df <- tidyr::gather(plot.df.1, "Stat", "Frac", -CellTyep, -Group)
plot.df$CellTyep <- factor(plot.df$CellTyep, levels = order.cells)
plot.df$Stat <- factor(plot.df$Stat, levels = c("G", "L"))

ggplot(plot.df, aes(x = CellTyep, y = Frac, fill = Stat)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_line(linetype = "dotted")) +
    xlab("") +
    ylab("Proportion") +
    scale_fill_manual(values = c("G" = alpha("red", 0.5), "L" = "grey"))
