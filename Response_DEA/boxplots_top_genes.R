shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))
shhh(library(ggplot2))
shhh(library(tidyr))
shhh(library(reshape2))

#Paths 
if (getwd()== "/home/mchacon"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
  
}else if(getwd()== "/home/mchacon/Desktop/cluster/Projects/scRNAseq/mchacon/"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
}else{
  basepath <- "/gpfs/projects/bsc83"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq"
}
chem_version <- "V2"
plots_path <- paste0(data_path, "/mchacon/1M-scBloodNL/plots/01_DEA_Age_Sex/", chem_version)
robjects <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects")
enrich_dir<- paste0(plots_path, "/Enrichments/")

#functions
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Functions.R"))
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Themes.R"))

cell_level <- "L2"; pheno <- "age"; cell <- "naive_CD8T"; combo <- "PA_24h.PA_3h"; stim <- "PA"; d <- "up"
df_all <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/nCells_per_donor/", 
                           chem_version, "/", cell_level, "/age_residuals.lmer_lm_ratio_log2cpm_dframe.rds"))

df_cell <- df_all[df_all$combo == combo & df_all$celltype == cell & df_all$adj.P.Val < 0.05, ]
df_cell <- df_cell %>% mutate(direction = ifelse(estimate > 0, "up", "down"))
df_cell <- df_cell[df_cell$direction == "up", ]

top_genes <- df_cell %>% mutate(abs_logFC = abs(logFC)) %>% slice_max(order_by = abs_logFC, n = 10) 

top_genes <- top_genes %>% mutate(direction = ifelse(estimate > 0, "up", "down"))

genes <- top_genes$ID
top_genes <- df_cell[df_cell$ID %in% genes, ]

rc1 <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/pseudobulk_dreamlet/sum/", chem_version, "/",
                       cell_level, "/", cell, "/Gender_age/min_prop.0.2/", stim, "_24h_pb.Rds"))

rc2 <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/pseudobulk_dreamlet/sum/", chem_version, "/",
                        cell_level, "/", cell, "/Gender_age/min_prop.0.2/",stim, "_3h_pb.Rds"))
rc3 <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/pseudobulk_dreamlet/sum/", chem_version, "/",
                      cell_level, "/", cell, "/Gender_age/min_prop.0.2/UT_0h_pb.Rds"))
##########################
expr1 <- assay(rc1, cell)[genes, ]
expr2 <- assay(rc2, cell)[genes, ]
expr_ut <- assay(rc3, cell)[genes, ]

# Find common donors across UT, 3h, and 24h
common_donors <- Reduce(intersect, list(colnames(expr1), colnames(expr2), colnames(expr_ut)))

#common_donors <- intersect(colnames(expr1), colnames(expr2))

common1 <- expr1[, common_donors]
common2 <- expr2[, common_donors]
common_ut <- expr_ut[, common_donors]

metadata1 <- colData(rc1)[common_donors, ]
metadata2 <- colData(rc2)[common_donors, ]  # both metadata should be the same
metadata_ut <- colData(rc3)[common_donors, ]  # Metadata for UT (same donors)

identical(rownames(common1), rownames(common2))
common1 <- as.data.frame(as.matrix(t(common1)))
common2 <- as.data.frame(as.matrix(t(common2)))
common_ut <- as.data.frame(as.matrix(t(common_ut)))

common1$Donor <- rownames(common1)
common2$Donor <- rownames(common2)
common_ut$Donor <- rownames(common_ut)
common1$Age <- metadata1$age
common2$Age <- metadata2$age # same metadata in all three
common_ut$Age <- metadata_ut$age

# Transformar a formato largo usando `melt`
expr_1_long <- melt(common1, id.vars = c("Donor", "Age"), variable.name = "ID", value.name = "Expression")
expr_1_long$Time <- "24h"

expr_2_long <- melt(common2, id.vars = c("Donor", "Age"), variable.name = "ID", value.name = "Expression")
expr_2_long$Time <- "3h"

expr_ut_long <- melt(common_ut, id.vars = c("Donor", "Age"), variable.name = "ID", value.name = "Expression")
expr_ut_long$Time <- "0h"  # Add unstimulated condition

# Combinar en un solo dataframe
expr_long <- rbind(expr_1_long, expr_2_long)
expr_long <- rbind(expr_1_long, expr_2_long, expr_ut_long)

expr_long$AgeGroup <- cut(expr_long$Age, breaks = c(-Inf, 40, 60, Inf), 
                            labels = c("Y", "40-60", "O"))

f_expr <- subset(expr_long, AgeGroup != "40-60")

#f_expr$Time <- factor(f_expr$Time, levels = c("3h", "24h"))
f_expr$Time <- factor(f_expr$Time, levels = c("0h", "3h", "24h"))

shhh(library(ggstatsplot))
library(ggplot2)
f_expr <- merge(f_expr, top_genes[, c("ID", "direction")], by.x = "ID", by.y = "ID", all.x = TRUE)
f_expr$Expression <- f_expr$Expression + 1
f_expr$Expression <- log2(f_expr$Expression)

p <- ggplot(f_expr, aes(x = Time, y = Expression, fill = AgeGroup, alpha = AgeGroup)) + geom_violin(width= 0.9)+  
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) + scale_alpha_manual(values= c(0.3, 0.8, 1)) +
  #scale_fill_manual(values = c("up" = "#bc4749", "down" = "#264653")) +
  scale_fill_manual(values = c("Y" = "#B5B8B7", "O" = "#555F6A")) +
  
  facet_grid(ID ~ AgeGroup, scales = "free_y") +  # Faceta por gen y tiempo
  labs(title = "Gene Expression at different times for young and old individuals",
       x = "Time", y = "log(Gene Expression)") +
  theme_paper + theme(
    axis.text.x = element_text(hjust = 1, size = 11),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 14),
    strip.background = element_rect(color = "black", fill = "white", linewidth =0.8)
    
    )
p
pdf(paste0(data_path, "/mchacon/1M-scBloodNL/Figures/", chem_version, 
           "/", "DEGS_FIG_supp_PA_", pheno, "_", cell_level, "_Plot.pdf"), width =10, height = 8)
p
dev.off()
