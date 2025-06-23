shhh <- suppressPackageStartupMessages
shhh(library(dplyr)); shhh(library(ggplot2)); shhh(library(tidyr)); shhh(library(ggpubr)); shhh(library(scater))
shhh(library(gridExtra)); shhh(library(SingleCellExperiment)); shhh(library(SeuratObject)); shhh(library(Seurat))
shhh(library(dreamlet)); shhh(library(MAST)); shhh(library(zenith)); shhh(library(Matrix)); shhh(library(lme4)); shhh(library(AnnotationDbi))
getwd()
#Paths 
if (getwd()== "/home/mchacon"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
  
}else if(getwd()== "/home/mchacon/Desktop/cluster/Projects/scRNAseq/mchacon/1M-scBloodNL/scripts/plots/log_ratio"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
}else{
  basepath <- "/gpfs/projects/bsc83"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq"
}

chem_version <- "V2"
plots_path <- paste0(data_path, "/mchacon/1M-scBloodNL/plots/02_log_ratio/", chem_version, "/", cell_level)
robjects <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects")
enrich_dir<- paste0(plots_path, "/Enrichments/")
figure_path <- paste0(data_path, "/mchacon/1M-scBloodNL/Figures/", chem_version, "/")

#functions
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Functions.R"))
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Themes.R"))

##################### UP/DOWN GENES ###########################
chem_version <- "V2"; cell_level <- "L2"; pheno <- "age"; method_1 <- "residuals.lmer_lm"
get_signif <- function(cell_level, chem_version, pheno, method){
  #read datasets for each pheno and cell_level 
  deg <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/nCells_per_donor/", 
                       chem_version, "/", cell_level, "/age__meta_", method, "_ratio_log2cpm_dframe.rds"))
  
  #modify some names for plotting
  deg$celltype <- gsub("_", " ", deg$celltype)
  deg$direction <- ifelse(deg$logFC > 0, "up", "down")
  #Select significant genes 
  deg_signif <- deg[deg$adj.P.Val < 0.05, ]
  return(deg_signif)
}
df <- get_signif("L2", "V2", "age", method_1)
df <- reorder_celltypes(df, cell_level = "L2", reverse = T)

df_final <- df %>%
  dplyr::group_by(celltype, stim, combo, direction) %>% dplyr::summarise(n = n()) %>% mutate(n = ifelse(direction == "down", -n, n)) %>%
  dplyr::mutate(vjust_label = ifelse(direction == "down", 1.5, -0.5))  # Adjust vjust based on direction

df_final$combo <- factor(df_final$combo, levels = c("CA_3h.UT_0h", "CA_24h.UT_0h", "CA_24h.CA_3h", 
                                          "PA_3h.UT_0h", "PA_24h.UT_0h", "PA_24h.PA_3h", 
                                          "MTB_3h.UT_0h", "MTB_24h.UT_0h", "MTB_24h.MTB_3h"))

df_final$stim <- factor(df_final$stim, levels = c("CA", "PA", "MTB"))

if(pheno == "Gender") {pheno <- "Sex"}
plot_cond <- ggplot(df_final, aes(y = combo, x = n, fill = combo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.8) + 
  geom_text(aes(label = abs(n)),
            position = position_dodge(width = 0.7), 
            hjust = ifelse(df_final$n > 0, -0.2, 1.2),
            size = 4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
  facet_wrap(~ celltype, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = pal_combo) +  
  xlab("Number of significant DEGs") + ylab(NULL) +
  ggtitle(paste0("Number of significant DEGs influenced by ", pheno, " (log ratio, ", method_1, ")")) + 
  theme_paper +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y.right = element_text(angle = 0, hjust = 0, size = 17),
        axis.text.x = element_text(size = 12),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

print(plot_cond)

pdf(paste0(figure_path, "num_degs_age_V2_", method_1, "_.pdf"), width = 10, height = 5.5)
plot_cond 
dev.off()

