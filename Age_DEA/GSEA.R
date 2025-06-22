shhh <- suppressPackageStartupMessages
shhh(library(optparse)); shhh(library(dplyr)); shhh(library(ggplot2)); shhh(library(tidyr)); shhh(library(clusterProfiler)); shhh(library(org.Hs.eg.db)) 
shhh(library(AnnotationDbi)); shhh(library(RColorBrewer)) ; shhh(library(rrvgo))
print(getwd())

library(dplyr)
library(tibble)
#Paths 
if (getwd()== "/home/mchacon"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
  
}else if(getwd()== "/home/mchacon/Desktop/cluster/Projects/scRNAseq/mchacon/1M-scBloodNL/scripts/plots/DEA"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
}else{
  basepath <- "/gpfs/projects/bsc83"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq"
}
chem_version <- "V2"

plots_path <- paste0(data_path, "/mchacon/1M-scBloodNL/plots/01_DEA_Age_Sex/", chem_version, "/Enrichments/")
robjects <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects")
enrich_dir<- paste0(plots_path, "/Enrichments/")

#functions
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Functions.R"))
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Themes.R"))

#main_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq/"
#main_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"


get_genes <- function(pheno, cell_level, chem_version){
  #read datasets for each pheno and cell_level 
  df_all <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/",
                           pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
  
  #modify some names for plotting
  #df_all <- df_all[df_all$fdr.log2cpm < 0.05, ]
  df_all$direction <- ifelse(df_all$estimate.log2cpm > 0, "up", "down")
  #Select significant genes 
  return(df_all)
}
genes <- get_genes("age", "L2", "V2")
timepoints<- c("UT0h", "CA3h", "CA24h", "PA3h", "PA24h", "MTB3h", "MTB24h")
for (tp in timepoints) {
  df <- genes[genes$timepoint == tp, ]
  
  # Step 1: Prepare ranked gene list
  # ---------------------------------
  gene_ranking <- df %>%
    filter(!is.na(estimate.log2cpm)) %>%
    arrange(desc(estimate.log2cpm)) %>%
    dplyr::select(ID, estimate.log2cpm) %>%
    deframe()
  
  # Step 2: Convert gene symbols to Entrez IDs
  # ------------------------------------------
  gene_df <- bitr(names(gene_ranking), fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Retain ranking values and attach Entrez IDs
  gene_ranking_named <- gene_ranking[gene_df$SYMBOL]
  names(gene_ranking_named) <- gene_df$ENTREZID
  
  # Step 3: Run GSEA using GO (can change to KEGG or MSigDB if needed)
  # -------------------------------------------------------------------
  gsea_result <- gseGO(geneList = gene_ranking_named,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",            
                       keyType = "ENTREZID",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = FALSE)
  
  # Step 4: Visualize Results
  # --------------------------
  library(enrichplot)
  dotplot <- dotplot(gsea_result, showCategory = 10) + ggtitle(tp) + theme_paper
  
  pdf(paste0(plots_path, "GSEA_age_L2_naive_CD8T_", tp, "_BP.pdf"), width = 9, height = 6)
  print(dotplot)
  dev.off() 
  
}
