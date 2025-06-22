shhh <- suppressPackageStartupMessages
shhh(library(optparse)); shhh(library(dplyr)); shhh(library(ggplot2)); shhh(library(tidyr)); shhh(library(clusterProfiler)); shhh(library(org.Hs.eg.db)) 
shhh(library(AnnotationDbi)); shhh(library(RColorBrewer)) ; shhh(library(rrvgo))
print(getwd())
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
cell_level <- "L2"
plots_path <- paste0(data_path, "/mchacon/1M-scBloodNL/plots/01_log2cpm_DEA_Age_Sex/", chem_version, "/", cell_level)
robjects <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects")
enrich_dir<- paste0(plots_path, "/Enrichments/")

#functions
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Functions.R"))
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Themes.R"))

get_degs <- function(pheno, cell_level, chem_version){
  #read datasets for each pheno and cell_level 
  df_all <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/",
                           pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
  
  #modify some names for plotting
  df_all <- df_all[df_all$fdr.log2cpm < 0.05, ]
  df_all$direction <- ifelse(df_all$estimate.log2cpm > 0, "up", "down")
  #Select significant genes 
  return(df_all)
}

GO_enrichment <- function(pheno, cell_level, cell_type, tp, ont){
  df <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/",
                       pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
  df_signif <- get_degs(pheno, cell_level, chem_version)
  
  if (tp == "shared") {
    df <- df[df$celltype == cell_type, ]
    
    tested_genes <- split(df$ID, df$timepoint)
    common_genes <- Reduce(intersect, tested_genes) # for shared genes between all conditions 
    
    #### DEGs
    df_signif <- df_signif[df_signif$celltype == cell_type,]
    signif_genes_up <- df_signif[df_signif$direction == "up",]
    signif_genes_down <- df_signif[df_signif$direction == "down",]
    
    tested_g_up <- split(signif_genes_up$ID, signif_genes_up$timepoint)
    genes_up <- Reduce(intersect, tested_g_up)
    tested_g_down <- split(signif_genes_down$ID, signif_genes_down$timepoint)
    genes_down <- Reduce(intersect, tested_g_down)
    
  } else {
    genes_signif <- df_signif[df_signif$celltype == cell_type & df_signif$timepoint == tp,]$ID
    signif_genes_up <- df_signif[df_signif$celltype == cell_type & df_signif$timepoint == tp & df_signif$direction == "up",]$ID
    signif_genes_down <- df_signif[df_signif$celltype == cell_type & df_signif$timepoint == tp & df_signif$direction == "down",]$ID
    
    common_genes <- df[df$celltype == cell_type & df$timepoint == tp,]$ID
  }
  
  if(length(signif_genes_up) > 1){
    print("More than 1 DE gene in up")
    enrich_up <- enrichGO(signif_genes_up, universe = common_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont= ont)
    enrich_up@result
    if (!is.null(enrich_up)) {
      if (length(which(enrich_up@result$p.adjust < 0.05)) > 0) {
              up_plot <-  dotplot(enrich_up, showCategory = 10, color = "p.adjust")+ggtitle(paste0( "GO Terms Upregulated for ", pheno, ", ", cell_type, "-", tp, " ", ont))+theme()
        #pdf(paste0(enrich_dir, "GO_", pheno, "_", cell_level, "_", cell_type, "_", tp, "_", ont,"_up.pdf"), width = 7, height = 6)
        print(up_plot)
        #dev.off()   
        saveRDS(enrich_up, paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, "/", cell_level, "/", 
        cell_type, "/Enrichments/", pheno, "_", tp, "_up_", ont,  "_GO_enrichments.rds"))
        print(".rds saved")
        }
    } else {
        print(paste0("No valid enrichment for: ", cell_type,"-up-", tp))
      }
  }
  
  if(length(signif_genes_down) > 1){
    print("More than 1 DE gene down")
    
    enrich_down <- enrichGO(signif_genes_down, universe = common_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont= ont)
    enrich_down@result
    
    if (!is.null(enrich_down)) {
      if (length(which(enrich_down@result$p.adjust < 0.05)) > 0) {            
          down_plot <-  dotplot(enrich_down, showCategory = 20, color = "p.adjust")+ggtitle(paste0( "GO Terms Downregulated for ", pheno, ", ",  cell_type, "-", tp, " ", ont))+theme()
          #pdf(paste0(enrich_dir, "GO_", pheno, "_", cell_level, "_", cell_type, "_", tp, "_", ont,"_down.pdf"), width = 7, height = 6)
          print(down_plot)
          #dev.off()   
          saveRDS(enrich_down, paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, "/", cell_level, "/",
                                 cell_type, "/Enrichments/", pheno, "_", tp, "_down_", ont,  "_GO_enrichments.rds"))
          print(".rds saved")
      }
        
    } else {
      print(paste0("No valid enrichment for: ", cell_type,"-down-", tp))
    }
  }
}


cell_types <- c("naive_CD8T")
timepoints <- c('CA3h', 'CA24h', 'PA3h', 'PA24h', 'MTB3h', 'MTB24h', 'UT0h') # or change it by "shared"
for (cell_type in cell_types) {
  for (tp in timepoints) {
    GO_enrichment("age", "L2", cell_type, tp, ont)
  }
}

pheno <- "age"; cell_level <- "L2"; ont <- "BP"; cell_type <- "naive_CD8T"

# Read every file for a celltype and all conditions and then merge all the result dataframes in only one. 

for (direction in c("up", "down")) {
  enrich_df <- data.frame()
  for(tp in timepoints) {
    file <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, 
                    "/", cell_level, "/", cell_type, "/Enrichments/", pheno, "_", tp, "_", 
                                                                 direction, "_", ont,  "_GO_enrichments.rds")
    if (file.exists(file)){
      print(paste0("Reading file: ", file))
      df <- readRDS(file)
      df <- df@result
      df$direction <- direction
      df$timepoint <- tp
      df <- df %>% arrange(desc(p.adjust))
      enrich_df <- rbind(enrich_df, df)
    }
  }
  saveRDS(enrich_df, paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, 
                "/", cell_level, "/", cell_type,"/Enrichments/All_df_", pheno, "_", direction, "_", ont,  "_GO_enrichments.rds"))
}  

# Get reducedTerms for the whole dataframe (CELL_TYPE AND ALL CONDITIONS DF)

for (direction in c("up", "down")) {
  #direction <- "down"; cell_level <- "L2"; cell_type <- "naive_CD8T"; pheno <- "age"; ont <- "BP"
  GO_all <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, "/", cell_level,
                          "/", cell_type, "/Enrichments/All_df_", pheno, "_", direction, "_", ont,  "_GO_enrichments.rds"))
  simMatrix <- calculateSimMatrix(GO_all$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
  scores <- setNames(-log10(GO_all$qvalue), GO_all$ID)
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.9, orgdb ="org.Hs.eg.db")
  saveRDS(reducedTerms, paste0(data_path, "/mchacon/1M-scBloodNL/robjects/01_log2cpm_DEA_Age_Sex/", chem_version, "/", cell_level,
                              "/", cell_type, "/Enrichments/", pheno, "_", direction, "_", ont,  "_0.9.reducedTerms.rds"))

}
