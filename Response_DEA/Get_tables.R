if (getwd()== "/home/mchacon"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
  
}else if(getwd()== "/home/mchacon/Desktop/cluster/Projects/scRNAseq/mchacon"){
  basepath <- "/home/mchacon/Desktop/cluster"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
}else{
  basepath <- "/gpfs/projects/bsc83"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq"
}
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Functions.R"))

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))
option_list = list(
  make_option(c("--chemistry"), action="store", default="V2", type="character",
              help="Chemistry version used"),
  make_option(c("--cell_level"), action="store", default="L2", type="character",
              help="Cell level of resolution")
  
)
opt = parse_args(OptionParser(option_list=option_list))
###########################################################################################################
chem_version <- opt$chemistry
cell_level <- opt$cell_level

#cell_types <- c("CD8T", "CD4T", "B", "NK", "DC", "monocyte")

cell_types <- c("mDC", "th2_CD4T", "naive_CD4T", "memory_CD8T", "NKdim", "mono_1", "mono_2", "NK", "naive_CD8T", "hemapoietic_stem", "th1_CD4T", 
 "mono_4", "unknown", "B", "NKbright", "megakaryocyte", "reg_CD4T", "pDC", "memory_CD4T", "plasma_B")

cell_types <- c("mono_2",  "naive_CD4T", "th1_CD4T", "mono_1", "memory_CD8T_left_and_naive_CD8T_right", "memory_CD8T", "NKdim","mDC",
  "megakaryocyte","NK","mono_4","naive_CD4T_transitioning_to_stim","mono_3","pDC","reg_CD4T","th2_CD4T","unknown","NKbright","B",
  "double_negative_T","T_helper","cyto_CD4T","naive_CD8T","hemapoietic_stem","plasma_B")

cell_types <- c("naive_CD8T")
cell <- "naive_CD8T"; chem_version <- "V3"; cell_level <- "L2"
phenotypes <- c("age")
combo <- c("CA_3h.UT_0h", "CA_24h.UT_0h", "PA_3h.UT_0h", "PA_24h.UT_0h", "MTB_3h.UT_0h", "MTB_24h.UT_0h", 
           "CA_24h.CA_3h", "PA_24h.PA_3h", "MTB_24h.MTB_3h")
nCells_per_donor <- TRUE; method <- "residuals.lmer_lm"
###########################################################################################################
pheno <- "age"
get_tables_pseudobulk <- function(cell_level, cell, c, pheno, method) {
  #c <- "MTB_24h.MTB_3h"
  if (nCells_per_donor == TRUE) {
    file <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/nCells_per_donor/", 
                   chem_version, "/", cell_level, "/", cell, "/", method, "/", pheno, "_", c, ".lmer_nearZeroVar.by_metric.rds")
  } else {
    file <- paste0(data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/", 
                   chem_version, "/", cell_level, "/", cell, "/", method, "/", pheno, "_", c, ".lmer_nearZeroVar.by_metric.rds")
  }
  
  print(paste0("Reading file:", file))
  if(file.exists(file)) {
    df <- readRDS(file)
    #if(pheno == "age") {
    #  toptable <- df$age.topTable
      
    #} else {
     # toptable <- df$Gender.topTable
    #}
    #df <- df[df$fdr.log2cpm < 0.05, ]
    df$logFC <- df$estimate
    df$adj.P.Val <- df$fdr.log2cpm
    df$celltype <- cell
    df$pheno <- pheno
    df$combo <- c
    df$stim <- sub("^(.*?)_.*", "\\1", c)
    
    #toptable$celltype <- cell
    #toptable$pheno <- pheno
    #toptable$combo <- c
    #toptable$stim <- sub("^(.*?)_.*", "\\1", c)
    #return(toptable)
    return(df)
  } else {
    print(paste0("File does not exist:", file))
    return(NULL)
  }
}

for (pheno in phenotypes) {
  temp_data_df <- list()
  
  for (cell in cell_types) {
  # Process each combination in combo
    deg_list <- lapply(combo, function(c) get_tables_pseudobulk(cell_level, cell, c, pheno, method))
    deg_list <- deg_list[!sapply(deg_list, is.null)]  # Remove null entries
  
    if (length(deg_list) > 0) {
      # Extract the listData from each DFrame and convert it to a data.frame
      #deg_list_extracted <- lapply(deg_list, function(x) as.data.frame(x@listData))
      deg_list_extracted <- deg_list
      # Combine all the extracted data.frames into one data frame
      deg_df <- do.call(rbind, deg_list_extracted)
      temp_data_df <- append(temp_data_df, list(deg_df))
      
    } else {
      print(paste0("No data available for ", pheno,",",  cell))
    }
  
   
  }
  if (length(temp_data_df) > 0) {
    print("Binding final df...")
    final_df <- do.call(rbind.data.frame, temp_data_df)
    # Save the final combined data frame as an RDS file
    saveRDS(final_df, paste0(data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/nCells_per_donor/", 
                           chem_version, "/", cell_level, "/", "age_residuals.lmer_lm_ratio_log2cpm_dframe.rds"))
    
    print(paste0("Saved file: ", data_path, "/mchacon/1M-scBloodNL/robjects/log_ratio/", 
                 chem_version, "/", cell_level, "/", "age_residuals.lmer_lm_ratio_log2cpm_dframe.rds"))
  } else {
    print(paste0("No data available for ", pheno))
  }
  
  
}

