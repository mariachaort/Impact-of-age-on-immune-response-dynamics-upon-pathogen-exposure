main_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"
main_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq/"

data_path <- "/gpfs/projects/bsc83/Data/scRNAseq/Oelen2022/"
source(paste0(main_path, "mchacon/1M-scBloodNL/scripts/Functions.R"))
source(paste0(main_path, "mchacon/1M-scBloodNL/scripts/Themes.R"))

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))
option_list = list(
  make_option(c("--chemistry"), action="store", default="V3", type="character",
              help="Chemistry version used"),
  make_option(c("--cell_level"), action="store", default="L2", type="character",
              help="Cell level of resolution")
)
opt = parse_args(OptionParser(option_list=option_list))
###########################################################################################################
chem_version <- opt$chemistry
cell_level <- opt$cell_level

cell_types <- c("mDC", "th2_CD4T", "naive_CD4T", "memory_CD8T", "NKdim", "mono_1", "mono_2", "NK", "naive_CD8T", "hemapoietic_stem", "th1_CD4T", 
                "mono_4", "unknown", "B", "NKbright", "megakaryocyte", "reg_CD4T", "pDC", "memory_CD4T", "plasma_B")

phenotypes <- c("age")
conditions <- c("UT", "PA", "MTB", "CA")
times <- c("0h", "3h", "24h")
###########################################################################################################
celltype <- "naive_CD8T"; pheno <- "Gender"; cond <- "MTB"; time <- "3h"
get_tables_pseudobulk <- function(cell_level, celltype, pheno, cond, time) {
  file <- paste0(main_path, "mchacon/1M-scBloodNL/pseudobulk_inrt_lmer/sum/nCells_per_donor/", chem_version, "/", cell_level, "/",
                 celltype, "/Gender_age/min_prop.0.2/", 
                 pheno, "_", cond, "_", time, ".lmer_nearZeroVar.by_metric.rds")
  
  print(paste0("Reading file:", file))
  if (file.exists(file)) {
    df <- readRDS(file)
    #df <- df[, c("fdr.log2cpm", "ID", "estimate.log2cpm")]     # Filter relevant columns
    df$celltype <- celltype
    df$pheno <- pheno
    df$condition <- cond
    df$time <- time
    df$timepoint <- paste0(cond, time, collapse = "_")
    df$adj.P.Val <- df$fdr.log2cpm
    df$logFC <- df$estimate.log2cpm
    #df <- df[df$fdr.log2cpm < 0.05, ]
    #df$direction <- ifelse(df$estimate.log2cpm > 0, "up", "down")
    return(df)
  } else {
    print(paste0("File does not exist: ", file))
    return(NULL)
  }
}

# To make sure dataframes have the same columns
standardize_columns <- function(df_list) {
  all_columns <- unique(unlist(lapply(df_list, colnames)))
  df_list_standardized <- lapply(df_list, function(df) {
    missing_cols <- setdiff(all_columns, colnames(df))
    df[missing_cols] <- NA  # Añadir las columnas faltantes con NA
    return(df)
  })
  return(df_list_standardized)
}

for (pheno in phenotypes) {
  temp_data_df <- list()
  
  for (cond in conditions) {
    for (time in times) {
      deg_list <- lapply(cell_types, function(celltype) get_tables_pseudobulk(cell_level, celltype, pheno, cond, time))
      
      deg_list <- deg_list[!sapply(deg_list, is.null)]
      
      if (length(deg_list) > 0) {
        deg_list <- standardize_columns(deg_list)
        deg_df <- do.call(rbind.data.frame, deg_list)
        temp_data_df <- append(temp_data_df, list(deg_df))
      } else {
        print(paste0("No data available for ", pheno, ", ", cond, ", ", time))
      }
    }
  }
  
  if (length(temp_data_df) > 0) {
    print("Binding final df...")
    final_df <- do.call(rbind.data.frame, temp_data_df)
    saveRDS(final_df, paste0(main_path, "mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/",
                             pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
    print(paste0("Saved file: ", pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
  } else {
    print(paste0("No data available for ", pheno))
  }
}
