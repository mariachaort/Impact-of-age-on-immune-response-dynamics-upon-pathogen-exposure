#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/mchacon/Desktop/cluster/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
  .libPaths(c(.libPaths(), "/gpfs/apps/MN5/GPP/R/4.3.0/GCC/lib64/R/library"))
  #.libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1/INTEL/lib64/R/library"))
  print(.libPaths())
} else if(file.exists(path_em)){
  setwd(paste(path_em))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--filtering"), action="store", default="meta", type="character",
              help="input genes to filter"),
  make_option(c("--chemistry"), action="store", default="V2", type="character",
              help="Chemistry version used"),
  make_option(c("--cell_level"), action="store", default='L2', type='character',
              help="L1 or L2"),
  make_option(c("--cell_type"), action="store", default="plasma_B", type='character',
              help="Cell type"),
  make_option(c("--combo"), action="store", default="PA_3h.UT_0h", type="character",
              help="Ratio to be modeled"),
  make_option(c("--method"), action="store", default="residuals.lmer_lm", type="character",
              help="date_combined, residuals.lmer_lm or residuals.lm_lmer"),
  make_option(c("--phenotype"), action="store", default="age", type='character',
              help="Gender/Age"),
  make_option(c("--nCells_per_donor"), action="store", default=TRUE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--in_dir"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/pseudobulk_dreamlet/sum/', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/robjects/log_ratio/', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(SingleCellExperiment))
shhh(library(edgeR))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyverse))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(lme4))
shhh(library(lmerTest))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(parallel))
shhh(library(caret))
shhh(library(DescTools))

################################## Set Variables and load Data ##################################
chem_version <- opt$chemistry
cell_level <- opt$cell_level
cell <- opt$cell_type
nCells_per_donor <- opt$nCells_per_donor
combo <- opt$combo
method <- opt$method
in_dir <- opt$in_dir
out_dir <- opt$out_dir

#combo <- "PA_3h.UT_0h"
split_combo <- strsplit(combo, "\\.")[[1]]

stim1 <- sub("_.*", "", split_combo[1])
t1 <- gsub("h", "", sub(".*_", "", split_combo[1]))

stim2 <- sub("_.*", "", split_combo[2])
t2 <- gsub("h", "", sub(".*_", "", split_combo[2]))

print("#Conditions#")
print(paste("Stim1:", stim1, "T1:", t1))
print(paste("Stim2:", stim2, "T2:", t2))
print(paste("Cell type:", cell))
print(paste("nCells_per_donor:", nCells_per_donor))
print(paste("Phenotype:", opt$phenotype))
print(paste("Method:", method))

in.dir <- paste0(path_cluster, in_dir)
if(opt$nCells_per_donor){in.dir <- paste0(in.dir, 'nCells_per_donor/')}
in.dir <- paste0(in.dir, chem_version, '/', cell_level, '/', cell, '/Gender_age/min_prop.0.2/')

# Input directory
print(paste0('Main input directory: ', in.dir))

# Output directory
if(opt$nCells_per_donor){out.dir <- paste0(path_cluster, out_dir, 'nCells_per_donor/')} else {out.dir <- paste0(path_cluster,out_dir)}
out.dir <- paste0(out.dir, chem_version, "/", cell_level, "/", cell, "/", method, "/")
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ', out.dir))

if (opt$phenotype == "age" | opt$phenotype == "Gender"){
  if (opt$filtering == "meta"){
    unique_sg <- readRDS(paste0(path_cluster, "Projects/scRNAseq/mchacon/1M-scBloodNL/scripts/analysis/significant_genes.union_meta_aida.rds"))
    print("Filtering for input genes is: genes from meta-analysis (aida's)")
    # this file already contains the union of: 
    #/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg3-sc_pseudobulk_DEA-consortium/I2_100.nDonors_100/pseudobulkDEA_inrt_lmer/L1/*/age/Intersect_3/common_recalculated.datasets_meta.list.rds
  } else if (is.null(opt$filtering)){
    genes_age <- readRDS(paste0(path_cluster, "Projects/scRNAseq/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/", 
                                opt$phenotype, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
    
    signif_age <- genes_age[genes_age$fdr.log2cpm < 0.05, ]
    
    unique_sg <- unique(signif_age$ID)
    print("Filtering for input genes is: age-DEGs from the age-DEA (with nCells_per_donor)")
    
  } else if (opt$filtering == "stim.cond_age"){
    genes_age <- readRDS(paste0(path_cluster, "Projects/scRNAseq/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/", 
                                opt$phenotype, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))
    
    signif_age <- genes_age[genes_age$fdr.log2cpm < 0.05, ]
    genes_stim <- readRDS(paste0(path_cluster, "Projects/scRNAseq/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/", 
                                 "Stim.cond_", chem_version, "_", cell_level, "_log2cpm_dframe.rds"))
    signif_stim <- genes_stim[genes_stim$fdr.log2cpm < 0.05, ]
    
    # Extract gene names
    genes_age_sig <- signif_age$ID
    genes_stim_sig <- signif_stim$ID
    
    # Get intersection
    sg <- intersect(genes_age_sig, genes_stim_sig)
    unique_sg <- unique(sg)

  }
  
}

################################## Functions ################################## 
# (main) Read input data + normalize
# in_dir = in.dir
process_data <- function(in_dir, cond, time){
  # Read input files
  pb_fn <- paste0(in_dir, cond, "_", time, 'h_pb.Rds')
  dreamlet_fn <- paste0(in_dir, cond, "_", time, 'h_dea_vp_topTable.rds')
  if(!file.exists(pb_fn) | !file.exists(dreamlet_fn)){
    stop('No pseudobulk or dreamlet output files.')
  }
  
  print(paste0('Reading pseudobulk file in: ', pb_fn))
  system.time(pb <- readRDS(pb_fn))
  
  print(paste0('Reading dreamlet output file in: ', dreamlet_fn))
  system.time(dreamlet_res <- readRDS(dreamlet_fn))
  res.proc <- dreamlet_res$processed #from processAssays()
  
  # Filter pb 
  ## View details of dropping samples
  details(res.proc)
  
  ## Keep samples & genes after processAssays()
  res_ct.proc <- res.proc[[1]]
  voomWithDreamWeights.mat <- res_ct.proc$E
  samples_kept <- colnames(voomWithDreamWeights.mat)
  genes_kept <- rownames(voomWithDreamWeights.mat)
  
  ## Check nSamples and nGenes tested
  genes_all <- rownames(pb)
  genes_tested <- rownames(voomWithDreamWeights.mat)
  genes_all.n <- nrow(pb)
  genes_tested.n <- nrow(voomWithDreamWeights.mat)
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(pb)
  samples_tested <- colnames(voomWithDreamWeights.mat)
  samples_all.n <- ncol(pb)
  samples_tested.n <- ncol(voomWithDreamWeights.mat)
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))
  
  ## Filter pb 
  pb.ge <- as.matrix(assays(pb)[[1]])
  pb_filt.ge <- pb.ge[rownames(pb.ge)%in%genes_kept, colnames(pb.ge)%in%samples_kept]
  
  # Get metadata
  pb.md <- as.data.frame(colData(pb))
  pb_filt.md <- pb.md[rownames(pb.md)%in%samples_kept,]
  
  # Get formula
  de_form <- res_ct.proc$formula
  
  # edgeR::cpm() transformation
  log2cpm.mat <- edgeR::cpm(pb_filt.ge, log = TRUE)
  log2cpm.mat <- log2cpm.mat[,match(rownames(pb_filt.md), colnames(log2cpm.mat))]
  identical(colnames(log2cpm.mat), rownames(pb_filt.md)) #check
  
  # variancePartition::voomWithDreamWeights() transformation
  #voomWithDreamWeights.mat <- voomWithDreamWeights.mat[,match(rownames(pb_filt.md), colnames(voomWithDreamWeights.mat))]
  #identical(colnames(voomWithDreamWeights.mat), rownames(pb_filt.md)) #check
  
  # Output
  ## check
  dim(pb_filt.ge)
  #dim(voomWithDreamWeights.mat)
  dim(log2cpm.mat)
  pb_filt.ge[1:3,1:3]
  #voomWithDreamWeights.mat[1:3,1:3]
  log2cpm.mat[1:3,1:3]
  
  ## return
  out <- list(expr_counts = pb_filt.ge, #raw
              #expr_voomWithDreamWeights = voomWithDreamWeights.mat, #voomWithDreamWeights
              expr_log2cpm = log2cpm.mat, #log2cpm
              donor_metadata = pb_filt.md,
              de_form = de_form)
  return(out)
}


# (acc) Inverse-normal rank transformation
# x <- gene_expr
rankTransform <- function(x){
  require(DescTools)
  x_norm <- x
  print(paste0('nSamples: ', length(x)))
  notNA <- which(!is.na(x))
  print(paste0('nSamples (not NA): ', length(x[notNA])))
  percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x)+1)
  # percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x[notNA])+1) #check with Maxime
  mean_level <- mean(Winsorize(x[notNA]))
  sd_level <- sd(Winsorize(x[notNA]))
  x[notNA] <- qnorm(percentile, mean_level, sd_level)
  
  # check normal distribution
  ## If the p-value > 0.05: Fail to reject the null hypothesis, meaning the data is likely normal.
  ## If the p-value ≤ 0.05: Reject the null hypothesis, meaning the data is not normal.
  # x_norm_inrt <- x
  # shapiro.test(x_norm)
  # shapiro.test(x_norm_inrt)
  # ks.test(x_norm, "pnorm", mean = mean(x_norm), sd = sd(x_norm))
  # ks.test(x_norm_inrt, "pnorm", mean = mean(x_norm_inrt), sd = sd(x_norm_inrt))
  
  return(x)
}

# 2. LMM using lme4::lmer()
# Function
# g <- genes_expressed[1]
# norm <- norm_methods[2]
# process_data_list = process_data.list
# phenotype = opt$phenotype
# phenotype_order = Group_order
# freqCut_nzv = 95/5
# uniqueCut_nzv = 10
lmer.lm_func <- function(g, norm, process_data_list = final_log_ratio, phenotype = opt$phenotype, phenotype_order = Group_order, freqCut_nzv = 95/5, uniqueCut_nzv = 10){
  print(paste0('# Normalize: ', norm))
  print(paste0('# Gene: ', g))
  print(paste0('# Method: ', method))
  
  ### prepare data ###
  # get expression data
  ## normalized
  norm_i <- paste0('expr_', norm)
  norm_data <- process_data_list[[norm_i]]
  gene_vec <- norm_data[g,]
  
  ## inverse-normal rank transformation
  donors <- names(gene_vec)
  gene_expr <- unname(gene_vec)
  gene_expr_inrt <- rankTransform(gene_expr)
  names(gene_expr_inrt) <- donors
  
  # get donor metadata
  donor_md <- process_data_list$donor_metadata
  donor_md$assignment <- rownames(donor_md)
  
  # dataframe
  gene_df <- data.frame(assignment=names(gene_expr_inrt), value=unname(gene_expr_inrt))
  df_i <- merge(gene_df, donor_md, by = 'assignment')
  rownames(df_i) <- df_i$assignment
  df_i <- df_i[,-1]
  df_i <- df_i %>% mutate_if(is.character, as.factor)
  if(is.factor(df_i[[phenotype]])){
    df_i[[phenotype]] <- factor(df_i[[phenotype]],
                                levels = phenotype_order)
  }

  # check covariates variance --> nearZeroVar()
  metrics <- 'value'
  nearZeroVar.metrics <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, names=TRUE)
  nearZeroVar.df <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, saveMetrics=TRUE)
  nearZeroVar.df$nObs <- nrow(df_i)
  if(length(nearZeroVar.metrics)>0){
    print(paste0('nearZeroVar: ', nearZeroVar.metrics))
  }

  # formula (testing)
  # covs <- c('Age', 'Gender', '(1|date)')
  # fmla_covs <- paste(covs, collapse = '+')
  
  form_tmp <- process_data_list$de_form
  fmla_covs <- deparse(form_tmp)
  fmla <- paste0('value',fmla_covs)
  form <- as.formula(fmla)
  print(paste0('Fitting model: ',fmla))

  # check variance of the fitted value
  if(var(df_i$value)!=0){
    if (method == "residuals.lm_lmer" | method == "date_combined") {
      mod <- lmerTest::lmer(form, data = df_i)
      tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
      tidy_mod <- as.data.frame(tidy_mod)
      
    } else if (method == "residuals.lmer_lm") {
      mod <- lm(form, data = df_i)
      tidy_mod <- broom::tidy(mod, conf.int = TRUE)
      tidy_mod <- as.data.frame(tidy_mod)
      
    } else {
      stop("Invalid method. Use 'residuals.lm_lmer' for mixed model or 'residuals.lmer_lm' for linear model.")
    }
  }else{
    print(paste0('There is no variance of the fitted value. All values equal to: ', as.character(unique(df_i$value))))
    tidy_mod <- data.frame(effect = rep('fixed',3),
                           term = c('(Intercept)', 'age', 'GenderF'))
    tidy_mod.na <- as.data.frame(matrix(NA, nrow = 3, ncol = 7))
    colnames(tidy_mod.na) <- c('estimate', 'std.error', 'statistic',
                               'df', 'p.value', 'conf.low', 'conf.high')
    tidy_mod <- cbind(tidy_mod, tidy_mod.na)
  }
  out <- list(model = tidy_mod,
              nearZeroVar = nearZeroVar.df)

  cat('\n')
  return(out)
}
################################## Analyses #################################### 
###### Read input data + normalize ######
system.time(process_data.list1 <- process_data(in.dir, stim1, t1))
system.time(process_data.list2 <- process_data(in.dir, stim2, t2))

log2cpm_m1 <- process_data.list1$expr_log2cpm
log2cpm_m2 <- process_data.list2$expr_log2cpm

print("FIltering by genes in the union...")
log2cpm_m1 <- log2cpm_m1[rownames(log2cpm_m1) %in% unique_sg, , drop = FALSE]
log2cpm_m2 <- log2cpm_m2[rownames(log2cpm_m2) %in% unique_sg, , drop = FALSE]

print("Number of genes in both matrices...")
nrow(log2cpm_m1)
nrow(log2cpm_m2)

colData1 <- process_data.list1$donor_metadata
colData1$assignment <- rownames(colData1)
colData2 <- process_data.list2$donor_metadata
colData2$assignment <- rownames(colData2)

print("Filtering by common donors and genes...")
common_donors <- intersect(colnames(log2cpm_m1), colnames(log2cpm_m2))
common_genes <- intersect(rownames(log2cpm_m1), rownames(log2cpm_m2))

log2cpm_common.1 <- log2cpm_m1[common_genes, common_donors, drop = FALSE]  

log2cpm_common.2 <- log2cpm_m2[common_genes, common_donors, drop = FALSE]  

# Verify both matrices have the same donors and genes in the same order
print("Checking if matrices have the same donors and genes...")
if ((identical(rownames(log2cpm_common.1), rownames(log2cpm_common.2)) == FALSE |identical(colnames(log2cpm_common.1), 
                                                                                           colnames(log2cpm_common.2)) == FALSE)){
  stop("Rownames or colnames in log2cpm matrices are not identical.")
  
}

genes <- common_genes
donors <- common_donors

################################## DATE COMBINED ##############################
### 1. Create a new variable for the batch effect ('date') --> maybe too much stratified
# Prepare data
md <- readRDS('/gpfs/projects/bsc83/Data/scRNAseq/Oelen2022/1M_full_metadata_modified.rds')
donor_comb_md <- unique(md[,c('assignment', 'Gender', 'age', 'time', 'condition', 'date', 'date_lane')])
donor_comb_md <- donor_comb_md[-which(is.na(donor_comb_md$assignment)),]
donor_comb_md_list <- split(donor_comb_md, donor_comb_md$assignment)
donor_comb_md_list <- lapply(donor_comb_md_list, function(x) x[order(x$date_lane),])
donor_comb_md.i <- droplevels(donor_comb_md[(donor_comb_md$condition == stim1 & donor_comb_md$time == t1) | (donor_comb_md$condition == stim2 & donor_comb_md$time == t2),])
donor_comb_md.i %>% group_by(assignment) %>% reframe(date_combined = paste(date,collapse='_')) %>% as.data.frame() -> donor_comb_md_date.i 
sort(table(donor_comb_md_date.i$date_combined),decreasing=T) #check 
donor_md <- droplevels(unique(donor_comb_md[,c('assignment', 'Gender', 'age')]))
donor_md_new <- merge(donor_md, donor_comb_md_date.i, by = 'assignment')
rownames(donor_md_new) <- donor_md_new$assignment
donor_md_new <- donor_md_new[,-which(colnames(donor_md_new)=='assignment')]
donor_md_new <- droplevels(donor_md_new)

library(stringi);library(stringr)
mapped_dates <- sapply(donor_md_new$date_combined, function(x) {
  sorted_dates <- paste(sort(unique(strsplit(x, "_")[[1]])), collapse = "_")
  return(sorted_dates)
})

donor_md_new$new_date_combined <- mapped_dates

donor_md_new_common <- donor_md_new[rownames(donor_md_new)%in%donors,] 

colData1$new_date_combined <- donor_md_new_common$new_date_combined[match(colData1$assignment, rownames(donor_md_new_common))]
colData2$new_date_combined <- donor_md_new_common$new_date_combined[match(colData2$assignment, rownames(donor_md_new_common))]

### 2. Extract residuals (2 different options) ###
# Prepare data
if (method == "residuals.lm_lmer" | method == "residuals.lmer_lm"){
  in_list1 <- list(ge_mat = log2cpm_common.1, md = colData1)
  in_list2 <- list(ge_mat = log2cpm_common.2, md = colData2)
  
  # Function
  # g <- genes[1]
  # in_list <- in_list1
  # in_list <- in_list2
  get_residuals <- function(g, in_list){
    # prepare DF
    ge_mat <- in_list[['ge_mat']]
    md <- in_list[['md']]
    g_vec <- ge_mat[rownames(ge_mat)==g,]
    g_df <- data.frame(assignment = names(g_vec), 
                       ge = unname(g_vec))
    g_df <- merge(g_df, md, by = 'assignment')
    g_df <- g_df[order(g_df$assignment),]
    g_df <- g_df[order(as.numeric(g_df$assignment)),]
    rownames(g_df) <- g_df$assignment
    g_df <- g_df[,-which(colnames(g_df)=='assignment')]
    g_df <- droplevels(g_df)
    sort(table(g_df$date),decreasing=T)
    
    # fit model
    if (method == "residuals.lmer_lm"){
      lmer_fit <- lmer(ge ~ Gender + nCells + (1|date), data=g_df)
      res <- residuals(lmer_fit)
    } else if (method == "residuals.lm_lmer") {
      lm_fit <- lm(ge ~ nCells, data=g_df)
      res <- residuals(lm_fit)
    }
    return(res)
  }
  # Apply function
  shhh(library(parallel))
  
  res_l1 <- mclapply(genes, function(i) get_residuals(i, in_list1))
  res_l2 <- mclapply(genes, function(i) get_residuals(i, in_list2))
  names(res_l1) <- genes
  names(res_l2) <- genes
  matrix1_common_res <- do.call("rbind",res_l1)
  matrix2_common_res <- do.call("rbind",res_l2)
  
  common_samples <- intersect(colnames(matrix1_common_res), colnames(matrix2_common_res))
  common_genes <- intersect(rownames(matrix1_common_res), rownames(matrix2_common_res))
  matrix1_common_res <- matrix1_common_res[common_genes, common_samples, drop = FALSE]
  matrix2_common_res <- matrix2_common_res[common_genes, common_samples, drop = FALSE]
  
  print("Checking if matrices for residuals option have the same donors and genes...")
  if ((identical(rownames(matrix1_common_res), rownames(matrix2_common_res)) == FALSE |identical(colnames(matrix1_common_res), 
                                                                                                 colnames(matrix2_common_res)) == FALSE)){
    stop("Rownames or colnames in log2cpm matrices are not identical.")
    
  }
  log_ratio <- matrix1_common_res - matrix2_common_res
  donor.metadata <- donor_md_new_common[common_samples, ]
  all(rownames(donor.metadata) == colnames(log_ratio))
  
} else if (method == "date_combined"){
  log_ratio <- log2cpm_common.1 - log2cpm_common.2
  donor.metadata <- donor_md_new_common[common_donors, ]
  all(rownames(donor.metadata) == colnames(log_ratio))
}

## Creating formula for next step ##
if (method == "residuals.lmer_lm"){
  form <- as.formula(~age)
} else if (method == "residuals.lm_lmer" | method == "date_combined") {
  form <- as.formula(~ age + Gender + (1|new_date_combined))
}
print("Preparing final log ratio...")
final_log_ratio <- list(expr_log2cpm = log_ratio,
                        donor_metadata = donor.metadata,
                        de_form = form) 

###### Fit lmer or lm ######
# Variables
genes_expressed <- rownames(log_ratio)
norm_methods <- c('log2cpm')

Group_order <- NULL
phenotype_term <- opt$phenotype
if(opt$phenotype%in%c('Gender')){
  Group_order <- c('M','F')
  phenotype_term <- paste0(phenotype_term, Group_order[2])
}

# Apply function (parallel)
## testing
# genes_expressed.all <- genes_expressed #testing
# set.seed(123) #testing
# genes_expressed <- sample(genes_expressed.all, 50) #testing
# system.time(lmer_out <- sapply(norm_methods, function(i)
#   sapply(genes_expressed, function(j) lmer_func(g = j, norm = i), simplify = FALSE), simplify = FALSE)) #to debug

print('Fitting lmer or lm...')  
system.time(lmer_out <- setNames(
    mclapply(norm_methods, function(i){
      setNames(
        mclapply(genes_expressed, function(j) lmer.lm_func(g = j, norm = i)), 
        genes_expressed  # Set names for inner mclapply output
      )
    }),
    norm_methods  # Set names for outer mclapply output
  )
) #faster

# Rearrange the data
items <- c('model', 'nearZeroVar')
lmer_out.items <- sapply(items, function(i) lapply(lmer_out, function(metric) 
  lapply(metric, function(x) x[[i]])), simplify = FALSE)

## nearZeroVar report
nearZeroVar.list <- lmer_out.items$nearZeroVar
nearZeroVar.bymetric.out <- lapply(nearZeroVar.list, function(x){
  df <- do.call("rbind", x)
  df$ID <- rownames(df)
  return(df)
})

# Save output
## LM output
lmer.by_metric <- lmer_out.items$model
lmer.by_metric.out <- lapply(lmer.by_metric, function(x) do.call("rbind", x))
if(!is.null(opt$interaction)){phenotype_term <- grep(':', unique(lmer.by_metric.out[[1]]$term), value = TRUE)}
lmer.by_metric.out <- lapply(lmer.by_metric.out, function(x) x[x$term==phenotype_term,])
lmer.by_metric.out <- lapply(lmer.by_metric.out, function(x){
  rownames(x) <- sub("\\.[^.]+$", "", rownames(x))
  x$ID <- rownames(x)
  return(x)
})
lmer.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", opt$filtering, "_", combo, '.lmer.by_metric.rds')
print(paste0('Saving output (tidy lmer) in: ', lmer.by_metric.fn))
saveRDS(lmer.by_metric.out, lmer.by_metric.fn)

## nearZeroVar output
nearZeroVar.bymetric.fn <- paste0(out.dir, opt$phenotype, "_", opt$filtering, "_", combo, '.nearZeroVar.by_metric.rds')
print(paste0('Saving output (nearZeroVar report) in: ', nearZeroVar.bymetric.fn))
saveRDS(nearZeroVar.bymetric.out, nearZeroVar.bymetric.fn)

# Report failed genes
check_list.norm_methods <- lapply(lmer.by_metric.out, function(x){xx <- is.na(unique(x[['estimate']])); names(xx) <- rownames(x); return(xx)})
check_failed <- function(i){
  check_list <- check_list.norm_methods[[i]]
  missing_cases <- length(check_list[check_list=='TRUE'])
  missing_genes <- names(check_list[check_list=='TRUE'])
  total_cases <- length(check_list)
  complete_cases <- total_cases - missing_cases
  print(paste0('There is no variability (', i, ') in: ', missing_cases, ' DEGs...'))
}
check_var <- lapply(names(check_list.norm_methods), function(i) check_failed(i))

# Arrange data for assigning the mode (x = voomWithDreamWeights and y = log2cpm)
print('Arranging data...')
## add nearZeroVar info
#voom_df <- merge(lmer.by_metric.out$voomWithDreamWeights, nearZeroVar.bymetric.out$voomWithDreamWeights, by = 'ID')
inrt_df <- merge(lmer.by_metric.out$log2cpm, nearZeroVar.bymetric.out$log2cpm, by = 'ID')
lmer.by_gene.df <- inrt_df 
#lmer.by_gene.df <- merge(voom_df, inrt_df, by = c("effect", "term", "ID"), suffixes = c('.voomWithDreamWeights', '.log2cpm'))

## calculate FDR and split by gene again
#lmer.by_gene.df$fdr.voomWithDreamWeights <- p.adjust(lmer.by_gene.df$p.value.voomWithDreamWeights, 'fdr')
lmer.by_gene.df$fdr.log2cpm <- p.adjust(lmer.by_gene.df$p.value, 'fdr')
lmer.by_gene <- split(lmer.by_gene.df, lmer.by_gene.df$ID)

## save
lmer_nearZeroVar.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", opt$filtering, "_", combo, '.lmer_nearZeroVar.by_metric.rds')
print(paste0('Saving joined output (tidy lmer + nearZeroVar report) in: ', lmer_nearZeroVar.by_metric.fn))
saveRDS(lmer.by_gene.df, lmer_nearZeroVar.by_metric.fn)

## check
print('FDR < 0.05...')
#table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.05, lmer.by_gene.df$estimate>0) #log2cpm

cat('\n')

print('FDR < 0.1...')
#table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.1, lmer.by_gene.df$estimate>0) #log2cpm

cat('\n')

print('p-nominal < 0.01...')
#table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.01, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value<0.01, lmer.by_gene.df$estimate>0) #log2cpm

cat('\n')

print('p-nominal < 0.05...')
#table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value<0.05, lmer.by_gene.df$estimate>0) #log2cpm

cat('\n')

print('p-nominal < 0.1...')
#table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value<0.1, lmer.by_gene.df$estimate>0) #log2cpm
