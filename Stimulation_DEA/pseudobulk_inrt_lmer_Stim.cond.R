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
  make_option(c("--chemistry"), action="store", default="V2", type="character",
              help="Chemistry version used"),
  make_option(c("--cell_level"), action="store", default='cell_type', type='character',
              help="From Azimuth: low resolution (cell_type_lowerres) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default="naive_CD8T", type='character',
              help="Cell types in low resolution (cell_type_lowerres) or high resolution (cell_type)"),
  make_option(c("--condition"), action="store", default="PA", type="character",
              help="Stimuli used: PA, CA, MTB"),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--phenotype"), action="store", default="Stim.cond", type='character',
              help="condition"),
  make_option(c("--nCells_per_donor"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--random"), action="store", default="date,assignment", type='character',
              help="date"),
  make_option(c("--min_prop"), action="store", default=0.2, type='character',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--in_dir"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/pseudobulk_dreamlet', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/pseudobulk_inrt_lmer', type='character',
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
shhh(library(emmeans))

################################## Set Variables and load Data ##################################
chem_version <- opt$chemistry
condition <- opt$condition

# Output directory
if(opt$cell_level == "cell_type") {cl <- "L2"} else {cl <- "L1"}
sdir <- paste0(opt$aggr_fun, '/')
if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'cell_permutations/')}

if(opt$nCells_per_donor){sdir <- paste0(sdir, 'nCells_per_donor/')}

sdir <- paste0(sdir, chem_version, '/', cl, '/', opt$cell_type, '/Stim.cond_age_Gender/', 'min_prop.', opt$min_prop, '/')

if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'it_', opt$permutation_idx, '/')}

# Input directory
in.dir <- paste0(opt$in_dir, '/', sdir)
print(paste0('Main input directory: ', in.dir))
if(!dir.exists(in.dir)){stop('Input directory missing.')}

# Output directory
out.dir <- paste0(opt$out_dir, '/sum/', chem_version, '/', cl, '/', opt$cell_type, '/Stim.cond/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ', out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Random: ', opt$random))

################################## Functions ################################## 
# (main) Read input data + normalize
# in_dir = in.dir
process_data <- function(in_dir = in.dir, condition){
  # Read input files
  pb_fn <- paste0(in_dir, condition, '_Alltimes_pb.Rds')
  dreamlet_fn <- paste0(in_dir, condition, '_Alltimes_dea_topTable.rds')
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
  
  # Subset gene expression data based on filtered metadata
  pb_filt.ge <- as.matrix(assays(pb)[[1]])
  pb_filt.ge <- pb_filt.ge[rownames(pb_filt.ge) %in% genes_kept, colnames(pb_filt.ge) %in% rownames(pb_filt.md)]
  
  # Get formula
  de_form <- res_ct.proc$formula
  
  # edgeR::cpm() transformation
  log2cpm.mat <- edgeR::cpm(pb_filt.ge, log = TRUE)
  log2cpm.mat <- log2cpm.mat[,match(rownames(pb_filt.md), colnames(log2cpm.mat))]
  # Ensure metadata and expression data are aligned
  stopifnot(identical(colnames(log2cpm.mat), rownames(pb_filt.md)))

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
  
  return(x)
}

# 2. LMM using lme4::lmer()
# Function
# g <- genes_expressed[1]
# g <- "GIMAP6"
# norm <- norm_methods[1]
# process_data_list = process_data.list
# phenotype = opt$phenotype
# phenotype_order = Group_order
# freqCut_nzv = 95/5
# uniqueCut_nzv = 10

lmer_func <- function(g, norm, process_data_list = process_data.list, phenotype = opt$phenotype, 
                      phenotype_order = Group_order, freqCut_nzv = 95/5, uniqueCut_nzv = 10){
  print(paste0('# Normalize: ', norm))
  print(paste0('# Gene: ', g))
  
  ### Prepare Data ###
  norm_i <- paste0('expr_', norm)
  norm_data <- process_data_list[[norm_i]]
  gene_vec <- norm_data[g,]
  
  # Apply inverse-normal rank transformation
  donors <- names(gene_vec)
  gene_expr <- unname(gene_vec)
  gene_expr_inrt <- rankTransform(gene_expr)
  names(gene_expr_inrt) <- donors
  
  # Get donor metadata
  donor_md <- process_data_list$donor_metadata
  donor_md$Donor <- rownames(donor_md)
  
  # Create Dataframe
  gene_df <- data.frame(Donor=names(gene_expr_inrt), value=unname(gene_expr_inrt))
  df_i <- merge(gene_df, donor_md, by = 'Donor')
  rownames(df_i) <- df_i$Donor
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
  
  # Check variance of expression
  if(var(df_i$value) != 0){
    form <- as.formula("value ~ Stim.cond + age + (1|assignment) + Gender + (1|date)")
    print(paste0('Fitting lmer: ', deparse(form)))
    
    # Fit the model
    mod <- lmer(form, data = df_i)
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    #tidy_mod <- as.data.frame(tidy_mod)
    anova(mod)
    
    # Compute pairwise contrasts: A model where you're studying the effects of Stim.cond on gene expression, adjusting for age, Gender (and random effects) 
    # The emmeans() function computes the estimated marginal means for the factor Stim.cond 
    # (i.e., the adjusted mean expression for each level of Stim.cond after considering the effects of other covariates).
    # When using emmeans() to estimate marginal means, you obtain the adjusted means for each level of a factor (or a combination of factors). 
    # These adjusted means are not raw means, but means that account for other variables in the model.
    emmeans_results <- emmeans(mod, ~ Stim.cond)
    contrast_df <- contrast(emmeans_results, method = "pairwise", adjust ="tukey")
    contrast_df <- as.data.frame(contrast_df)
    print(contrast_df)

  } else {
    print(paste0('No variance in fitted value. All values are: ', as.character(unique(df_i$value))))
    
    # Create empty results if there's no variance
    tidy_mod <- data.frame(effect = rep('fixed', 4),
                           term = c('(Intercept)', 'GenderF', 'condition', 'age'))
    tidy_mod.na <- as.data.frame(matrix(NA, nrow = 4, ncol = 7))
    colnames(tidy_mod.na) <- c('estimate', 'std.error', 'statistic',
                               'df', 'p.value', 'conf.low', 'conf.high')
    tidy_mod <- cbind(tidy_mod, tidy_mod.na)
    
    contrast_df <- data.frame(contrast = c("PA24h - PA3h", "PA24h - UT0h", "PA3h - UT0h"),
                              estimate = NA, std.error = NA, df = NA, t.ratio = NA, p.value = NA)
  }
  
  out <- list(model = contrast_df, nearZeroVar = nearZeroVar(df_i[, "value"], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, saveMetrics = TRUE))
  
  cat('\n')
  return(out)
}


################################## Analyses #################################### 
###### Read input data + normalize ######
process_data.list <- process_data(in.dir, condition)

###### Fit lmer ######
# Variables
genes_expressed <- rownames(process_data.list$expr_log2cpm)
norm_methods <- c('log2cpm')
stim_levels <- list(
  CA = c('CA24h', 'CA3h', 'UT0h'),
  PA = c('PA24h', 'PA3h', 'UT0h'),
  MTB = c('MTB24h', 'MTB3h', 'UT0h')
)
Group_order <- NULL
phenotype_term <- opt$phenotype
if(opt$phenotype%in%c('Stim.cond')){
  Group_order <- stim_levels[[opt$condition]]
}

# Apply function (parallel)
# testing
# genes_expressed.all <- genes_expressed #testing
# set.seed(123) #testing
# genes_expressed <- sample(genes_expressed.all, 50) #testing
system.time(lmer_out <- sapply(norm_methods, function(i)
  sapply(genes_expressed, function(j) lmer_func(g = j, norm = i), simplify = FALSE), simplify = FALSE)) #to debug

# print('Fitting lmer...')  
# system.time(lmer_out <- setNames(
#     mclapply(norm_methods, function(i){
#       setNames(
#         mclapply(genes_expressed, function(j) lmer_func(g = j, norm = i)), 
#         genes_expressed  # Set names for inner mclapply output
#       )
#     }),
#     norm_methods  # Set names for outer mclapply output
#   )
# ) #faster
  
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
lmer_results <- lmer.by_metric$log2cpm

final_df_list <- lapply(names(lmer_results), function(gene_name) {
  # Get the summary_emm object for each gene
  gene_results <- lmer_results[[gene_name]]  # Extract the dataframe for this gene
  # Add the gene name as a new column
  gene_results$ID <- gene_name
  
  # Adjust p-values using FDR
  gene_results$adj.p.value <- p.adjust(gene_results$p.value, method = "fdr")
  
  # Select the columns we want in the final dataframe
  gene_results <- gene_results[, c("ID", "contrast", "estimate", "SE", "df", "t.ratio", "p.value", "adj.p.value")]
  
  return(gene_results)
})
lmer.by_metric.out <- do.call(rbind, final_df_list)

lmer.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", condition, '.lmer.by_metric.rds')
print(paste0('Saving output (tidy lmer) in: ', lmer.by_metric.fn))
saveRDS(lmer.by_metric.out, lmer.by_metric.fn)

## nearZeroVar output
nearZeroVar.bymetric.fn <- paste0(out.dir, opt$phenotype, "_", condition, '.nearZeroVar.by_metric.rds')
print(paste0('Saving output (nearZeroVar report) in: ', nearZeroVar.bymetric.fn))
saveRDS(nearZeroVar.bymetric.out, nearZeroVar.bymetric.fn)

# Report failed genes
# check_list.norm_methods <- lapply(lmer.by_metric.out, function(x){xx <- is.na(unique(x[['estimate']])); names(xx) <- rownames(x); return(xx)})
# check_failed <- function(i){
#   check_list <- check_list.norm_methods[[i]]
#   missing_cases <- length(check_list[check_list=='TRUE'])
#   missing_genes <- names(check_list[check_list=='TRUE'])
#   total_cases <- length(check_list)
#   complete_cases <- total_cases - missing_cases
#   print(paste0('There is no variability (', i, ') in: ', missing_cases, ' DEGs...'))
# }
# check_var <- lapply(names(check_list.norm_methods), function(i) check_failed(i))

# Arrange data
print('Arranging data...')
## add nearZeroVar info
inrt_df <- merge(lmer.by_metric.out, nearZeroVar.bymetric.out$log2cpm, by = 'ID')
lmer.by_gene.df <- inrt_df 

## save
lmer_nearZeroVar.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", condition, '.lmer_nearZeroVar.by_metric.rds')
print(paste0('Saving joined output (tidy lmer + nearZeroVar report) in: ', lmer_nearZeroVar.by_metric.fn))
saveRDS(lmer.by_gene.df, lmer_nearZeroVar.by_metric.fn)

## check
print('FDR < 0.05...')
table(lmer.by_gene.df$adj.p.value<0.05, lmer.by_gene.df$estimate>0)

cat('\n')

print('FDR < 0.1...')
table(lmer.by_gene.df$adj.p.value<0.1, lmer.by_gene.df$estimate>0)

cat('\n')

print('p-nominal < 0.01...')
table(lmer.by_gene.df$adj.p.value<0.01, lmer.by_gene.df$estimate>0) 

cat('\n')

print('p-nominal < 0.05...')
table(lmer.by_gene.df$adj.p.value<0.05, lmer.by_gene.df$estimate>0)

cat('\n')

print('p-nominal < 0.1...')
table(lmer.by_gene.df$adj.p.value<0.1, lmer.by_gene.df$estimate>0)
