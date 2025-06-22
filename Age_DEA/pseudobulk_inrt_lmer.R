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
#f <- readRDS("Projects/scRNAseq/mchacon/1M-scBloodNL/pseudobulk_inrt_lmer/sum/L2/naive_CD8T/Gender_age/min_prop.0.2/age.lmer_nearZeroVar.by_metric.rds")
#genes <- f[f$fdr.log2cpm < 0.05 & f$estimate.log2cpm >0, ]

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--chemistry"), action="store", default="V3", type="character",
              help="Chemistry version used"),
  make_option(c("--cell_level"), action="store", default='cell_type', type='character',
              help="From Azimuth: low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--cell_type"), action="store", default="naive_CD8T", type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--condition"), action="store", default="PA", type="character",
              help="Stimuli used: PA, CA, MTB"),
  make_option(c("--time"), action="store", default="3h", type="character",
              help="3h or 24h (run it separately for UT_0h)"),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--phenotype"), action="store", default="age", type='character',
              help="Gender/Age"),
  make_option(c("--nCells_per_donor"), action="store", default=TRUE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--permutation_idx"), action="store", default=NULL, type='character',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--phenotypes"), action="store", default="Gender_age", type='character',
              help="Gender_age"),
  make_option(c("--random"), action="store", default="date", type='character',
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

################################## Set Variables and load Data ##################################
### Testing ###
# Dreamlet v1.1.9 (https://github.com/GabrielHoffman/dreamlet/blob/devel/NEWS.md, Nov 5 2024) --> module load R/4.3.0; dreamlet in MN5
## Main DEA
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab

## Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$nCells_per_donor <- TRUE

## Cell-level permutations + Add nCells per donors in the model
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$nCells_per_donor <- TRUE
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab

## Main DEA w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Cell-level permutations w/ aggregateToPseudoBulk(fun="mean", ...)
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab

## Cell-level permutations w/ aggregateToPseudoBulk(fun="mean", ...) & Add nCells per donors in the model
# opt$cell_type <- 'CD8_TEM' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.selected.tab
# opt$permutation_idx <- 1 #for array file = pseudobulk_dreamlet.permutation_idx_20times.tab
# opt$aggr_fun <- 'mean' #for array file = pseudobulk_dreamlet.aggr_fun.tab
# opt$nCells_per_donor <- TRUE

#  Dreamlet v1.4.1 (https://github.com/GabrielHoffman/dreamlet/blob/devel/NEWS.md, Nov 5 2024) --> module load miniconda; conda activate dreamlet_1.4.0 (/home/bsc/bsc083616/bin/miniconda3/envs/dreamlet_1.4.0)
## Main DEA
# opt$out_dir <- 'Projects/scRNAseq/aripol1/OneK1K_Age/pseudobulk_dreamlet_1.4.1'
# opt$cell_type <- 'MAIT' #for array file = Azimuth_l2.cell_type.filt.tab

### End Testing ###

# aggr func 
valid_aggr_funcs <- c('sum', 'mean', 'median', 'prop.detected', 'num.detected', 'sem', 'number')
if (!(opt$aggr_fun %in% valid_aggr_funcs)) {
  stop(paste(
    "Invalid value for --aggr_fun:", opt$aggr_fun, 
    "\nValid options are:", paste(valid_aggr_funcs, collapse = ", ")
  ))
}

chem_version <- opt$chemistry
# Output directory
if(opt$cell_level == "cell_type") {cl <- "L2"} else {cl <- "L1"}
sdir <- paste0(opt$aggr_fun, '/')
if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'cell_permutations/')}

if(opt$nCells_per_donor){sdir <- paste0(sdir, 'nCells_per_donor/')}

sdir <- paste0(sdir, chem_version, '/', cl, '/', opt$cell_type, '/', opt$phenotypes, '/', 'min_prop.', opt$min_prop, '/')

if(!is.null(opt$permutation_idx)){sdir <- paste0(sdir, 'it_', opt$permutation_idx, '/')}

# Input directory
in.dir <- paste0(opt$in_dir, '/', sdir)
print(paste0('Main input directory: ', in.dir))
if(!dir.exists(in.dir)){stop('Input directory missing.')}

# Output directory
out.dir <- paste0(opt$out_dir, '/', sdir)
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main output directory: ', out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Cell-level permutations (iteration): ',  opt$permutation_idx))
print(paste0('Phenotypes: ', opt$phenotypes))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Random: ', opt$random))

################################## Functions ################################## 
# (main) Read input data + normalize
# in_dir = in.dir
cond <- opt$condition
t <- opt$time
process_data <- function(in_dir = in.dir){
  # Read input files
  pb_fn <- paste0(in_dir, cond, "_", t, '_pb.Rds')
  dreamlet_fn <- paste0(in_dir, cond, "_", t, '_dea_vp_topTable.rds')
  if(!file.exists(pb_fn) | !file.exists(pb_fn)){
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
  voomWithDreamWeights.mat <- voomWithDreamWeights.mat[,match(rownames(pb_filt.md), colnames(voomWithDreamWeights.mat))]
  identical(colnames(voomWithDreamWeights.mat), rownames(pb_filt.md)) #check

  # Output
  ## check
  dim(pb_filt.ge)
  dim(voomWithDreamWeights.mat)
  dim(log2cpm.mat)
  pb_filt.ge[1:3,1:3]
  voomWithDreamWeights.mat[1:3,1:3]
  log2cpm.mat[1:3,1:3]
  
  ## return
  out <- list(expr_counts = pb_filt.ge, #raw
              expr_voomWithDreamWeights = voomWithDreamWeights.mat, #voomWithDreamWeights
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
lmer_func <- function(g, norm, process_data_list = process_data.list, phenotype = opt$phenotype, phenotype_order = Group_order, freqCut_nzv = 95/5, uniqueCut_nzv = 10){
  print(paste0('# Normalize: ', norm))
  print(paste0('# Gene: ', g))
  
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
  print(paste0('Fitting lmer: ',fmla))

  # check variance of the fitted value
  if(var(df_i$value)!=0){
    # lmer
    mod <-  lmerTest::lmer(form, data = df_i)
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    tidy_mod <- as.data.frame(tidy_mod)
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
system.time(process_data.list <- process_data())

###### Fit lmer ######
# Variables
genes_expressed <- rownames(process_data.list$expr_counts)
norm_methods <- c('voomWithDreamWeights', 'log2cpm')

Group_order <- NULL
phenotype_term <- opt$phenotype
if(opt$phenotype%in%c('Gender','Age_cat', 'Age_cat_all')){
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

print('Fitting lmer...')  
system.time(lmer_out <- setNames(
    mclapply(norm_methods, function(i){
      setNames(
        mclapply(genes_expressed, function(j) lmer_func(g = j, norm = i)), 
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
lmer.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", cond, "_", t, '.lmer.by_metric.rds')
print(paste0('Saving output (tidy lmer) in: ', lmer.by_metric.fn))
saveRDS(lmer.by_metric.out, lmer.by_metric.fn)

## nearZeroVar output
nearZeroVar.bymetric.fn <- paste0(out.dir, opt$phenotype, "_", cond, "_", t, '.nearZeroVar.by_metric.rds')
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
voom_df <- merge(lmer.by_metric.out$voomWithDreamWeights, nearZeroVar.bymetric.out$voomWithDreamWeights, by = 'ID')
inrt_df <- merge(lmer.by_metric.out$log2cpm, nearZeroVar.bymetric.out$log2cpm, by = 'ID')
lmer.by_gene.df <- merge(voom_df, inrt_df, by = c("effect", "term", "ID"), suffixes = c('.voomWithDreamWeights', '.log2cpm'))

## calculate FDR and split by gene again
lmer.by_gene.df$fdr.voomWithDreamWeights <- p.adjust(lmer.by_gene.df$p.value.voomWithDreamWeights, 'fdr')
lmer.by_gene.df$fdr.log2cpm <- p.adjust(lmer.by_gene.df$p.value.log2cpm, 'fdr')
lmer.by_gene <- split(lmer.by_gene.df, lmer.by_gene.df$ID)

## save
lmer_nearZeroVar.by_metric.fn <- paste0(out.dir, opt$phenotype, "_", cond, "_", t, '.lmer_nearZeroVar.by_metric.rds')
print(paste0('Saving joined output (tidy lmer + nearZeroVar report) in: ', lmer_nearZeroVar.by_metric.fn))
saveRDS(lmer.by_gene.df, lmer_nearZeroVar.by_metric.fn)

## check
print('FDR < 0.05...')
table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.05, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('FDR < 0.1...')
table(lmer.by_gene.df$fdr.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$fdr.log2cpm<0.1, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.01...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.01, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.01, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.05...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.05, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.05, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm

cat('\n')

print('p-nominal < 0.1...')
table(lmer.by_gene.df$p.value.voomWithDreamWeights<0.1, lmer.by_gene.df$estimate.voomWithDreamWeights>0) #voomWithDreamWeights
table(lmer.by_gene.df$p.value.log2cpm<0.1, lmer.by_gene.df$estimate.log2cpm>0) #log2cpm
