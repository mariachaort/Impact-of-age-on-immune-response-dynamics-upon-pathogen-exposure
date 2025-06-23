#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/home/mchacon/Desktop/cluster/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
  .libPaths(c(.libPaths(), "/gpfs/apps/MN5/GPP/R/4.3.0/GCC/lib64/R/library"))
  #.libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1/INTEL/lib64/R/library"))
} else if(file.exists(path_em)){
  setwd(paste(path_em))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--chemistry"), action="store", default="V2", type="character",
              help="Chemistry version used"),
  make_option(c("--condition"), action="store", default="PA", type="character",
              help="Stimuli used: PA, CA, MTB or AllStim"),
  make_option(c("--cell_level"), action="store", default='cell_type', type='character',
              help="L2 --> cell_type or L1 -> cell_type_lowerres"),
  make_option(c("--cell_type"), action="store", default="memory_CD8T", type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--aggr_fun"), action="store", default="sum", type='character',
              help="Aggregation function."),
  make_option(c("--nCells_per_donor"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE"),
  make_option(c("--permutation_idx"), action="store", default=NULL, type='integer',
              help="Iteration index, set a seed() with this value."),
  make_option(c("--phenotypes"), action="store", default="Projects/scRNAseq/mchacon/1M-scBloodNL/scripts/analysis/pseudobulk_dreamlet.phenotypes.tab", type='character',
              help="Gender/age/Stim.cond"),
  make_option(c("--covs"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/scripts/analysis/pseudobulk_dreamlet.covariates.tab',  type='character',
              help="Covariates file."),
  make_option(c("--random"), action="store", default="assignment,date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--min_prop"), action="store", default=0.2, type='double',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained. If sex !is.null, then 0.4."),
  make_option(c("--vp_reduced"), action="store", default=FALSE, type='logical',
              help="Not estimate the effect of the batch factor in the VariancePartition analysis."),
  make_option(c("--in_dir"), action="store", default='Data/scRNAseq/Oelen2022', type='character',
              help="Input directory"),
  make_option(c("--out_dir"), action="store", default='Projects/scRNAseq/mchacon/1M-scBloodNL/pseudobulk_dreamlet', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(SingleCellExperiment))
shhh(library(dreamlet))
shhh(library(pkgload))
shhh(library(zenith))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyverse))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(Seurat))
shhh(library(SeuratObject))

################################## Set Variables and load Data #################################
# Input sc-file
in.dir <- opt$in_dir
out.dir <- opt$out_dir
cond <- opt$condition
chem_version <- opt$chemistry
ctype <- opt$cell_type
cell_level <- opt$cell_level
if (cell_level == "cell_type") {
  cl <- "L2"
} else {
  cl <- "L1"
}

in.fn <- paste0(in.dir, "/", chem_version , "/", cond, "/", ctype, "_", cl, "_", cond, "_Alltimes_so.rds")


# Phenotypes
phenotypes_fn <- opt$phenotypes
phenotypes_fn <- paste0(path_cluster, phenotypes_fn)
phenotypes <- read.table(phenotypes_fn)$V1
phenotypes_tag <- paste(phenotypes, collapse='_')

# Covariates
covs_fn <- opt$covs
covs_fn <- paste0(path_cluster, covs_fn)
covs_df <- read.table(covs_fn, header=TRUE)
print("Reading info...")
# Output directory
out.dir <- paste0(out.dir, '/', opt$aggr_fun, '/')
if(!is.null(opt$permutation_idx)){out.dir <- paste0(out.dir, 'cell_permutations/')}

if(opt$nCells_per_donor){out.dir <- paste0(out.dir, 'nCells_per_donor/')}

out.dir <- paste0(out.dir, chem_version, '/', 
                  cl, '/', opt$cell_type, '/',  
                  phenotypes_tag, '/', 
                  'min_prop.', as.character(opt$min_prop), '/')

if(opt$vp_reduced){out.dir <- paste0(out.dir, 'vp_reduced/')}

if(!is.null(opt$permutation_idx)){out.dir <- paste0(out.dir, 'it_', as.character(opt$permutation_idx), '/')}
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
out.dir <- paste0(path_cluster, out.dir)
print(paste0('Main output directory: ',out.dir))

# Report
print(paste0('Cell level: ', opt$cell_level, ",", cl))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Aggregate function: ', opt$aggr_fun))
print(paste0('Cell-level permutations (iteration): ',  as.character(opt$permutation_idx)))
print(paste0('Random: ', opt$random))
print(paste0('Covariates: ', opt$covs))

################################## Functions ################################## 
# Acessory functions
# 1. Read input data
# in_fn <- in.fn
# permutation_idx = opt$permutation_idx
read_data <- function(in_fn = in.fn, permutation_idx = opt$permutation_idx){
  print(paste0('Reading SCE file in: ', in_fn))
  system.time(seurat.obj <- readRDS(in_fn))
  sce <- as.SingleCellExperiment(seurat.obj)
  if(ncol(sce)>0){
    print(paste0('Number of donors: ', length(unique(sce$assignment))))
    print(paste0('Number of cells: ', ncol(sce)))
    print('Number of cells per donor: ')
    print(summary(as.vector(table(sce$assignment))))
  }else{
    print('No cells remains.')
    sce <- NULL
  }
  return(sce)
}
# 2. Define formula (VP or DEA)
#gt <- "DEA"
#df <- covs_df
#vp <- opt$vp_reduced
define_form <- function(gt, df, vp){
  # forms
  print(gt)
  random_var_dea <- df[df$DEA=='random',]$covariate
  if(vp){
    if(gt=='VP'){
      df[df$DEA=='random',]$VP <- NA
    }
  }
  cnames <- c('covariate',gt)
  df_i <- df[,cnames]
  colnames(df_i)[2] <- 'type'
  fixed_var <- df_i[!is.na(df_i$type) & df_i$type=='fixed',]$covariate
  fixed.fmla <- NULL
  if(length(fixed_var)>0){
    fixed.fmla <- paste(fixed_var,collapse='+')
  }
  
  random_var <- df_i[!is.na(df_i$type) & df_i$type=='random',]$covariate
  random.fmla <- NULL
  if(length(random_var)>0){
      random.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
  }
  form_vars <- paste(c(fixed.fmla,random.fmla), collapse='+')
  form_vars <- paste0('~',form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  form <- as.formula(form_vars)
  
  # specificy colors
  model_vars <- c(fixed_var, random_var)
  Gender.hex <- brewer.pal(9, 'Greens')[7]
  age.hex <- brewer.pal(9, 'Blues')[7]
  #nCells.hex <- brewer.pal(9, 'Reds')[7]
  random.hex <- brewer.pal(9, 'Greys')[seq_len(length(random_var_dea))]
  Residuals.hex <- brewer.pal(9, 'Greys')[3]
  cols_vars <- c(Gender.hex, age.hex, random.hex, Residuals.hex)
  names(cols_vars) <- c('Gender', 'age', random_var_dea, 'Residuals')
  cols_vars.in <- c(model_vars, 'Residuals')
  cols_vars <- cols_vars[names(cols_vars)%in%cols_vars.in]
  
  # output
  out <- list(form = form,
              cols = cols_vars)
  
  return(out)
}
# 3. DEA extract and plots
# i <- phe
# dea_res <- dea_res
# contrast_var <- contrast_coefName
# vp_res <- vp_res
# cols <- cols
# phe_dir <- out_sdir
extract_plots <- function(i, dea_res, contrast_var, phe_dir){
  # Extract results (topTable --> DEGs)
  ### Each entry in res.dl stores a model fit by dream(), and results can be extracted using topTable() as in limma by specifying the coefficient of interest. 
  ### The results shows the gene name, log fold change, average expression, t-statistic, p-value, FDR (i.e. adj.P.Val).
  genes <- rownames(dea_res[[1]]$residuals)
  all_topTables <- list()
  for (coef in contrast_var) {
    print(paste("Using contrast:", coef))
    
    topTable.res <- topTable(dea_res, coef = coef, number = length(genes))
    all_topTables[[coef]] <- topTable.res
    
    degs <- topTable.res[topTable.res$adj.P.Val <= 0.05, ]$ID
    
    plotVolcano.p <- plotVolcano(dea_res, coef = coef)
    plotVolcano.fn <- paste0(phe_dir, coef, '_Alltimes_plotVolcano.png')
    print(paste0('Saving plotVolcano in: ', plotVolcano.fn))
    ggsave(plotVolcano.fn, plotVolcano.p)
    
    if (length(degs) > 0) {
      print(paste0('# of DEGs: ', length(degs)))
      plotGeneHeatmap.p <- plotGeneHeatmap(dea_res, coef = coef, genes = degs)
      plotGeneHeatmap.fn <- paste0(phe_dir, coef, '_Alltimes_plotGeneHeatmap.png')
      print(paste0('Saving plotGeneHeatmap in: ', plotGeneHeatmap.fn))
      ggsave(plotGeneHeatmap.fn, plotGeneHeatmap.p)
    }
  }
  
  return(all_topTables)  # returns list with all toptables 
  
} 

# 4. DEA + VP extract and plots (by phenotype)
# phe <- names(contrast_list)[3]
# dea_res <- res.dl
# vp_res <- vp.lst
# c_list <- contrast_list
# cols <- cols_vars
# o_dir <- out_dir
extract_plots_by_phe <- function(phe, dea_res, c_list, o_dir){
  print(phe)
  
  # create output dir
  out_sdir <- paste0(o_dir, '/', phe, '/')
  if(!dir.exists(out_sdir)){dir.create(out_sdir, recursive = T)}
  
  # pick contrast variable
  contrast_coefName <- c_list[[phe]]
  
  # extract and plots
  extract_plots.res <- extract_plots(i = phe,
                                     dea_res = dea_res,
                                     contrast_var = contrast_coefName,
                                     #vp_res = vp_res,
                                     #cols = cols,
                                     phe_dir = out_sdir)
  return(extract_plots.res)
  
}
# Main function
# 5. dreamlet
# ge_dge = pb
#covariates = covariates_df
# contrast_list = contrast_coefName.list #for DEA
#gene_test = c('VP','DEA')
#vp_reduced = opt$vp_reduced
# out_dir = out.dir
dreamlet.func <- function(ge_dge, covariates, contrast_list, gene_test = c('VP','DEA'), vp_reduced = opt$vp_reduced, out_dir = out.dir){
  ### Defining the VP/DEA formulas ###
  print('Defining the VP/DEA formulas...')
  gene_test.forms <- sapply(gene_test, function(i) define_form(i, covariates, vp_reduced), simplify = FALSE)
  
  #### Normalize and apply voom/voomWithDreamWeights ####
  # Run processAssays()
  form <- gene_test.forms$DEA$form
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(sceObj = ge_dge, 
                                        formula = form,
                                        min.cells = 5,
                                        min.count = 5,
                                        min.samples = 4,
                                        min.prop=opt$min_prop))

  # View details of dropping samples
  details(res.proc)

  # Check nSamples and nGenes tested
  genes_all <- rownames(ge_dge)
  genes_tested <- rownames(as.data.frame(res.proc))
  genes_all.n <- nrow(ge_dge)
  genes_tested.n <- nrow(as.data.frame(res.proc))
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(ge_dge)
  samples_tested <- colnames(as.data.frame(res.proc))
  samples_all.n <- ncol(ge_dge)
  samples_tested.n <- ncol(as.data.frame(res.proc))
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))

  # Show voom plot for each cell clusters
  ## Here the mean-variance trend from voom is shown for each cell type. Cell types with sufficient number of cells and reads show a clear mean-variance trend. 
  ## While in rare cell types like megakaryocytes, fewer genes have sufficient reads and the trend is less apparent.
  plotVoom.p <- plotVoom(res.proc)
  plotVoom.fn <- paste0(out_dir, cond, '_Alltimes_plotVoom.png')
  ggsave(plotVoom.fn, plotVoom.p)

  print('Running DEA...')
  system.time(res.dl <- dreamlet(res.proc, form))

  ### Extract results and plots (DEA and VP) ###
  # phe <- names(contrast_list)[3]
  # dea_res <- res.dl
  # vp_res <- vp.lst
  # c_list <- contrast_list
  # cols <- cols_vars
  # o_dir <- out_dir
  extract_plots_by_phe.res <- sapply(names(contrast_list),
                                     function(i) extract_plots_by_phe(phe = i,
                                                                      dea_res = res.dl,
                                                                      #vp_res = vp.lst,
                                                                      c_list = contrast_list,
                                                                      #cols = cols_vars,
                                                                      o_dir = out_dir), simplify = FALSE)

  ### Save outputs ###
  out <- list(processed = res.proc,
              dea = res.dl,
              topTable = extract_plots_by_phe.res)
  out_fn <- paste0(out_dir, cond, '_Alltimes_dea_topTable.rds')
  print(paste0('Saving dreamlet results: ',out_fn))
  saveRDS(out, out_fn)
  
  return(out)
}

################################## Analyses #################################### 
# Read input data
system.time(sce <- read_data(in.fn))
sce$Stim.cond <- paste0(sce$condition, sce$time, "h")
sce$Donor <- paste0(sce$Stim.cond, "_", sce$assignment)
# Check
if(is.null(sce)){stop('No cells in the SCE object.')}

## Add nCells per donor in the model 
#################### I had to change a little this part of the code, as Aida seems to have in her dfs some variables that I don't (bare_barcode_lane)
################### so I adapted this part according to Oelen's SCEs
if(opt$nCells_per_donor){
  print('Calculating the nCells per donor...')
  ct_metadata <- as.data.frame(colData(sce))
  
  ct_metadata %>% 
    group_by(assignment) %>% 
    count(name = 'nCells') %>% 
    arrange(desc(nCells)) %>% as.data.frame() -> ncells_per_donor.df
  ct_metadata.added <- left_join(ct_metadata, ncells_per_donor.df, by = 'assignment')
  
  rownames(ct_metadata.added) <- ct_metadata.added$barcode
  print(summary(ct_metadata.added$nCells))
  
  print('Adding nCells per donor in the SCE ctcell metadata (colData(sce))...')
  # Check that the rownames of `ct_metadata.added` match the colnames of the SCE object
  if (!all(rownames(ct_metadata.added) %in% colnames(sce))) {
    print("Row names of ct_metadata.added do not match SCE colnames --> Reordering...")
    
    # Ensure the order matches the SCE's column names
    ct_metadata.added <- ct_metadata.added[match(colnames(sce), rownames(ct_metadata.added)), ]
  }
  
  # Add the new metadata to the SCE object's colData
  colData(sce)$nCells <- ct_metadata.added$nCells
}
stim_levels <- list(
  CA = c('UT0h', 'CA3h', 'CA24h'),
  PA = c('UT0h', 'PA3h', 'PA24h'),
  MTB = c('UT0h', 'MTB3h', 'MTB24h')
)
if (!opt$condition %in% names(stim_levels)) stop("Pathogen not recognized.")
## Contrast order
Group_order.vec <- list(
  age = 'age',
  Gender = c('M', 'F'),
  Stim.cond = stim_levels[[opt$condition]]
)

if('Gender'%in%phenotypes){
    colData(sce)$Gender <- as.factor(colData(sce)$Gender)
    colData(sce)[['Gender']] <- factor(colData(sce)[['Gender']],
                                       levels = Group_order.vec[['Gender']])
}
colnames(colData(sce))
colData(sce)$Stim.cond <- factor(colData(sce)$Stim.cond, levels = Group_order.vec[['Stim.cond']])
colnames(colData(sce))[colnames(colData(sce))==opt$cell_level] <- 'celltype'
cnames <- c('Donor', 'celltype', 'assignment', 'date', 'Gender', 'age', 'Stim.cond')
if(opt$nCells_per_donor){cnames <- c(cnames, 'nCells')}
cnames <- unique(c(cnames, phenotypes))
colData(sce) <- colData(sce)[,cnames]

# Aggregate to pseudobulk
## Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
## aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
print('Performing pseudobulk...')
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "celltype", 
                                        sample_id = "Donor",
                                        fun = opt$aggr_fun,
                                        verbose = FALSE))
pb_raw <- pb
pb_fn <- paste0(out.dir, cond, '_Alltimes_pb.Rds')
saveRDS(pb, pb_fn)

##########################################################
# Define contrasts by phenotype
model_vars <- c(phenotypes, "assignment", "date")
if(opt$nCells_per_donor){model_vars <- c(model_vars, 'nCells')}

covariates_df <- covs_df[covs_df$covariate%in%model_vars,]
contrast_coefName.list <- sapply(names(Group_order.vec), function(i){
  contrast_var <- Group_order.vec[[i]]
  if(length(contrast_var) > 1){
    # Instead of taking only the second level, take all except the first (reference)
    contrast_var <- paste0(i, contrast_var[-1]) # I changed this to make the comparisons of 3h vs 0h and 24h vs 0h 
  }
  return(contrast_var)
}, simplify = FALSE)

contrast_coefName.list <- contrast_coefName.list[names(contrast_coefName.list)%in%phenotypes]

# Run dreamlet (DEA, VP and extract results/plots)
print('Running dreamlet...')
# ngenes <- 50
# ngenes <- ifelse(nrow(pb)>=ngenes, ngenes, nrow(pb))
# pb <- pb[1:ngenes,]
system.time(dreamlet.res <- dreamlet.func(ge_dge = pb,
                                          covariates = covariates_df,
                                          contrast_list = contrast_coefName.list))

# Check sessionInfo()
print(sessionInfo())
cat('\n')
print('#####################################################################')
cat('\n')

# Check bias
# phe <- phenotypes[2]
check_bias <- function(phe){
  print(paste0('################### ', phe, ' ###################'))
  table_fdr0.05 <- table(dreamlet.res$topTable[[phe]]$adj.P.Val<0.05, dreamlet.res$topTable[[phe]]$logFC>0)
  table_fdr0.1 <- table(dreamlet.res$topTable[[phe]]$adj.P.Val<0.1, dreamlet.res$topTable[[phe]]$logFC>0)
  table_pval0.01 <- table(dreamlet.res$topTable[[phe]]$P.Value<0.01, dreamlet.res$topTable[[phe]]$logFC>0)
  table_pval0.05 <- table(dreamlet.res$topTable[[phe]]$P.Value<0.05, dreamlet.res$topTable[[phe]]$logFC>0)
  
  table_list <- list(fdr0.05 = table_fdr0.05, 
                     fdr0.1 = table_fdr0.1,
                     pval0.01 = table_pval0.01,
                     pval0.05 = table_pval0.05)
  # i <- names(table_list)[3]
  check_res <- sapply(names(table_list), function(i){
    print(paste0('### ', i, ' ####'))
    table.ss_logFC <- table_list[[i]]
    print(table.ss_logFC)
    cat('\n')
    if('TRUE'%in%rownames(table.ss_logFC)){
      ss_down.bias <- binom.test(unname(table.ss_logFC[2,]))
      print('Sign (down):')
      print(ss_down.bias)
      cat('\n')
      print('Sign (up):')
      ss_up.bias <- binom.test(rev(unname(table.ss_logFC[2,])))
      print(ss_up.bias)
    }
    cat('\n')
    cat('\n')
  }, simplify = FALSE)
}
check_bias.list <- lapply(phenotypes, function(i) check_bias(i))
