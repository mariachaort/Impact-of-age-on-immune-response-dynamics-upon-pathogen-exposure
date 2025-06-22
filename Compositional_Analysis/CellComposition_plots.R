# Cell type composition changes - Figures 

library(ggplot2); library(dplyr); library(RColorBrewer);library(ggside);library(ggbeeswarm);library(ggpubr);library(ggrepel)


#paths 
if (getwd()== "/home/mchacon/"){
  basepath <- "/home/mchacon/Desktop/cluster/"
  data_path <- "/home/mchacon/Desktop/cluster/Projects/scRNAseq"
}else{
  basepath <- "/gpfs/projects/bsc83/"
}
datapath <- paste0(basepath, "Data/scRNAseq/Oelen2022/") 
robejctspath <- paste0(basepath,"Projects/scRNAseq/msopena/04_Stimulation_Age_Sex/robjects/")

#functions
source(paste0(basepath, "Projects/scRNAseq/msopena/04_Stimulation_Age_Sex/scripts/utils.R"))
source(paste0(data_path, "/mchacon/1M-scBloodNL/scripts/Themes.R"))

# get order of cells 
metadata <- readRDS(paste0(datapath, "1M_full_metadata_modified.rds")) %>% dplyr::group_by(cell_type, chem) %>% dplyr::count()
order_v2 <- metadata %>% filter(chem=="V2") %>% arrange(-n) %>% pull(cell_type)
order_v3 <- metadata %>% filter(chem=="V3") %>% arrange(-n) %>% pull(cell_type)

# Set parameters ---
stim <- c("CA", "MTB", "PA", "UT")
time <- c("0h","3h", "24h")
chem <- c("V2", "V3")

# Create a data frame of all combinations of parameters
param_combinations <- expand.grid( stim = stim, time = time, chem = chem)

# Filter the combinations based on the conditions:
# For stim == "UT", remove rows where time is 3h or 24h
param_combinations <- param_combinations[!(param_combinations$stim == "UT" & param_combinations$time %in% c("3h", "24h")), ]
# For stim in "MTB", "CA", "PA", remove rows where time is "0h"
param_combinations <- param_combinations[!(param_combinations$stim %in% c("MTB", "CA", "PA") & param_combinations$time == "0h"), ]

coda <- do.call(rbind, lapply(1:nrow(param_combinations), function(i) {
  # Extract parameters for this combination
  stim_val <- param_combinations$stim[i]
  time_val <- param_combinations$time[i]
  chem_val <- param_combinations$chem[i]
  
  # Set file path 
  file_name <- paste0(robejctspath, "/02_CellComposition/", chem_val,"/CellType_Composition_", stim_val,"_", time_val, ".rds")
  
  # Read file and add the parameters as columns
  df <- tryCatch({
    df <- readRDS(file_name)
    df$stim <- stim_val
    df$time <- time_val
    df$chem <- chem_val
    return(df)
  }, error = function(e) {
    message("Error reading file: ", file_name)
    return(NULL) # Return NULL if there's an error
  })
  
  return(df)
}))
alpha_vec <- c("ss" = 1, "ns" = 0.2)

plot_coda <- function(chemistry, order_cells, pheno){
  coda_chem <- coda[coda$chem == chemistry ,]
  coda_chem <- coda_chem[coda_chem$term == pheno, ]
  coda_chem$celltype <- factor(coda_chem$celltype, levels=rev(order_cells))
  coda_chem$significance <- ifelse(coda_chem$p.value< 0.05, "ss","ns")
  # Replace _ with space for nicer labels
  celltype_labels <- gsub("_", " ", levels(coda_chem$celltype))
  names(celltype_labels) <- levels(coda_chem$celltype)
  coda_chem[coda_chem$significance== "ss", ]
  coda_chem <- coda_chem[coda_chem$celltype != "unknown", ]
  coda_chem[coda_chem$stim == "UT", ]$time <-"3h"
  coda_chem$stim <- factor(coda_chem$stim, levels=c("UT", "CA", "MTB", "PA"))
  plot <- ggplot(coda_chem, aes(x=estimate, y=celltype)) + 
    geom_point(aes(alpha=significance,color=stim), size=4) + xlab("CoDA estimate")+ylab(NULL)+
    geom_pointrange(aes(xmin=conf.low, xmax=conf.high, alpha=significance, color=stim), fatten = .1) +
    geom_vline(xintercept=0, linetype = "dashed")  + 
    theme_paper + theme(axis.text.x = element_text(size = 14), legend.position="top", legend.title=element_blank(),legend.key.size = unit(0.5, "lines") ) +
    scale_fill_manual(values = pal_stimulation)+scale_color_manual(values=pal_stimulation)+
    scale_alpha_manual(values=alpha_vec) +facet_grid(time~stim, space="free", scale="free") + scale_y_discrete(labels = celltype_labels)
  return(plot)
}

plot <- plot_coda("V2", order_v3, "age")                      
plot_coda("V2", order_v2, "GenderF")   

fig_dir <- paste0(data_path, "/mchacon/1M-scBloodNL/Figures/V2/")

pdf(paste0(fig_dir, "CoDA_Sex_V2.pdf"), width = 17, height = 14)
plot
dev.off()
