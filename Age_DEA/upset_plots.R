library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)
library(scales)

pheno <- "age"; cell_level <- "L2"
df <- readRDS(paste0(data_path, "/mchacon/1M-scBloodNL/robjects/final_dataframes/", chem_version, "/",
                     pheno, "_", chem_version, "_", cell_level, "_nCells_log2cpm_dframe.rds"))

df <- df[df$adj.P.Val < 0.05, ]
df$direction <- ifelse(df$logFC > 0, "up", "down")
# Filtrar por tipo celular
cell <- "naive_CD8T"; d <- "up"
df_cell <- df[df$celltype == cell & df$direction == d, ]
df_cell$timepoint <- paste0(df_cell$time, df_cell$condition)
# Create new df
new_df <- df_cell %>%
  dplyr::select(timepoint, ID) %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = timepoint, values_from = Presence, values_fill = list(Presence = 0))


# Specify conditions
conditions <- colnames(new_df)[2:8]

upset_plot <- upset(
  new_df, conditions, min_size = 4, wrap = TRUE, name = "conditions",
  base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE)
  ),
  queries = list(
    upset_query(intersect = c('0hUT'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['0hUT'], color = pal_tp['0hUT']),
    upset_query(intersect = c('3hCA'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['3hCA'], color = pal_tp['3hCA']),
    upset_query(intersect = c('24hCA'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['24hCA'], color = pal_tp['24hCA']),
    upset_query(intersect = c('3hPA'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['3hPA'], color = pal_tp['3hPA']),
    upset_query(intersect = c('24hPA'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['24hPA'], color = pal_tp['24hPA']),
    upset_query(intersect = c('3hMTB'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['3hMTB'], color = pal_tp['3hMTB']),
    upset_query(intersect = c('24hMTB'), only_components = c("intersections_matrix", "Intersection size"), 
                fill = pal_tp['24hMTB'], color = pal_tp['24hMTB'])
    )
)

print(upset_plot)

pdf(paste0(data_path, "/mchacon/1M-scBloodNL/Figures/", chem_version, 
           "/", "upset_FIG_S2_", d, "_", pheno, "_", cell_level, ".pdf"), width =17, height = 8)
upset_plot
dev.off()


