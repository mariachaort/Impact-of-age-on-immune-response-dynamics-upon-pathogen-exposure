####################### FOR L2 ###########################
for (pheno in phenotypes) {
  pheno <- "age"
  degs <- get_degs(pheno, cell_level, chem_version)
  degs <- degs  %>% dplyr::filter(!celltype %in% c("megakaryocyte", "unknown", "hemapoietic stem", "th1 CD4T", "th2 CD4T"))
  degs$stim.con <- paste(degs$time, degs$condition, sep ="")
  unique(degs$celltype)
  library(dplyr)
  degs <- degs %>% filter(!is.na(ID))
  
  df_final <- degs

  df_final <- df_final %>%
    dplyr::group_by(celltype, direction, condition) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n) * 100) %>% mutate(n = ifelse(direction == "down", -n, n)) %>%
    dplyr::mutate(vjust_label = ifelse(direction == "down", 1.5, -0.5))  # Adjust vjust based on direction
  
  #########################################################################################
  ############### BINOMIAL TEST FOR UP/DOWN BIAS #######################
  
  result_df <- data.frame(celltype = character(0), stim.cond = character(0), binom.test.p = numeric(0))
  stim.conds <- c("0hUT", "3hCA", "24hCA", "3hPA", "24hPA", "3hMTB", "24hMTB")
  
  for (cell in unique(df_final$celltype)) {
    for (sc in stim.conds) {
      degs_filtered <- degs[degs$celltype == cell & degs$stim.cond == sc, ]
      # Calculate the number of down DEGs and total DEGs
      x <- length(degs_filtered[degs_filtered$direction == "down", ]$ID) ; n <- length(degs_filtered$ID)  # Total number of DEGs
      
      # Perform the binomial test only if n > 0 and n >= x
      if (n > 0 && n >= x) {
        binom <- binom.test(x, n)$p.value
      } else {
        binom <- NA  # If not valid, assign NA
      }
      result_df <- rbind(result_df, data.frame(celltype = cell, stim.cond = sc, binom.test.p = binom))
    }
  }
  result_df <- semi_join(result_df, df_final, by = c("celltype", "stim.cond"))
  df_final <- df_final %>% left_join(result_df, by = c("celltype", "stim.cond"))
  
  df_final$signif <- ifelse(df_final$binom.test.p< 0.01, "*", "")
  
  for (cell in unique(df_final$celltype)) {
    for (sc in stim.conds) {
      df_subset <- df_final[df_final$celltype == cell & df_final$stim.cond == sc, ] # check if the subset exists
      
      if (nrow(df_subset) > 0 && length(df_subset$celltype) == 2 && 
          any(df_subset$direction == "down" & df_subset$signif == "*")) {
        
        df_cell <- df_final[df_final$celltype == cell & df_final$stim.cond == sc, ]
        
        if (abs(df_cell[df_cell$direction == "up", ]$n) < abs(df_cell[df_cell$direction == "down", ]$n)) {
          df_final[df_final$celltype == cell & df_final$stim.cond == sc & df_final$direction == "up", ]$signif <- ""
        } else {
          df_final[df_final$celltype == cell & df_final$stim.cond == sc & df_final$direction == "down", ]$signif <- ""
        }
      }
    }
  }
  df_final$n_signif <- paste(abs(df_final$n), df_final$signif, sep = " ")
  df_final <- reorder_celltypes(df_final, cell_level ="L2", reverse = T) 
  df_final$condition <- factor(df_final$condition, levels = c('UT', 'CA', 'PA', 'MTB'))
  df_final$stim.cond <- factor(df_final$stim.cond, levels = c("0hUT", "3hCA", "24hCA", "3hPA", "24hPA", "3hMTB", "24hMTB"))
  
  #pheno <- if (pheno == "Gender" || pheno == "Sex") "Sex" else "Age"
  plot_cond <- ggplot(df_final, aes(x = n, y = stim.cond, fill = stim.cond)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = abs(n), hjust = ifelse(n > 0, -0.1, 1.1)), size = 4.5) +
    facet_wrap(~ celltype, ncol = 1, strip.position = "right") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = contrast_palette) +
    labs(x = "Number of significant DEGs", y = NULL) +
    theme_paper +
    theme(
      strip.text.y.right = element_text(angle = 0, hjust = 0, size = 17),
      strip.placement = "outside",
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.spacing = unit(0.4, "lines"),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(-max(abs(df_final$n)) * 1.1, max(abs(df_final$n)) * 1.1))  # Add margin on x-axis
  
  print(plot_cond)
  
  
  filename <- paste0(data_path, "mchacon/1M-scBloodNL/Figures/", chem_version, 
                     "/", "DEGS_S4_", pheno, "_", cell_level, "_Plot.png")
  ggsave(filename, plot = plot_cond, bg = "white", width = 20, height = 9, dpi = 250)
  
  
}
