# Suppl_5b_contingency_chisquare_genetics

``` r
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
```

Load genetic information and cluster assignments, then merge for
analysis

``` r
all_genetic_info <- read_table("data/ALL_selected_genetic_information_used_for_heatmap_03APRILE.txt")
c1 <-  read.csv("../Fig4/data/cluster_motif_1_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", sep = "\t", header = F) %>% mutate(cluster="C1")
c2 <- read.csv("../Fig4/data/cluster_motif_2_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", sep = "\t", header = F) %>% mutate(cluster="C2")

df_cluster <- rbind(c1,c2)

df_contingency <- merge(df_cluster,all_genetic_info,by.x = "V1", by.y = "Patient_ID") %>% 
    mutate(plusplus15_INFO=case_when(plusplus15_INFO=="plusplus15_YES"  ~ 1, TRUE ~ 0)) %>%
    mutate(MYC_STR=case_when(MYC_STR=="MYC_STR_YES" ~ 1, TRUE ~ 0)) %>%
  mutate(PR_index=case_when(PR_index=="PR_Q3_YES" ~ 1, TRUE ~ 0)) 
```

Split the merged dataframe by cluster for separate analyses

``` r
df_contingency_C1 <- df_contingency %>% filter(cluster=="C1") 
df_contingency_C2 <-   df_contingency %>% filter(cluster=="C2")
```

Define the genetic features to analyze and define a function to analyze
each genetic feature: build contingency table, perform chi-square test,
and create balloon plot

``` r
# Define the vector of column names to check
column_vector <- c("traslocationNSD2_CALL", "traslocationCCND3_CALL", "traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL",
                   "traslocationCCND2_CALL", "traslocationMAF_CALL", "traslocationMAFB_CALL", "mutationTP53", "mutationFGFR3", "mutationDIS3", 
                   "mutationATM", "mutationKRAS", "mutationBRAF", "mutationHIST1H1E", "mutationLTB")

plot_list <- list()  # For the balloon plots
summary_df <- data.frame(Column = character(), ChiSquare = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

# Function to compute contingency table and plot for one column
analyze_column <- function(column_name, df_contingency_C1, df_contingency_C2) {
  
  if (!(column_name %in% colnames(df_contingency_C1)) | !(column_name %in% colnames(df_contingency_C2))) {
    message(paste("Column", column_name, "not found in both dataframes. Skipping."))
    return(NULL)
  }
  # Build contingency table
  contingency_table <- matrix(c(
    sum(df_contingency_C1[[column_name]] == 1, na.rm = TRUE), sum(df_contingency_C1[[column_name]] == 0, na.rm = TRUE),
    sum(df_contingency_C2[[column_name]] == 1, na.rm = TRUE), sum(df_contingency_C2[[column_name]] == 0, na.rm = TRUE)
  ), nrow = 2, byrow = TRUE)
  colnames(contingency_table) <- c("Present", "Absent")
  rownames(contingency_table) <- c(paste0("C1_", column_name), paste0("C2_", column_name))
  #print(paste("Contingency table for", column_name))
  #print(contingency_table)

  # Perform chi-square test
  chisq_result <- chisq.test(contingency_table)
  #print(paste("Chi-square test p-value for", column_name, ":", chisq_result$p.value))
  
    # Save column results to summary dataframe
  summary_df <<- rbind(summary_df, data.frame(
    Column = column_name,
    ChiSquare = chisq_result$statistic,
    PValue = chisq_result$p.value
  ))
  
  # Plot using ggballoonplot
  plot <- ggballoonplot(contingency_table, fill = "value", color = "lightgray", 
                        size = "value", show.label = TRUE) +
    gradient_fill(c("#1C90BF", "white", "#EB8A3D")) +
    border() +
    ggtitle(paste0(column_name, ": p-value = ", round(chisq_result$p.value, 3))) +
    scale_x_discrete(labels = c("Present", "Absent")) + 
    theme(legend.position = "none")
    plot_list[[column_name]] <<- plot
}

# Loop through the vector of columns and run the analysis
results <- lapply(column_vector, analyze_column, 
                  df_contingency_C1 = df_contingency_C1, 
                  df_contingency_C2 = df_contingency_C2)

plots_grid <- grid.arrange(grobs = plot_list, ncol = 4)  
```

![](Suppl_5b_contingency_chisquare_genetics_files/figure-markdown_github/unnamed-chunk-4-1.png)
![Suppl_5b_contingency_chisquare_genetics](https://github.com/cleliacort/NRF1_paper/blob/main/Fig5/figures/Suppl/contingency_table_between_C1_and_C2_for_all_genetics__09APRILE25.png)
