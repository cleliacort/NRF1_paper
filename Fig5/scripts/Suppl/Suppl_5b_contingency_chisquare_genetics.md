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

# Define the genetic features to analyze and define a function to analyze each genetic feature: build contingency table, perform chi-square test, and create balloon plot

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
  print(paste("Contingency table for", column_name))
  print(contingency_table)

  # Perform chi-square test
  chisq_result <- chisq.test(contingency_table)
  print(paste("Chi-square test p-value for", column_name, ":", chisq_result$p.value))
  
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
```

    ## [1] "Contingency table for traslocationNSD2_CALL"
    ##                          Present Absent
    ## C1_traslocationNSD2_CALL      59    320
    ## C2_traslocationNSD2_CALL      23    236
    ## [1] "Chi-square test p-value for traslocationNSD2_CALL : 0.0183776779985601"
    ## [1] "Contingency table for traslocationCCND3_CALL"
    ##                           Present Absent
    ## C1_traslocationCCND3_CALL       3    376
    ## C2_traslocationCCND3_CALL       5    254
    ## [1] "Chi-square test p-value for traslocationCCND3_CALL : 0.364224480467375"
    ## [1] "Contingency table for traslocationMYC_CALL"
    ##                         Present Absent
    ## C1_traslocationMYC_CALL      38    341
    ## C2_traslocationMYC_CALL      54    205
    ## [1] "Chi-square test p-value for traslocationMYC_CALL : 0.00020989839290892"
    ## [1] "Contingency table for traslocationMAFA_CALL"
    ##                          Present Absent
    ## C1_traslocationMAFA_CALL       3    376
    ## C2_traslocationMAFA_CALL       2    257
    ## [1] "Chi-square test p-value for traslocationMAFA_CALL : 0.999999999999997"
    ## [1] "Contingency table for traslocationCCND1_CALL"
    ##                           Present Absent
    ## C1_traslocationCCND1_CALL      81    298
    ## C2_traslocationCCND1_CALL      49    210
    ## [1] "Chi-square test p-value for traslocationCCND1_CALL : 0.512240049781944"
    ## [1] "Contingency table for traslocationCCND2_CALL"
    ##                           Present Absent
    ## C1_traslocationCCND2_CALL       5    374
    ## C2_traslocationCCND2_CALL       5    254
    ## [1] "Chi-square test p-value for traslocationCCND2_CALL : 0.774978312613853"
    ## [1] "Contingency table for traslocationMAF_CALL"
    ##                         Present Absent
    ## C1_traslocationMAF_CALL      11    368
    ## C2_traslocationMAF_CALL      15    244
    ## [1] "Chi-square test p-value for traslocationMAF_CALL : 0.107692797729835"
    ## [1] "Contingency table for traslocationMAFB_CALL"
    ##                          Present Absent
    ## C1_traslocationMAFB_CALL       6    373
    ## C2_traslocationMAFB_CALL       5    254
    ## [1] "Chi-square test p-value for traslocationMAFB_CALL : 0.982961142590471"
    ## [1] "Contingency table for mutationTP53"
    ##                 Present Absent
    ## C1_mutationTP53      13    366
    ## C2_mutationTP53      33    226
    ## [1] "Chi-square test p-value for mutationTP53 : 1.6368265000644e-05"
    ## [1] "Contingency table for mutationFGFR3"
    ##                  Present Absent
    ## C1_mutationFGFR3      15    364
    ## C2_mutationFGFR3       6    253
    ## [1] "Chi-square test p-value for mutationFGFR3 : 0.360158401307681"
    ## [1] "Contingency table for mutationDIS3"
    ##                 Present Absent
    ## C1_mutationDIS3      42    337
    ## C2_mutationDIS3      19    240
    ## [1] "Chi-square test p-value for mutationDIS3 : 0.14901823521252"
    ## [1] "Contingency table for mutationATM"
    ##                Present Absent
    ## C1_mutationATM      17    362
    ## C2_mutationATM       4    255
    ## [1] "Chi-square test p-value for mutationATM : 0.068942781353768"
    ## [1] "Contingency table for mutationKRAS"
    ##                 Present Absent
    ## C1_mutationKRAS      95    284
    ## C2_mutationKRAS      65    194
    ## [1] "Chi-square test p-value for mutationKRAS : 1"
    ## [1] "Contingency table for mutationBRAF"
    ##                 Present Absent
    ## C1_mutationBRAF      34    345
    ## C2_mutationBRAF      13    246
    ## [1] "Chi-square test p-value for mutationBRAF : 0.0850587205980446"
    ## [1] "Contingency table for mutationHIST1H1E"
    ##                     Present Absent
    ## C1_mutationHIST1H1E      15    364
    ## C2_mutationHIST1H1E      12    247
    ## [1] "Chi-square test p-value for mutationHIST1H1E : 0.829048581603587"
    ## [1] "Contingency table for mutationLTB"
    ##                Present Absent
    ## C1_mutationLTB       9    370
    ## C2_mutationLTB      12    247
    ## [1] "Chi-square test p-value for mutationLTB : 0.178862769467262"

``` r
plots_grid <- grid.arrange(grobs = plot_list, ncol = 4)  # Adjust ncol for your desired layout
```

![](Suppl_5b_contingency_chisquare_genetics_files/figure-markdown_github/unnamed-chunk-4-1.png)
![Suppl_5b_contingency_chisquare_genetics](https://github.com/cleliacort/NRF1_paper/blob/main/Fig5/figures/Suppl/contingency_table_between_C1_and_C2_for_all_genetics__09APRILE25.png)
