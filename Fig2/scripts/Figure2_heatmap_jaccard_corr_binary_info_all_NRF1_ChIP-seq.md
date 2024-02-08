# Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP

Upload the needed packages.

``` r
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(paletteer)
library(ComplexHeatmap)
library(circlize)
```

Define the function to compute the jaccard coefficient (JC).

``` r
compute_jaccard_similarity_df <- function(df) {
  num_rows <- nrow(df)
  similarity_matrix <- matrix(0, nrow = num_rows, ncol = num_rows,
                              dimnames = list(rownames(df), rownames(df)))
  for (i in 1:num_rows) {
    for (j in 1:num_rows) {
      set1 <- df[i, ]
      set2 <- df[j, ]
      intersection <- sum(set1 & set2)
      union_val <- sum(set1 | set2)
      similarity_matrix[i, j] <- intersection / union_val
    }
  }
  return(similarity_matrix)
}
```

Read a binary matrix that is populated based on the presence/absence of
the CISTROME ChIP-seq signal and the ENCODE ChIP-seq signal at NRF1
binding sites listed in the consensus master list.

``` r
file_chip_consensus="data/matrix_binary_filled_chip_KMS27_MM196_MM217_KMS18_our_and_cistrome_and_encode_chip_1123.txt"

binary_matrix <- read_delim(file_chip_consensus,delim = "\t", escape_double = FALSE,trim_ws = TRUE) %>% distinct() %>% 
  mutate(coordinates=paste0(chromosome,"-",start,"-",end)) %>% select(-chromosome, -start,-end) %>% column_to_rownames(var="coordinates") %>% select(-"GSM1891655_NRF1_CHIP_HMEC_C-DB",-"GSM1891656_NRF1_CHIP_HMEC_C-DB")

# Transpone the matrix
t_matrix <- t(binary_matrix)
```

Compute the JC and customize the results.

``` r
jaccard_result_t <- as.data.frame(compute_jaccard_similarity_df(t_matrix))
sub_jaccard <- jaccard_result_t %>% dplyr::select(contains("IRE") | contains("kms"))
t_sub_jaccard <- t(sub_jaccard)
```

Customize the cell type labels.

``` r
metadata <- data.frame(sample_id=colnames(t_sub_jaccard)) %>% 
    mutate(type_disease=case_when(grepl("K562",sample_id) ~ "leukemia",
                                grepl("HMEC",sample_id) ~ "normal_breast",
                                grepl("MCF-7",sample_id) ~ "breast_cancer",
                                grepl("T47D",sample_id) ~ "breast_cancer",
                                grepl("HCC1954",sample_id) ~ "breast_cancer",
                                grepl("GM12878",sample_id) ~ "lymphoblastoid",
                                grepl("SK-N-SH",sample_id) ~ "neuroblastoma",
                                grepl("HepG2",sample_id) ~ "hepatocellular_carcinoma",
                                TRUE ~ "MM")) 
```

Define the column annotation for the heatmap.

``` r
colAnn <- HeatmapAnnotation(
  type_disease = metadata$type_disease,
  
  col = list(type_disease=c("leukemia"= "#FF69B4",
                            "MM"="#BA55D3" ,
                            "breast_cancer"="#FFD700",
                            "neuroblastoma"="#98FB98",
                            "hepatocellular_carcinoma"="#D3D3D3",
                            "normal_breast"="#ADD8E6",
                            "lymphoblastoid"="#87CEFA")),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = c(TRUE),
  show_annotation_name = FALSE
)
```

Create the heatmap with JC information.

``` r
set.seed(123)

# Set color palette for the heatmap
mycols <- colorRamp2(breaks = c(0,0.6,0.8,1), colors = c("#33608CFF","#BBD6EAFF","#FDECED","#9A68A4FF"))

#Set the parameters for the heatmap 
distance_rows <- "manhattan"
method_rows <- "ward.D2"
distance_col <- "manhattan"
method_col <- "ward.D2"

h_corr <- Heatmap(t_sub_jaccard,
  column_title = paste0("master_list_CHIP_NRF1_CONSENSUS_KMS27_MM196_MM217_KMS18 (n=",length(rownames(binary_matrix)),")"),
  name = "Jaccard Coefficient", 
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  
  bottom_annotation = colAnn,

  clustering_distance_rows = distance_rows,
  clustering_method_rows = method_rows,
  
  clustering_distance_columns = distance_col,
  clustering_method_columns = method_col,
  
  col =  mycols,
  na_col = "white",
  show_row_names=T,
  show_column_names=T,
  rect_gp = gpar(col = "black", lwd = 1),
  border = TRUE,
  cluster_columns = T,
  cluster_rows=T,
  
  show_row_dend = F,
  show_column_dend = F,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", t_sub_jaccard[i, j]), x, y, gp = gpar(fontsize = 12))}

)
#show(h_corr)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/heatmap_jaccardCoeff_of_OUR_AND_ENCODE_AND_CISTROME_CHIP_NRF1_on_CONSENSUS_masterlist_KMS27_MM196_MM217_KMS18_1123.png"
alt="Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq" />
<figcaption
aria-hidden="true">Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq</figcaption>
</figure>
