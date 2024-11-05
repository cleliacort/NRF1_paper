Upload the needed packages.

``` r
library(ggpubr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library("edgeR")
library(paletteer)
library(scales)
```

Read matrix input file containing enrichment information derived from
ATAC-seq experiments conducted on the in-house MM cohort (n=55) and MGUS
cohort (n=11), specifically evaluating 699 NRF1-tumor-specific sites and
300 random sites shared with MGUS.

``` r
file_chip_consensus="data/matrix_on_11_master_list_consensus_our_chip_1123_master_list_mgus_ANNOTATED_TSS_300random_plus_10_master_list_consensus_our_chip_1123.txt"

reads_matrix <- read.delim(file_chip_consensus) %>% 
  mutate(coordinates=paste(chromosome,start, end, sep="_")) %>% select(-chromosome,-start, -end) %>% 
  column_to_rownames(var="coordinates")
```

Customize the input matrix.

``` r
reads_matrix <- as.data.frame(subset(reads_matrix, rowMeans(reads_matrix) > 1))
sorted_columns <- sort(names(reads_matrix))
reads_matrix_sorted <- reads_matrix[, sorted_columns]
```

Retrieve the clinical and genetic information of the in-house MM cohort,
and incorporate this data into the design dataframe. Customize the
design dataframe by also including the absence of genetic and clinical
information from the MGUS.

``` r
# Retrive the ID and the groups to which the samples belong 
pheno <- data_frame(sample_id=colnames(reads_matrix_sorted)) %>% arrange(sample_id) %>% 
  mutate(type=case_when(str_detect( sample_id, "tumour")~ "tumour", TRUE ~ "mgus"))

samplesheet_orig <- read.delim("data/samplesheet_with_clinical_data_customized_cytogenetic_pc_info_mm59bis_4kpeaks_necessary_info_03september2024.tsv",header = T) %>% distinct() 

df_sub1 <- samplesheet_orig %>% select(patient_id,state,perc_PC_pat)
df_sub2 <-   samplesheet_orig %>% select(official_labelling,state,perc_PC_pat) %>% rename(patient_id=official_labelling)
df_samplesheet <- rbind(df_sub1,df_sub2) %>% filter(!str_detect(patient_id, "MGUS")) %>% 
  filter(!str_detect(patient_id, "MM")) %>% rename(sample_id=patient_id)

pheno_status <- merge(pheno, df_samplesheet,by="sample_id",all.x = TRUE) %>% mutate(state = if_else(is.na(state), "MGUS", state)) %>%
  arrange(sample_id)
```

Normalize the matrix.

``` r
groups_sample <- pheno$type
e <- as.matrix(reads_matrix_sorted)
regions <- DGEList(counts=e, group=groups_sample)
regions <- calcNormFactors(regions, method="TMM") 
normalized_counts <- cpm(regions, normalized.lib.sizes = TRUE)
```

Scale the matrix.

``` r
t_normalized_count <-  as.data.frame(t(normalized_counts)) 
scaled_matrix <- t(scale(t_normalized_count)) 
```

Define the columns annotation for the heatmap.

``` r
color_gradient <- colorRamp2(c(0,20,40,60,80,100), c("white","#EAC0BDFF", "#F57667FF","#D8392CFF","#C21417FF","#9C0824FF"))

colAnn <- HeatmapAnnotation(na_col = "#434343FF",
  group = pheno_status$type,
  status = pheno_status$state,
  perc_PC = pheno_status$perc_PC_pat,
  
  col = list(group = c("tumour"="#9A68A4FF","mgus" = "#7DBE60FF"),
             status = c('ONSET'='#FED9A6', 'TREATED'='#DECBE4',"MGUS"= "white"),
             perc_PC = color_gradient),
  simple_anno_size = unit(0.3, "cm"),
  annotation_label = c("Disease Type","Disease status","Perc. of malignant plasmacell"),
  show_legend = c(TRUE),
  annotation_name_gp= gpar(fontsize = 7)
)
```

Define the rows annotation for the heatmap.

``` r
annot_cluster_sites="data/11_master_list_consensus_our_chip_1123_master_list_mgus_ANNOTATED_TSS_300random_plus_10_master_list_consensus_our_chip_TYPE_1123.bed"
annot_cluster_sites_df <- read.delim(annot_cluster_sites,header = F) %>%
  mutate(coordinates=paste(V1,V2, V3, sep="_")) %>% select(-V1,-V2, -V3) %>%
  rename(type=V4) %>% arrange(coordinates)
sub_sites <- data.frame(coordinates=rownames(scaled_matrix))
annot_site_selected <- merge(annot_cluster_sites_df,sub_sites,by="coordinates")%>% arrange(type)

# Order the matrix based on the type of sites (only or common)
scaled_matrix <- scaled_matrix[annot_site_selected$coordinates, ]

# Create the rowannotation
rowAnn <- rowAnnotation(
  type_of_sites = annot_site_selected$type,
  col = list(type_of_sites=c("only"="#225188FF","common"="#ADD8E6")),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = c(TRUE),
  show_annotation_name = FALSE
)
```

Create the heatmap

``` r
set.seed(123)
num_col_split <- 2

distance_col <- "euclidean"
method_col <- "ward.D"

mycols <- colorRamp2(breaks = c(-4,-2,0,2,4), colors = c("#08519C","#4292C6","#FFFFFF","#D7301F","#B30000"))#roma


h <- ComplexHeatmap::Heatmap(scaled_matrix,

             name = "Z-score Cpm-TMM expression", 
             column_names_gp = gpar(fontsize = 3),

             col = mycols,
             show_row_names=F,
             show_column_names=T,
             
             top_annotation=colAnn,
             right_annotation=rowAnn,

             clustering_distance_columns = distance_col,
             clustering_method_columns = method_col,
             
             na_col = "white",

             column_split = num_col_split,
             
             cluster_columns = T,
             cluster_rows=F,

             border = TRUE
)
#show(h)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/heatmap_tot_cluster_2_normalized_scaled_1111_chip_nrf1_kms27_chip_nrf1_mm196_chip_nrf1_mm217_chip_nrf1_kms18_1223_cluster_col_euclidean_ward.D_rows_no_clustering_1223.png"
alt="Fig3_heatmap_tumour_mgus" />
<figcaption aria-hidden="true">Fig3_heatmap_tumour_mgus</figcaption>
</figure>
