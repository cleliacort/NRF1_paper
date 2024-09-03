# Fig1_dataset_description

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
```

Reading the clinical information of our tumour cohort and ATAC-seq peaks
related to each patient.

``` r
samplesheet_orig <- read_delim("data/samplesheet_with_clinical_data_customized_cytogenetic_pc_info_mm59bis_4kpeaks_03september2024.tsv",delim = "\t", col_names = T) 

peaks <- read.delim("data/number_of_peaks_per_sample_atac_tumour_mgus_2023.txt",sep = "\t", header = T) %>% arrange(NumLine)

samplesheet <- merge(samplesheet_orig,peaks,by.x="official_labelling",by.y="Filename")
```

Make the sub-heatmap with patient status information.

``` r
sub_state <- samplesheet %>% select("official_labelling","state","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine)
df_state <- as.data.frame(t(sub_state))
df_matrix_state <- as.matrix(df_state)

h_state <- Heatmap(df_matrix_state,
                 name = "Disease Status", 
                 row_names_gp = gpar(fontsize = 6),
                 col =  c('ONSET'='#FED9A6', 'TREATED'='#DECBE4'),
                 na_col = "white",
                 show_row_names=T,
                 show_column_names=,
                 rect_gp = gpar(col = "black", lwd = 1),
                 border = TRUE,
                 cluster_columns = F,
                 cluster_rows=F
                 )
#print(h_state)
```

Make the sub-heatmap with cytogenic information.

``` r
sub_cit <- samplesheet %>% select("official_labelling","t(4;14)","t(11;14)","1qgain","del13q","IGH_rearr","delp53","trisomia11","trisomia4",
                                  "iperdiploidia","negativa","unknown","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine)
df_cit <- as.data.frame(t(sub_cit))
df_matrix_cit <- as.matrix(df_cit)


h_cit <- Heatmap(df_matrix_cit,
                 name = "Presence/Absence", 
                 row_names_gp = gpar(fontsize = 6),
                 col =  c('1'='#FDAE6B', '0'='white'),
                 na_col = "white",
                 show_row_names=T,
                 show_column_names=,
                 rect_gp = gpar(col = "black", lwd = 1),
                 border = TRUE,
                 cluster_columns = F,
                 cluster_rows=F
                 )
# print(h_cit)
```

Make the sub-heatmap with sex information.

``` r
sub_sex <- samplesheet %>% select("official_labelling","sex","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine)
df_sex <- as.data.frame(t(sub_sex))
df_matrix_sex <- as.matrix(df_sex)

h_sex <- Heatmap(df_matrix_sex,
                 name = "sex", 
                 column_names_gp = gpar(fontsize = 0.1),
                 row_names_gp = gpar(fontsize =6),
                 col =  c('M'='#80B1D3', 'F'='#FCCDE5'),
                 na_col = "white",
                 show_row_names=T,
                 show_column_names=F,
                 rect_gp = gpar(col = "black", lwd = 1),
                 border = TRUE,
                 cluster_columns = F,
                 cluster_rows=F
)

# show(h_sex)
```

Make the barlot with number of peaks per patient information.

``` r
number_peaks_barplot <- samplesheet %>% select("official_labelling","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) 
```

Make the sub-heatmap with the number of tumour plasma cell information.

``` r
pc_perc <- samplesheet %>% select("official_labelling","PC% MIDOLLARI (PM)  esordio","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine) %>% 
  separate("PC% MIDOLLARI (PM)  esordio",c("perc",NA),sep="%") %>% 
  mutate(perc=as.numeric(perc))
df_pc <- as.data.frame(t(pc_perc))
df_matrix_pc <- as.matrix(df_pc)

sub_pc <- samplesheet %>% select("official_labelling","perc_PC_pat","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine)
df_pc <- as.data.frame(t(sub_pc))
df_matrix_pc <- as.matrix(df_pc)


purples_palette <- brewer.pal(n = 8, name = "Purples")
selected_colors_p <- purples_palette[3:8]

h_pc <- Heatmap(df_matrix_pc,
                 name = "% malignant PC", 
                 column_names_gp = gpar(fontsize = 0.1),
                row_names_gp = gpar(fontsize = 6),
                 col =  selected_colors_p,
                 na_col = "white",
                 show_row_names=T,
                 show_column_names=F,
                 rect_gp = gpar(col = "black", lwd = 1),
                 border = TRUE,
                 cluster_columns = F,
                 cluster_rows=F
)
#show(h_pc)
```

Make the sub-heatmap with age information.

``` r
sub_age <- samplesheet %>% select("official_labelling","age","NumLine") %>% 
  column_to_rownames(var="official_labelling") %>% arrange(NumLine) %>% select(-NumLine) 
df_age <- as.data.frame(t(sub_age))
df_matrix_age <- as.matrix(df_age)

greens_palette <- brewer.pal(n = 8, name = "Greens")
selected_colors_g <- greens_palette[4:8]

h_age <- Heatmap(df_matrix_age,
                 name = "age", 
                 column_names_gp = gpar(fontsize = 6),
                 row_names_gp = gpar(fontsize = 6),
                 col = selected_colors_g,
                 na_col = "white",
                 show_row_names=T,
                 show_column_names=T,
                 rect_gp = gpar(col = "black", lwd = 1),
                 border = TRUE,
                 cluster_columns = F,
                 cluster_rows=F
)
# show(h_age)
```

Create the multi layers heatmap.

``` r
ha = HeatmapAnnotation(Num_of_peaks = anno_barplot(number_peaks_barplot$NumLine, height = unit(2, "cm")))
ht_list <- ha %v% h_state %v% h_cit %v% h_pc %v% h_sex %v% h_age
show(ht_list)
```

![](Figure1_dataset_pheno_description_files/figure-markdown_github/unnamed-chunk-8-1.png)

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/heatmap_clinical_and_peaks_information_with_disease_status_030924.png"
alt="Fig1_dataset_description" />
<figcaption aria-hidden="true">Fig1_dataset_description</figcaption>
</figure>
