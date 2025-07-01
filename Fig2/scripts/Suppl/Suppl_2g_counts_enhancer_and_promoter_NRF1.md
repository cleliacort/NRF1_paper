# Suppl_2g_counts_enhancer_and_promoter_NRF1
Upload the needed packages.

``` r
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(paletteer)
library(edgeR)
```

Read the raw read counts

``` r
count_tumour_mgus="data/Suppl//matrix_multicov_tumout_and_MGUS_on_tumour_master_list_0124.txt"

count_matrix_file <- read.delim(count_tumour_mgus, header = TRUE, sep = "\t") 
count_matrix_cust <- count_matrix_file %>%   
  rename(chr=colnames(count_matrix_file)[1],start=colnames(count_matrix_file)[2],end=colnames(count_matrix_file)[3]) %>% 
  mutate(coordinates = paste0(chr,"_",start,"_",end)) %>% 
  select(-chr,-start,-end) %>% 
  column_to_rownames(var="coordinates")
```

Filter the raw count matrix

``` r
e=data.matrix(count_matrix_cust)
e=subset(e, rowMeans(e) > 5)
e=subset(e, rowMeans(e) < 5000)

groups_pheno <- data.frame(ID=colnames(count_matrix_cust)) %>% 
  mutate(phenotype=case_when(str_detect( ID, "tumour")~ "tumour", TRUE ~ "MGUS")) 
```

Normalize the raw counts

``` r
regions <- DGEList(counts=e, group=groups_pheno$phenotype)
regions <- estimateCommonDisp(regions)
regions <- estimateTagwiseDisp(regions)
regions <- calcNormFactors(regions, method="TMM")
cpm_tmm <- cpm(regions, normalized.lib.size=T)
```

Customize the normalized reads count matrix.

``` r
cpm_tmm_cust <- as.data.frame(cpm_tmm) 
t_cpm_tmm_cust <- as.data.frame(t(cpm_tmm_cust)) %>% rownames_to_column(var="ID")
```

Add the state to the normalized reads count matrix.

``` r
samplesheet_orig <- read_delim("../Fig1/data/samplesheet_with_clinical_data_customized_cytogenetic_pc_info_mm59bis_4kpeaks_necessary_info_03september2024.tsv",delim = "\t", col_names = T) %>% distinct() %>% 
  filter(state!="RELAPSE") %>% filter(state!="REMISSION") %>% filter(state!="MGUS")

peaks <- read.delim("../Fig1/data/Suppl/number_of_peaks_per_sample_atac_tumour_mgus_2023.txt",sep = "\t", header = T) %>% arrange(NumLine)

samplesheet <- merge(samplesheet_orig,peaks,by.x="official_labelling",by.y="Filename")

df_plot_with_stage <- merge(t_cpm_tmm_cust,samplesheet,by.x = "ID", by.y = "official_labelling", all.x = TRUE)
```

Make the plot related to the NRF1 enhancer.

``` r
region_corodinates <- "7_129419085_129420959" #ENHANCER NRF1
df_plot_enhancer <-  df_plot_with_stage %>% dplyr::select(ID,!!region_corodinates,state) %>% 
  mutate(state= case_when(is.na(state) ~ "MGUS", state=="ONSET" ~ "ONSET", TRUE ~ "TREATED")) %>% 
  rename(cpm_tmm=region_corodinates)

comparisons <- list(c("MGUS", "ONSET"), c("MGUS", "TREATED"), c("TREATED", "ONSET"))

p <- ggboxplot(df_plot_enhancer,x = "state", y = "cpm_tmm", fill = "state", color = "#696969",
             width = 0.6,outlier.shape = NA,palette = c("#8B9DAC","#FED9A6","#DECBE4"),
             lwd=0.5,xlab="state",ylab="cpm_tmm" ,
             bxp.errorbar =TRUE,
        add = "jitter")+
stat_compare_means(aes(group = state)) #+ 
#stat_compare_means(aes(group = state), comparisons = comparisons, method = "wilcox.test")
    
p <- ggpar(p,ylim=c(0,200))+border()
# show(p) 
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/Suppl/NRF1_enrichment_ONSET_TREATED_0724_cpm_tmm_correct_ENHANCER.png"
alt="Suppl_2g_NRF1_enrichment_ONSET_TREATED_enhancer" />
<figcaption
aria-hidden="true">Suppl_2g_NRF1_enrichment_ONSET_TREATED_enhancer</figcaption>
</figure>

Make the plot related to the NRF1 promoter.

``` r
region_corodinates <- "7_129251005_129252566" #PROMOTER NRF1

df_plot_promoter <-  df_plot_with_stage %>% dplyr::select(ID,!!region_corodinates,state) %>% 
  mutate(state= case_when(is.na(state) ~ "MGUS", state=="ONSET" ~ "ONSET", TRUE ~ "TREATED")) %>% 
  rename(cpm_tmm=region_corodinates)

comparisons <- list(c("MGUS", "ONSET"), c("MGUS", "TREATED"), c("TREATED", "ONSET"))

p <- ggboxplot(df_plot_promoter,x = "state", y = "cpm_tmm", fill = "state", color = "#696969",
             width = 0.6,outlier.shape = NA,palette = c("#8B9DAC","#FED9A6","#DECBE4"),
             lwd=0.5,xlab="state",ylab="cpm_tmm" ,
             bxp.errorbar =TRUE,
        add = "jitter")+
stat_compare_means(aes(group = state)) + 
stat_compare_means(aes(group = state), comparisons = comparisons, method = "wilcox.test")
    
p <- ggpar(p,ylim=c(0,200))+border()

# show(p)  
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/Suppl/NRF1_enrichment_ONSET_TREATED_0724_cpm_tmm_correct_PROMOTER.png"
alt="Suppl_2g_NRF1_enrichment_ONSET_TREATED_promoter" />
<figcaption
aria-hidden="true">Suppl_2g_NRF1_enrichment_ONSET_TREATED_promoter</figcaption>
</figure>
