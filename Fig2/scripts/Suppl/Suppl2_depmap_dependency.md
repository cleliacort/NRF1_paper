Load required libraries

``` r
library(ggpubr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(paletteer)
library(scales)
library(HGNChelper)
```

Load and prepare CRISPR DepMap data for NRF1 gene analysis.

``` r
crisp <- read.csv("data/Suppl/CRISPR_DepMap_21Q3_Public_Score_Chronos.csv", header = T)

crisp_nrf1 <- crisp %>% select(cell_line_display_name,lineage_1,lineage_2,NRF1) %>% arrange(NRF1)
crisp_nrf1 <- crisp_nrf1 %>% dplyr::mutate(rank=seq(1:nrow(crisp_nrf1)))

crisp_nrf1_mm <- crisp_nrf1 %>% dplyr::filter(lineage_2=='Multiple Myeloma') %>% arrange(NRF1) 
crisp_nrf1_mm <- crisp_nrf1_mm 
```

Load motif cluster data and transform gene names to HGNC-compliant
format.

Generate a scatter plot of TF fitness values for cluster 2 genes,
highlighting multiple myeloma cell lines in pink.

``` r
df_plot_long_cells <- merge(df_plot_long,crisp_nrf1,by.x = "cell_id",by.y = "cell_line_display_name")
labels_mm <- df_plot_long_cells %>% filter(lineage_2=='Multiple Myeloma')

ggplot_c2 <- ggplot(data=df_plot_long, aes(x= genes , y=Valore) )+
  geom_point(size = 3,show.legend = FALSE,colour="gray74")+
  geom_point(data = labels_mm,aes(x=genes, y=Valore),size=3,colour="#D57AA2")+
  #geom_hline(yintercept=0)+
  xlab("TF from cluster 2")+ylab("CRISP DepMap Fitness")+
  ylim(-2,1)+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

show(ggplot_c2)
```

![](Suppl2_depmap_dependency_files/figure-markdown_github/unnamed-chunk-7-1.png)

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/Suppl/depmap_dependency_cluster_motif_2_0724.png"
alt="Fig2_NRF1_median_expression_survival" />
<figcaption
aria-hidden="true">Fig2_NRF1_median_expression_survival</figcaption>
</figure>
