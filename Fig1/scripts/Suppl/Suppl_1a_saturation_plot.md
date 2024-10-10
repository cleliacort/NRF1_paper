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
```

``` r
name_file_RG <- "data/Suppl/saturation_plot_all_TUMOUR_REGULATORY_REGIONS_0224.txt"
name_file_TSS <- "data/Suppl/saturation_plot_all_TUMOUR_TSS_0224.txt"
name_file_all <- "data/Suppl/saturation_plot_all_tumour_all_peaks_0224.txt"

df_RG <- read.delim(name_file_RG, sep = "\t", header = T) %>% mutate(group="RG")
df_TSS <- read.delim(name_file_TSS, sep = "\t", header = T) %>% mutate(group="TSS")
df_all <- read.delim(name_file_all, sep = "\t", header = T) %>% mutate(group="RG+TSS")

df_plot <- rbind(df_RG,df_TSS,df_all)
```

``` r
plot_line <- ggline(df_plot, "Number.of.features", "Median",
  linetype = "group", color = "group",
  palette = c("#696969","#00AFBB","#E7B800"),
  xlab = "SI",
  ylab = "Median number of covered nucleotides") +
  scale_x_continuous(breaks = seq(min(df_plot$Number.of.features), max(df_plot$Number.of.features),by = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + labs(linetype ="Region type",color = "Region type")+
  scale_y_continuous(labels = scales::comma) 

#show(plot_line)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/prova_saturation_plot_diveded_by_promoter_and_enhancer_plus_all_0724.png"
alt="Fig1_dataset_description" />
<figcaption aria-hidden="true">Fig1_dataset_description</figcaption>
</figure>
