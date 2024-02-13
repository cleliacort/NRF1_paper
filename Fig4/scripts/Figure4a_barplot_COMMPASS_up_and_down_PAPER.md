# Figure4a_barplot_COMMPASS_up_and_down

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

Reading and customizing information about which genes significantly
increase(UP) or decrease(DOWN) their expression during disease
progression.

``` r
file_name="data/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_05_wilcox_OR_ISS_pv1_3.txt"
gene_list <- read.delim(file_name,sep = "\t", header  = T) 
counts_genes_diff <- gene_list %>% group_by(EXPRESSION) %>% count() %>% mutate(tot_tss=length(gene_list$GENE_NAME)) %>% 
  mutate(perc=round((n*100)/tot_tss,2))
```

Create a bar plot showing the number of genes classified as UP, DOWN,
and NA.

``` r
cust_palette <- c("#E79E85","#A0A0A0","#b3bdcc")

counts_genes_diff$EXPRESSION <- factor(counts_genes_diff$EXPRESSION, levels = c("UP","NAN","DOWN"))

gg_perc <- ggplot(counts_genes_diff,aes(fill= EXPRESSION, y= perc, x= reorder(EXPRESSION,perc))) + 
  geom_bar(stat = "identity",aes(fill=EXPRESSION))+
  coord_flip()+
  scale_fill_manual(values = cust_palette)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Percentage (%)")+
  xlab("")+
  labs(title = "Percentage of site annotation types") +
  guides(fill = "none")+
  geom_text(aes(label = paste(perc, "%", sep = ""), angle = 0),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4,
            hjust = -1)+
      labs(title = "Num. of genes UP/DOWN in COMMPASS coming from TSS +-2kb (1123)",size = "small")#+
#show(gg_perc)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/barplot_COMMPASS_number_of_genes_up_down_na_LABELS_0124.png"
alt="Figure4a_barplot_COMMPASS_up_and_down" />
<figcaption
aria-hidden="true">Figure4a_barplot_COMMPASS_up_and_down</figcaption>
</figure>
