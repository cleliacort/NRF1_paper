Upload the needed packages.

``` r
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(kableExtra)
```

Read input file containing the regulatory regions and the related ABC
score.

``` r
CRE_intersect_ABC <- read_delim("data/non_coding_regions_intersected_EnhancerPredictions_greater2kb_0524.bed", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)

CRE_intersect_ABC_mod <- CRE_intersect_ABC %>% select(X7,X8) %>% separate(X7,c("genes"),sep="_")
```

Create the plot.

``` r
CRE_intersect_ABC_01 <- CRE_intersect_ABC_mod %>% filter(X8>0.1)

CRE_intersect_ABC_mod1 <- CRE_intersect_ABC_mod %>% mutate(genes = if_else(near(X8, 0.03002231), "C8orf76_bis", genes))


plot <- ggplot(CRE_intersect_ABC_mod1,aes(x=reorder(genes,X8),y=X8)) +
  geom_point(aes(size = 2),colour="gray74")+ #set features of dots
  geom_point(data = CRE_intersect_ABC_01,aes(x= reorder(genes,X8),y=X8,size=2),colour="indianred3")+
  #geom_text_repel(data = comm, aes(label = genes, x=reorder(genes,X8),y=X8),min.segment.length = 0.7,box.padding = 0.7) +
  scale_fill_identity() +
  theme_bw()+ 
  theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())+ 
  ylab("ABC-score")+xlab("Gene target of selected regulatory regions")+
  guides(size = FALSE)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey30") + # add horizontal line
  scale_x_discrete(expand = expansion(mult = c(0.04,0.04))) #+
  #geom_point(data = CRE_intersect_ABC_01,shape = 1,size = 5,stroke = 0.6,colour="grey30")


# show(plot)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/dotplot_our_regulatory_regions_intersected_with_abc_score_NO_LABELS_0624.png"
alt="Fig3_ABC_score_CRE_non_coding" />
<figcaption
aria-hidden="true">Fig3_ABC_score_CRE_non_coding</figcaption>
</figure>
