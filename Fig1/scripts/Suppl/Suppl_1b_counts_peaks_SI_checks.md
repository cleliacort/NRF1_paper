``` r
library(argparse)
library(ggpubr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
```

``` r
masterlist_SI <- read_delim("data/Suppl/master_list_atac_PEACKS_tumour_with_SI_0524.bed", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)
```

``` r
df_count <- masterlist_SI %>% group_by(X4) %>% count() 
df_count$five_labels <- c(0, rep(1:(nrow(df_count)-1)%/%5))

df_count_five <- df_count %>% group_by(five_labels) %>% summarise(sum_five = sum(n))
df_count_five <- df_count_five %>% mutate(type="five") %>% mutate(
    color_category = case_when(
      sum_five <= 3000 ~ "Low",
      sum_five > 3000 & sum_five <= 15000 ~ "Medium",
      sum_five > 15000 ~ "High"
    )
  )
```

``` r
custom_ticks <- c(seq(0,10,by=1))
custom_labels <- c("SI_1-SI_5", "SI_6-SI_10","SI_11-SI_15","SI_16-SI_20","SI_21-SI_25","SI_26-SI_30","SI_31-SI_35","SI_36-SI_40","SI_41-SI_45","SI_46-SI_50","SI_51-SI_55")

# Plot with customized x-axis labels
gg_bar <- ggplot(df_count_five, aes(x = five_labels, y = sum_five, fill = color_category)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_x_reverse(
    breaks = custom_ticks,  # Custom tick positions
    labels = custom_labels  # Custom tick labels
  ) +
  scale_fill_manual(values = c("Low" = "#FDD0A2", "Medium" = "#FDAE6B", "High" = "#AC532C")) +
  ylab("") +
  xlab("") +
  geom_text(aes(label = paste0("n=", sum_five), angle = 0),
            position = position_nudge(y = 0.5),
            color = "black", size = 4,
            hjust = -0.1) +
  scale_y_continuous(limits = c(0, 200000)) +
  coord_flip()
#show(gg_bar)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/barplot_number_per_group_5_peaks_SI_0724.png"
alt="Suppl_1b_counts_peaks_SI_checks" />
<figcaption
aria-hidden="true">Suppl_1b_counts_peaks_SI_checks</figcaption>
</figure>
