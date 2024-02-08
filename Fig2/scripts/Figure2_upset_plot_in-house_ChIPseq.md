# Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq

Upload the needed packages.

``` r
library("UpSetR")
library(ggpubr)
```

Retrieve the count of peak intersections among NRF1 ChIP-seq experiments
conducted in-house. This information is sourced from the ‘upset’ module
of the previously used Intervene tools.

``` r
expressionInput <- c('chip_nrf1_kms18'=45,'chip_nrf1_mm217'=4224,'chip_nrf1_mm217&chip_nrf1_kms18'=3,'chip_nrf1_mm196'=5854,'chip_nrf1_mm196&chip_nrf1_kms18'=14,'chip_nrf1_mm196&chip_nrf1_mm217'=1499,'chip_nrf1_mm196&chip_nrf1_mm217&chip_nrf1_kms18'=48,'chip_nrf1_kms27'=13651,'chip_nrf1_kms27&chip_nrf1_kms18'=493,'chip_nrf1_kms27&chip_nrf1_mm217'=1325,'chip_nrf1_kms27&chip_nrf1_mm217&chip_nrf1_kms18'=356,'chip_nrf1_kms27&chip_nrf1_mm196'=513,'chip_nrf1_kms27&chip_nrf1_mm196&chip_nrf1_kms18'=125,'chip_nrf1_kms27&chip_nrf1_mm196&chip_nrf1_mm217'=2270,'chip_nrf1_kms27&chip_nrf1_mm196&chip_nrf1_mm217&chip_nrf1_kms18'=5749)
```

Create the upset plot.

``` r
plot <- upset(fromExpression(expressionInput), nsets=6, nintersects=30, show.numbers="yes", main.bar.color="#969696", sets.bar.color="#FDAE6B", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size",text.scale = 1.5)
#show(plot)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/Intervene_upset_KMS27_MM196_MM217_KMS18_1123.png"
alt="Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq" />
<figcaption
aria-hidden="true">Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq</figcaption>
</figure>
