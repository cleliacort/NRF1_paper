# Fig2g_and_Suppl2_emission_and_transition_heatmap_chromHMM

``` r
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(paletteer)
```

Uploading emission file.

``` r
emissions_6 <- read_delim("data/emissions_6.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE) %>% 
  gather(key = "modifications", value = "Value", c(-`State (Emission order)`)) %>% 
  mutate(modifications=as.factor(modifications))

emissions_6$modifications <- factor(emissions_6$modifications, levels = c("H3k4me3", "H3k27Ac", "H3k4mono","H3k27me3","H3k9me3"))
```

Making the heatmap with emission information.

``` r
heatmap_emi <- ggplot(emissions_6, aes(x = modifications, y = as.factor(`State (Emission order)`), fill = Value)) +
  geom_tile(color = "#595959", size = 0.5) +
  scale_fill_gradient(low = "white",high = "#08519C") + 
  geom_text(aes(label = round(Value, 3)), color = "#292929", size = 4) +  # Add text labels
  theme_minimal()
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
#show(heatmap_emi)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/heatmap_emission_made_0924.png"
alt="Fig2_NRF1_median_expression_survival" />
<figcaption
aria-hidden="true">Fig2_NRF1_median_expression_survival</figcaption>
</figure>

Uploading transition file.

``` r
transitions_6 <- read_delim("data/transitions_6.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE) %>% 
  gather(key = "transitions", value = "Value", c(-`State (from\\to) (Emission order)`)) %>% 
  mutate(transitions=as.factor(transitions))

#show(transitions_6)
```

Making the heatmap with transition information.

``` r
heatmap_transitions <- ggplot(transitions_6, aes(x = transitions, y = as.factor(`State (from\\to) (Emission order)`), fill = Value)) +
  geom_tile(color = "#595959", size = 0.5) +
  scale_fill_gradient(low = "white",high = "darkred") + 
  geom_text(aes(label = round(Value, 3)), color = "#292929", size = 4) + 
  xlab("")+ylab("")+labs(title="Transition")+# Add text labels
  theme_minimal()
#show(heatmap_transitions)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/Suppl/heatmap_transitions_made_0924.png"
alt="Fig2g_emission_and_transition_heatmap_chromHMM" />
<figcaption
aria-hidden="true">Fig2g_emission_and_transition_heatmap_chromHMM</figcaption>
</figure>
