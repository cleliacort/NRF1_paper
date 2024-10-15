Reading the necessary files data.

``` r
patients_class <- "data/patient_classification_based_on_median_expression_of_NRF1.txt"
patients_class_df <- read.delim(patients_class, sep = "\t", header = T) #%>% distinct() 

iss <- read.delim("../COMMPASS_IA17/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% dplyr::select(PUBLIC_ID,D_PT_iss)

patient_class_iss_df <- merge(patients_class_df,iss,by= "PUBLIC_ID")
```

Making the plot.

``` r
patient_class_iss_df$status <- factor(patient_class_iss_df$status, levels = c("low", "high"))
my_comparisons <- list(c("low", "high") )

p <- ggboxplot(patient_class_iss_df,x = "D_PT_iss", y = "NRF1", fill = "status", color = "#696969",
               palette = c("#8B9DAC","#efbf9a"),width = 0.6,outlier.shape = NA,
               lwd=0.5,xlab="ISS",ylab="Normalized reads counts(tpm)" ,
               bxp.errorbar =TRUE)+
  stat_compare_means(aes(group = status), label = "p.signif") #wilcox
```

    ## Warning in (function (mapping = NULL, data = NULL, geom = "boxplot", position =
    ## "dodge2", : Ignoring unknown aesthetics: fill

``` r
p <- ggpar(p,ylim=c(0,15))+border()
#show(p)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/TPM_per_ISS_in_population_divided_per_NRF1_expression_low_high_201123.png"
alt="Fig2_NRF1_median_expression_survival" />
<figcaption
aria-hidden="true">Fig2_NRF1_median_expression_survival</figcaption>
</figure>
