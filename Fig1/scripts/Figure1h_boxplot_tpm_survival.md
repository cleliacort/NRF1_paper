# Figure1h_boxplot_tpm_survival

Upload the needed packages.

Read patient classification based on NRF1 expression and ISS
information, and merge them.

``` r
patients_class <- "data/patient_classification_based_on_median_expression_of_NRF1.txt"
patients_class_df <- read.delim(patients_class, sep = "\t", header = T) 
iss <- read.delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% dplyr::select(PUBLIC_ID,D_PT_iss)
patient_class_iss_df <- merge(patients_class_df,iss,by= "PUBLIC_ID")
```

Create boxplot showing the NRF1 enrichment per ISS in high and low NRF1
expression level populations.

``` r
patient_class_iss_df$status <- factor(patient_class_iss_df$status, levels = c("low", "high"))
my_comparisons <- list(c("low", "high") )

p <- ggboxplot(patient_class_iss_df,x = "D_PT_iss", y = "NRF1", fill = "status", color = "#696969",
               palette = c("#8B9DAC","#efbf9a"),width = 0.6,outlier.shape = NA,
               lwd=0.5,xlab="ISS",ylab="Normalized reads counts(tpm)" ,
               bxp.errorbar =TRUE)+
  stat_compare_means(aes(group = status), label = "p.signif") #wilcox
    
p <- ggpar(p,ylim=c(0,15))+border()
#show(p)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/TPM_per_ISS_in_population_divided_per_NRF1_expression_low_high_201123.png"
alt="Figure1h_boxplot_tpm_survival" />
<figcaption
aria-hidden="true">Figure1h_boxplot_tpm_survival</figcaption>
</figure>
