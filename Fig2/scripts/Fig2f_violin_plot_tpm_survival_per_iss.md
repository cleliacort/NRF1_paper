# Fig2f_violin_plot_survival_per_ISS

First merge: keep only patients analyzed in the paper, i.e., those with
full clinical info. Second merge: add ISS info.

``` r
patients_class <- "/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig2/data/patient_classification_based_on_median_expression_of_NRF1.txt"
patients_class_df <- read.delim(patients_class, sep = "\t", header = T) 
iss <- read.delim("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig4/data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% dplyr::select(PUBLIC_ID,D_PT_iss)
nrf1_patient_sel_iss <- merge(patients_class_df,iss,by="PUBLIC_ID")
```

For each ISS group (1, 2, 3), calculate the median NRF1 expression,then
classify patients as “HIGH” or “LOW” within each ISS group

``` r
iss1 <- nrf1_patient_sel_iss %>% filter(D_PT_iss==1)
median1 <- median(iss1$NRF1)
iss1 <- iss1 %>% mutate(expr_level=case_when(NRF1>=median1 ~ "HIGH", TRUE ~ "LOW"))

iss2 <- nrf1_patient_sel_iss %>% filter(D_PT_iss==2)
median2 <- median(iss2$NRF1)
iss2 <- iss2 %>% mutate(expr_level=case_when(NRF1>=median2 ~ "HIGH", TRUE ~ "LOW"))

iss3 <- nrf1_patient_sel_iss %>% filter(D_PT_iss==3)
median3 <- median(iss3$NRF1)
iss3 <- iss3 %>% mutate(expr_level=case_when(NRF1>=median3 ~ "HIGH", TRUE ~ "LOW"))
```

Combine all ISS groups back together

``` r
iss_df <- rbind(iss1,iss2,iss3) 
all_patient_df <- iss_df %>% mutate(x_name=paste0(expr_level,"_",D_PT_iss))
```

LOW expression group analysis

``` r
df_low <- all_patient_df %>% filter(expr_level=="LOW")
df_low$D_PT_iss <- ordered(df_low$D_PT_iss, levels = c("1", "2", "3"))
j_test_low <- jonckheereTest(NRF1 ~ D_PT_iss, data = df_low, alternative = "greater")
df_low <- all_patient_df %>% filter(expr_level=="LOW")

df_low$x_name <- factor(df_low$x_name, levels = c("LOW_1","LOW_2","LOW_3"),
                                 labels = c("Low_ISS1","Low_ISS2","Low_ISS3"))

pv <- paste0("Jonckheere-test, p=",round(j_test_low$p.value,4))

p <- ggviolin(df_low,x = "x_name", y = "NRF1", fill = "expr_level", color = "#696969",
               palette = c("#8B9DAC"),width = 0.6,outlier.shape = NA,
               lwd=0.5,xlab="ISS",ylab="Normalized reads counts(tpm)" ,
               alpha = 0.5,
         bxp.errorbar =TRUE,
         add = c("boxplot"),add.params = list(size = 0.5, width = 0.2, alpha = 0.5))+
  scale_y_continuous(limits = c(0, 10),breaks = seq(0, 10, by = 2))+
      theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) + border()+
  xlab("")+  labs(fill = "Expression Groups")  +
    annotate("text", x = "Low_ISS2", y = 8, label = pv, size = 4, color = "black")
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/TPM_per_ISS_in_population_divided_per_NRF1_expression_vertical_boxplot_per_iss_with_statistics_ONLY_LOW_WITH_jonckheereTest_INCRESING_WITH_ISS_25feb25.png"
alt="Fig2f_violin_plot_survival_per_ISS_group_LOW" />
<figcaption
aria-hidden="true">Fig2f_violin_plot_survival_per_ISS_group_LOW</figcaption>
</figure>

HIGH expression group analysis

``` r
df_high <- all_patient_df %>% filter(expr_level=="HIGH")
df_high$D_PT_iss <- ordered(df_high$D_PT_iss, levels = c("1", "2", "3"))
j_test_high <- jonckheereTest(NRF1 ~ D_PT_iss, data = df_high, alternative = "greater")
df_high <- all_patient_df %>% filter(expr_level=="HIGH")

df_high$x_name <- factor(df_high$x_name, levels = c("HIGH_1","HIGH_2","HIGH_3"),
                                labels = c("High_ISS1","High_ISS2","High_ISS3"))

pv <- paste0("Jonckheere-test, p=",round(j_test_high$p.value,4))

p <- ggviolin(df_high,x = "x_name", y = "NRF1", fill = "expr_level", color = "#696969",
              palette = c("#EFBF9A"),width = 0.6,outlier.shape = NA,
               lwd=0.5,xlab="ISS",ylab="Normalized reads counts(tpm)" ,
               alpha = 0.5,
         bxp.errorbar =TRUE,
         add = c("boxplot"),add.params = list(size = 0.5, width = 0.2, alpha = 0.5))+
scale_y_continuous(limits = c(0, 18),breaks = seq(0, 18, by = 2))+
      theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) + border()+
  xlab("")+  labs(fill = "Expression Groups")  +
  annotate("text", x = "Low_ISS2", y = 17, label = pv, size = 4, color = "black")
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/TPM_per_ISS_in_population_divided_per_NRF1_expression_vertical_boxplot_per_iss_with_statistics_ONLY_HIGH_WITH_jonckheereTest_INCRESING_WITH_ISS_25feb25.png"
alt="Fig2f_violin_plot_survival_per_ISS_group_HIGH" />
<figcaption
aria-hidden="true">Fig2f_violin_plot_survival_per_ISS_group_HIGH</figcaption>
</figure>
