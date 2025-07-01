# Suppl_4a_bubbleplot_iss_genetics

Upload the needed packages.

Retrieve the clinical (disease stage) and survival information

``` r
# Upload the survival data 
survival <- read.csv("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig4/data/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv", sep = "\t", header = T)
survival <- survival %>%dplyr::select(PUBLIC_ID,oscdy,censos) 

# Retrieve the ISS information
stage <- read.delim("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig4/data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv" , sep = "\t", header = T) %>% 
  dplyr::select(Specimen_ID,PUBLIC_ID) %>% distinct()
iss <- read.delim("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig4/data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% 
  dplyr::select(PUBLIC_ID,D_PT_iss)
stage_iss <- merge(stage,iss,by="PUBLIC_ID") %>% mutate(Patient_ID=paste(Specimen_ID,"CD138pos",sep="_"))

# Merge iss with survival data to select only those having survival info
surv_stage_iss <- merge(survival,stage_iss,by="PUBLIC_ID", all.x = TRUE)
```

Retrieve the genetic (mutation and translocation) information

``` r
mut_and_trasloc <- read.delim("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig4/data/MMRF_CoMMpass_A17_merged_locally_mutations_and_traslocation.txt", sep = "\t", header = T) %>% 
  dplyr::select(SAMPLE,
                traslocationNSD2_CALL,
                traslocationCCND3_CALL,
                traslocationMYC_CALL,
                traslocationMAFA_CALL,
                traslocationCCND1_CALL,
                traslocationCCND2_CALL,
                traslocationMAF_CALL,
                traslocationMAFB_CALL,
                mutationTP53,mutationFGFR3,
                mutationDIS3,mutationATM,
                mutationKRAS,mutationBRAF,
                mutationHIST1H1E, mutationLTB)
```

Merge the clinical, survival, and genetic information into a single
dataframe

``` r
surv_stage_iss_mut_trasloc <- merge(surv_stage_iss,mut_and_trasloc,by.x = "Patient_ID", by.y = "SAMPLE")
surv_stage_iss_mut_trasloc <- surv_stage_iss_mut_trasloc %>% filter(!is.na(D_PT_iss)) %>% separate(Specimen_ID,c("Visit_ID"),sep="_BM")
```

Add new genetic features (if available) and merge with the main
phenotype data

``` r
genetic_NEW_INFO <- read_delim("data/genetic_NEW_INFO_added_to_those_used_in_our_paper_03APRILE.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE) %>% select(-traslocationNSD2_CALL,-traslocationCCND3_CALL, -traslocationMYC_CALL,-traslocationMAFA_CALL,-traslocationCCND1_CALL,-traslocationCCND2_CALL,-traslocationMAF_CALL,-traslocationMAFB_CALL)

pheno_order <- merge(surv_stage_iss_mut_trasloc,genetic_NEW_INFO,by="Visit_ID")
```

Load cluster assignments for two groups (C1 and C2), label them, and
combine

``` r
C1 <- read_delim("data/cluster_motif_1_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE) %>% mutate(cluster_group="C1")

C2 <- read_delim("data/cluster_motif_2_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE) %>% mutate(cluster_group="C2")

df_results <- rbind(C1,C2) %>% rename(public_id=X1)
```

Merge cluster assignments with phenotype/genetic data for downstream
analysis

``` r
df_bubble_plot <- merge(pheno_order,df_results,by.x = "Patient_ID", by.y = "public_id")
df_bubble_plot_count <- df_bubble_plot %>% group_by(cluster_group,D_PT_iss,traslocationNSD2_CALL) %>% count() %>% filter(traslocationNSD2_CALL==1)
column_vector <- c("traslocationNSD2_CALL", "traslocationCCND3_CALL", "traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL",
                   "traslocationCCND2_CALL", "traslocationMAF_CALL", "traslocationMAFB_CALL", "mutationTP53", "mutationFGFR3", "mutationDIS3", 
                   "mutationATM", "mutationKRAS", "mutationBRAF", "mutationHIST1H1E", "mutationLTB")#,
```

Prepare a list of all genetic features to analyze and compute
statistical analysis for each genetic feature, test if its frequency
differs between ISS1 and ISS3

``` r
column_vector <- c("traslocationNSD2_CALL", "traslocationCCND3_CALL", "traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL",
                   "traslocationCCND2_CALL", "traslocationMAF_CALL", "traslocationMAFB_CALL", "mutationTP53", "mutationFGFR3", "mutationDIS3", 
                   "mutationATM", "mutationKRAS", "mutationBRAF", "mutationHIST1H1E", "mutationLTB")#,

df_bubble_grouped <- data.frame()
df_sign_all_variable <- data.frame()  

for (i in seq_along(column_vector)) {
  
  print(column_vector[i])
  
  df_bubble_plot_count <- df_bubble_plot %>% group_by(D_PT_iss,!!sym(column_vector[i])) %>% 
    count() %>% filter(!!sym(column_vector[i])==1) %>% mutate(genetic_char=column_vector[i]) %>% 
    rename(genetic_status=!!sym(column_vector[i]))
  
  df_bubble_grouped <- rbind(df_bubble_grouped,df_bubble_plot_count)

  stage1_3 <- df_bubble_plot_count[df_bubble_plot_count$D_PT_iss %in% c(1, 3), ] %>%
    mutate(tot_pop_size=length(df_bubble_plot$Patient_ID))
  
  df_sign <- prop.test(stage1_3$n, stage1_3$tot_pop_size)
  pv <- round(df_sign$p.value,2)
  
  df_sign_small <- data.frame(genetic_char=column_vector[i],pvalue=pv,prop_iss1=df_sign$estimate[1],prop_iss3=df_sign$estimate[2]) %>% mutate(direction=case_when(prop_iss3 >     prop_iss1 ~ "increase", prop_iss3<prop_iss1 ~ "decrease",TRUE ~ "NAN" ))
  
  df_sign_all_variable <- rbind(df_sign_all_variable,df_sign_small)
}
```

    ## [1] "traslocationNSD2_CALL"
    ## [1] "traslocationCCND3_CALL"
    ## [1] "traslocationMYC_CALL"
    ## [1] "traslocationMAFA_CALL"
    ## [1] "traslocationCCND1_CALL"
    ## [1] "traslocationCCND2_CALL"
    ## [1] "traslocationMAF_CALL"
    ## [1] "traslocationMAFB_CALL"
    ## [1] "mutationTP53"
    ## [1] "mutationFGFR3"
    ## [1] "mutationDIS3"
    ## [1] "mutationATM"
    ## [1] "mutationKRAS"
    ## [1] "mutationBRAF"
    ## [1] "mutationHIST1H1E"
    ## [1] "mutationLTB"

``` r
df_overall_plot <- merge(df_sign_all_variable,df_bubble_grouped,by="genetic_char")
```

Create the bubble plot

``` r
df_final_plot <- df_overall_plot %>% filter(D_PT_iss==1) %>%  mutate(D_PT_iss=as.factor(D_PT_iss)) %>% 
   mutate(pvalue_safe = ifelse(pvalue == 0, 1e-3, pvalue),
         minuslog2pv = -log2(pvalue_safe))

ggballon <- ggplot(df_final_plot,aes(x=D_PT_iss, y=genetic_char)) + 
  geom_point(aes(size = minuslog2pv, fill=direction),shape = 21, color = "#595959")+
     scale_size(name="PValue(proportion test)", range = c(1,10), breaks=c(0,3,4,8)) +#,labels=c(">=0",">=1",">=10",">=100",">=1000",">=1000000"),guide="legend")
  scale_fill_manual(name="Population change direction",
                    values = c("increase" = "#D46780FF", "decrease" = "#798234FF", "NAN"="grey"),  
                    labels = c("increase" = "increase_C1_to_C3", "decrease" = "decrease_C1_to_C3"))+
  geom_text_repel(aes(label = pvalue, 3), size = 3, max.overlaps = Inf) +
  ylab("")+xlab("")+
  theme_bw()+
      theme(
    axis.text.x = element_text(color = "white"),
    panel.grid.major = element_line(color = "white", linewidth = 0.25), 
    axis.ticks.x = element_line(color = "white", linewidth = 0.5)
  )
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/ISS1_VS_ISS3_statistic_proportional_test_verticale_15APRILE25_v2.png"
alt="Suppl_4a_bubbleplot_iss_genetics" />
<figcaption
aria-hidden="true">Suppl_4a_bubbleplot_iss_genetics</figcaption>
</figure>
