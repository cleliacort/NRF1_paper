# Multivariate Analysis

Retrieve and customize clinical information.

``` r
survival <- read.csv("../Fig4/data/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv", sep = "\t", header = T)
survival <- survival %>%dplyr::select(PUBLIC_ID,oscdy,censos) 
ids <- read.delim("../Fig4/data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv" , sep = "\t", header = T) %>% 
  dplyr::select(Specimen_ID,PUBLIC_ID,VISITDY) %>% distinct()
iss <- read.delim("../Fig4/data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% 
  dplyr::select(PUBLIC_ID,D_PT_iss,D_PT_age, D_PT_gender,D_PT_therclass,D_PT_lvisit,D_PT_trtstdy,sctflag,ecog)
```

Retrieve and customize mutation and translocation information.

``` r
mut_and_trasloc <- read.delim("../Fig4/data/MMRF_CoMMpass_A17_merged_locally_mutations_and_traslocation.txt", sep = "\t", header = T) %>% 
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
                mutationHIST1H1E, mutationLTB,
                mutationMAFA,
                mutationMAFB,
                mutationMAF,
                mutationCCND1,
                mutationSP140,
                mutationCYLD,mutationHUWE1,
                mutationEGR1,mutationPRKD2,
                mutationTRAF3,
                mutationFAT3,mutationIGLL5,
                mutationATM,mutationNRAS,
                mutationACTG1,
                mutationDUSP2,mutationMAX,
                mutationATM,mutationFAM46C
)
```

Upload the patient cluster classification

``` r
ref_heatmap_class <- paste0("/Users/cleliacortile/Desktop/Lavoro/00_PAPER_NRF1/Fig3/data/Suppl/","patient_ID_with_cluster_motif_and_iss_info_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_1023.txt")
ref_heatmap_class_file <- read.table(ref_heatmap_class, header = T) %>% select(-D_PT_iss)
colnames(ref_heatmap_class_file) <- c("Patient_ID","classification_ref")
```

Merge clinical and genetic information in a unique dataframe.

``` r
ids_iss <- merge(ids,iss,by="PUBLIC_ID") %>% mutate(Patient_ID=paste(Specimen_ID,"CD138pos",sep="_"))
surv_ids_iss <- merge(survival,ids_iss,by="PUBLIC_ID", all.x = TRUE)
mut_trasloc_surv_stage_iss <- merge(surv_ids_iss,mut_and_trasloc,by.x = "Patient_ID", by.y = "SAMPLE")
mut_trasloc_surv_stage_iss <- mut_trasloc_surv_stage_iss %>% filter(!is.na(D_PT_iss)) #%>% filter(oscdy>=0)
class_mut_trasloc_surv_stage_iss <- merge(mut_trasloc_surv_stage_iss,ref_heatmap_class_file,by="Patient_ID")
```

``` r
SeqFISH <- read_delim("data/SeqFISHFiles_MMRF_CoMMpass_IA17_genome_gatk_cna_seqFISH.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE) %>%  
  select(SAMPLE,matches("1p22|1q21|13q14|17p13|Hyper"))%>% 
  select(-matches("percent"))%>% 
  separate(SAMPLE,c("Visit_ID"),sep="_CD138pos") %>% 
   mutate(del_1p22_manual_Call= case_when(SeqWGS_Cp_1p22 <= -0.3 ~ 1 , TRUE ~ 0),
          gain_1q21_manual_Call= case_when(SeqWGS_Cp_1q21 >= 0.3 ~ 1 , TRUE ~ 0),
          del_13q14_manual_Call= case_when(SeqWGS_Cp_13q14 <= -0.3 ~ 1 , TRUE ~ 0),
          del_17p13_manual_Call= case_when(SeqWGS_Cp_17p13 >= 0.3 ~ 1 , TRUE ~ 0))

surv_stage_iss_mut_trasloc <- merge(class_mut_trasloc_surv_stage_iss,SeqFISH,by.x = "Specimen_ID", by.y = "Visit_ID", all.x = TRUE )
```

Rename clnical features for consistency

``` r
info_multi <- surv_stage_iss_mut_trasloc %>% rename(ECOG_old=ecog , 
                                                    Age=D_PT_age,
                                                    ASCT=sctflag, 
                                                    Gender=D_PT_gender,
                                                    ISS=D_PT_iss,
                                                    Therapy=D_PT_therclass              
                                                    ) %>% 
  mutate(ECOG=case_when(ECOG_old==0 ~ "full",ECOG_old==1  ~ "moderate", ECOG_old>=2 ~ "non-active")) 
```

Exclude the population receiving rare treatment, as they are too small
to provide any significant hazard ratio (HR).

``` r
info_multi_clean <- info_multi %>% filter(Therapy!="Carfilzomib-based") %>% 
  filter( Therapy!="combined bortezomib/carfilzomib-based") %>% 
  filter(Therapy!="combined daratumumab/IMIDs/carfilzomib-based")
```

Define the levels of analyzed variables.

``` r
info_multi_clean <- within(info_multi_clean, {
  classification_ref <- factor(classification_ref,levels = c("cluster1", "cluster2"), labels = c("low","high"))
  classification_ref = relevel(classification_ref, ref = "low")
  
  ASCT <- factor(ASCT,levels = c(0, 1), labels = c("No","Yes"))
  ASCT = relevel(ASCT, ref = "No")
  
  Gender <- factor(Gender, labels = c("M","F"))
  Gender = relevel(Gender, ref = "M")
  
  ISS <- factor(ISS,levels = c(1,2,3), labels = c("1","2","3"))
  
  mutationTP53 <- factor(mutationTP53,levels = c(0, 1), labels = c("No","Yes"))
  mutationFGFR3 <- factor(mutationFGFR3,levels = c(0, 1), labels = c("No","Yes"))
  mutationDIS3 <- factor(mutationDIS3,levels = c(0, 1), labels = c("No","Yes"))
  mutationATM <- factor(mutationATM,levels = c(0, 1), labels = c("No","Yes"))
  mutationKRAS <- factor(mutationKRAS,levels = c(0, 1), labels = c("No","Yes"))
  mutationBRAF <- factor(mutationBRAF,levels = c(0, 1), labels = c("No","Yes"))
  mutationHIST1H1E <- factor(mutationHIST1H1E,levels = c(0, 1), labels = c("No","Yes"))
  mutationLTB <- factor(mutationLTB, levels = c(0, 1),labels = c("No","Yes"))
  mutationMAFA <- factor(mutationMAFA, levels = c(0, 1),labels = c("No","Yes"))
  mutationMAF <- factor(mutationMAF, levels = c(0, 1),labels = c("No","Yes"))
  mutationMAFB <- factor(mutationMAFB, levels = c(0, 1),labels = c("No","Yes"))
  mutationCCND1 <- factor(mutationCCND1, levels = c(0, 1),labels = c("No","Yes"))
  
  mutationSP140 <- factor(mutationSP140, levels = c(0, 1),labels = c("No","Yes"))
  mutationCYLD <- factor(mutationCYLD, levels = c(0, 1),labels = c("No","Yes"))
  mutationHUWE1 <- factor(mutationHUWE1, levels = c(0, 1),labels = c("No","Yes"))
  mutationEGR1 <- factor(mutationEGR1, levels = c(0, 1),labels = c("No","Yes"))
  mutationPRKD2 <- factor(mutationPRKD2, levels = c(0, 1),labels = c("No","Yes"))
  mutationTRAF3 <- factor(mutationTRAF3, levels = c(0, 1),labels = c("No","Yes"))
  mutationFAT3 <- factor(mutationFAT3, levels = c(0, 1),labels = c("No","Yes"))
  mutationIGLL5 <- factor(mutationIGLL5, levels = c(0, 1),labels = c("No","Yes"))
  mutationNRAS <- factor(mutationNRAS, levels = c(0, 1),labels = c("No","Yes"))
  mutationACTG1 <- factor(mutationACTG1, levels = c(0, 1),labels = c("No","Yes"))
  mutationDUSP2 <- factor(mutationDUSP2, levels = c(0, 1),labels = c("No","Yes"))
  mutationMAX <- factor(mutationMAX, levels = c(0, 1),labels = c("No","Yes"))
  mutationFAM46C <- factor(mutationFAM46C, levels = c(0, 1),labels = c("No","Yes"))


  traslocationNSD2_CALL <- factor(traslocationNSD2_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND3_CALL <- factor(traslocationCCND3_CALL, levels = c(0, 1),labels = c("No","Yes"))
  traslocationMYC_CALL<- factor(traslocationMYC_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationMAFA_CALL<- factor(traslocationMAFA_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND1_CALL<- factor(traslocationCCND1_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND2_CALL<- factor(traslocationCCND2_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationMAF_CALL<- factor(traslocationMAF_CALL, levels = c(0, 1),labels = c("No","Yes"))
  traslocationMAFB_CALL<- factor(traslocationMAFB_CALL, levels = c(0, 1),labels = c("No","Yes"))

  del_1p22_manual_Call <- factor(del_1p22_manual_Call,levels = c(0, 1), labels = c("No","Yes"))
  gain_1q21_manual_Call <- factor(gain_1q21_manual_Call, levels = c(0, 1),labels = c("No","Yes"))
  del_13q14_manual_Call<- factor(del_13q14_manual_Call,levels = c(0, 1), labels = c("No","Yes"))
  del_17p13_manual_Call<- factor(del_17p13_manual_Call,levels = c(0, 1), labels = c("No","Yes"))
  SeqWGS_Cp_Hyperdiploid_Call<- factor(SeqWGS_Cp_Hyperdiploid_Call,levels = c(0, 1), labels = c("No","Yes"))

  
  Therapy <- factor(Therapy, levels = c("Bortezomib-based","combined bortezomib/IMIDs-based",
                                        "combined bortezomib/IMIDs/carfilzomib-based","combined IMIDs/carfilzomib-based",
                                        "IMIDs-based" ))
  Therapy=relevel(Therapy, ref = "Bortezomib-based")
  Therapy <- factor(Therapy, labels = c("V-based","combo K/IMIDs-based",
                                        "combo V/IMIDs/K-based","combo IMIDs/K-based",
                                        "IMIDs-based"))
})
```

Compute the univariate of selected variables

``` r
options(scipen = 0, digits = 2)
covariates <- c( "ISS","Age","Gender","Therapy","ASCT","ECOG","classification_ref", "traslocationNSD2_CALL","traslocationCCND3_CALL",
                  "traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL","traslocationCCND2_CALL",
                  "traslocationMAF_CALL","traslocationMAFB_CALL","mutationTP53","mutationFGFR3",
                  "mutationDIS3","mutationKRAS","mutationBRAF",
                  "mutationHIST1H1E","mutationLTB","mutationMAFA","mutationMAFB","mutationMAF","mutationCCND1",
                  "mutationSP140",
                  "mutationCYLD","mutationHUWE1",
                  "mutationEGR1","mutationPRKD2",
                  "mutationTRAF3",
                  "mutationFAT3","mutationIGLL5",
                  "mutationATM","mutationNRAS",
                  "mutationACTG1",
                  "mutationDUSP2","mutationMAX",
                  "mutationATM","mutationFAM46C",
                  "del_1p22_manual_Call","gain_1q21_manual_Call",
                  "del_13q14_manual_Call","del_17p13_manual_Call","SeqWGS_Cp_Hyperdiploid_Call")

df_results <- data.frame(variable="remove_it",beta=0,HR=0,HR_lower=0,HR_upper=0,HR_inter=0,wd_pvalue=0,log_rank_pvalue=0)

for (i in 1:length(covariates)){
  print(covariates[i])
  # run the coxph univariate analysis for each variable 
  formula_cox <- as.formula(paste0("Surv(oscdy, censos) ~ ", covariates[i]))
  os_multi<-coxph(formula_cox,data = info_multi_clean,na.action = na.exclude)
  # save the statistical results of univariate analysis
  sum <- summary(os_multi)
  beta <- sum$coefficients[,1]
  # save beta value
  beta <- as.data.frame(t(t(beta))) %>% rownames_to_column(var="variable") %>% rename(beta=V1)
  if("1" %in% beta$variable){beta$variable=covariates[i]}
  #save hazard ratio
  HR <- signif(sum$coefficients[,2],2)
  HR<- as.data.frame(t(t(HR))) %>% rownames_to_column(var="variable") %>% rename(HR=V1)
  if("1" %in% HR$variable){HR$variable=covariates[i]}
  #save the lower limit for the hazard ratio 
  HR_lower <- signif(sum$conf.int[,"lower .95"],2)
  HR_lower<- as.data.frame(t(t(HR_lower))) %>% rownames_to_column(var="variable")%>% rename(HR_lower=V1)
  if("1" %in% HR_lower$variable){HR_lower$variable=covariates[i]}
  #save the upper limit for the hazard ratio 
  HR_upper <- signif(sum$conf.int[,"upper .95"],2)
  HR_upper<- as.data.frame(t(t(HR_upper))) %>% rownames_to_column(var="variable")%>% rename(HR_upper=V1)
  if("1" %in% HR_upper$variable ){HR_upper$variable=covariates[i]}
  #merge the info in a unique df 
  HR <- merge(HR,HR_lower,by="variable")
  HR <- merge(HR,HR_upper,by="variable")
  HR <- HR %>% mutate(HR_inter=paste0("(",HR_lower,"-",HR_upper,")",sep=""))
  # retrieve the significance
  wd <- signif(sum$waldtest["pvalue"],2)
  pvalue <- signif(sum$sctest["pvalue"],2)
  # add the significance to the final df
  HR <- HR %>% mutate(wd_pvalue=wd, log_rank_pvalue=pvalue)
  HR <- merge(beta,HR,by="variable")
  df_results <- rbind(df_results,HR)
}
```

    ## [1] "ISS"
    ## [1] "Age"
    ## [1] "Gender"
    ## [1] "Therapy"
    ## [1] "ASCT"
    ## [1] "ECOG"
    ## [1] "classification_ref"
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
    ## [1] "mutationKRAS"
    ## [1] "mutationBRAF"
    ## [1] "mutationHIST1H1E"
    ## [1] "mutationLTB"
    ## [1] "mutationMAFA"
    ## [1] "mutationMAFB"
    ## [1] "mutationMAF"
    ## [1] "mutationCCND1"
    ## [1] "mutationSP140"
    ## [1] "mutationCYLD"
    ## [1] "mutationHUWE1"
    ## [1] "mutationEGR1"
    ## [1] "mutationPRKD2"
    ## [1] "mutationTRAF3"
    ## [1] "mutationFAT3"
    ## [1] "mutationIGLL5"
    ## [1] "mutationATM"
    ## [1] "mutationNRAS"
    ## [1] "mutationACTG1"
    ## [1] "mutationDUSP2"
    ## [1] "mutationMAX"
    ## [1] "mutationATM"
    ## [1] "mutationFAM46C"
    ## [1] "del_1p22_manual_Call"
    ## [1] "gain_1q21_manual_Call"
    ## [1] "del_13q14_manual_Call"
    ## [1] "del_17p13_manual_Call"
    ## [1] "SeqWGS_Cp_Hyperdiploid_Call"

Create a table with univariate results

``` r
options(scipen = 100, digits = 4)
name_table=paste0("figures/","univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_REVISION_v2_0325.pdf")

name_table=paste0("figures/","univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_REVISION_v2_0625.pdf")

custom_results <- df_results %>% filter(variable!="remove_it") %>% arrange(log_rank_pvalue) %>% 
  mutate(status=case_when(log_rank_pvalue<=0.05 ~ "signif", TRUE ~ "NO_signif")) %>% 
  mutate(log_rank_pvalue = formatC(log_rank_pvalue, format = "e", digits = 2) ) %>%
  mutate(wd_pvalue = formatC(wd_pvalue, format = "e", digits = 2) ) %>%
  mutate(variable = case_when(
    variable == "traslocationNSD2_CALL" ~ "t(4;14)",  
    variable == "traslocationCCND3_CALL" ~ "t(6;14)", 
        variable == "traslocationMYC_CALL" ~ "t(MYC)", 
    variable == "traslocationMAFA_CALL" ~ "t(8;14)", 
        variable == "traslocationCCND1_CALL" ~ "t(11;14)", 
    variable == "traslocationMAF_CALL" ~ "t(14;16)", 
        variable == "traslocationMAFB_CALL" ~ "t(14;20)",  
    TRUE ~ variable  
  ))
```

``` r
custom_results %>% kbl() %>% kable_styling() 
custom_results %>% kbl() %>% kable_styling() %>% save_kable(name_table)
write.xlsx(custom_results, "data/Suppl/univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_0124.xlsx",colNames = TRUE)
write.xlsx(custom_results, "figures/univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_REVISION_v2_0625.xlsx",colNames = TRUE)
```

Count the number of patient for each group and the varianza of censos
info (just for check)

``` r
options(scipen = 0, digits = 2)
covariates <- c(  "traslocationNSD2_CALL","traslocationCCND3_CALL",
                  "traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL","traslocationCCND2_CALL",
                  "traslocationMAF_CALL","traslocationMAFB_CALL","mutationTP53","mutationFGFR3",
                  "mutationDIS3","mutationKRAS","mutationBRAF",
                  "mutationHIST1H1E","mutationLTB","mutationMAFA","mutationMAFB","mutationMAF","mutationCCND1",
                  "mutationSP140",
                  "mutationCYLD","mutationHUWE1",
                  "mutationEGR1","mutationPRKD2",
                  "mutationTRAF3",
                  "mutationFAT3","mutationIGLL5",
                  "mutationATM","mutationNRAS",
                  "mutationACTG1",
                  "mutationDUSP2","mutationMAX","mutationFAM46C",
                  "del_1p22_manual_Call","gain_1q21_manual_Call",
                  "del_13q14_manual_Call","del_17p13_manual_Call","SeqWGS_Cp_Hyperdiploid_Call")

df_results_population <- data.frame(variable="",samples_yes="",samples_no="",variability_censos_YES="",variability_censos_NO="")

# Customized way to organize the results
for (i in 1:length(covariates)){
        print(covariates[i])
        var_name <- covariates[i]
        sample_size_yes <- sum(info_multi_clean[[var_name]] == "Yes", na.rm = TRUE)
        sample_size_no <- sum(info_multi_clean[[var_name]] == "No", na.rm = TRUE)
        df_yes <- info_multi_clean %>% select(!!var_name,censos) %>% filter(info_multi_clean[[var_name]]=="Yes")
        df_no <- info_multi_clean %>% select(!!var_name,censos) %>% filter(info_multi_clean[[var_name]]=="No")
        
        # Compute proportion of "Yes" (1)
        p1 <- mean(df_yes$censos)
        p2 <- mean(df_no$censos)
        
        df_sub <- data.frame(variable=var_name,samples_yes=sample_size_yes,samples_no=sample_size_no, variability_censos_YES=p1,variability_censos_NO=p2)
        
        df_results_population <- rbind(df_results_population,df_sub)
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
    ## [1] "mutationKRAS"
    ## [1] "mutationBRAF"
    ## [1] "mutationHIST1H1E"
    ## [1] "mutationLTB"
    ## [1] "mutationMAFA"
    ## [1] "mutationMAFB"
    ## [1] "mutationMAF"
    ## [1] "mutationCCND1"
    ## [1] "mutationSP140"
    ## [1] "mutationCYLD"
    ## [1] "mutationHUWE1"
    ## [1] "mutationEGR1"
    ## [1] "mutationPRKD2"
    ## [1] "mutationTRAF3"
    ## [1] "mutationFAT3"
    ## [1] "mutationIGLL5"
    ## [1] "mutationATM"
    ## [1] "mutationNRAS"
    ## [1] "mutationACTG1"
    ## [1] "mutationDUSP2"
    ## [1] "mutationMAX"
    ## [1] "mutationFAM46C"
    ## [1] "del_1p22_manual_Call"
    ## [1] "gain_1q21_manual_Call"
    ## [1] "del_13q14_manual_Call"
    ## [1] "del_17p13_manual_Call"
    ## [1] "SeqWGS_Cp_Hyperdiploid_Call"

``` r
df_big <- merge(custom_results,df_results_population,by="variable",all.x = TRUE) %>% arrange(status)
```

``` r
name_table=paste0("figures/","POPOLAZIONE_MUT_TRASL_GAIN_LOSS_univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_REVISION_v2_0325.pdf")
df_big %>% kbl() %>% kable_styling() %>% save_kable(name_table)
```

Calculate the multivariate using variables significant at univariate
analysis.

``` r
os_multi<-coxph(Surv(oscdy, censos) ~ classification_ref + 
                  ISS + Age + Gender+ ASCT + Therapy + 
                  
                  traslocationNSD2_CALL+
                  traslocationCCND1_CALL+
                  traslocationMAF_CALL+
                  traslocationMYC_CALL+
                  
                  SeqWGS_Cp_Hyperdiploid_Call+
                  
                  del_1p22_manual_Call+
                  gain_1q21_manual_Call+
                  del_13q14_manual_Call+
                  del_17p13_manual_Call,
                  data = info_multi_clean, na.action = na.exclude)
```

Make a plot showing the results of the multivariate analysis

``` r
gg_multivar <- ggforest(
  os_multi,
  data = info_multi_clean,
  main = "OS HR selected variables 4 + therapy",
  cpositions = c(0.02, 0.22, 0.39),
  fontsize = 1,
  noDigits = 2
)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig5/figures/Suppl/multivariata_ALL_SIGN_classification_ref_MUT_AND_TRASLC_OF_INTEREST_WITH_MYC_15APRIL25.png"
alt="Suppl3_Multivariate_analysis" />
<figcaption aria-hidden="true">Suppl3_Multivariate_analysis</figcaption>
</figure>
