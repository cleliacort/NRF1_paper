Retrieve and customize clinical information.

``` r
survival <- read.csv("data/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv", sep = "\t", header = T)
survival <- survival %>%dplyr::select(PUBLIC_ID,oscdy,censos) 

ids <- read.delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv" , sep = "\t", header = T) %>% 
  dplyr::select(Specimen_ID,PUBLIC_ID,VISITDY) %>% distinct()

iss <- read.delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% 
  dplyr::select(PUBLIC_ID,D_PT_iss,D_PT_age, D_PT_gender,D_PT_therclass,D_PT_lvisit,D_PT_trtstdy,sctflag,ecog)
```

Retrieve and customize mutation and translocation information.

``` r
mut_and_transloc_to_select <- c("SAMPLE","mutationFGFR3","mutationKRAS","mutationTP53","mutationLTB","mutationSP140","mutationDIS3","mutationCYLD","mutationHUWE1","mutationEGR1",
                                "mutationPRKD2","mutationTRAF3","mutationBRAF","mutationFAT3","mutationIGLL5","mutationATM","mutationNRAS","mutationACTG1","mutationHIST1H1E","mutationDUSP2",
                                "mutationMAX","mutationATM","mutationFAM46C",
                                "traslocationNSD2_CALL","traslocationCCND3_CALL","traslocationMYC_CALL","traslocationMAFA_CALL","traslocationCCND1_CALL","traslocationCCND2_CALL",
                                "traslocationMAF_CALL","traslocationMAFB_CALL")
mut_and_transloc <- read.delim("data/MMRF_CoMMpass_A17_merged_locally_mutations_and_traslocation.txt",sep = "\t") %>% select(!!mut_and_transloc_to_select) 
```

Upload the patient cluster classification

``` r
ref_heatmap_class <- paste0("data/Suppl/","patient_ID_with_cluster_motif_and_iss_info_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_1023.txt")

ref_heatmap_class_file <- read.table(ref_heatmap_class, header = T) %>% select(-D_PT_iss)
colnames(ref_heatmap_class_file) <- c("Patient_ID","classification_ref")
```

Merge clinical and genetic information in a unique dataframe.

``` r
# Merge 
ids_iss <- merge(ids,iss,by="PUBLIC_ID") %>% mutate(Patient_ID=paste(Specimen_ID,"CD138pos",sep="_"))

# Merge iss with survival data to select only those having survival info
surv_ids_iss <- merge(survival,ids_iss,by="PUBLIC_ID", all.x = TRUE)

# Merge the survival-iss information with th mut-traslocation info 
mut_trasloc_surv_stage_iss <- merge(surv_ids_iss,mut_and_transloc,by.x = "Patient_ID", by.y = "SAMPLE")

# Filter out patient that do not have ISS info
mut_trasloc_surv_stage_iss <- mut_trasloc_surv_stage_iss %>% filter(!is.na(D_PT_iss)) %>% filter(oscdy>=0)
class_mut_trasloc_surv_stage_iss <- merge(mut_trasloc_surv_stage_iss,ref_heatmap_class_file,by="Patient_ID")
```

Rename clnical features for consistency

``` r
info_multi <- class_mut_trasloc_surv_stage_iss %>% rename(ECOG_old=ecog , 
                                                    Age=D_PT_age,
                                                    ASCT=sctflag, 
                                                    Gender=D_PT_gender,
                                                    ISS=D_PT_iss,
                                                    Therapy=D_PT_therclass) %>% 
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
  
  mutationFGFR3 <- factor(mutationFGFR3,levels = c(0, 1), labels = c("No","Yes"))
  mutationDIS3 <- factor(mutationDIS3,levels = c(0, 1), labels = c("No","Yes"))
  mutationTP53 <- factor(mutationTP53,levels = c(0, 1), labels = c("No","Yes"))
  mutationATM <- factor(mutationATM,levels = c(0, 1), labels = c("No","Yes"))
  mutationKRAS <- factor(mutationKRAS,levels = c(0, 1), labels = c("No","Yes"))
  mutationBRAF <- factor(mutationBRAF,levels = c(0, 1), labels = c("No","Yes"))
  mutationHIST1H1E <- factor(mutationHIST1H1E,levels = c(0, 1), labels = c("No","Yes"))
  mutationLTB <- factor(mutationLTB, levels = c(0, 1),labels = c("No","Yes"))

  mutationSP140 <- factor(mutationSP140, levels = c(0, 1),labels = c("No","Yes"))
  mutationCYLD <- factor(mutationCYLD, levels = c(0, 1),labels = c("No","Yes"))
  mutationHUWE1 <- factor(mutationHUWE1, levels = c(0, 1),labels = c("No","Yes"))
  mutationEGR1 <- factor(mutationEGR1, levels = c(0, 1),labels = c("No","Yes"))
  mutationPRKD2 <- factor(mutationPRKD2, levels = c(0, 1),labels = c("No","Yes"))
  mutationTRAF3 <- factor(mutationTRAF3, levels = c(0, 1),labels = c("No","Yes"))
  mutationFAT3 <- factor(mutationFAT3, levels = c(0, 1),labels = c("No","Yes"))
  mutationIGLL5 <- factor(mutationIGLL5, levels = c(0, 1),labels = c("No","Yes"))
  mutationNRAS <- factor(mutationNRAS, levels = c(0, 1),labels = c("No","Yes"))
  mutationDUSP2 <- factor(mutationDUSP2, levels = c(0, 1),labels = c("No","Yes"))
  mutationMAX <- factor(mutationMAX, levels = c(0, 1),labels = c("No","Yes"))
  mutationFAM46C <- factor(mutationFAM46C, levels = c(0, 1),labels = c("No","Yes"))
  mutationACTG1 <- factor(mutationACTG1, levels = c(0, 1),labels = c("No","Yes"))

  traslocationNSD2_CALL <- factor(traslocationNSD2_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND3_CALL <- factor(traslocationCCND3_CALL, levels = c(0, 1),labels = c("No","Yes"))
  traslocationMYC_CALL<- factor(traslocationMYC_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationMAFA_CALL<- factor(traslocationMAFA_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND1_CALL<- factor(traslocationCCND1_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationCCND2_CALL<- factor(traslocationCCND2_CALL,levels = c(0, 1), labels = c("No","Yes"))
  traslocationMAF_CALL<- factor(traslocationMAF_CALL, levels = c(0, 1),labels = c("No","Yes"))
  traslocationMAFB_CALL<- factor(traslocationMAFB_CALL, levels = c(0, 1),labels = c("No","Yes"))

  Therapy <- factor(Therapy, levels = c("Bortezomib-based","combined bortezomib/IMIDs-based",
                                        "combined bortezomib/IMIDs/carfilzomib-based","combined IMIDs/carfilzomib-based",
                                        "IMIDs-based" ))
  Therapy=relevel(Therapy, ref = "Bortezomib-based")
  Therapy <- factor(Therapy, labels = c("V-based","combo K/IMIDs-based",
                                        "combo V/IMIDs/K-based","combo IMIDs/K-based",
                                        "IMIDs-based"))
})
```

Compute the univariate of selected variables.

``` r
options(scipen = 0, digits = 2)

# Choosing variable
covariates <- c( "ISS","Age","Gender","Therapy","ASCT","ECOG","classification_ref", mut_and_transloc_to_select[-1])

df_results <- data.frame(variable="remove_it",beta=0,HR=0,HR_lower=0,HR_upper=0,HR_inter=0,wd_pvalue=0,log_rank_pvalue=0)

# Customized way to organize the results
for (i in 1:length(covariates)){
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

Create a table with univariate results

``` r_table
name_table=paste0("figures/","univariate_table_patient_ID_based_cluster_heatmap_COMMPASS_derived_from_high_correlated_cluster1_signature_paper_features_selection_0124.pdf")

custom_results <- df_results %>% filter(variable!="remove_it") %>% arrange(log_rank_pvalue) %>% 
  mutate(status=case_when(log_rank_pvalue<0.05 ~ "signif",TRUE ~ "NO_signif")) %>% 
  mutate(log_rank_pvalue = formatC(log_rank_pvalue, format = "e", digits = 2) ) %>%
  mutate(wd_pvalue = formatC(wd_pvalue, format = "e", digits = 2) ) #%>%

custom_results %>% kbl() %>% kable_styling() 

custom_results %>% kbl() %>% kable_styling() %>% save_kable(name_table)
```

Calculate the multivariate using variables significant at univariate
analysis.

``` r
os_multi<-coxph(Surv(oscdy, censos) ~ classification_ref + ISS + Age + Gender+ ASCT + Therapy + 
                  mutationTP53+
                  mutationFGFR3+
                  mutationTRAF3+
                  mutationACTG1+
                  mutationPRKD2+
                  mutationLTB,
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
  #refLabel = "reference",
  noDigits = 2
)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/Suppl/multivariata_ALL_SIGN_classification_ref_ISS_Age_Gender_ASCT_Therapy_mutationTP53_mutationFGFR3_mutationLTB.png"
alt="Suppl3_Multivariate_analysis" />
<figcaption aria-hidden="true">Suppl3_Multivariate_analysis</figcaption>
</figure>
