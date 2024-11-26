Upload the needed packages.

``` r
library(ggpubr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library("edgeR")
library(paletteer)
library(scales)
library(survival)
library(survminer)
```

## Heatmap based on gene signature using COMMPASS

Reading the COMMPASS IA17 tpm file.

``` r
prefix="TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124"

tpm <- "data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv"
matrix_orig <- read.table(tpm, header = TRUE, row.names = 1) 
matrix_mod <- matrix_orig %>% distinct() %>% rownames_to_column(var="ensembl_gene_id")
```

Upload the gene-ensmblid correspondence.

``` r
mart37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = "grch37.ensembl.org")
genes37 <-  biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart=mart37)
```

Retrieve the Ensembl ID of the gene passed as input

``` r
matrix_mod_gene <- merge(matrix_mod,genes37,by="ensembl_gene_id")
matrix_mod_gene <- matrix_mod_gene %>% select(-ensembl_gene_id)
```

Upload the gene list of interest and use it to subset the original
COMMPASS TPM dataframe.

``` r
genes <- "data/cluster_motif_1_heatmap_correlation_compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_0124.txt"
gene_list <- read.table(genes, header = F) %>% rename(gene=V1)
matrix_sub <- merge(matrix_mod_gene,gene_list,by.x = "hgnc_symbol", by.y="gene")
```

Retrieve the clinical(disease stage) and survival information.

``` r
# Upload the survival data 
survival <- read.csv("data/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv", sep = "\t", header = T)
survival <- survival %>%dplyr::select(PUBLIC_ID,oscdy,censos) 

# Retrieve the ISS information
stage <- read.delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv" , sep = "\t", header = T) %>% 
  dplyr::select(Specimen_ID,PUBLIC_ID) %>% distinct()
iss <- read.delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", sep = "\t", header = T) %>% 
  dplyr::select(PUBLIC_ID,D_PT_iss)
stage_iss <- merge(stage,iss,by="PUBLIC_ID") %>% mutate(Patient_ID=paste(Specimen_ID,"CD138pos",sep="_"))

# Merge iss with survival data to select only those having survival info
surv_stage_iss <- merge(survival,stage_iss,by="PUBLIC_ID", all.x = TRUE)
```

Retrieve the genetic information.

``` r
# Select only the mutation and translocation information needed 
mut_and_trasloc <- read.delim("data/MMRF_CoMMpass_A17_merged_locally_mutations_and_traslocation.txt", sep = "\t", header = T) %>% 
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

Merge the clinical, survival and genetic information in a unique
dataframe.

``` r
# Merge the survival-iss information with th mut-traslocation info 
surv_stage_iss_mut_trasloc <- merge(surv_stage_iss,mut_and_trasloc,by.x = "Patient_ID", by.y = "SAMPLE")

# Filter out patient that do not have ISS info
surv_stage_iss_mut_trasloc <- surv_stage_iss_mut_trasloc %>% filter(!is.na(D_PT_iss))
```

Define patient groups based on the expression of the NRF1 gene and add
this information to the dataframe created earlier.

``` r
matrix_nrf1 <- matrix_mod_gene %>% filter(hgnc_symbol=="NRF1") 
matrix_nrf1 <- as.data.frame(t(matrix_nrf1)) 
matrix_nrf1 <- matrix_nrf1 %>% filter(V1!="NRF1") %>% mutate(median_gene=median(as.numeric(V1))) 

threshold_5 <- quantile(as.numeric(matrix_nrf1$V1), probs = c(0.20,0.40,0.60,0.80,1))
matrix_nrf1 <- matrix_nrf1 %>%  mutate(status_5thr=case_when( as.numeric(V1) <= threshold_5[1] ~ "low",
                                                              as.numeric(V1) > threshold_5[1] & as.numeric(V1) <= threshold_5[2] ~ "low_medium",
                                                              as.numeric(V1) > threshold_5[2] & as.numeric(V1) <= threshold_5[3] ~ "medium",
                                                              as.numeric(V1) > threshold_5[3] & as.numeric(V1) <= threshold_5[4] ~ "high_medium",
                                                              as.numeric(V1) > threshold_5[4] ~ "high",
                                                              TRUE ~ "no")) %>% rownames_to_column(var="Patient_ID")

surv_stage_iss_mut_trasloc <- merge(surv_stage_iss_mut_trasloc,matrix_nrf1,by="Patient_ID") 
```

Subset the original matrix to include only patients with clinical,
survival, and mutational information.

``` r
# Customize the reads matrix 
matrix_sub_mod <- matrix_sub %>% column_to_rownames(var="hgnc_symbol")
matrix_sub_mod_t <- as.data.frame(t(matrix_sub_mod)) %>% rownames_to_column(var="id")

# Merge the matrix with the metadata information 
matrix_prova <- merge(matrix_sub_mod_t,surv_stage_iss_mut_trasloc, by.x = "id", by.y = "Patient_ID") %>%
  dplyr::select(-one_of(names(surv_stage_iss_mut_trasloc)[-1])) #%>% arrange(id)

# Customize the filtered matrix for the heatmap creation
matrix_prova_mod <- matrix_prova %>% column_to_rownames(var="id")
matrix_prova_scale <- t(scale(matrix_prova_mod))
```

Define the columns annotation for the heatmap.

``` r
metadata_df <- data_frame(Patient_ID=colnames(matrix_prova_scale))
pheno_order <- merge(metadata_df,surv_stage_iss_mut_trasloc,by="Patient_ID")

pheno_order$status_5thr <- factor(pheno_order$status_5thr, levels = c("high","high_medium","medium","low_medium","low"))

colAnn <- HeatmapAnnotation(
  iss_stage = pheno_order$D_PT_iss,

  #MUT
  TP53_mut = pheno_order$mutationTP53,
  BRAF_mut = pheno_order$mutationBRAF,
  FGFR3_mut = pheno_order$mutationFGFR3,
  DIS3_mut = pheno_order$mutationDIS3,
  ATM_mut = pheno_order$mutationATM,
  KRAS_mut = pheno_order$mutationKRAS,
  HIST1H1E_mut = pheno_order$mutationHIST1H1E,
  LTB_mut = pheno_order$mutationLTB,
  
  #TRASLOCATION
  trasl_myc =pheno_order$traslocationMYC_CALL,
  trasl_ccnd3=pheno_order$traslocationCCND3_CALL,
  trasl_mafa=pheno_order$traslocationMAFA_CALL,
  trasl_ccnd1=pheno_order$traslocationCCND1_CALL,
  trasl_maf=pheno_order$traslocationMAF_CALL,

  NRF1_expr_median=pheno_order$status_5thr,
  
  col = list(
    iss_stage = c("1"="#EACEB4","2" = "#E79E85", "3"="#BB5A5A"),
    #MUT
    TP53_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    BRAF_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    FGFR3_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    DIS3_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    ATM_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    KRAS_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    HIST1H1E_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    LTB_mut = c("1"="#555555FF", "0"="#E2E2E2FF"),
    #TRASLOCATION
    trasl_myc = c("1"="#555555FF", "0"="#E2E2E2FF"),
    trasl_ccnd3 = c("1"="#555555FF", "0"="#E2E2E2FF"),
    trasl_mafa = c("1"="#555555FF", "0"="#E2E2E2FF"),
    trasl_ccnd1 = c("1"="#555555FF", "0"="#E2E2E2FF"),
    trasl_maf = c("1"="#555555FF", "0"="#E2E2E2FF"),
    
    NRF1_expr_median=c("high"="#3e6f64","high_medium"="#658C83","medium"="#8BA9A2","low_medium"="#B2C5C1","low"="#D8E2E0")
    ),  

  simple_anno_size = unit(0.3, "cm"),
  annotation_label = c(
    "iss_stage",
    
    "TP53_mut",
    "BRAF_mut",
    "FGFR3_mut",
    "DIS3_mut",
    "ATM_mut",
    "KRAS_mut",
    "HIST1H1E_mut",
    "LTB_mut",
    
    "trasl_myc",
    "trasl_ccnd3",
    "trasl_mafa",
    "trasl_ccnd1",
    "trasl_maf",
    "NRF1_expr_median"),
  show_legend = c(TRUE, F,F,F,F,F,F,F,F,F,F,F,F,F,T),
  annotation_name_gp= gpar(fontsize = 7)
)
```

Create the heatmap

``` r
mycols <- colorRamp2(c(-2,-1,0,1,2) ,c("#27374D","#b3bdcc","#FFFFFF","#e9b597","#D36C2F"))

num_col_split <- 2

distance_rows <- "manhattan"
method_rows <- "ward.D2"

distance_col <- "canberra"
method_col <- "ward.D2"

h <- Heatmap(matrix_prova_scale,
             name = "Z-score Tpm expression", 
             column_names_gp = gpar(fontsize = 1),
             col = mycols,
             show_row_names=F,
             show_column_names=F,
             
             top_annotation=colAnn,
             
             clustering_distance_rows = distance_rows,
             clustering_method_rows = method_rows,
             
             clustering_distance_columns = distance_col,
             clustering_method_columns = method_col,
             
             column_dend_side = "top",
             na_col = "white",
             
             column_split = num_col_split,

             border = TRUE,
             cluster_columns = T,
             cluster_rows=T
)
#print(h)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/heatmap_cluster_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_NRF1_expr_thr5_0124.png"
alt="Fig4d_heatmap_COMMPASS_highly_enriched_NRF1-regualted_signature" />
<figcaption
aria-hidden="true">Fig4d_heatmap_COMMPASS_highly_enriched_NRF1-regualted_signature</figcaption>
</figure>

Store the cluster information to which each patient belongs.

``` r
output_dir <- "figures/"
ht = draw(h)
```

![](Figure4_correct_scale_heatmap_highly_correlated_TSS_and_survivals_files/figure-markdown_github/column_cluster_save-1.png)

``` r
col_cluster <- column_order(ht)

for (i in 1:num_col_split) {
  cat(paste0("Cluster ",i,sep=""))
  cat("\n")
  c <- matrix_prova_scale[,col_cluster[[i]]]
  c_name <- data_frame(public_id=colnames(c)) 
  file_path <- file.path(output_dir, paste0("cluster_motif_", i,"_CLUSTERS_TOT_",num_col_split,"_",prefix,"_cluster_col_",distance_col,"_",method_col,"_rows_",distance_rows,"_",method_rows,"_0124.txt"))
  #write.table(c_name,file_path, sep="\t", quote=F, row.names=F,col.names = FALSE)
}
```

    ## Cluster 1
    ## Cluster 2

## Visualizing the number of patients by ISS for each cluster

Count the number of patients belonging to each disease stage for each
cluster of the heatmap.

``` r
df_plot_iss <- data.frame()

for (i in 1:num_col_split) {
  cat(paste0("Cluster ",i,sep=""))
  cat("\n")
  c <- matrix_prova_scale[,col_cluster[[i]]]
  c_name <- data_frame(public_id=colnames(c)) %>% mutate(type=paste0("cluster",i))
  df_plot_iss <- rbind(df_plot_iss,c_name)
}
```

    ## Cluster 1
    ## Cluster 2

Count the number of patients at each disease stage for each cluster in
the heatmap.

``` r
# Subset the suv info only to the ISS
surv_stage_iss_mut_trasloc_sub <- surv_stage_iss_mut_trasloc %>%  select(Patient_ID,D_PT_iss) %>% rename(public_id=Patient_ID)

# Merge the cluster ID info with the ISS
df_plot_iss_info <- left_join(df_plot_iss,surv_stage_iss_mut_trasloc_sub)

# Retrieve the total number of sample for each cluster
df_tot <- df_plot_iss %>% group_by(type) %>% count()

# Count the number of ISS for each cluster 
df_plot_count <- df_plot_iss_info %>% group_by(type,D_PT_iss) %>% count() %>% rename(num=n)

# Merge the total info and the ISS and compute the perc
df_plot_count <- left_join(df_plot_count,df_tot) %>% rename(tot=n) %>% 
  mutate(perc=round((num*100)/tot,2))
```

Create donut plots.

``` r
# cluster_sel="cluster1"
cluster_sel="cluster2"

df_plot_count$D_PT_iss <- factor(df_plot_count$D_PT_iss, levels = c("3","2","1"))
df_plot_count_C <- df_plot_count %>% filter(type==cluster_sel)
palette <- c("#BB5A5A","#E79E85","#EACEB4")

# Compute percentages
df <- df_plot_count_C %>%
  group_by(type) %>%
  mutate(perc = perc / sum(perc) * 100)

# Add cumulative percentage for ordering the segments in the plot
df <- df %>%
  arrange(D_PT_iss) %>%
  mutate(ymax = cumsum(perc),
         ymin = lag(ymax, default = 0))

# Calculate label position
df <- df %>%
  mutate(labelPosition = (ymax + ymin) / 2,
         label = paste0("ISS-", D_PT_iss, "\n", round(perc, 1), "%"))

# Create the donut chart
donut_chart <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = D_PT_iss)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  #geom_text(aes(x = 3.5, y = labelPosition, label = label), size = 5, color = "white") +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_manual(values =  c("#BB5A5A","#E79E85","#EACEB4")) +
  labs(title = cluster_sel)

# Display the plot
# print(donut_chart)
```

![Fig4_ISS_per_cluster1](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/circularPlot_ISS_cluster1_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_1023.png)
![Fig4_ISS_per_cluster2](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/circularPlot_ISS_cluster2_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_1023.png)

## Survival analysis based on patient subdivision from the heatmap

Obtain and customize the survival information for the patient in the
clusters of the heatmap.

``` r
#Retrieve the needed survival info 
surv_info <- surv_stage_iss_mut_trasloc %>% select(Patient_ID,oscdy,censos) %>% rename(public_id=Patient_ID)

# Merge the information with the cluster info per patient
df_surv <- merge(df_plot_iss,surv_info, by="public_id")

# Convert days into months
df_surv <- df_surv %>% mutate(oscdy=round(oscdy/30,0))

# Filter out patients based on months
df_surv <- df_surv %>% filter(oscdy<=65)
```

Fit the survival model

``` r
fit <- survfit(Surv(oscdy, censos) ~ type, data = df_surv)
```

Make the survival plot.

``` r
blank_theme <- theme_bw()+ theme(panel.grid = element_blank())

gg_surv <- ggsurvplot(fit,
                 risk.table = TRUE, 
                 xlim = c(0, 59),
                 break.x.by = 6,  break.y.by = .1, axes.offset = FALSE,
                 pval = T, 
                 surv.median.line = "hv",
                 conf.int = F,
                 size = 1.5,
                 palette=c("#d2d2d2","#698d85"),
                 xlab = "Months", ylab = "Probability of Survival",
                 legend = c(0.85,0.9),
                 legend.title = "",
                 font.legend = 14,
                 font.xlab = c(14, "bold"),
                 font.ylab = c(14, "bold"),
                 risk.table.title = "N. patients at risk", risk.table.y.text = F, ggtheme = theme_pubr())

#show(gg_surv)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/survival_cluster_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_1023.png"
alt="Fig4f_survival_by_NRF1_gene_regulated_signature" />
<figcaption
aria-hidden="true">Fig4f_survival_by_NRF1_gene_regulated_signature</figcaption>
</figure>

## Survival analysis at each ISS based on patient subdivision from the heatmap

Conduct a survival analysis and create a plot based on the patient
subdivisions derived from the heatmap for each disease stage.

``` r
# List to store survival plots
splots <- list()
# Customize surv information
surv_stage_iss_mut_trasloc_cust <- surv_stage_iss_mut_trasloc %>% select(-oscdy,-censos)
df_surv_iss <- merge(df_surv,surv_stage_iss_mut_trasloc_cust,by.x = "public_id", by.y="Patient_ID")

# Loop through three different values of iss
for (iss in 1:3) {
  # Work with survival data for the current iss
  df_surv_single_iss <- df_surv_iss %>% filter(D_PT_iss == iss) %>% select(public_id, type, oscdy, censos, D_PT_iss) 
  
  # Fit the survival
  fit <- survfit(Surv(oscdy, censos) ~ type, data = df_surv_single_iss)
  
  # Make the survival plot
  sp <- ggsurvplot(fit,
                   risk.table = TRUE, 
                   xlim = c(0, 59),
                   break.x.by = 6,  break.y.by = .1, axes.offset = FALSE,
                   pval = T, 
                   surv.median.line = "hv",
                   conf.int = F,
                   size = 1.5,
                   palette=c("#d2d2d2","#698d85"),
                   xlab = "Months", ylab = "Probability of Survival",
                   #legend = c(0.7,0.9),
                   legend.title = "",
                   legend = "none",  # Disable the legend
                   font.legend = 14,
                   font.xlab = c(14, "bold"),
                   font.ylab = c(14, "bold"),
                   risk.table.title = "N. patients at risk", risk.table.y.text = F, ggtheme = theme_pubr()) +
    ggtitle(paste0("Survival with population at ISS ", iss))
 
   ggexport(sp,res=300,filename=paste0(output_dir,"/survival_COMMPASS_signature_high_corr_divided_per_ISS",iss,"_0124.png"),width = 3000, height = 2000)

  # Add the survival plot to the list
  splots[[iss]] <- sp
}
re <- arrange_ggsurvplots(splots, print = TRUE,ncol = 3, nrow = 1, risk.table.height = 0.4)
```

![](Figure4_correct_scale_heatmap_highly_correlated_TSS_and_survivals_files/figure-markdown_github/multi_surv_plot-1.png)

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/survival_COMMPASS_signature_high_corr_divided_per_ISS_0124.png"
alt="Fig4_survival_by_iss" />
<figcaption aria-hidden="true">Fig4_survival_by_iss</figcaption>
</figure>
