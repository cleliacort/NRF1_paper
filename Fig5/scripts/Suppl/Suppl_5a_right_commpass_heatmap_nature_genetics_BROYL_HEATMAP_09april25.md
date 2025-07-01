# Suppl_5a_right_commpass_heatmap_nature_BROYL_signatures

Upload the list of patients used for the paper.

``` r
fig4_c1 <- read.delim("data/cluster_motif_1_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", sep = "\t", header = F)  

fig4_c2 <- read.delim("data/cluster_motif_2_CLUSTERS_TOT_2_TSS_gene_highly_correlated_among_increase_disease_COMMPASS_only_tumour_with_NRF1_consensus_4cell_lines_0124_cluster_col_canberra_ward.D2_rows_manhattan_ward.D2_0124.txt", sep = "\t", header = F) 
  
patient_sel <- rbind(fig4_c1,fig4_c2) %>% separate(V1,c("Visit_ID"),sep="_BM_CD138pos")
```

Define the ++15=tretrasomy info.

``` r
plusplus15 <- read_table("data/CopyNumberEstimates_MMRF_CoMMpass_IA17_genome_ichorCNA.seg")

df_plusplus15 <- plusplus15 %>% filter(chr=="chr15") %>% 
  mutate(plusplus15_INFO=case_when(chr=="chr15" & copy.number == 4 ~ "plusplus15_YES", TRUE ~ "plusplus15_NO")) %>%
  group_by(SAMPLE) %>%   # Replace with the actual column representing the patient
  summarise(plusplus15_INFO = ifelse(any(plusplus15_INFO == "plusplus15_YES"), "plusplus15_YES", "plusplus15_NO"), .groups = "drop") %>% select(SAMPLE,plusplus15_INFO) %>% 
  filter(!str_detect(SAMPLE, "Plasma")) %>% 
  separate(SAMPLE,c("Visit_ID"),sep="_BM_CD138pos")

df_all_info <- merge(patient_sel,df_plusplus15,by="Visit_ID")
```

Define what the Nat.Genet. call MYC_STR. This stepp caused loww of 45
patient because the info is not present, i decided to keep them and put
no-info.

``` r
myc_str_delly <- read_table("data/StructuralEventFiles_MMRF_CoMMpass_IA17_genome_delly.tsv")

# Define MYC genomic location
MYC_CHROM <- "chr8"
MYC_START <- 128748315
MYC_END <- 128753680
MYC_CENTROMERIC_LIMIT <- MYC_START - 3000000  # 3Mb centromeric to MYC

# Define MYC STR flag
myc_str_delly_INFO <- myc_str_delly %>%
  mutate(
    MYC_STR = case_when(
      CHROM == MYC_CHROM & CHR2 != MYC_CHROM ~ "MYC_STR_YES",
      # MYC intrachromosomal deletion (DEL on chr8 with breakpoints within 3Mb centromeric of MYC)
      SVTYPE == "DEL" & CHROM == MYC_CHROM & POS < MYC_START & ENDPOSSV > MYC_CENTROMERIC_LIMIT ~ "MYC_STR_YES",
      TRUE ~ "MYC_STR_NO"
    )
  ) %>% separate(SAMPLE,c("Visit_ID"),sep="_BM_CD138pos") %>% 
  group_by(Visit_ID) %>%   # Replace with the actual column representing the patient
  summarise(MYC_STR = ifelse(any(MYC_STR == "MYC_STR_YES"), "MYC_STR_YES", "MYC_STR_NO"), .groups = "drop")

df_all_info <- merge(df_all_info,myc_str_delly_INFO,by="Visit_ID",all.x = TRUE)
```

Define the proliferative index.

``` r
mmrf_patient_info <- read_delim("data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)

PR_Q3 <- quantile(mmrf_patient_info$Prolif_Index, 0.75, na.rm = TRUE)

mmrf_patient_info_mod <- mmrf_patient_info %>% select(PUBLIC_ID,Prolif_Index) %>%  #drop_na() %>% 
    mutate(PR_index=case_when(Prolif_Index>=PR_Q3 ~ "PR_Q3_YES", is.na(Prolif_Index) ~ "PR_Q3_NO", TRUE ~ "PR_Q3_NO"))  

df_all_info <- df_all_info %>% mutate(Visit_ID_orig=Visit_ID) %>% separate(Visit_ID,c("mmrf","id","visit"),sep="_") %>% 
  mutate(patient_id=paste0(mmrf,"_",id)) %>% select(-mmrf,-id,-visit)

df_all_info <- merge(df_all_info,mmrf_patient_info_mod, by.x="patient_id", by.y = "PUBLIC_ID",all.x = TRUE)
```

Retrieve the gain/del info.

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
          del_17p13_manual_Call= case_when(SeqWGS_Cp_17p13 >= 0.3 ~ 1 , TRUE ~ 0)
          ) %>% select(Visit_ID,del_1p22_manual_Call,gain_1q21_manual_Call,del_13q14_manual_Call,del_17p13_manual_Call,SeqWGS_Cp_Hyperdiploid_Call) %>% 
  separate(Visit_ID,c("Visit_ID"),sep="_BM")

df_all_info <- df_all_info %>% separate(Visit_ID_orig,c("Visit_ID"),sep="_CD138pos")
df_all_info <- merge(df_all_info,SeqFISH, by = "Visit_ID",all.x = TRUE)
```

Retrieve traslocations info.

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
                traslocationMAFB_CALL
                ) %>% separate(SAMPLE,c("Visit_ID"),sep="_BM_CD138pos")

df_all_info <- merge(df_all_info,mut_and_trasloc, by = "Visit_ID",all.x = TRUE)
```

Upload the tpm matrix and select the patients.

``` r
tpm <- "../Fig4/data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv"
matrix_orig <- read.table(tpm, header = TRUE, row.names = 1) 
matrix_mod <- matrix_orig %>% distinct()# %>% rownames_to_column(var="ensembl_gene_id")

t_df_matrix <- as.data.frame(t(matrix_mod)) %>% rownames_to_column(var="Patient_ID") %>% 
  separate(Patient_ID,c("Visit_ID"),sep="_BM_CD138pos") %>% 
    separate(Visit_ID,c("Visit_ID"),sep="_CD138pos")

df_all_info_name <- df_all_info %>% select(Visit_ID)

df_expr <- merge(df_all_info_name,t_df_matrix,by = "Visit_ID") %>% column_to_rownames(var="Visit_ID")

t_df_expr <- as.data.frame(t(df_expr))
```

Define genes to select.

``` r
log2_tpm_matrix <- log2(t_df_expr + 1)

gcv <- function(x) {
    mean_x <- mean(x, na.rm = TRUE)
    var_x <- var(x, na.rm = TRUE)
    exp_sd <- sqrt(exp(var_x) - 1)
    exp_mean <- exp(mean_x)
    gcv_value <- exp_sd / exp_mean
    return(gcv_value)
}

# Apply to each gene (row)
gcv_values <- apply(log2_tpm_matrix, 1, gcv)

# Filter to only keep genes with GCV >= 1
keep_genes <- gcv_values >= 0.4
filtered_matrix <- log2_tpm_matrix[keep_genes, ]
```

Subset matrix

``` r
Broyl_UP <- read_delim("data/Broyl_et_2010_gene_signatures_08april25_UP.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)

Broyl_DOWN <- read_delim("data/Broyl_et_2010_gene_signatures_08april25_DOWN.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)

df_broyl_all <- rbind(Broyl_UP,Broyl_DOWN) 
colnames(df_broyl_all) <- c("gene","signature_group","expr_type")

mart37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = "grch37.ensembl.org")
genes37 <-  biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart=mart37)

df_broyl_all_adj <- merge(df_broyl_all,genes37,by.x = "gene", by.y="hgnc_symbol") %>% select(ensembl_gene_id)
filtered_matrix_mod <- filtered_matrix %>% rownames_to_column(var="ensembl_gene_id")
df_sub <- merge(filtered_matrix_mod,df_broyl_all_adj,by="ensembl_gene_id") %>% distinct() %>% column_to_rownames(var="ensembl_gene_id")
```

Scale the matrix with selected gene.

``` r
matrix_prova_scale <- as.data.frame(t(scale(t(df_sub)))) %>% rownames_to_column(var="ensembl_gene_id")
matrix_prova_scale <- merge(matrix_prova_scale,genes37,by="ensembl_gene_id") %>% select(-ensembl_gene_id) %>% column_to_rownames(var="hgnc_symbol")
matrix_prova_scale <- as.matrix(matrix_prova_scale)
```

Define the annotation.

``` r
metadata_df <- data_frame(Visit_ID=colnames(matrix_prova_scale))
pheno_order <- merge(metadata_df,df_all_info,by="Visit_ID") %>% distinct()

colAnn <- HeatmapAnnotation(
  HRD=pheno_order$SeqWGS_Cp_Hyperdiploid_Call,
  gain1q=pheno_order$gain_1q21_manual_Call,
  plusplus15_INFO=pheno_order$plusplus15_INFO,

  #TRASLOCATION
  trasl_ccnd1_t11_14 = pheno_order$traslocationCCND1_CALL,
  trasl_ccnd2 = pheno_order$traslocationCCND2_CALL,
  trasl_ccnd3 = pheno_order$traslocationCCND3_CALL,
  trasl_maf_t14_16 = pheno_order$traslocationMAF_CALL,
  trasl_mafa_t8_14 = pheno_order$traslocationMAFA_CALL,
  trasl_mafb_t14_20 = pheno_order$traslocationMAFB_CALL,
  trasl_nsd2_t4_14 = pheno_order$traslocationNSD2_CALL,
  trasl_myc = pheno_order$traslocationMYC_CALL,
  
  MYC_STR=pheno_order$MYC_STR,
  PR_index=pheno_order$PR_index,
  foo = anno_block(gp = gpar(fill = "grey80")),

  col = list(
    PR_index = c("PR_Q3_YES"="#6CA167FF","PR_Q3_NO"="#E2E2E2FF"),
    MYC_STR=  c("MYC_STR_YES"="#145A76FF", "MYC_STR_NO"="#E2E2E2FF"),
    plusplus15_INFO = c("plusplus15_YES"="#950404FF", "plusplus15_NO"="#E2E2E2FF"),
    gain1q = c("1"="#950404FF", "0"="#E2E2E2FF"),
    trasl_ccnd1_t11_14 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_ccnd2 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_ccnd3 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_nsd2_t4_14 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_maf_t14_16 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_mafa_t8_14 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_mafb_t14_20 = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    trasl_myc = c("1"="#145A76FF", "0"="#E2E2E2FF"),
    HRD = c("1"="#FED105FF", "0"="#E2E2E2FF")),  

  simple_anno_size = unit(0.3, "cm"),
  annotation_name_gp = gpar(fontsize = 6),

  show_legend = c(F,F,F,F,F,F, F,F,F,F,F,F,F,F)
)
```

``` r
metadata_df_row <- data_frame(gene_names=rownames(matrix_prova_scale))
pheno_order_row <- merge(metadata_df_row,df_broyl_all,by.x = "gene_names", by.y="gene") %>% distinct()

df_wide <- pheno_order_row %>%
  pivot_wider(
    names_from = signature_group,
    values_from = expr_type
  ) %>%
  mutate(across(-gene_names, ~ case_when(
    . == "UP" ~ 1,
    . == "DOWN" ~ -1,
    is.na(.) ~ 0,
    TRUE ~ 0  # catch-all just in case
  ))) %>% 
    arrange(factor(gene_names, levels = metadata_df_row$gene_names))

# Create the rowannotation
rowAnn <- rowAnnotation(
    MF = df_wide$MF,
    CM = df_wide$CLUSTER_MYELOID,
    CD2 = df_wide$CD2,
    MS = df_wide$MS,
    HY = df_wide$HY,
    PR = df_wide$PR,
    NFkB = df_wide$NFkB,
    CD1 = df_wide$CD1,
    PRL3 = df_wide$PRL3,
    CTA = df_wide$CTA,

  col = list(
    MF =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    CM = c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    CD2 =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    MS =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    HY =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    PR = c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    NFkB =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    CD1 =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    PRL3 =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF"),
    CTA =  c("-1"="#1BB6AFFF", "0"="#E2E2E2FF","1"="#CB87B4FF")),
  
  annotation_label = c(
    "MF","CM","CD2","MS","HY","PR","NFkB","CD1","PRL3","CTA"),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = c(TRUE,F,F,F,F,F,F,F,F,F),
  show_annotation_name = TRUE
)
```

Make the clustering and the heatmap avoiding the displaying of heatmap.

``` r
distance_col <- "canberra"
method_col <- "ward.D2"

distance_row <- "canberra"
method_row <- "ward.D2"

# CLUSTERING WORK ON THE ROW, PUT ON THE ROW WHAT YOU WANT TO CLUSTER
t_matrix_prova_scale <- t(matrix_prova_scale)
dd <- dist(t_matrix_prova_scale, method = distance_col)
hc <- hclust(dd, method = method_col)

dd_rows <- dist(matrix_prova_scale, method = distance_row)
hc_rows <- hclust(dd_rows, method = method_row)

mycols <- colorRamp2(c(-2,-1,0,1,2) ,c("#678096FF","#ACC2CFFF","#FFFFFF","#e9b597","#A12A19FF"))

h <- Heatmap(matrix_prova_scale,
        cluster_columns = hc,
        col = mycols,
        top_annotation = colAnn,
        right_annotation=rowAnn,
        column_split = 10,
        border = TRUE,
        show_row_names=T,
        show_column_names=F,
        row_names_gp = gpar(fontsize = 4)
)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig5/figures/Suppl/heatmap_genetic_on_BROYL_GENE_SIGNATURE_09APRIL25_COLUMN_BORDER_BLOCK_ANNO.png"
alt="Suppl_5a_right_commpass_heatmap_nature_BROYL_signatures" />
<figcaption
aria-hidden="true">Suppl_5a_right_commpass_heatmap_nature_BROYL_signatures</figcaption>
</figure>
