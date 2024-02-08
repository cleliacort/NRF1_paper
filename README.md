# README

This repository houses the code used for the analysis and creation of figures in "TITLE".

Considering to run all the script from a project located in the paper directory.

**Figure 1**

[Fig1b_dataset_description](scripts_github/Figure1_dataset_pheno_description.html): 

[Fig1c_tf_enrichment_heatmap](scripts_github/Figure1c_motif_matrix_heatmap.R): The R script was executed using the following parameters

```r
Rscript Fig1/scripts_github/Figure1c_motif_matrix_heatmap.R -i Fig1/data/matrix_motif_atac_tumour_mgus_0423_groupv2.txt -o Fig1/ -p heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2 -c 3 -store_rc TRUE
```

[Fig1d_fp_correlation_with_Nilson_dataset](scripts_github/Figure1d_correlation_with_Nilson.html):

[Fig1g_commpass_survival_gene_NRF1](Figure1g_commpass_survival_gene_auto.R): The R script was executed using the following parameters

```r
Rscript Fig1/scripts_github/Figure1g_commpass_survival_gene_auto.R -i "NRF1" -r COMMPASS_IA17/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o Fig1/figures/ -p survival_commpass_NRF1_median -surv COMMPASS_IA17/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv -c "CD138pos" -t "median”
```

[Fig1h_boxplot_tpm_survival](scripts_github/Figure1h_boxplot_tpm_survival.html):

**Figure 2**

[Fig2_upset_plot_in-house_ChIPseq](Fig2/scripts_github/Figure2_upset_plot_in-house_ChIPseq.html):

[Fig2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq](Fig2/scripts_github/Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq.html):

**Suppl1**

[PCA_analysis]( scripts_github/Supp_scripts_github/PCA_analysis_auto_0124.R):

```r
Rscript Fig1/scripts_github/Supp_scripts_github/PCA_analysis_auto_0124.R -i Fig1/data/Suppl/matrix_with_multicov_atac_mgus_MAXIMUM_VALUE_0423.txt -o Fig1/figures/Suppl/ -p "PCA_master_list_tumour_of_tumour_samples.png" -pheno Fig1/data/Suppl/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv
```

[correlation_with_hint_ATAC_Roma](scripts_github/Supp_scripts_github/correlation_with_hint_ATAC_Roma.html):