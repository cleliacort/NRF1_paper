# README

This repository houses the code used for the analysis and creation of figures in "TITLE".

https://github.com/cleliacort/NRF1_paper/blob/main/

Considering to run all the script from a project located in the paper directory.

**Figure 1**

[Fig1b_dataset_description](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1_dataset_pheno_description.md): 

[Fig1c_tf_enrichment_heatmap](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1c_motif_matrix_heatmap.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2_manhattan_ward.D2.png)

```r
Rscript Fig1/scripts/Figure1c_motif_matrix_heatmap.R -i Fig1/data/matrix_motif_atac_tumour_mgus_0423_groupv2.txt -o Fig1/ -p heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2 -c 3 -store_rc TRUE
```

[Fig1d_fp_correlation_with_Nilson_dataset](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1d_correlation_with_Nilson.md):

[Fig1g_commpass_survival_gene_NRF1](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1g_commpass_survival_gene_auto.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/survival_cluster_survival_commpass_NRF1_median_1123.png)

```r
Rscript Fig1/scripts/Figure1g_commpass_survival_gene_auto.R -i "NRF1" -r COMMPASS_IA17/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o Fig1/figures/ -p survival_commpass_NRF1_median -surv COMMPASS_IA17/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv -c "CD138pos" -t "median”
```

[Fig1h_boxplot_tpm_survival](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1h_boxplot_tpm_survival.md):

**Figure 2**

[Fig2_upset_plot_in-house_ChIPseq](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Figure2_upset_plot_in-house_ChIPseq.md):

[Fig2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Figure2_heatmap_jaccard_corr_binary_info_all_NRF1_ChIP-seq.md):

**Figure 3**

[Fig3_heatmap_tumour_mgus](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/Fig3_heatmap_tumour_mgus.md):

**Suppl1**

[PCA_analysis](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/PCA_analysis_auto_0124.R):

```r
Rscript Fig1/scripts/Suppl/PCA_analysis_auto_0124.R -i Fig1/data/Suppl/matrix_with_multicov_atac_mgus_MAXIMUM_VALUE_0423.txt -o Fig1/figures/Suppl/ -p "PCA_master_list_tumour_of_tumour_samples.png" -pheno Fig1/data/Suppl/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv
```

Click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/PCA_master_list_tumour_of_tumour_samples.png) to view the resulting PCA plot.

[correlation_with_hint_ATAC_Roma](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/correlation_with_hint_ATAC_Roma.md):