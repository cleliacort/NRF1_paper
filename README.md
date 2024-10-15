# README

This repository houses the code used for the analysis and creation of figures in "TITLE".

https://github.com/cleliacort/NRF1_paper/blob/main/

Considering to run all the script from a project located in the paper directory.

**Figure 1. Epigenetic profiling of NRF1 reveals NRF1 as a binder of highly penetrant loci**

[Fig1b_dataset_description](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1b_dataset_pheno_description.md): 

[Fig1c_tf_enrichment_heatmap](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1c_motif_matrix_heatmap.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2_manhattan_ward.D2.png)

```r
Rscript ./scripts/Figure1c_motif_matrix_heatmap.R -i ./data/matrix_motif_atac_tumour_mgus_0423_groupv2.txt -o ./figures -p heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2 -c 3 -store_rc TRUE
```

[Figure1f_correlation_with_Nilson.md](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1f_correlation_with_Nilson.md)

**Figure 2**

[Fig2d_NRF1_median_expression_survival](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Fig2e_boxplot_tpm_survival.md)

[Fig2e_boxplots_showing_NRF1_enrichment_vs_ISS](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/commpass_survival_gene_auto.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/TPM_per_ISS_in_population_divided_per_NRF1_expression_low_high_201123.png)

```bash
Rscript ./commpass_survival_gene_auto.R -i "NRF1" -r ../COMMPASS_IA17/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o ./figures -p survival_commpass_NRF1_median_1123 -surv ../COMMPASS_IA17/CoMMpass_IA17_FlatFiles_0323/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv -c "CD138pos" -t "median"
```

[Fig2g_emission_and_transition_heatmap_chromHMM](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Fig2g_and_Suppl2_emission_and_transition_heatmap_chromHMM_0924.md)

**Figure 3**

[Fig3_heatmap_tumour_mgus](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/Fig3_heatmap_tumour_mgus.md):

[Fig3_heatmap_cell_lines_mgus_atac](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/deeptool_heatmap_atac_1023_PAPER.sh):

Click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/deptools_plotprofile_NRF1_consensus_TUMOUR_OF_MGUS_U266_RPMI_KMS27_ON_10_master_list_consensus_our_chip_1123_2000_5000_5000_0_2_kmeans1_colorList7.png) to view the resulting heatmap plot.

[Fig3_heatmap_chip_mgus_esordio](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/deeptool_heatmap_chip_mgus_esordio_0923_PAPER.sh):

Click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/deptools_heatmap_plotprofile_master_list_consensus_nrf1_cleaned_from_mgsu_on_CHIP_ESORDIO_ands_MGUS_1123_2000_5000_5000_0_20_kmeans1_colorList7.png) to view the resulting heatmap plot.

**Figure 4**

[Fig4a_barplot_COMMPASS_up_and_down](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Figure4a_barplot_COMMPASS_up_and_down_PAPER.md):

[Fig4b_commpass_for_signature_pv1_3](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/commpass_for_signature_auto_ensID_2_pv1_3_PAPER.R):

The R script was run using specific parameters. This generated two plots. The first plot, which can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_pv1_3.png), illustrates the enrichment of genes that have notably increased their expression during the progression of the disease. 

```bash
Rscript Fig4/scripts/commpass_for_signature_auto_ensID_2_pv1_3.R -v data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv -g data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv -i data/10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_GENE_NAME.bed -r data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o figures -p compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -t compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -c "CD138pos" -s TRUE -S "UP" -sub T
```

The second plot, accessible [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_DOWN_pv1_3.png), displays the enrichment of genes whose expression significantly decreased during the same disease progression.

```bash
Rscript Fig4/scripts/commpass_for_signature_auto_ensID_2_pv1_3.R -v data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv -g data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv -i data/10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_GENE_NAME.bed -r data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o figures -p compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -t compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -c "CD138pos" -s TRUE -S "DOWN" -sub T
```

[Fig4c_commpass_search_for_correlated_genes](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Figure4d_correct_scale_heatmap_highly_correlated_TSS_PAPER.Rmd):

The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/heatmap_correlation_compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_0124_manhattan_ward.D_manhattan_ward.D.png)

```bash
Rscript Fig4/scripts/commpass_search_for_correlated_genes_auto_PAPER.R -r data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -g data/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_GENE_NAME.txt -o figures -p heatmap_correlation_compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_0124 -R 2 -C 2 -store_rc T -sub T -t UP -title "Correlation of genes increasing-UP during disease progression matching 10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3"
```

[Figure4def](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Figure4_correct_scale_heatmap_highly_correlated_TSS_PAPER.md):

**Suppl1**

[Suppl_1_saturation_plot](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/Suppl_1a_saturation_plot.md)

[Suppl_1_counts_peaks_SI_checks](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/Suppl_1b_counts_peaks_SI_checks.md)

[Suppl_1_make_ggbarplot](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/Suppl_1b_make_ggbarplot.R): the R script was executed with the following parameters and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/peaks_for_SI_atac_tumour_mgus_0423.png)

```jsx
Rscript make_ggbarplot.r -i data/Suppl/number_of_peaks_per_sample_atac_tumour_mgus_2023.txt -o figures/Suppl -p peaks_for_SI_atac_tumour_mgus_0423.png -x "Filename" -y "NumLines" -pheno data/Suppl/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv -y_lab "Number of peaks" -x_lab ""
```

[PCA_analysis](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/PCA_analysis_auto_0124.R):

```r
Rscript Fig1/scripts/Suppl/PCA_analysis_auto_0124.R -i Fig1/data/Suppl/matrix_with_multicov_atac_mgus_MAXIMUM_VALUE_0423.txt -o Fig1/figures/Suppl/ -p "PCA_master_list_tumour_of_tumour_samples.png" -pheno Fig1/data/Suppl/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv
```

Click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/PCA_master_list_tumour_of_tumour_samples.png) to view the resulting PCA plot.

[correlation_with_hint_ATAC_Roma](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/correlation_with_hint_ATAC_Roma.md):

**Suppl4**

[Multivariate_analysis](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Suppl/multivariate_analysis_PAPER.md):