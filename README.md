# README

[![DOI](https://zenodo.org/badge/754646203.svg)](https://doi.org/10.5281/zenodo.14330213)

This repository contains the code for analyzing data and generating figures in "**Nuclear Respiratory Factor 1 (NRF1) sustains the development, progression and resistance of Multiple Myeloma**."

The scripts for reproducing the work are organized by figure. The raw data used to generate the plots will soon be available on the GEO website. Please email me if you need access to specific processed datasets. 

## **Figure 1**

[Fig1b_dataset_description](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1b_dataset_pheno_description.md): 

[Fig1c_tf_enrichment_heatmap](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1c_motif_matrix_heatmap.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2_manhattan_ward.D2.png)

```r
Rscript ./scripts/Figure1c_motif_matrix_heatmap.R -i ./data/matrix_motif_atac_tumour_mgus_0423_groupv2.txt -o ./figures -p heatmap_motifs_score_obs_exp_atac_tumour_mgus_0423_groupv2 -c 3 -store_rc TRUE
```

[Figure1f_correlation_with_Nilson.md](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Figure1f_correlation_with_Nilson.md)

## **Figure 2**

[Fig2e_NRF1_median_expression_survival](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/commpass_survival_gene_auto.R): The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/figures/survival_cluster_survival_commpass_NRF1_median_1123.png)

```bash
Rscript ./commpass_survival_gene_auto.R -i "NRF1" -r ../COMMPASS_IA17/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o ./figures -p survival_commpass_NRF1_median_1123 -surv ../COMMPASS_IA17/CoMMpass_IA17_FlatFiles_0323/MMRF_CoMMpass_IA17_STAND_ALONE_SURVIVAL_V2.tsv -c "CD138pos" -t "median"
```

[Fig2f_violin_plot_tpm_survival_per_iss](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Fig2f_violin_plot_tpm_survival_per_iss.md)

[Fig2h_emission_and_transition_heatmap_chromHMM](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Fig2g_and_Suppl2_emission_and_transition_heatmap_chromHMM_0924.md)

## **Figure 3**

[Fig3d_heatmap_atac_tumour_mgus](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/Fig3d_heatmap_tumour_mgus.md)

[Fig3f_heatmap_cell_lines_mgus_atac](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/deeptool_heatmap_atac_1023_PAPER.sh): click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/deptools_plotprofile_NRF1_consensus_TUMOUR_OF_MGUS_U266_RPMI_KMS27_MM196new_ON_10_master_list_consensus_our_chip_1123_2000_5000_5000_0_2_kmeans1_colorList7.png) to view the resulting heatmap plot.

[Fig3g_heatmap_chip_mgus_esordio](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/deeptool_heatmap_chip_mgus_esordio_0923_PAPER.sh): click [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/deptools_heatmap_plotprofile_master_list_consensus_nrf1_cleaned_from_mgsu_on_CHIP_ESORDIO_ands_MGUS_1123_2000_5000_5000_0_20_kmeans1_colorList7.png) to view the resulting heatmap plot.

## **Figure 4**

[Fig4b_genes_significantly_modulated_during_disease_progression_commpass_pv1_3](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/commpass_for_signature_auto_ensID_2_pv1_3_PAPER.R):

The R script was executed with the following parameters, and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_pv1_3.png)

```bash
Rscript Fig4/scripts/commpass_for_signature_auto_ensID_2_pv1_3.R -v data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv -g data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv -i data/10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_GENE_NAME.bed -r data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -o figures -p compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -t compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb -c "CD138pos" -s TRUE -S "UP" -sub T
```

[Fig4e_Fig4f_Fig4g_Fig4h_NRF1_gene_regualted_signature_plus_survival_commpass](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Figure4_correct_scale_heatmap_highly_correlated_TSS_and_survivals.md)

## **Suppl1**

[Suppl_1a_saturation_plot](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/Suppl_1a_saturation_plot.md)

[Suppl_1b_top_make_ggbarplot](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/make_ggbarplot.R): the R script was executed with the following parameters and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/peaks_for_SI_atac_tumour_mgus_0423.png)

```jsx
Rscript make_ggbarplot.r -i data/Suppl/number_of_peaks_per_sample_atac_tumour_mgus_2023.txt -o figures/Suppl -p peaks_for_SI_atac_tumour_mgus_0423.png -x "Filename" -y "NumLines" -pheno data/Suppl/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv -y_lab "Number of peaks" -x_lab ""
```

[Suppl_1b_bottom_counts_peaks_SI_checks](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/Suppl_1b_counts_peaks_SI_checks.md)

## **Suppl2**

[Suppl_2a_CRISP_depmap_dependency](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Suppl/Suppl_2a_depmap_dependency.md) 

[Suppl_2g_NRF1_enhancer_promoter_enrichment_evaluation](https://github.com/cleliacort/NRF1_paper/blob/main/Fig2/scripts/Suppl/Suppl_2g_counts_enhancer_and_promoter_NRF1.md) 

## **Suppl3**

[Suppl_3a_PCA_analysis_tumour_plus_mgus](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/Suppl/PCA_analysis_auto_0124.R): the R script was executed with the following parameters and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/Suppl/PCA_mgus_plus_tumour_treated_plus_tumour_onset_1224.png)

```r
Rscript PCA_analysis_auto_0124.R -i data/Suppl/matrix_multicov_tumout_and_MGUS_on_tumour_master_list_0124.txt -o figures/Suppl -p "PCA_mgus_plus_tumour_treated_plus_tumour_onset_1224.png" -pheno data/Suppl/sample_sheet_clinical_PHENOTYPE_tumour_and_mgus_0124.txt -labels FALSE 
```

[Suppl_3b_make_ggbarplot_MGUS_peaks_per_patient](https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/scripts/Suppl/make_ggbarplot.R): the R script was executed with the following parameters and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/figures/Suppl/barplot_num_peaks_MGUS_2023_sel.png)

```r
Rscript make_ggbarplot.r -i data/Suppl/number_of_peaks_per_sample_atac_tumour_mgus_2023_MGUS.txt -o figures -p barplot_num_peaks_MGUS_2023_sel -x "Filename" -y "NumLines" -pheno data/Suppl/sample_sheet_MGUS_PHENOTYPE.txt -y_lab "Number of peaks" -x_lab ""
```

## **Suppl4**

[Suppl_4a_bubbleplot_iss_genetics](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Suppl/Suppl_4a_bubbleplot_iss_genetics.Rmd)

[Suppl_4b_detect_correlated_target_gene_from_COMMPASS](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/scripts/Suppl/commpass_search_for_correlated_genes_auto_PAPER.R): the R script was executed with the following parameters and the resulting plot can be found [here](https://github.com/cleliacort/NRF1_paper/blob/main/Fig4/figures/Suppl/heatmap_correlation_compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_0124_manhattan_ward.D_manhattan_ward.D.png)

```r
Rscript commpass_search_for_correlated_genes_auto_PAPER.R -r data/Suppl/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv -g data/compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_GENE_NAME.txt -o figures -p heatmap_correlation_compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3_0124 -R 2 -C 2 -store_rc T -sub T -t UP -title "Correlation of genes increasing-UP during disease progression matching 10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_UP_05_wilcox_OR_ISS_pv1_3"
```

[Suppl3_Multivariate_analysis](https://github.com/cleliacort/NRF1_paper/blob/main/Fig3/scripts/Suppl/Suppl3_multiv