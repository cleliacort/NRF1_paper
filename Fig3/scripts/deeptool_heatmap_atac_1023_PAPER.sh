#!/bin/bash 
source /data/anaconda/etc/profile.d/conda.sh
conda activate ATAC_seq


# CHUNCK1
source /data/anaconda/etc/profile.d/conda.sh
conda activate ATAC_seq

dataDir="data"

bed_regions="data/10_master_list_consensus_our_chip_1123.bed" 


MGUS="data/atac_seq_cell_lines/merged_bigwig_MGUS_0923_mean_31_08_2023.bw"
U266="data/atac_seq_cell_lines/sorted_IRE-U266-ATAC-MM-131123_R1.fastq___IRE-U266-ATAC-MM-131123_R2.fastq_treat_pileup_S.bw"
RPMI="data/atac_seq_cell_lines/sorted_IRE-RPMI-ATAC-MM-131123_R1.fastq___IRE-RPMI-ATAC-MM-131123_R2.fastq_treat_pileup_S.bw"
KMS27="data/atac_seq_cell_lines/sorted_IRE-KMS27-ATAC-MM-131123_R1.fastq___IRE-KMS27-ATAC-MM-131123_R2.fastq_treat_pileup_S.bw"
MM196="data/atac_seq_cell_lines/mm_196_REP1.mLb.clN.bigWig"

outputDir="/figures/KMS27_MM196_MM217_KMS18" 
prefix="deptools_plotprofile_NRF1_consensus_TUMOUR_OF_MGUS_U266_RPMI_KMS27_MM196new_ON_10_master_list_consensus_our_chip_1123"

len=2000
after=5000
before=5000
min_z=0
max_z=2
kmeans_info=1
colorList7="white,#c14406"

name_output=$prefix"_"$len"_"$after"_"$before"_"$min_z"_"$max_z"_kmeans"$kmeans_info"_colorList7"
file_cluster=$outputDir"/cluster_information_plotHeatmap_"$name_output".bed"

echo "computing the matrix neat"
computeMatrix scale-regions -R $bed_regions -S $U266 $RPMI $KMS27 $MM196 $MGUS -o $outputDir/$name_output".gz" --verbose -a $after -b $before --regionBodyLength $len --skipZeros --missingDataAsZero

echo "making heatMap"
plotHeatmap -m $outputDir/$name_output".gz" -out $outputDir/$name_output".png"  --sortRegions descend \
--sortUsing max --colorList $colorList7 --zMin $min_z --zMax $max_z --dpi 600 --kmeans $kmeans_info --outFileSortedRegions $file_cluster \
--samplesLabel U266 RPMI KMS27 MM196 MGUS --plotTitle $prefix
