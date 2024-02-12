#!/bin/bash 
source /data/anaconda/etc/profile.d/conda.sh
conda activate ATAC_seq


# CORRELATION HEATMAP CELL LINES 
bw_files=" "
id_file=" "

outputDir="data/Figures"

#ESORDIO-CHIP-NRF1
bw_es="data/chip_seq_NRF1_ES_0823/macs2"
for filenames in `ls $bw_es | grep "bw"`; #using sort the filenmaes is alwasy R1
do
    file=$bw_es/$filenames
    id=$( echo "$filenames" | cut -d "_" -f 2 | cut -d "-" -f 1,2 )"_ES"

    bw_files=$bw_files$file" "
    id_file=$id_file$id" "
done 
echo $bw_files;
echo;

#MGUS-CHIP-NRF1
bw_mgus="data/chip_seq_NRF1_MGUS_072023/macs2"
for filenames in `ls $bw_mgus | grep "fastq_treat_pileup_S.bw" | grep -v "MM257" `; #using sort the filenmaes is alwasy R1
do
    file=$bw_mgus/$filenames
    id=$( echo "$filenames" | cut -d "_" -f 2 | cut -d "-" -f 1,2,3 )"_MGUS"

    bw_files=$bw_files$file" "
    id_file=$id_file$id" "
done 
echo $bw_files
echo;
echo $id_file



prefix="deptools_heatmap_plotprofile_master_list_consensus_nrf1_cleaned_from_mgsu_on_CHIP_ESORDIO_ands_MGUS_1123"
len=2000
after=5000
before=5000
min_z=0
max_z=20
kmeans_info=1
colorList7="white,#9c3e94"

name_output=$prefix"_"$len"_"$after"_"$before"_"$min_z"_"$max_z"_kmeans"$kmeans_info"_colorList7"
file_cluster=$outputDir"/cluster_information_plotHeatmap_"$name_output".bed"

master_list_only_tumour_conseuns_all_chip_nrf1="data/10_master_list_consensus_our_chip_1123.bed"

echo "computing the matrix neat"
computeMatrix scale-regions -R $master_list_only_tumour_conseuns_all_chip_nrf1 -S $bw_files -o $outputDir/$name_output".gz" \
--verbose -a $after -b $before --regionBodyLength $len --skipZeros --missingDataAsZero


echo "making heatMap"
plotHeatmap -m $outputDir/$name_output".gz" -out $outputDir/$name_output".png"  --sortRegions descend \
--sortUsing max --colorList $colorList7 --zMin $min_z --zMax $max_z --dpi 600 --kmeans $kmeans_info --outFileSortedRegions $file_cluster \
--samplesLabel $id_file --plotTitle $prefix