# DESCRIPTION: compute a box plot showing the enrichment of a passed set of genes in COMMPASS data set;
# NOTA0: The script is compatible with the A17 version of COMMPASS.
# NOTA1: If you wish to use the script with a different version of COMMPASS, you need to adjust the name check, specifically the cell name (e.g., CD138). By default, the script considers the Specimen_ID, which appears like "MM138_1_BM". However, if this doesn't match the name in the reads count file, you will need to modify it. For instance, in the AI17 case, we had to add the "CD138" suffix to the Specimen_ID to ensure a match between the reads and metadata files.
# NOTA2: You can use the argument "significant TRUE" to only consider genes that are significantly modulated.
# NOTA3: The "selection" argument allows you to retain significant genes that either increase (UP) or decrease (DOWN) during disease progression.
#--------------------
# Upload packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(HGNChelper))
library(biomaRt)
library(dbplyr)
library(BiocManager)
# define the command-line arguments
parser <- ArgumentParser(description = "Create a barplot from a tab-separated file.")
parser$add_argument("--visit_file", "-v", type = "character", help = "Input the path of the file with visit number info.")
parser$add_argument("--stage_file", "-g", type = "character", help = "Input the path of the file with stage info.")
parser$add_argument("--genes_file", "-i", type = "character", help = "Input the path of the file with the list of genes target.")
parser$add_argument("--reads_file", "-r", type = "character", help = "Input the path of with the normalized reads counts.")
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.")
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name.")
#parser$add_argument("--colour_sel", "-c", type = "character", help = "Color of highlighted cell lines.",default="#659CC6")
parser$add_argument("--x_label","-x_lab", type = "character", help = "Label for x-axis.",default="ISS")
parser$add_argument("--y_label","-y_lab", type = "character", help = "Label for y-axis.",default="Normalized reads counts(tpm)")
parser$add_argument("--title","-t", type = "character", help = "Title of the figure.",default="")
parser$add_argument("--cell","-c", type = "character", help = "Cell type detected.",default="")
parser$add_argument("--significant","-s", type = "logical", help = "Do you want to consider only the significant genes among those passed?",default=TRUE,required = TRUE)
parser$add_argument("--selection","-S", type = "character", help = "Do you want to consider genes that increase(UP) or decrease(DOWN) with the disease progression?",default="ALL",required = TRUE)
parser$add_argument("--subset_matrix","-sub", type = "logical", help = "Do you want to subset the COMMPASS matrix based on the target genes?.", default=FALSE )


# parse the command-line arguments
args <- parser$parse_args()

# # check if the input file, output directory, plot file, colors, x column, and y column were specified
# if (is.null(args$input_file) || is.null(args$plot_name) || is.null(args$output_dir)) {
#   stop("Usage: Rscript script.R --input_file path/input_file.txt --prefix_name file_name --output_dir output_dir --colour_bar color_name --plot_name plot_name.png --x_lab label_x-axis --y_lab label_y-axis")
# }

############################
# PRE-PROCESSING
############################
#-----------
cat("#Set parameters generic names")
cat("\n")
want_significant <- args$significant
which_selection <- args$selection
cell <- args$cell
# cell="CD138pos"
label_x_cust <- args$x_label
labe_y_cust <- args$y_label
visit_file <- args$visit_file
stage_file <- args$stage_file
gene_file <- args$genes_file
reads_file <- args$reads_file
subset_matrix <- args$subset_matrix
output_dir <- args$output_dir
plot_name <- args$plot_name
title_plot <- args$title
# # -------------------
# MANUAL EXAMPLE
cell="CD138pos"
want_significant=TRUE
which_selection="UP"
label_x_cust <- args$x_label
labe_y_cust <- args$y_label
visit_file <- "data/MMRF_CoMMpass_IA17_PER_PATIENT_VISIT_V2.tsv"
stage_file <- "data/MMRF_CoMMpass_IA17_PER_PATIENT_V2.tsv"
gene_file <- "data/10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_GENE_NAME.bed"
reads_file <- "data/MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv"
subset_matrix <- args$subset_matrix
output_dir <- args$output_dir
plot_name <- "compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_100624"
title_plot <-"compass_progression_10_master_list_consensus_our_chip_1123_ANNOTATED_selected_TSS_minus_plus_2kb_100624" 
labe_y_cust <- "Normalized reads counts(tpm)"
label_x_cust <- "ISS"
output_dir <- "figures"
#-----------
cat("# Reads COMMPASS file with name with visit info")
cat("\n")
stage <- read.delim(visit_file , sep = "\t", header = T) %>% 
  dplyr::select(Specimen_ID,PUBLIC_ID) %>% distinct()

#-----------
cat("# Reads COMMPASS file with ISS info")
cat("\n")
iss <- read.delim(stage_file, sep = "\t", header = T) %>% 
  dplyr::select(PUBLIC_ID,D_PT_iss)

#-----------
cat("# Merge the name and ISS info")
cat("\n")
stage_iss <- merge(stage,iss,by="PUBLIC_ID")

#-----------
cat("# Reads the list of selected genes to examine")
cat("\n")
list_gene <- read.delim(gene_file, sep = "\t", header = F) 
# Change the colnames 
colnames(list_gene) <- c("gene_name")
list_gene <- list_gene %>% filter(!is.na(gene_name)) %>% distinct()

#-----------
cat("# Upload the COMMPASS file containing normalized reads counts of all transcripts examined")
cat("\n")
norm_reads_counts <- read.delim(reads_file, sep = "\t", header = T)

#-----------
cat("Retrieve biomart...")
cat("\n")
mart37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = "grch37.ensembl.org")
genes37 <-  biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart=mart37)
# Retrieve the Ensembl ID of the gene passed as input
ens_gene <- merge(list_gene,genes37,by.x="gene_name",by.y="hgnc_symbol") 

#-----------
cat("# Subset the normalized reads counts to those selected genes")
cat("\n")
reads_sub <- merge(norm_reads_counts,ens_gene,by.x="Gene",by.y="ensembl_gene_id")  %>% distinct()

if(subset_matrix == TRUE | subset_matrix==T ){
  #-----------
  cat("# Storing the subsetted matrix with reads counts")
  cat("\n")
  file_path <- file.path(output_dir, paste0("subset_COMMPASS_matrix_on_target_genes_",plot_name,"_pv1_3.txt"))
  write.table(reads_sub,file_path, sep="\t", quote=F, row.names=T,col.names = T)
}

#-----------
cat("#Customize the subset normalized reads count matrix and transpose it")
cat("\n")
reads_sub_mod <- reads_sub %>% mutate(ens_hgnc=paste(gene_name,Gene,sep="-")) %>%  dplyr::select(-gene_name,-Gene) %>% column_to_rownames(var="ens_hgnc")
t_reads_sub_mod <- as.data.frame(t(reads_sub_mod)) %>% rownames_to_column(var = "Specimen_ID") 

# #------------- MANUAL CHECK FOR OUTLIERS
# #block the scientific annotation
# options(scipen = 999)
# reads_sub_mod_custom <- reads_sub_mod
# #compute the median of expression for each gene
# reads_sub_mod_custom$median <-  apply(reads_sub_mod_custom, 1, median)
# #reads_sub_mod$mean <-  apply(reads_sub_mod, 1, mean)
# reads_sub_mod_custom <- reads_sub_mod_custom %>% select(median)
# #count the 3rd quantile
# third_quartile <- quantile(reads_sub_mod_custom$median, 0.75)
# #filter the genes with median greater than 3rd quantile
# reads_sub_mod_custom1 <- reads_sub_mod_custom %>% filter(median>third_quartile) %>% rownames_to_column(var="genes")
# 
# df_prova <- data.frame(genes=df_info)
# prova <- merge(reads_sub_mod_custom1,df_prova,by="genes")
# #-------------

#-----------
if (!is.null(cell)) {  
  stage_mod <- stage_iss %>% mutate(Specimen_ID=paste(Specimen_ID,cell,sep="_"))
       cat("# Customize the name of patient sample") 
       cat("\n")###############
}

#-----------
cat("# Select the patient ID for each ISS stage")
cat("\n")
t1 <- stage_mod %>% dplyr::filter(D_PT_iss==1) %>% dplyr::select(Specimen_ID) %>% distinct()
t2 <- stage_mod %>% filter(D_PT_iss==2) %>% dplyr::select(Specimen_ID) %>% distinct()
t3 <- stage_mod %>% filter(D_PT_iss==3) %>% dplyr::select(Specimen_ID) %>% distinct()
# initialize a datafram in r
df_long_results <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_long_results) <- c("Specimen_ID","type","Genes","Valore")

#-----------
cat("# Making the dataframe for the plot...")
cat("\n")
if(want_significant == FALSE){ 
  cat("# NO significant genes were selected..")
  
  #–------------------
  #this will take all the genes passed as input WITHOUT CHECKING FOR THE SIGNIFICANT MODULATION ACCROSS DISEASE PROGRESSIO
  #-------------------
  # Select the reads for each ISS stage 
  reads_t1 <- merge(t_reads_sub_mod,t1,by="Specimen_ID") %>% mutate(type="t1")
  reads_t2 <- merge(t_reads_sub_mod,t2,by="Specimen_ID") %>% mutate(type="t2")
  reads_t3 <- merge(t_reads_sub_mod,t3,by="Specimen_ID") %>% mutate(type="t3")
  
  # Put together the reads count with ISS information
  df_results <- rbind(reads_t1,reads_t2,reads_t3)
  
  # Transform the data frame in long form 
  df_long_results <- gather(df_results, key = "Genes", value = "Valore", c(-Specimen_ID,-type))
  
  # #----------- MANUAL CUSTOMIZATION
  # vector_rem <- c("RPL10", "TXNDC5")
  # df_long_results <- df_long_results %>% filter(!(Genes %in% vector_rem))
  # #----------------
} else {
  #–------------------
  # SELECT ONLY THE SIGNIFICANT GENES FROM THE INPUTED LIST
  #-------------------
  cat("# Selecting only the significant genes...")
  cat("\n")
  p_thr=0.05
  list_TRANSCRIPTS=c()
  list_INFO=c()

  for (i in 2:length(colnames(t_reads_sub_mod))) {
    ID_transcript=colnames(t_reads_sub_mod[i])
    df <- t_reads_sub_mod %>% dplyr::select(Specimen_ID,ID_transcript)
    
    cpm_st1 <- merge(df,t1,by="Specimen_ID")
    cpm_st2 <- merge(df,t2,by="Specimen_ID")
    cpm_st3 <- merge(df,t3,by="Specimen_ID")
    
    pv1_2 <- wilcox.test(cpm_st1[,2] ,cpm_st2[,2])$p.value
    pv1_3 <- wilcox.test(cpm_st1[,2] ,cpm_st3[,2])$p.value
    pv2_3 <- wilcox.test(cpm_st2[,2] ,cpm_st3[,2])$p.value
    info=paste(ID_transcript,pv1_2,pv1_3,pv2_3,median(cpm_st1[,2]),median(cpm_st2[,2]),median(cpm_st3[,2]),"NAN",sep = "_")
    
    if(!is.na(pv1_2) & !is.na(pv1_3) & !is.na(pv2_3)){
      #info=paste(ID_transcript,pv1_2,pv1_3,pv2_3,median(cpm_st1[,2]),median(cpm_st2[,2]),median(cpm_st3[,2]),"NAN",sep = "_")
      
      if( pv1_3<=p_thr ){
        list_TRANSCRIPTS <- append(list_TRANSCRIPTS,ID_transcript)
        
        
        if (median(cpm_st1[,2])!=0 & median(cpm_st2[,2])!=0  & median(cpm_st3[,2])!=0 ){
          if (median(cpm_st1[,2])>median(cpm_st3[,2])){
            info=paste(ID_transcript,pv1_2,pv1_3,pv2_3,median(cpm_st1[,2]),median(cpm_st2[,2]),median(cpm_st3[,2]),"DOWN",sep = "_")
          }
          else{
            info=paste(ID_transcript,pv1_2,pv1_3,pv2_3,median(cpm_st1[,2]),median(cpm_st2[,2]),median(cpm_st3[,2]),"UP",sep = "_")
          }
      }
        # list_INFO <- append(list_INFO,info)
        
        print(ID_transcript)
    }
    list_INFO <- append(list_INFO,info)
      
    }
  }
  
  # Make a dataframe with the list of information for each transcript
  df_info <- data.frame(d=list_INFO)
  df_info <- df_info %>% separate(d,c("GENE_NAME","pv1_2","pv1_3","pv2_3","median_iss1","median_iss2","median_iss3","EXPRESSION"),sep="_")
  
  # Save significant genes with wilcox information
  #name_file_pv <- paste0(output_dir, "/", plot_name,"_05_wilcox_OR_ISS_pv1_3.txt")
  #write.table(df_info,name_file_pv, sep="\t", quote=F, col.names=T, row.names=F)

  # Select the transcript for the behavior in the disease progression
  if(which_selection == "DOWN"){
    df_info <- df_info %>% filter(EXPRESSION=="DOWN") %>% dplyr::select(GENE_NAME)
    df_info <- c(df_info$GENE_NAME)
  } else if(which_selection == "UP") {
    df_info <- df_info %>% filter(EXPRESSION=="UP") %>% dplyr::select(GENE_NAME)
    df_info <- c(df_info$GENE_NAME)
  }else{
    df_info <- c(df_info$GENE_NAME)
  }
  num_diff_genes <- length(df_info)
  
    #Save the selected gene based on the behavior during disease progression
  #name_file_pv_filt <- paste0(output_dir, "/", plot_name,"_",which_selection,"_05_wilcox_OR_ISS_pv1_3.txt")
  #write.table(df_info,name_file_pv_filt, sep="\t", quote=F, col.names=F, row.names=F)

  # Subset the dataframe with the reads count of each selected transcript of all sample
  significant_stage_transcripts <- t_reads_sub_mod %>% dplyr::select(1,all_of(df_info))
  
  # Select the reads for each ISS stage 
  reads_t1 <- merge(significant_stage_transcripts,t1,by="Specimen_ID") %>% mutate(type="t1")
  reads_t2 <- merge(significant_stage_transcripts,t2,by="Specimen_ID") %>% mutate(type="t2")
  reads_t3 <- merge(significant_stage_transcripts,t3,by="Specimen_ID") %>% mutate(type="t3")
  
  # Put together the reads count with ISS information
  df_results <- rbind(reads_t1,reads_t2,reads_t3)
  
  # Transform the data frame in long form 
  df_long_results <- gather(df_results, key = "Genes", value = "Valore", c(-Specimen_ID,-type))
  head(df_long_results)
  # #----------- MANUAL CUSTOMIZATION
  # vector_rem <- c("PTMA", "TXNDC5")
  # df_long_results <- df_long_results %>% filter(!(Genes %in% vector_rem))
  # #----------------
  
}


#-----------
cat("# Making the plot")
cat("\n")
my_comparisons <- list( c("t1", "t2"),c("t1","t3"),c("t2","t3"))

p <- ggboxplot(df_long_results,x = "type", y = "Valore", fill = "type", color = "#696969",
               palette = c("#EACEB4","#E79E85","#BB5A5A"),width = 0.7,outlier.shape = NA,
               lwd=0.6,
               xlab=label_x_cust,ylab=labe_y_cust ,
               bxp.errorbar =TRUE,
               title=paste0(title_plot,"(n=",num_diff_genes,")")
               )+
  #stat_compare_means(comparisons = my_comparisons,label = "p.format")+
  border()#+
  #geom_jitter(position = position_jitter(width = 0.2), size = 1, color = "#4D4D4D", alpha = 0.5)

show(p)

#-----------
cat("# Saving the plot")
cat("\n")
ggexport(p,res=300,filename=paste0(output_dir, "/", plot_name,"_",which_selection,"_STATISTICS_pv1_3.png"),width = 2000, height = 2500)
# Save with border
max_y=60
p <- ggpar(p,ylim=c(0,max_y))+scale_y_continuous(breaks = seq(0, max_y, by = 10)) +
  theme(plot.title = element_text(size = 5)) +border()

# Save the plot to a file
ggexport(p,res=300,filename=paste0(output_dir, "/", plot_name,"_",which_selection,"_pv1_3.png"),width = 1500, height = 2500)






