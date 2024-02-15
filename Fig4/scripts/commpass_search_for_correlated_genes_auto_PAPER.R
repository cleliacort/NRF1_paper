# DESCRIPTION: This tool accepts a list of genes as input. It extrapolates a submatrix related to these genes from the COMMPASS IA17 dataset, computes the correlation among each gene, and displays the correlation through a heatmap.
# INPUT: It accepts a list of genes, identified by either gene name or Ensembl ID.
#--------------------
# Upload packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# define the command-line arguments
parser <- ArgumentParser(description = "Create a barplot from a tab-separated file.")
parser$add_argument("--reads_file", "-r", type = "character", help = "Input the path of with the normalized reads counts.")
parser$add_argument("--genes_file", "-g", type = "character", help = "Input the path of the file with the list of genes target.")
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.")
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name.")
parser$add_argument("--type_expression","-t", type = "character", help = "Do you want to select those genes that increase=UP or decrease=DOWN with the disease progression?",default=TRUE,required = TRUE)
parser$add_argument("--cluster_rows", "-R", type = "numeric", help = "Number of rows cluster for the heatmap.",default=3)
parser$add_argument("--cluster_columns", "-C", type = "numeric", help = "Number of columns cluster for the heatmap.",default=3)
parser$add_argument("--store_rows_cluster","-store_rc", type = "logical", help = "Store row_cluster names.", default=FALSE , required = FALSE)
parser$add_argument("--subset_matrix","-sub", type = "logical", help = "Do you want to subset the COMMPASS matrix based on the target genes?.", default=FALSE )
parser$add_argument("--plot_title","-title", type = "character", help = "Indicate the plot name.", default="" )

# parse the command-line arguments
args <- parser$parse_args()

############################
# PRE-PROCESSING
############################
#-----------
cat("#Set parameters generic names")
cat("\n")
reads_file <- args$reads_file
genes_file <- args$genes_file
output_dir <- args$output_dir
plot_name <- args$plot_name
type_expr <- args$type_expression
store_rows_cluster <- args$store_rows_cluster
cluster_rows <- args$cluster_rows
cluster_columns <- args$cluster_columns
subset_matrix <- args$subset_matrix
plot_title <- args$plot_title
#-----------
cat("# Reads the COMMPASS file containing normalized reads counts of all transcripts examined")
cat("\n")
norm_reads_counts <- read.delim(reads_file, sep = "\t", header = T)

#-----------
cat("# Reads the list of genes to examine")
cat("\n")
list_gene <- read.delim(genes_file, sep = "\t", header = F)

#-----------
cat("# Check if the list contains the ensembl ID or hgnc symbol")
cat("\n")
exists_in_dataframe <- any(sapply(list_gene, function(column) any(grepl("ENSG", column))))

ens_gene <- data.frame()
mart37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = "grch37.ensembl.org")
genes37 <-  biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart=mart37)

if (!exists_in_dataframe) {
  #-----------
  cat("# Adding the corresponding ensembl ID")
  cat("\n")
  colnames(list_gene) <- c("gene_name")
  ens_gene <- merge(list_gene,genes37,by.x="gene_name",by.y="hgnc_symbol")
} else {
  #-----------
  cat("# Adding the corresponding hgnc symbol")
  cat("\n")
  colnames(list_gene) <- c("ensembl_gene_id")
  ens_gene <- merge(list_gene,genes37,by="ensembl_gene_id")
  ens_gene <- ens_gene %>% rename(gene_name="hgnc_symbol")
}


#-----------
cat("# Subset the matrix with reads counts")
cat("\n")
reads_sub <- merge(norm_reads_counts,ens_gene,by.x="Gene",by.y="ensembl_gene_id") %>% distinct() 
reads_sub_mod <- reads_sub %>% dplyr::select(-Gene) %>% column_to_rownames(var="gene_name")

if(subset_matrix == TRUE | subset_matrix==T ){
  #-----------
  cat("# Storing the subsetted matrix with reads counts")
  cat("\n")
  file_path <- file.path(output_dir, paste0("subset_COMMPASS_matrix_on_target_genes_",plot_name,".txt"))
  write.table(reads_sub,file_path, sep="\t", quote=F, row.names=T,col.names = T)
}


#-----------
cat("# Compute the correlation matrix of subsetted matrix")
cat("\n")
reads_sub_mod_t <- t(reads_sub_mod)
correlation_matrix <- cor(reads_sub_mod_t)

file_path <- file.path(output_dir, paste0("correlation_matrix_COMMPASS_of_target_genes_",plot_name,".txt"))
write.table(correlation_matrix,file_path, sep="\t", quote=F, row.names=T,col.names = T)


#-----------
#
# CORRELATION HEATMAP
#
#-------------
cat("# Make the correlation heatmap")
cat("\n")
color_ramp_function <- colorRamp2(c(0,0.2,0.4,0.6,0.8, 1) ,c("#3e6f64","#94AFA9","#f6e1d5","#e9b597","#db8958","#D36C2F"))
clustering_distance_rows_set = "manhattan"
clustering_method_rows_set = "ward.D"

clustering_distance_columns_set = "manhattan"
clustering_method_columns_set = "ward.D"


h <- Heatmap(correlation_matrix,
             name = "Pearson Correlation", 
             column_names_gp = gpar(fontsize = 6),
             col = color_ramp_function,
             show_row_names=F,
             show_column_names=T,
             column_names_side = "top",
             
             clustering_distance_rows = clustering_distance_rows_set,
             clustering_method_rows = clustering_method_rows_set,
             
             clustering_distance_columns = clustering_distance_columns_set,
             # ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
             clustering_method_columns = clustering_method_columns_set,
             # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid"
             
             column_dend_side = "bottom",
             na_col = "white",
             
             column_split = cluster_columns, 
             row_split = cluster_rows, 
             
             border = TRUE,
             cluster_columns = T,
             column_title = plot_title
)


#-----------
cat("# Saving the plot")
cat("\n")
complete_name_plot_png=paste0(output_dir, "/",plot_name,"_",clustering_distance_rows_set,"_",clustering_method_rows_set,"_",
                              clustering_distance_columns_set,"_",clustering_method_columns_set,".png",sep="")
ggexport(h,res=300,filename=complete_name_plot_png,width = 4500, height = 2000)


#-----------
if(store_rows_cluster == TRUE | store_rows_cluster==T ){
  cat("# Storing the cluster rownames...")
  cat("\n")
  ht = draw(h)
  
  r.dend <- row_dend(ht)
  
  rcl.list <- row_order(ht)
  
  for (i in 1:cluster_rows) {
    cat(paste0("Cluster ",i,sep=""))
    cat("\n")
    c <- correlation_matrix[rcl.list[[i]],] %>% as.data.frame() %>% rownames_to_column() %>% select("rowname")
    file_path <- file.path(output_dir, paste0("cluster_motif_", i,"_",plot_name,".txt"))
    write.table(c,file_path, sep="\t", quote=F, row.names=F,col.names = FALSE)
  }
}