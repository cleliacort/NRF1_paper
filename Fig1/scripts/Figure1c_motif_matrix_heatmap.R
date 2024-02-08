suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggrepel))

########################
#
#     UNSUPERVISED
#
########################
# define the command-line arguments
parser <- ArgumentParser(description = "Create a heatmap from a tab-separated file.")
parser$add_argument("--input_file", "-i", type = "character", help = "Input the path for the file matrix to use.")
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.")
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name.")
parser$add_argument("--cluster_rows", "-c", type = "numeric", help = "Number of rows cluster for the heatmap.",default=3)
parser$add_argument("--store_rows_cluster","-store_rc", type = "logical", help = "Store row_cluster names.", default=FALSE , required = FALSE)

# parse the command-line arguments
args <- parser$parse_args()

# check if the input file, output directory, plot file, colors, x column, and y column were specified
if (is.null(args$input_file) || is.null(args$plot_name)) {
  stop("Usage: Rscript script.R --input_file path/input_file.txt --prefix_name file_name --output_dir output_dir --colour_bar color_name --plot_name plot_name.png --x_lab label_x-axis --y_lab label_y-axis")
}

set.seed(123)

############################
# PRE-PROCESSING
############################
#-----------
cat("#Set parameters generic names")
cat("\n")
input_file <- args$input_file
output_dir <- args$output_dir
plot_name <- args$plot_name
cluster_rows <- args$cluster_rows
store_rows_cluster <- args$store_rows_cluster 

#-----------
cat("# Reading the matrix of motif")
cat("\n")
matrix <- read.delim(input_file, sep = "\t", header = T)
matrix_cust <- matrix %>% rename(motifs = !!names(matrix)[1]) %>% column_to_rownames(var="motifs")

cat("# Converting the matrix values in numeric")
cat("\n")
matrix_mod <- as.data.frame(lapply(matrix_cust, as.numeric))
rownames(matrix_mod) <- rownames(matrix_cust)

#-----------
cat("# Scaling the matrix of motif and order it based on SI.")
cat("\n")
matrix_mod_t <- as.data.frame(scale(t(matrix_mod))) %>% rownames_to_column(var="SI") %>% mutate(num_SI = parse_number(SI)) %>% 
  select(-SI) %>% arrange(as.numeric(num_SI)) %>% column_to_rownames(var="num_SI") #%>% 

#-----------
cat("# Removing motifs with all zero.")
cat("\n")
motifs_score_matrix <- as.data.frame(t(matrix_mod_t)) %>% mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% #substitute the NA with zero  
  mutate_all(as.numeric) %>% filter(rowSums(.) != 0)
matrix_heatmap <- as.matrix(motifs_score_matrix)

##-----------------------------------------#
#   COLUMN ANNOTATION
#-----------------------------------------#
#-----------
cat("# Retrive the column annotation.")
cat("\n")
metadata_col <- data.frame(SI=as.numeric(colnames(matrix_heatmap)))
# Define annotation for the col (SI)
mycols_col <- colorRamp2(breaks = c(0,15,30,45,55,70), colors = c( "#E9F6F2", "#D3EDE5", "#BDE4D8", "#A7DCCB", "#91D3BE", "#7BCAB1"))

colAnn <- HeatmapAnnotation(df = metadata_col,
                            which = 'col',
                            col = list(SI=mycols_col),
                            annotation_label = c("SI"))


#-----------------------------------------#
#   HEATMAP
#-----------------------------------------#
#-----------
cat("# Making the heatmap.")
cat("\n")
mycols <- colorRamp2(breaks = c(-4,-2,0,2,4), colors = c("#08519C","#4292C6","#FFFFFF","#D7301F","#B30000"))#roma
#mycols <- colorRamp2(breaks = c(-4,-2,0,2,4), colors = c("#238B45","#41AB5D","#FFFFFF","#D7301F","#B30000"))#nilson
#mycols <- colorRamp2(breaks = c(-4,-2,0,2,4), colors = c("#525252", "#BDBDBD", "#FFFFFF", "#D7301F","#B30000"))#roma-hint
distance_row="manhattan"
method_row="ward.D2"

h <- Heatmap(matrix_heatmap,
            column_title = plot_name,
            name = "Score(obs/exp)", 
            column_names_gp = gpar(fontsize = 1),
            col = mycols,
            show_row_names=F,
            top_annotation=colAnn,
            clustering_distance_rows = distance_row,
            clustering_method_rows = method_row,
            row_split = cluster_rows, 
            border = TRUE,
            cluster_columns = F,
            cluster_rows=T
)

#-----------
cat("# Save the heatmap.")
cat("\n")  
name_file_png=paste0(output_dir, "/",plot_name,"_",distance_row,"_",method_row,".png",sep="")
ggexport(h,res=300,filename=name_file_png,width = 3000, height = 2000)
  
name_file_svg=paste0(output_dir, "/",plot_name,"_",distance_row,"_",method_row,".svg",sep="")
svg(name_file_svg, width = 12, height = 9)
print(h)
dev.off()

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
    c <- motifs_score_matrix[rcl.list[[i]],] %>% as.data.frame() %>% rownames_to_column() %>% select("rowname")
    file_path <- file.path(output_dir, paste0("cluster_motif_", i,"_",plot_name,".txt"))
    write.table(c,file_path, sep="\t", quote=F, row.names=F,col.names = FALSE)
  }
}
  