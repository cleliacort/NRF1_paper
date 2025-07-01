# DESCRIPTION:
#   Perform a Principal Component Analysis (PCA) using a data matrix and a phenotype file.
#
# INPUT:
#   - Data matrix: The first three rows contain 'chromosome', 'start', and 'end' information.
#                  Each subsequent column represents a sample.
#
#   - Phenotype file: Must contain a column named 'ID' (matching sample names in the matrix)
#                     and a column named 'phenotype' with phenotype information for each sample.
#-----------------------------
# Upload packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))

# define the command-line arguments
parser <- ArgumentParser(description = "Create a barplot from a tab-separated file.")
parser$add_argument("--input_file", "-i", type = "character", help = "Input the gene counts matrix of your samples.")
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.")
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name with the extension.")
parser$add_argument("--phenotype", "-pheno", type = "character", help = "File containing the groups to which each sample belongs.")
parser$add_argument("--labels", "-l", type = "character", default="FALSE", help = "Indicate whether you want the label to be TRUE or FALSE.")

# parse the command-line arguments
args <- parser$parse_args()


#-----------
cat("#Set parameters generic names")
cat("\n")
count_matrix <- args$input_file
output_dir <- args$output_dir
plot_name <- args$plot_name
phenotype <- args$phenotype
labels <- args$labels

# count_matrix <- "data/Suppl/matrix_multicov_tumout_and_MGUS_on_tumour_master_list_0124.txt"
# plot_name <- "PCA_mgus_plus_tumour_treated_plus_tumour_onset_1224.png"
# phenotype <- "data/Suppl/sample_sheet_clinical_PHENOTYPE_tumour_and_mgus_0124.txt"
# labels <- "FALSE"
# output_dir <- "figures/Suppl"

#-----------
cat("# Upload the raw count matrix")
cat("\n")
count_matrix_file <- read.delim(count_matrix, header = TRUE, sep = "\t") 
count_matrix_cust <- count_matrix_file %>%   
  dplyr::rename(chr=colnames(count_matrix_file)[1],start=colnames(count_matrix_file)[2],end=colnames(count_matrix_file)[3]) %>% 
  mutate(coordinates = paste0(chr,"_",start,"_",end)) %>% 
  select(-chr,-start,-end) %>% 
  column_to_rownames(var="coordinates")

#-----------
cat("# Filter the raw count matrix")
cat("\n")
e=data.matrix(count_matrix_cust)
e=subset(e, rowMeans(e) > 5)
e=subset(e, rowMeans(e) < 5000)



#-----------
cat("# Retrieve the phenotype info")
cat("\n")
pheno_file <- read.delim(phenotype, header = TRUE, sep = "\t") 
names <- data.frame(ID=colnames(count_matrix_cust))
groups_pheno <- merge(names,pheno_file,by="ID") 


#-----------
cat("# Normalize the raw counts")
cat("\n")
regions <- DGEList(counts=e, group=groups_pheno$phenotype)
regions <- estimateCommonDisp(regions)
regions <- estimateTagwiseDisp(regions)
regions <- calcNormFactors(regions, method="TMM")
cpm_tmm <- cpm(regions, normalized.lib.size=T)

#-----------
cat("# Compute PCA analysis")
cat("\n")
t_cpm <- as.data.frame(t(cpm_tmm))
t_cpm_pca <- prcomp(t_cpm,scale. = F)

#-----------
cat("# Customize the PCA results")
cat("\n")
var_explained=round(100*t_cpm_pca$sdev^2/sum(t_cpm_pca$sdev^2),1)
df_pca <-  data.frame(PC1 = t_cpm_pca$x[,1], PC2 = t_cpm_pca$x[,2]) %>% rownames_to_column(var="ID")
df_pca_pheno <- merge(df_pca,pheno_file,by="ID")

#-----------
cat("# Defining the color palette")
cat("\n")
pheno_mod <- pheno_file %>% group_by(phenotype) %>% summarise()
color_num <- nrow(pheno_mod)

cust_palette <- c("#8B9DAC","#FED9A6","#DECBE4")
# if(color_num>2){
#   cust_palette <- brewer.pal(n = color_num, name = "Pastel1")
# }

#-----------
cat("# Visualize the PCA results")
cat("\n")
pca_plot <- ggplot(df_pca_pheno, aes(PC1,PC2, color = type))+
 geom_point(aes(colour=phenotype,size=3))+
  scale_color_manual(values = cust_palette)+
  labs(x=paste0("PC1(",round(var_explained[1],1),"%)"),
       y=paste0("PC2(",round(var_explained[2],1),"%)"))+
  theme_classic() +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


if(labels=="TRUE"){
  text <- geom_text_repel(aes(label = ID), color="#525252", size=2, vjust=0.8)
  pca_plot <- pca_plot+text 
}  


# show(pca_plot)

#-----------
cat("# Saving the plot")
cat("\n")
name_file_plot=file.path(output_dir,plot_name) 
ggexport(pca_plot,res=300,filename=name_file_plot,width = 3000, height = 2000)


