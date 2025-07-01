# Upload packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(HGNChelper))

# define the command-line arguments
parser <- ArgumentParser(description = "Create a barplot from a tab-separated file.")
parser$add_argument("--input_file", "-i", type = "character", help = "Input the path of the file with sample and number information.",required = TRUE)
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.",required = TRUE)
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name.",default="Number_of_lines.txt")
#parser$add_argument("--colour_sel", "-c", type = "character", help = "Color of highlighted cell lines.",default="#659CC6")
parser$add_argument("--x_label","-x_lab", type = "character", help = "Label for x-axis.")
parser$add_argument("--y_label","-y_lab", type = "character", help = "Label for y-axis.")
parser$add_argument("--x_axis","-x", type = "character", help = "Name of the column data to show on x-axis.",required = TRUE)
parser$add_argument("--y_axis","-y", type = "character", help = "Name of the column data to show on x-axis.",required = TRUE)
parser$add_argument("--phenotype","-pheno", type = "character", help = "A file with phenotype information about the sample.",required = TRUE)

# parse the command-line arguments
args <- parser$parse_args()

############################
# PRE-PROCESSING
############################
#-----------
cat("#Set parameters generic names")
cat("\n")
input_file <- args$input_file
output_dir <- args$output_dir
plot_name <- args$plot_name
x_label <- args$x_label
y_label <- args$y_label
x_axis <- args$x_axis
y_axis <- args$y_axis
phenotype <- args$phenotype
# input_file <- "/Users/cleliacortile/Desktop/Lavoro/atac_tumour_mgus_2023/Analysis_august23_check/number_of_peaks_per_sample_atac_tumour_mgus_august23.txt"
# phenotype <- "/Users/cleliacortile/Desktop/Lavoro/atac_tumour_mgus_2023/Analysis/samplesheet/sample_sheet_official_clinical_2023_subsetted_PHENOTYPE.csv"
# x_axis="Filename"

#-----------
cat("# Reading the list of files")
cat("\n")
list_file <- read.delim(input_file, sep = "\t", header = T) %>% 
  rename( ID = .data[[x_axis]]) %>% 
  rename( NumLines = .data[[y_axis]])
#-----------
cat("# Reading the phenotype file")
cat("\n")
pheno <- read.delim(phenotype, sep = "\t", header = T)
colnames(pheno) <- c("ID","phenotype")
#-----------
cat("# Defining the color palette")
cat("\n")
pheno_mod <- pheno %>% group_by(phenotype) %>% summarise()
color_num <- nrow(pheno_mod)
#-----------
cat("# Customizing the dataframe")
cat("\n")
df_plot <- merge(list_file,pheno, by = "ID")
#-----------
cat("# Creating the plot")
cat("\n")
cust_palette <- c("#FED9A6","#DECBE4")
if(color_num>2){
  cust_palette <- brewer.pal(n = color_num, name = "Pastel1")
}

gg <- ggbarplot(df_plot, x = "ID", y = "NumLines",
               fill = "phenotype",
               color = "black",            # Set bar border colors to white
               palette = cust_palette,            # jco journal color palett. see ?ggpar
               sort.val = "asc",           # Sort the value in dscending order
               sort.by.groups = TRUE,      # Sort inside each group
               x.text.angle = 90)
# show(p)
gg <- ggpar(gg,xlab = x_label,ylab = y_label)

#-----------
cat("# Saving the plot")
cat("\n")
complete_name_plot_png=paste0(output_dir, "/",plot_name,".png",sep="")
ggexport(gg,res=300,filename=complete_name_plot_png,width = 4500, height = 2000)

complete_name_plot_svg=paste0(output_dir, "/",plot_name,".svg",sep="")
ggsave(complete_name_plot_svg, plot = gg)








