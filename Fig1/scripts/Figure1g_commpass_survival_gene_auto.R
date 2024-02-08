# DESCRIPTION: Compute a survival analysis with COMPASS data by considering specified gene as input and subset the population based on the median or 3rd quantile of the gene expression.
# NOTA0: The script is compatible with the A17 version of COMMPASS.
# NOTA1: If you wish to use the script with a different version of COMMPASS, you need to adjust the name check, specifically the cell name (e.g., CD138). By default, the script considers the Specimen_ID, which appears like "MM138_1_BM". However, if this doesn't match the name in the reads count file, you will need to modify it. For instance, in the AI17 case, we had to add the "CD138" suffix to the Specimen_ID to ensure a match between the reads and metadata files.
#-----------------------------
# Upload packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(HGNChelper))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

# Define the command-line arguments
parser <- ArgumentParser(description = "Create a barplot from a tab-separated file.")
parser$add_argument("--gene_name", "-i", type = "character", help = "Input the gene name to look for.")
parser$add_argument("--reads_file", "-r", type = "character", help = "Input the path of with the normalized reads counts.")
parser$add_argument("--survival_file", "-surv", type = "character", help = "Input the path of the file with survival info.")
parser$add_argument("--output_dir", "-o", type = "character", help = "Output directory.")
parser$add_argument("--plot_name", "-p", type = "character", help = "Plot name.")
parser$add_argument("--cell","-c", type = "character", help = "Cell type detected.",default="")
parser$add_argument("--division_type", "-t", type = "character", help = "Define the method for diving the population into two subgroups(median or 3rd quantile)", default="median")
# parse the command-line arguments
args <- parser$parse_args()
############################
# INTRODUCTIVE PROCEDURES
############################

#-----------
cat("#Set parameters generic names")
cat("\n")
cell <- args$cell
division_type=args$division_type
reads_file <- args$reads_file
survival_file <- args$survival_file

#-----------
cat("# Reads the genes name to examine")
cat("\n")
gene_to_look <- args$gene_name

#-----------
cat("# Upload the COMMPASS file containing normalized reads counts of all transcripts examined")
cat("\n")
norm_reads_counts <- read.delim(reads_file, sep = "\t", header = T)

#-----------
cat("# Retrieve biomart...")
mart37 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")#,host = "grch37.ensembl.org")
genes37 <-  biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart=mart37)
# Retrieve the Ensembl ID of the gene passed as input
ens_gene <- genes37 %>% filter(hgnc_symbol==gene_to_look) 

#-----------
cat("# Subset the normalized reads counts to those selected genes")
cat("\n")
reads_sub <- merge(norm_reads_counts,ens_gene,by.x="Gene",by.y="ensembl_gene_id") %>% distinct() %>% rename(gene_name=hgnc_symbol)

#--------------------
cat("# Check if the selected gene is present in COMMPASS..")
cat("\n")
if (nrow(reads_sub) == 0) {
  cat("# Sorry, the gene is not present in COMMPASS!!")
} else {
  #----------
  cat("# Great, the gene is present in COMMPASS!!")
  cat("\n")
  cat("# Customize the subset normalized reads count matrix and transpose it")
  cat("\n")
  reads_sub_mod <- reads_sub %>% select(-Gene) %>%  column_to_rownames(var="gene_name")
  t_reads_sub_mod <- as.data.frame(t(reads_sub_mod)) %>% rownames_to_column(var = "Specimen_ID") 

  #-----------------------------------------------#
  # Assign the low/high level based on the median or on the quantile
  #-----------------------------------------------#
  if(division_type == "median"){ 
  
    t_reads_sub_mod <- t_reads_sub_mod %>% 
      mutate(median_gene=median(t_reads_sub_mod[,2])) %>% 
      mutate(status=case_when(.data[[gene_to_look]] <= median_gene ~ "low", TRUE ~ "high")) %>% 
      separate(Specimen_ID,c("project","id","visit","space"),sep = "_") %>% filter(visit==1) %>% filter_all(all_vars(!is.na(.))) %>% 
      mutate(PUBLIC_ID=paste0(project,"_",id)) %>% 
      select(-project,-id)
  } else if(division_type == "quantile") {
    t_reads_sub_mod <- t_reads_sub_mod %>% 
      mutate(quantile_gene=quantile(t_reads_sub_mod[,2],0.75)) %>% 
      mutate(status=case_when(.data[[gene_to_look]] <= quantile_gene ~ "low", TRUE ~ "high")) %>% 
      separate(Specimen_ID,c("project","id","visit"),sep = "_") %>% filter(visit==1) %>% filter_all(all_vars(!is.na(.))) %>% 
      mutate(PUBLIC_ID=paste0(project,"_",id)) %>% 
      select(-project,-id)
  }else{
    cat("# Error method selection: you can only choose between median and quantile!!")
    cat("\n")
  }
  
  
  #-----------
  cat("# Reads COMMPASS file with survival info")
  cat("\n")
  survival <- read.csv(survival_file, sep = "\t", header = T)
  survival <- survival %>%dplyr::select(PUBLIC_ID,oscdy,censos) 
  
  #-----------
  cat("# Convert days into months")
  cat("\n")
  df_surv <- merge(t_reads_sub_mod, survival, by="PUBLIC_ID")
  df_surv <- df_surv %>% mutate(oscdy=round(oscdy/30,0))
  
  #-----------
  cat("# Filter out patients based on months")
  cat("\n")
  df_surv <- df_surv %>% filter(oscdy<=65)
  
  cat("# Storing the classification and surv information")
  cat("\n")
  file_path <- file.path(args$output_dir, paste0("/patient_classification_based_on_",division_type,"_expression_of_",gene_to_look,".txt"))
  write.table(df_surv,file_path, sep="\t", quote=F, row.names=F,col.names = T)
  
  #-----------
  cat("# Fit the survival")
  cat("\n")
  fit <- survfit(Surv(oscdy, censos) ~ status, data = df_surv)
  
  #-----------
  cat("# Compute the median suvival time")
  cat("\n")
  time <- surv_median(fit)
  
  #-----------
  cat("# Make the survival plot")
  cat("\n")
  
  blank_theme <- theme_bw()+ theme(panel.grid = element_blank())
  
  gg <- ggsurvplot(fit,
                   risk.table = TRUE, 
                   xlim = c(0, 59),
                   break.x.by = 6,  break.y.by = .1, axes.offset = FALSE,
                   pval = T, 
                   surv.median.line = "hv",
                   conf.int = F,
                   size = 1.5,
                   palette=c("#efbf9a","#8B9DAC"),
                   xlab = "Months", ylab = "Probability of Survival",
                   legend = c(0.85,0.9),
                   legend.title = "",
                   font.legend = 14,
                   font.xlab = c(14, "bold"),
                   font.ylab = c(14, "bold"),
                   risk.table.title = "N. patients at risk", risk.table.y.text = F, ggtheme = theme_pubr()) +
    labs(title = paste0(gene_to_look,"_",division_type))
  
  print(gg)
  #-----------
  cat("# Saving the plot")
  cat("\n")
  ggexport(gg,res=300,filename=paste0(args$output_dir, "/", "survival_cluster_",args$plot_name,".png"),width = 2500, height = 2000)
  

}

