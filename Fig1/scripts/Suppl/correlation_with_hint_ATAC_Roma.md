# Correlation_with_hint_ATAC_Roma

Upload the needed packages.

``` r
library(tidyverse)
library(here)
library(ggpubr)
library(ggplot2)
library(ggrepel)
```

Reading the footprinting percentage enrichment data from Roma dataset
computed with wellington algorithm and hint-ATAC algorithm.

``` r
df_roma <- read_delim(here("data","master_list_MOTIF_LOGMINU20_GROUPED_SI_HG19_ROMA_max_AGOSTO23.txt"),delim = "\t", col_names = T) %>% column_to_rownames(var="motif") 

df_hint_atac_roma <- read_delim(here("data","master_list_MOTIF_GREATER15_GROUPED_SI_HG19_ROMA_max_HINT_ATAC_AGOSTO23.txt"),delim = "\t", col_names = T) %>% column_to_rownames(var="motif")
```

Customize and merge footprinting percentage enrichment from Roma and
hint_atac_roma dataset.

``` r
matrix_roma <- as.matrix(df_roma)
matrix_hint_atac_roma <- as.matrix(df_hint_atac_roma)

tf_names_hint_atac_roma <- data.frame(tf_names=rownames(matrix_hint_atac_roma))
tf_names_roma <- data.frame(tf_names=rownames(matrix_roma))
tf_names_comm <- merge(tf_names_roma,tf_names_hint_atac_roma,by="tf_names")
```

Compute the correlation coefficient and store information into a
customized dataframe.

``` r
# Define the name of the rows to compare
row_names <- c(tf_names_comm$tf_names)

# Initialize a vector to store the correlations
correlations <- vector(length=length(row_names))

# initialize a datafram in r
df_results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_results) <- c("TF_name","correlation_spearman")

# Loop over each row name and calculate the correlation
for (i in 1:length(row_names)) {
  # Find the index of the row in each matrix
  index1 <- match(row_names[i], rownames(matrix_roma))
  index2 <- match(row_names[i], rownames(matrix_hint_atac_roma))
  
  # Extract the rows from each matrix and calculate the correlation
  row1 <- as.numeric(matrix_roma[index1,])
  row2 <- as.numeric(matrix_hint_atac_roma[index2,])
  correlation <- cor(row1, row2,method = "spearman")
  
  df_tmp <- data.frame(TF_name= row_names[i],correlation_spearman=correlation )
  df_results <- rbind(df_results,df_tmp)
}
```

Construct a dataframe for data labeling.

``` r
df_results <- na.omit(df_results)
df_results <- df_results %>% arrange(correlation_spearman) %>%  mutate(rank = rank(correlation_spearman)) 
df_label_color <- df_results  %>% filter(correlation_spearman>=0.55) %>% separate(TF_name,c("TF_name",NA),sep="/") 

# Define the dataframe for the labels names
head_label <- df_label_color %>% arrange(-rank) %>% slice_head(n=5)
nomi <- c(head_label$TF_name,"NRF(NRF)", "NRF1(NRF)","IRF4(IRF)","IRF8(IRF)") 

# Retrieve the information for the plotting related to TF of interest
df_label <- df_results  %>% filter(correlation_spearman>=0.55) %>% separate(TF_name,c("TF_name",NA),sep="/") %>% filter(TF_name %in% nomi)
```

Create the plot with correlation information.

``` r
col_cust="#807DBA" #hint

plot <- ggplot(df_results,aes(y= correlation_spearman, x=rank)) + 
  geom_point(aes(size = 8),colour="gray74")+ #set features of dots
  geom_point(data = df_label_color,aes(x= rank, y=correlation_spearman),size=5,colour=col_cust)+
 geom_text_repel(data = df_label, aes(label = TF_name, x = rank, y = correlation_spearman),min.segment.length = 0.5,box.padding = 0.5) +
  scale_fill_identity() +
  xlim(0,300)+
  theme_bw()+ 
  theme(panel.grid = element_blank())+ 
  ylab("Spearman Correlation")+xlab("Ranking of Spearman correlation")+
  ggtitle("Correlation between Roma WG and hint_atac_roma WG footprinting")+
  guides(size = FALSE) 
#show(plot)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/Suppl/spearman_correlation_ROMA_wg_footprinting_ROMA_hint_atac_roma_footprinting_percentage_labels.png"
alt="correlation_with_hint_ATAC_Roma" />
<figcaption
aria-hidden="true">correlation_with_hint_ATAC_Roma</figcaption>
</figure>
