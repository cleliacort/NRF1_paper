# Figure1d_correlation_with_Nilson_dataset

Upload the needed packages.

``` r
library(tidyverse)
library(here)
library(ggpubr)
library(ggplot2)
library(ggrepel)
```

Reading the footprinting percentage enrichment data from Roma and Nilson
dataset.

``` r
df_roma <- read_delim(here("data","master_list_MOTIF_LOGMINU20_GROUPED_SI_HG19_ROMA_max_AGOSTO23.txt"),delim = "\t", col_names = T) %>% column_to_rownames(var="motif") 

df_nilson <- read_delim(here("data","master_list_MOTIF_LOGMINU20_GROUPED_SI_HG19_NILSON_max_AGOSTO23.txt"),delim = "\t", col_names = T) %>% column_to_rownames(var="motif")
```

Customize and merge footprinting percentage enrichment from Roma and
Nilson dataset.

``` r
matrix_roma <- as.matrix(df_roma)
matrix_nilson <- as.matrix(df_nilson)

tf_names_nilson <- data.frame(tf_names=rownames(matrix_nilson))
tf_names_roma <- data.frame(tf_names=rownames(matrix_roma))
tf_names_comm <- merge(tf_names_roma,tf_names_nilson,by="tf_names")
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
  index2 <- match(row_names[i], rownames(matrix_nilson))
  
  # Extract the rows from each matrix and calculate the correlation
  row1 <- as.numeric(matrix_roma[index1,])
  row2 <- as.numeric(matrix_nilson[index2,])
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
col_cust="#41AB5D" #nilson

plot <- ggplot(df_results,aes(y= correlation_spearman, x=rank)) + 
  geom_point(aes(size = 8),colour="gray74")+ #set features of dots
  geom_point(data = df_label_color,aes(x= rank, y=correlation_spearman),size=5,colour=col_cust)+
 geom_text_repel(data = df_label, aes(label = TF_name, x = rank, y = correlation_spearman),min.segment.length = 0.5,box.padding = 0.5) +
  scale_fill_identity() +
  xlim(0,250)+
  theme_bw()+ 
  theme(panel.grid = element_blank())+ 
  ylab("Spearman Correlation")+xlab("Ranking of Spearman correlation")+
  ggtitle("Correlation between Roma WG and Nilson WG footprinting")+
  guides(size = FALSE) 
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
#show(plot)
```

<figure>
<img
src="https://github.com/cleliacort/NRF1_paper/blob/main/Fig1/figures/spearman_correlation_ROMA_vs_NILSON_footprinting_wg_percentage_labels.png"
alt="Figure1d_correlation_with_Nilson_dataset" />
<figcaption
aria-hidden="true">Figure1d_correlation_with_Nilson_dataset</figcaption>
</figure>
