set.seed(1221)
library(Seurat)
library(ggplot2)
library(dplyr)
setwd("/N/project/kim_lab/Mustaklem_Richard/IsolationCollab")
file <- readRDS("20231206_snIsolation_collab.rds")

View(file)
# Get unique technologies
unique(file@meta.data[["tech"]])

length(file@meta.data$tech)

sum(file@meta.data$tech=="10x_kit")

sum(file@meta.data$tech=="Singulator")

sum(file@meta.data$tech=="Sucrose")

summary(factor(file@meta.data$tech))
 
# Subset the dataset 

TENX_kit_file <- subset(file, subset = tech == "10x_kit")

Singulator_file <- subset(file, subset = tech == "Singulator")

Sucrose_file <- subset(file, subset = tech == "Sucrose")

Idents(TENX_kit_file) <- "annotated"

Idents(Singulator_file) <- "annotated"

Idents(Sucrose_file) <- "annotated"
 
 
# Change the file depending on which technology we are using

markers <- FindAllMarkers(Sucrose_file, only.pos = TRUE, logfc.threshold = 0.25)

tech_used <- unique(Sucrose_file@meta.data$tech)
 
# Read in PanglaoDB

# Need to download this file into working directory

data=read.delim("PanglaoDB_markers_27_Mar_2020.tsv",sep = "\t",header = T) 

brain_data=subset(data,organ=="Brain") # Subset only brain markers

summary(factor(brain_data$cell.type)) # What cell types exist in brain markers
 
# How to look up a gene 

hmf=read.csv("20230911_HumanMouseFly_gene_nomenclature.csv") #Read in mouse/human/fly nomenclature file

na.omit(hmf[hmf$human_symbol.x=="MAPT",]) 

na.omit(hmf[hmf$mouse_symbol=="Mapt",])

na.omit(hmf[hmf$fruit.fly_symbol=="tau",])
 
# Creating SAHA loops

# Use to look up any cluster. 

saha_lookup_cluster <- function(cluster) {

  i=markers[markers$cluster==cluster,"gene"] # Need to convert i to human genes for pangloa

  j=unique(na.omit(hmf[hmf$mouse_symbol%in%i,"human_symbol.x"]))# Convert genes to human orthologs

  return(summary(factor(brain_data[brain_data$official.gene.symbol %in% j, "cell.type"], levels(factor(brain_data$cell.type)))))  

}
 
# Write a function to graph the cell types quickly using ggplot

saha_graph <- function(dataframe){

  sumdata_df <- as.data.frame(dataframe)

  colnames(sumdata_df) <- c("Count")

  ggplot(data = sumdata_df, aes(x = rownames(sumdata_df), y = Count)) +

    geom_bar(stat = "identity", fill = "Black") +

    labs(title = "Cell Type Counts", x = "", y = "Count") +

    theme(axis.text.x = element_text(angle = 90, hjust = 1))

}
 
# Making interactive graph for respective cluster

# ggplotly(saha_graph(saha_lookup_cluster("Microglia")))
 
# Making our Master Data Frame

master_df <- data.frame(summary(factor(brain_data$cell.type)))

master_df$cluster="REF"

colnames(master_df)[1]<-"total_marker" 

master_df$celltype <- rownames(master_df)
 
# Loop through each cluster, find their markers, bind to masterdf 

for (i in c(1:length(unique(markers$cluster)))) {

  x=master_df[master_df$cluster=="REF",]

  x$cluster = unique(markers$cluster)[i]

  x$total_marker=saha_lookup_cluster(unique(markers$cluster)[i])

  x$celltype = rownames(x)

  master_df=rbind(master_df,x)

}
 
# Find proportion

master_df$prop="NA"

for (j in c(1:length(unique(markers$cluster)))) {

  master_df[master_df$cluster==unique(markers$cluster)[j],"prop"]=master_df[master_df$cluster==unique(markers$cluster)[j],"total_marker"]/master_df[master_df$cluster=="REF","total_marker"]  

}
 
# Find p value

master_df$pvalue = 1 

for (j in c(1:length(unique(markers$cluster)))) {

  # Size of possible markers in a given cluster

  A = length(markers[markers$cluster == unique(markers$cluster)[j], "gene"]) 

  # Size of signature testing (number of genes in pangloa cell type)

  for (k in c(1:length(unique(master_df$celltype)))) {

    B = master_df[master_df$cluster=="REF",][k,"total_marker"] 

    # Overlap 

    t = master_df[master_df$cluster==unique(markers$cluster)[j],][k,"total_marker"] 

    # Length of all possible cluster markers (not just the one testing)

    n = length(unique(rownames(markers))) 

    master_df[master_df$cluster==unique(markers$cluster)[j],][k,"pvalue"]=sum(stats::dhyper(t:B, A, n - A, B))

  }

}
 
# Default everything to F, not significant

master_df$sig = "F"

master_df[master_df$pvalue <= 0.05,"sig"]="T" # Label only significant groups
 
# Convert the 'cluster' variable to a factor with custom levels in ascending order

master_df$cluster <- factor(master_df$cluster, levels = unique(master_df$cluster))
 
# Create dataframe that only finds minimum pvalue for each cluster

master_df <- master_df %>%

  group_by(cluster) %>%

  mutate(top_sig = ifelse(pvalue == min(pvalue) & sig == "T", "T", "F"),

         top_sig_prop = ifelse(top_sig == "F", 0, as.numeric(prop))) %>%

  ungroup()
 
master_df <- data.frame(master_df)
 
# Ordering the Y Axis

unique_clusters <- unique(master_df$cluster)[-1] # Remove REF 

# Split the clusters into alphabetical and numeric parts

cluster_parts <- lapply(unique_clusters, function(x) {

  numeric_part <- gsub("\\D", "", x)#Any non-digit character.

  # Sometimes the cluster has no number, use if statement to check

  if (numeric_part == "") {

    return(list(x, NA))

  } else {

    return(list(gsub("\\d", "", x), as.numeric(numeric_part)))

  }

})

# Order clusters by both alphabetical and numeric parts

ordered_clusters <- unique_clusters[order(sapply(cluster_parts, function(x) sprintf("%s%03d", x[[1]], x[[2]])))]

# Set the levels for the factor variable

master_df$cluster <- factor(master_df$cluster, levels = ordered_clusters)
 
 
# Create the ggplot plot with the ordered 'cluster' variable -- ALL

ggplot(data = subset(master_df, cluster != "REF"), aes(x = celltype, y = cluster, size = as.numeric(prop), alpha = as.numeric(-log(pvalue,10)))) +

  geom_point(aes(color = sig)) +

  labs(title = paste(tech_used," - Every Cell Type by Cluster"), x = " ", y = "Cluster") +

  scale_color_manual(values = c("black", "red")) +

  theme_bw() +

  theme(axis.text.x = element_text(angle = 90, hjust = 1),

        plot.title = element_text(size = 10))
 
 
#10xkit, singulator, sucrose
 
 
