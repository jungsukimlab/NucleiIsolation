set.seed(1221)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

library(UCell)
library(Seurat)
library(ggplot2)
library(scMCA)
#load seurat object
isolation <- readRDS("20231206_snIsolation_collab.rds")

#Cluster markers
cluster.markers <- FindAllMarkers(isolation, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% 
  group_by(cluster) %>%
  slice_max(n=10, order_by = avg_log2FC) %>%
  print(n=400)
write.csv(cluster.markers, "markers.csv")

#scMCA
p <- isolation@assays$RNA@counts
#mca result - take 2-3 hrs for large datasets
mca_result <- scMCA(scdata = p, numbers_plot = 3)
x=mca_result$scMCA_probility
saveRDS(x,file = "isolation_scMCA.rds")

#Annotations
cluster.short <- c("ExN", "AS", "ExN", "ExN", "ExN", "ExN", "O", "ExN", "ExN", "O", "InN","InN", "InN", "MG", "InN","ExN", "OPC", "ExN", "InN", "ExN", "ExN", "AS", "ExN", "ExN", "AS", "MG", "ExN", "EC", "FB", "O", "ExN", "InN", "AS", "RBC", "O", "ExN")
names(cluster.short) <- levels(isolation)
isolation <- RenameIdents(isolation, cluster.short)
isolation[["Cluster_Short"]] <-Idents(isolation)
Idents(isolation) <- "Cluster_Short"
DimPlot(isolation)

cluster.long <- c("Excitatory Neuron", "Astrocyte", "Excitatory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Oligodendrocyte", "Excitatory Neuron", "Excitatory Neuron", "Oligodendrocyte", "Inhibitory Neuron","Inhibitory Neuron", "Inhibitory Neuron", "Macrophage", "Inhibitory Neuron","Excitatory Neuron", "Oligodendrocyte Progenitor Cell", "Excitatory Neuron", "Inhibitory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Astrocyte", "Excitatory Neuron", "Excitatory Neuron", "Astrocyte", "Macrophage", "Excitatory Neuron", "Endothelial Cell", "Fibroblast", "Oligodendrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Red Blood Cell", "Oligodendrocyte", "Excitatory Neuron")
names(cluster.long) <- levels(isolation)
isolation <- RenameIdents(isolation, cluster.long)
isolation[["Cluster_Long"]] <-Idents(isolation)
Idents(isolation) <- "Cluster_Long"
DimPlot(isolation)

isolation$Cluster_Long <- as.factor(isolation$Cluster_Long)
isolation$Cluster_Long <- factor(isolation$Cluster_Long, levels = c("Astrocyte", "Endothelial Cell", "Fibroblast", "Macrophage", "Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte", "Oligodendrocyte Progenitor Cell", "Red Blood Cell"))

cluster.simple <- c("Excitatory Neuron", "Astrocyte", "Excitatory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Oligo", "Excitatory Neuron", "Excitatory Neuron", "Oligo", "Inhibitory Neuron","Inhibitory Neuron", "Inhibitory Neuron", "Microglia", "Inhibitory Neuron","Excitatory Neuron", "Oligo", "Excitatory Neuron", "Inhibitory Neuron", "Excitatory Neuron", "Excitatory Neuron", "Astrocyte", "Excitatory Neuron", "Excitatory Neuron", "Astrocyte", "Other", "Excitatory Neuron", "Other", "Other", "Oligo", "Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Other", "Oligo", "Excitatory Neuron")
Idents(isolation) <- "seurat_clusters"
names(cluster.simple) <- levels(isolation)
isolation <- RenameIdents(isolation, cluster.simple)
isolation[["Cluster_Simple"]] <-Idents(isolation)
Idents(isolation) <- "Cluster_Simple"
DimPlot(isolation, reduction = "umap")

isolation$Cluster_Simple <- factor(isolation$Cluster_Simple, levels = c("Astrocyte", "Microglia", "Excitatory Neuron", "Inhibitory Neuron", "Oligo", "Other"))
Idents(isolation) <- "Cluster_Simple"
DimPlot(isolation, reduction="umap")

Idents(isolation) <- "seurat_clusters"
cluster.shorter <- c("eN", "A", "eN", "eN", "eN", "eN", "O", "eN", "eN", "O", "iN","iN", "iN", "M", "iN","eN", "O", "eN", "iN", "eN", "eN", "A", "eN", "eN", "A", "Ot", "eN", "Ot", "Ot", "O", "eN", "iN", "A", "Ot", "O", "eN")
names(cluster.shorter) <- levels(isolation)
isolation <- RenameIdents(isolation, cluster.shorter)
isolation[["Cluster_Shorter"]] <-Idents(isolation)
Idents(isolation) <- "Cluster_Shorter"
DimPlot(isolation)

isolation$Num_Shorter <- paste0(isolation$seurat_clusters, "_",isolation$Cluster_Shorter)
Idents(isolation) <- "Num_Shorter"
DimPlot(isolation)
Idents(isolation) <- "seurat_clusters"
DimPlot(isolation)

isolation$Num_Shorter <- factor(isolation$Num_Shorter, levels = c("0_eN", "1_A", "2_eN", "3_eN", "4_eN", "5_eN", "6_O", "7_eN", "8_eN", "9_O", "10_iN","11_iN", "12_iN", "13_M", "14_iN","15_eN", "16_O", "17_eN", "18_iN", "19_eN", "20_eN", "21_A", "22_eN", "23_eN", "24_A", "25_Ot", "26_eN", "27_Ot", "28_Ot", "29_O", "30_eN", "31_iN", "32_A", "33_Ot", "34_O", "35_eN"))
Idents(isolation) <- "Num_Shorter"

DimPlot(isolation)+
  theme_bw()+
  theme(axis.text = element_text(size=8), axis.title = element_text(size = 8))+
  theme(legend.text = element_text(size = 6))+
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size=4), ncol=3))
