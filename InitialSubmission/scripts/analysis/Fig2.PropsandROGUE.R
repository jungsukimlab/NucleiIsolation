set.seed(1221)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

library(UCell)
library(Seurat)
library(ggplot2)
library(propeller)
library(ROGUE)

#load seurat object
isolation <- readRDS("20231206_snIsolation_collab.rds")

#Proportions
x = propeller(clusters = isolation$Cluster_Simple, 
              sample = isolation$orig.ident, 
              group = isolation$tech)
x
write.csv(x, "Isolation_Simple_by_Tech_Props.csv")
y = plotCellTypeProps(clusters = isolation$Cluster_Simple, 
                      sample = isolation$tech)
y
key = c(levels(factor(isolation$orig.ident)))
key= data.frame(key)
key$Tech <- c("Kit", "Kit", "Machine", "Machine", "Custom", "Custom")
key
colnames(key) <- c("Samples", "Tech")
key
z= y$data
z
z$Samples <- factor(z$Samples, levels = c("Sucrose", "10x_kit", "Singulator"))
z$Clusters <- factor(z$Clusters, levels = c("Astrocyte", "Microglia", "Excitatory Neuron", "Inhibitory Neuron", "Oligo", "Other"))
z$Percent <- z$Proportions*100
z$short <- c("A", "A", "A", "M", "M", "M", "eN","eN","eN","iN","iN","iN","O","O","O","Ot","Ot","Ot")
z$short <- factor(z$short, levels = c("A", "M", "eN", "iN", "O", "Ot"))
write.csv(z, "Cell Type Proportions.csv")

ggplot(z, aes(x=Samples, y=Percent, fill = Samples))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  facet_grid(~short, scales = "free")+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("darkslategray3", "palevioletred", "darkgoldenrod1"))+
  #ggtitle("Proportion")+
  theme(legend.position = "top")+
  xlab("")+
  ylab("Percent (%)")+
  ggtitle("Proportion of Cell Types")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  #theme(axis.text.y = element_blank())+
  theme(legend.text = element_text(size = 8))+
  theme(axis.title.y = element_text(margin = margin(r=20)))+
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  NoLegend()

#ROGUE
                   
Idents(isolation) <- "tech"
summary(factor(isolation$tech))

iso.10x <- subset(isolation, idents = "10x_kit")
iso.sing <- subset(isolation, idents = "Singulator")
iso.suc <- subset(isolation, idents = "Sucrose")

Idents(iso.10x) <- "Cluster_Simple"
Idents(iso.sing) <- "Cluster_Simple"
Idents(iso.suc) <- "Cluster_Simple"

rouge.res.type.10x <- rogue(iso.10x@assays$RNA@counts, labels = iso.10x$Cluster_Simple, samples = iso.10x$Cluster_Simple, platform = "UMI")
rogue.boxplot(rouge.res.type.10x)

rouge.res.type.suc <- rogue(iso.suc@assays$RNA@counts, labels = iso.suc$Cluster_Simple, samples = iso.suc$Cluster_Simple, platform = "UMI")
rogue.boxplot(rouge.res.type.suc)

rouge.res.type.sing <- rogue(iso.sing@assays$RNA@counts, labels = iso.sing$Cluster_Simple, samples = iso.sing$Cluster_Simple, platform = "UMI")
rogue.boxplot(rouge.res.type.sing)

df1 <- data.frame(matrix(nrow = 18, ncol = 3))
colnames(df1) <- c("Cell_Type", "Rogue_Value", "Tech")
df1$Cell_Type <- c("Astrocyte", "Microglia", "Excitatory Neuron","Inhibitory Neuron", "Oligo","Other", "Astrocyte", "Microglia", "Excitatory Neuron","Inhibitory Neuron", "Oligo","Other","Astrocyte", "Microglia", "Excitatory Neuron","Inhibitory Neuron", "Oligo","Other")
df1$Tech <- c("Kit", "Kit","Kit","Kit","Kit","Kit", "Custom", "Custom", "Custom", "Custom", "Custom", "Custom", "Machine","Machine","Machine","Machine","Machine","Machine")
df1$Rogue_Value <- c(0.7459964, 0.6664026, 0.7806038, 0.7001894, 0.7933353, 0.762212, 0.7656841, 0.5855848, 0.5865414, 0.4733671, 0.687059, 0.7086275, 0.7202627, 0.7047332, 0.7367978, 0.6130273, 0.7095755, 0.7364299)

df1$Cell_Type <- factor(df1$Cell_Type, levels = c("Astrocyte", "Microglia", "Excitatory Neuron", "Inhibitory Neuron", "Oligo", "Other"))
df1$Tech <- factor(df1$Tech, levels = c("Custom", "Kit", "Machine"))

ggplot(df1, aes(x = Tech, y=Rogue_Value, fill= Tech))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  facet_grid(cols = vars(Cell_Type))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90))+
  scale_fill_manual(values = c("darkslategray3", 'palevioletred', "darkgoldenrod1"))+
  ylim(0, 1)+
  ylab("ROGUE Value")+
  ggtitle("ROGUE")+
  theme(axis.text.y = element_text(size=8), axis.title.y = element_text(size = 8))+
  #theme(axis.text.x = element_blank())+
  #theme(axis.tick.x = element_blank())+
  #theme(legend.text = element_text(size = 8))+
  #theme(legend.title = element_blank())+
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))+
  scale_x_discrete(labels=c("Astrocyte" = "A", "Microglia" = "M", "Excitatory Neuron" = "eN", "Inhibitory Neuron" = "iN", "Oligodendrocyte" = "O", "Other" = "O"))
  


#Figures were edited in Adobe illustrator to change style components, i.e., colors, titles 
