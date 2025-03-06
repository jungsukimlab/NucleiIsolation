set.seed(1221)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

library(Seurat)
library(SoupX)
library(DoubletFinder)
library(ggplot2)
library(gridExtra)
#load seurat object
isolation <- readRDS("20231206_snIsolation_collab.rds")

#Rho values
#Found in pre-processing steps with SoupX. refer to preprocess code
df2 <- data.frame(matrix(nrow = 6, ncol = 3))
colnames(df2) <- c("Sample", "Rho", "Tech")
df2$Sample <- c("Custom_1", "Custom_2", "Kit_1", "Kit_2", "Machine_1", "Machine_2")
df2$Rho <- c(0.174, 0.067, 0.523, 0.25, 0.136, 0.143)
df2$Tech <- c("Custom", "Custom", "Kit", "Kit", "Machine", "Machine")
#Manuscript figures were created in GraphPad but can visualize with ggplot2
ggplot(df2, aes(x = Sample, y = Rho, fill = Tech))+
  geom_bar(stat = "identity")+
  theme_bw()+
  facet_grid(cols = vars(Tech), scales = "free")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  #geom_hline(yintercept=0.9,linetype=2)+
  scale_fill_manual(values = c("darkslategray3", 'palevioletred', "darkgoldenrod1"))+
  #scale_shape_manual(values = c(15, 16, 17, 18,8,4))+
  ylab("Rho Value")+
  ggtitle("Contamination Fraction")+
  theme(axis.text = element_text(size=8), axis.title = element_text(size = 8))+
  theme(legend.text = element_text(size = 8))+
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

  #DoubletFinder 
  
