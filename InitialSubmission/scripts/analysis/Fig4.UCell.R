set.seed(1221)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

library(UCell)
library(Seurat)
library(ggplot2)
library(gridExtra)
#load seurat object
isolation <- readRDS("20231206_snIsolation_collab.rds")

#subset object based on cell type
Idents(isolation) <- "simple_anno"
isolation.microglia <- subset(isolation, idents="Microglia")
isolation.astrocyte <- subset(isolation, idents="Astrocyte")
isolation.oligo <- subset(isolation, idents="Oligo")

#define marker genes for astrocytes               
markers.as <- list()
markers.as$HM_AS = c("Kcnj10","Glul","Slc1a2","Slc1a3","Slc6a11","Slc16a1","Ldha","Srebf1")
markers.as$DA_AS = c("Hmgb1","Hmgb3","Hmgb2","Id3","Gfap","S100b","Cd81","Cebpa","Ptprg")
#define marker genes for microglia
markers.mg <- list()
markers.mg$HM_MG = c("Tmem119","P2ry12","P2ry13","Cx3cr1", "C1qa", "Csf1r", "Hexb")
markers.mg$DA_MG= c("Apoe", "Cst7", "Trem2", "Itgax", "B2m", "Cst7")
#define marker genes for oligos
markers.og <- list()
markers.og$HM_OL = c("Olig1","Olig2","Mog","Mbp","Mobp","Plp1","Sox10","Gpr37","Mag","Cnp","Myrf")

#set identity to isolation method
Idents(isolation.microglia) <- "tech"
#Compute UCell scores
isolation.microglia <- AddModuleScore_UCell(isolation.microglia, features = markers.mg)
signature.names.mg <- paste0(names(markers.mg), "_UCell")
P <- VlnPlot(isolation.microglia, features = signature.names.mg, pt.size = FALSE)
write.csv(P[[1]][['data']], "ucell.HM_MG.csv") #save homeostatic microglia scores 
write.csv(P[[2]][['data']], "ucell.DA_MG.csv") #save disease-associated scores 


Idents(isolation.astrocyte) <- "tech"
isolation.astrocyte <- AddModuleScore_UCell(isolation.astrocyte, features = markers.as)
signature.names.as <- paste0(names(markers.as), "_UCell")
Q <- VlnPlot(isolation.astrocyte, features = signature.names.as, pt.size = FALSE)
write.csv(Q[[1]][['data']], "ucell.HM_AS.csv") #save homeostatic microglia scores
write.csv(Q[[2]][['data']], "ucell.DA_AS.csv") #save disease-associated microglia scores


Idents(isolation.oligo) <- "tech"
isolation.oligo <- AddModuleScore_UCell(isolation.oligo, features = markers.og)
signature.names.og <- paste0(names(markers.og), "_UCell")
R <- VlnPlot(isolation.oligo, features = signature.names.og, pt.size = FALSE)
write.csv(R[[1]][['data']], "ucell.HM_OG.csv") #save homeostatic oligo scores

#For visualization purposes 
p1 <- VlnPlot(isolation.oligo, features = "HM_OL_UCell", pt.size = FALSE)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")+
  ylim(0,0.7)+
  theme_classic()+NoLegend()

p2 <- VlnPlot(isolation.microglia, features = "DA_MG_UCell", pt.size = FALSE)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")+
  ylim(0,0.7)+
  theme_classic()+NoLegend()
p3 <- VlnPlot(isolation.microglia, features = "HM_MG_UCell", pt.size = FALSE)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")+
  ylim(0,0.7)+
  theme_classic()+NoLegend()

p4 <- VlnPlot(isolation.astrocyte, features = "DA_AS_UCell", pt.size = FALSE)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")+
  ylim(0,0.7)+
  theme_classic()+NoLegend()
p5 <- VlnPlot(isolation.astrocyte, features = "HM_AS_UCell", pt.size = FALSE)+
  stat_summary(fun.y = median, geom='point', size = 2, colour = "black")+
  ylim(0,0.7)+
  theme_classic()+NoLegend()

pdf(file = "Ucell.pdf")
grid.arrange(p5,p4,p3,p2,p1, nrow=1,ncol=5)
dev.off()
#Plots were edited in Adobe Illustrator for stylistic components, i.e., titles, colors... 
