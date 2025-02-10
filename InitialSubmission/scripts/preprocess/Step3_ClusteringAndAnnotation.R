setwd("/N/project/kim_lab/Acri_Dom/Isolation_Exp/analysis_2023/")
#fromLuke load all just in case, #out all packages not using
library(dplyr)
library(Seurat)#rgeos module must be loaded prior to RStudio on RED
library(patchwork)
library(SoupX)
library(cowplot)
library(multtest)
library(metap)
library(ggplot2)
library(paletteer)
library(scales)
#library(monocle3)
library(Matrix)
library(SeuratWrappers)
library(scMCA)
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)
#library(scran)
#library(scWGCNA)
library(sctransform)
library(DoubletFinder)

set.seed(46202)

####everything in this is now post mt- removal. For mt- issue please refer to clustering done in IsolationPreProcess_DFandSX.R after line 535. Last edited 20221111
all.combined.sct=readRDS("20221207_all.combined.sct.rds")
ElbowPlot(all.combined.sct, ndims = 50)
all.combined.sct = RunUMAP(all.combined.sct, reduction = "pca", dims = 3:46, verbose = F)#selected 46; removed first 2 PC's based on mt- bias
all.combined.sct = FindNeighbors(all.combined.sct, reduction = "pca", dims = 3:46)
#choose cluster resolution with ClustTree (https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects)
#library(clustree)
#or just set 3 different ones
all.combined.sct = FindClusters(all.combined.sct, resolution = c(0.4,0.8,1.2))
levels(all.combined.sct$integrated_snn_res.1.2)
#FIX: all.combined.sct <- SetAllIdent(seuratobject, id='0.4')
all.combined.sct<-SetIdent(all.combined.sct,value = 'integrated_snn_res.0.8')
DimPlot(all.combined.sct, reduction = "umap",label = T)
FeaturePlot(all.combined.sct, features = c("Olig1","Gfap","Tmem119"))
VlnPlot(all.combined.sct, features = "percent.mt",split.by = "tech")
saveRDS(all.combined.sct, file = "20221106_all.combined.sct_unannotated.rds")
#start point
all.combined.sct=readRDS("20221106_all.combined.sct_unannotated.rds")

##RES SELECTION AND ANNOTATION
Idents(all.combined.sct)<-"integrated_snn_res.1.2"
DimPlot(all.combined.sct, reduction = "umap",label = T)
library(clustree)
#prefix: integrated_snn_res
clustree(all.combined.sct, prefix = "integrated_snn_res.")
#if you name to x can investigate structure
#annotate 0.4, 0.8, 1.2; scMCA -> PangalaoDB
library(scMCA)
Idents(all.combined.sct)<-"integrated_snn_res.0.8"
DimPlot(all.combined.sct, reduction = "umap",label = T)
mca_result <- scMCA(scdata = all.combined.sct@assays$SCT@counts, numbers_plot = 3)#SCT counts, not integrated since it only did var features

#create way to get information of sc_MCA for a cluster, sort by score, take MODE of top 10% of scores; maybe manually look at length of summary of top 10% if there is conflicting
x=mca_result$scMCA_probility
saveRDS(x,file = "20221214_all.combined.sct_mca_result_scMCA_probility.rds")
#annotate
y=data.frame(all.combined.sct$integrated_snn_res.0.8)
#check clustering relationship
z=merge(x,y,by.x=1,by.y=0)
Idents(all.combined.sct)<-"integrated_snn_res.0.8"
DimPlot(all.combined.sct,reduction = "umap",label = T)
j=subset(z,all.combined.sct.integrated_snn_res.0.8 == "1")
head(j)
k=j
k=j[j$Score >= quantile(j$Score)[4],]
colnames(k)[2]="cell_type"
#need histogram of cell type
ggplot(k,aes(x=as.factor(cell_type)))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90))

for (i in c(1:length(levels(all.combined.sct$integrated_snn_res.0.8)))) {
  j=subset(z,all.combined.sct.integrated_snn_res.0.8 == levels(all.combined.sct$integrated_snn_res.0.8)[i])
  head(j)
  k=j[j$Score >= quantile(j$Score)[4],]
  colnames(k)[2]="cell_type"
  #need histogram of cell type
  pdf(paste0("./scMCA_annotation_20221214/","scMCAplot_clust",i-1,".pdf"))
  print(ggplot(k,aes(x=as.factor(cell_type)))+
          geom_bar()+
          theme(axis.text.x = element_text(angle = 90)))
  dev.off()
}
#now get cluster markers
Idents(all.combined.sct)<-"integrated_snn_res.0.8"
all.combined.sct.markers <- FindAllMarkers(all.combined.sct,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.combined.sct.markers,"20221214_all.combined.sct_res0.8_36cluster_markers.csv")



Idents(all.combined.sct)
all.combined.sct<-RenameIdents(all.combined.sct,#scMCA->followed by PanglaoDB refinement
                               `0`="N_0",#cell type; Rgs16 neuron -> Neuron (N)
                               `1`="O_1",#cell type; Oligo Lineage/Myelinating Oligo -> Oligo(O)
                               `2`="A_2",#cell type; Astrocyte (Me8 high, Atp1b2 high) -> Astrocyte(A)
                               `3`="EC_3",#cell type; Slc17a6/Rgs6 Neuron -> Endothelial Cell (EC)
                               `4`="N_4",#cell type; N -> N
                               `5`="N_5",#cell type; N -> N
                               `6`="N_6",#cell type; N -> N
                               `7`="N_7",#cell type; N -> N (mt/Rp high)
                               `8`="N_8",#cell type; N -> N
                               `9`="N_9",#cell type; N -> N (but higher in germ cells, and RGC)
                               `10`="N_10",#cell type; N -> N
                               `11`="N_11",#cell type; N -> N (some Oligo signature, some Fibroblast, unclear)
                               `12`="N_12",#cell type; N -> N
                               `13`="Mx_13",#cell type; N -> Unknown/Fibro/Macrophage/T (Mx?)
                               `14`="N_14",#cell type; N -> N (unclear though... why is it so far?)
                               `15`="M_15",#cell type; N -> microglia (M)
                               `16`="N_16",#cell type; N -> N
                               `17`="N_17",#cell type; N -> N
                               `18`="N_18",#cell type; N -> N
                               `19`="EC_19",#cell type; scMCA -> Endothelial cell close second to neuron (EC)
                               `20`="Mx_20",#cell type; N -> Unkown/Macropage/Germ/EC/Neuron (Mx)
                               `21`="N_21",#cell type; N -> N
                               `22`="Mx_22",#cell type; N -> Unknown/Macrophage/Neutrophil/Neuron/Ec (Mx)
                               `23`="O_23",#cell type; OPC -> OPC (O)
                               `24`="N_24",#cell type; N -> N
                               `25`="N_25",#cell type; N -> N
                               `26`="Mx_26",#cell type; N -> Macrophage/Oligo/T/Dendritic/Neuron (Mx)
                               `27`="N_27",#cell type; N -> N
                               `28`="N_28",#cell type; N -> N
                               `29`="N_29",#cell type; N -> N (but some Oligo, RGC, fibroblast)
                               `30`="EC_30",#cell type; scMCA -> Alot of unknown, endothelial cell second (EC)
                               `31`="N_31",#cell type; N -> N
                               `32`="N_32",#cell type; N -> N
                               `33`="EC_33",#cell type; scMCA -> No return, EC by consensus
                               `34`="M_34",#cell type; Microglia
                               `35`="O_35"#cell type; OPC -> Oligo
)
all.combined.sct$annotated <- Idents(all.combined.sct)
saveRDS(all.combined.sct, file = "20221214_all.combined.sct_annotated.rds")


###########LOAD IN##########################
all.combined.sct=readRDS(file = "20221214_all.combined.sct_annotated.rds")
Idents(all.combined.sct)<-"annotated"
DimPlot(all.combined.sct, reduction = "umap", label = T)
DimPlot(all.combined.sct, reduction = "umap", split.by = "tech",label = T)
all.combined.sct
#note these are still pretty confusing, Marker genes calculated by integrated slot, redo with RNA or SCT? Would that help some of these?

#simple annotation
all.combined.sct$simple_annotation="N"
all.combined.sct$simple_annotation[all.combined.sct$annotated == "A_2"]="A"
all.combined.sct$simple_annotation[all.combined.sct$annotated == "EC_3"|all.combined.sct$annotated == "EC_19"|all.combined.sct$annotated == "EC_30"|all.combined.sct$annotated == "EC_33"]="EC"
all.combined.sct$simple_annotation[all.combined.sct$annotated == "Mx_13"|all.combined.sct$annotated == "Mx_20"|all.combined.sct$annotated == "Mx_22"|all.combined.sct$annotated == "Mx_26"]="Mx"
all.combined.sct$simple_annotation[all.combined.sct$annotated == "M_15"|all.combined.sct$annotated == "M_34"]="M"
all.combined.sct$simple_annotation[all.combined.sct$annotated == "O_1"|all.combined.sct$annotated == "O_23"|all.combined.sct$annotated == "O_35"]="O"

Idents(all.combined.sct)<-"simple_annotation"
DimPlot(all.combined.sct, reduction = "umap", label = T)
DimPlot(all.combined.sct, reduction = "umap", split.by = "tech",label = T)
#do simple markers
all.combined.sct.simple.markers <- FindAllMarkers(all.combined.sct,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.combined.sct.simple.markers,"20230103_v2_SCT_all.combined.sct_simple_annotation.csv")



#redo marker genes on SCT slot (was previously done on integration markers)
DefaultAssay(all.combined.sct)<-"SCT"
all.combined.sct
#no need to redo scMCA because it was done on SCA already
Idents(all.combined.sct)<-"integrated_snn_res.0.8"
all.combined.sct.markers <- FindAllMarkers(all.combined.sct,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.combined.sct.markers,"20221214_v2_SCT_all.combined.sct_res0.8_36cluster_markers.csv")

#save additional Seurat files for collaborators (done 20230106-1500)
neuron<-subset(all.combined.sct, simple_annotation == "N")
saveRDS(neuron, file = "20230106_neuron.rds")
rm(neuron)
microglia<-subset(all.combined.sct, simple_annotation == "M")
saveRDS(microglia, file = "20230106_microglia.rds")
rm(microglia)
astrocyte<-subset(all.combined.sct, simple_annotation == "A")
saveRDS(astrocyte, file = "20230106_astrocyte.rds")
rm(astrocyte)
oligo<-subset(all.combined.sct, simple_annotation == "O")
saveRDS(oligo, file = "20230106_oligo.rds")
rm(oligo)
endo<-subset(all.combined.sct, simple_annotation == "EC")
saveRDS(endo, file = "20230106_endo.rds")
rm(endo)
mixed<-subset(all.combined.sct, simple_annotation == "Mx")
saveRDS(mixed, file = "20230106_mixed.rds")
rm(mixed)

#clean for new analysis
nuc_iso = all.combined.sct
nuc_iso$DF.classifications_0.25_0.19_4444 <- NULL
nuc_iso$DF.classifications_0.25_0.13_7412 <- NULL
nuc_iso$DF.classifications_0.25_0.07_1869 <- NULL
nuc_iso$DF.classifications_0.25_0.1_2921 <- NULL
nuc_iso$DF.classifications_0.25_0.02_1807 <- NULL
nuc_iso$DF.classifications_0.25_0.02_306<- NULL
nuc_iso$pANN_0.25_0.19_4444 <- NULL
nuc_iso$pANN_0.25_0.13_7412 <- NULL
nuc_iso$pANN_0.25_0.07_1869 <- NULL
nuc_iso$pANN_0.25_0.1_2921 <- NULL
nuc_iso$pANN_0.25_0.02_1807 <- NULL
nuc_iso$pANN_0.25_0.02_306 <- NULL
nuc_iso$seurat_clusters <- NULL
factor(nuc_iso$integrated_snn_res.0.8)#selected for resolution!!!
nuc_iso$integrated_snn_res.0.4 <- NULL
nuc_iso$integrated_snn_res.1.2<- NULL
saveRDS(nuc_iso,file = "20230321_nuc_iso_object.rds")
nuc_iso=readRDS("./20230321_nuc_iso_object.rds")
nuc_iso

Idents(nuc_iso)<-"tech"
DimPlot(nuc_iso,reduction = "umap",label = T,repel = T)#5x7.5
VizDimLoadings(nuc_iso,dims = 1:4)

factor(nuc_iso$tech)
FeaturePlot(nuc_iso,features = "percent.mt")#5x5
FeaturePlot(nuc_iso,features = "percent.mt",split.by = "tech")#5x15
Idents(nuc_iso)<-"sample"
DefaultAssay(nuc_iso)<-"SCT"
VlnPlot(nuc_iso, features = c("percent.mt"))


###########FUTURE?###############################

#resculter with dims 1:2 to show more drastic effects?
ElbowPlot(all.combined.sct, ndims = 50)
nuc_iso_mito = RunUMAP(nuc_iso, reduction = "pca", dims = 1:46, verbose = F)#selected 46; removed first 2 PC's based on mt- bias
nuc_iso_mito = FindNeighbors(nuc_iso_mito, reduction = "pca", dims = 1:46)

Ident(nuc_iso_mito)
DimPlot(nuc_iso_mito,reduction = "umap")
