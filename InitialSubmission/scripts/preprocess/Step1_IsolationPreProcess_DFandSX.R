setwd("/N/project/kim_lab/Acri_Dom/Isolation_Exp/")
#IsolationPreProcess_DFandSX
library(Seurat)
library(SoupX)
library(DoubletFinder)
#fromLuke load all just in case
library(dplyr)
library(patchwork)
library(cowplot)
library(multtest)
library(metap)
library(ggplot2)
library(paletteer)
library(scales)
library(Matrix)
library(SeuratWrappers)
library(fields)
library(KernSmooth)
library(ROCR)

set.seed(46202)

######################################
####10x_s1############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
    #rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/10x_1/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "10x_s1", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "A"
Seurat[["tech"]] <- "10x_kit"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221031_10x_s1_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: , 23502)|B(10x_s2:)
pK = 0.19
mult_rate=.1891
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.19_4444))
Seurat$DF=Seurat$DF.classifications_0.25_0.19_4444
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
saveRDS(Seurat2,file="20221031_10x_s1_postSXandDF.rds")

######################################
####10x_s2############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
#rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/10x_2/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:0.250, 30351)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "10x_s2", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "B"
Seurat[["tech"]] <- "10x_kit"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:0.13)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221102_10x_s2_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: 18.91, 23502)|B(10x_s2:24.42, 30351)
pK = 0.13
mult_rate=.2442
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.13_7412))
Seurat$DF=Seurat$DF.classifications_0.25_0.13_7412
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:24.42)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
#already messed up, saved singlet object only. would have to go back through from postSX and redo DF to get doublets from each..
saveRDS(Seurat2,file="20221102_10x_s2_postSXandDF.rds")

######################################
####S2_s1############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
#rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/S2_1/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:0.250, 30351)|C(S2_s1: 0.14, 15241)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "S2_s1", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "C"
Seurat[["tech"]] <- "Singulator"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:0.13)|C(S2_s2:0.07)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221102_S2_s1_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: 18.91, 23502)|B(10x_s2:24.42, 30351)|C(S2_s1:12.26,15241)
pK = 0.07
mult_rate=.1226
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.07_1869))
Seurat$DF=Seurat$DF.classifications_0.25_0.07_1869
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:24.42)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
#already messed up, saved singlet object only. would have to go back through from postSX and redo DF to get doublets from each..
saveRDS(Seurat2,file="20221103_S2_s1_postSXandDF.rds")


######################################
####S2_s2############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
#rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/S2_2/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:0.250, 30351)|C(S2_s1: 0.14, 15241)|D(S2_s2:0.14, 19053)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "S2_s2", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "D"
Seurat[["tech"]] <- "Singulator"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:0.13)|C(S2_s1:0.07)|D(S2_s2:0.1)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221102_S2_s2_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: 18.91, 23502)|B(10x_s2:24.42, 30351)|C(S2_s1:12.26,15241)|D(S2_s1:15.33,19053)
pK = 0.1
mult_rate=.1533
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.1_2921))
Seurat$DF=Seurat$DF.classifications_0.25_0.1_2921
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:24.42)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
#already messed up, saved singlet object only. would have to go back through from postSX and redo DF to get doublets from each..
saveRDS(Seurat2,file="20221103_S2_s2_postSXandDF.rds")

######################################
####Suc_s1############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
#rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/Suc_1/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:0.250, 30351)|C(S2_s1: 0.14, 15241)|D(S2_s2:0.14, 19053) |E(Suc_s1:0.17)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "Suc_s1", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "E"
Seurat[["tech"]] <- "Sucrose"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:0.13)|C(S2_s1:0.07)|D(S2_s2:0.1)|E(Suc_s1:0.02)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221102_S2_s2_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: 18.91, 23502)|B(10x_s2:24.42, 30351)|C(S2_s1:12.26,15241)|D(S2_s1:15.33,19053)|E(Suc_s1:12.06,14985)
pK = 0.02
mult_rate=.1206
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.02_1807))
Seurat$DF=Seurat$DF.classifications_0.25_0.02_1807
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:24.42)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
#already messed up, saved singlet object only. would have to go back through from postSX and redo DF to get doublets from each..
saveRDS(Seurat2,file="20221103_Suc_s1_postSXandDF.rds")

######################################
####Suc_s2############################
######################################

# Clean ambient RNA using SoupX
#SP1_219Soup<-load10X("./output/Suc_2/")
#SP1_219Soup=setClusters(SP1_219Soup,clusters = SP1_219Soup$metaData$nUMIs)
#rename ./analysis/clustering/*_graphclust/ to ./graphclust/
# Estimate rho
#SP1_219Soup = autoEstCont(SP1_219Soup)
#if too high, forceAccept=T
#SP1_219Soup = autoEstCont(SP1_219Soup,forceAccept = T)
# Estimated global rho of 0.074
# Clean the data
# 14079 cells
#SP1_219Soupout = adjustCounts(SP1_219Soup)
#SP1_219Seurat<-CreateSeuratObject(counts = SP1_219Soupout, project = "SP1_219", min.cells = 3, min.features = 200)
#rm(SP1_219Soup)
#rm(SP1_219Soupout)
#VlnPlot(SP1_219Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#SP1_219Seurat

set.seed(46202)

# Clean ambient RNA using SoupX
Soup<-load10X("./output/Suc_2/")
# Estimate rho
Soup = autoEstCont(Soup)
#if too high, forceAccept=T
Soup = autoEstCont(Soup,forceAccept = T)
# Clean the data
# Sample(name: rho, cells)|A(10x_s1: 0.523, 23502)|B(10x_s2:0.250, 30351)|C(S2_s1: 0.14, 15241)|D(S2_s2:0.14, 19053) |E(Suc_s1:0.17)|F(Suc_s2:0.07)
Soupout = adjustCounts(Soup)
Seurat<-CreateSeuratObject(counts = Soupout, project = "Suc_s2", min.cells = 3, min.features = 200)
rm(Soup)
rm(Soupout)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^mt-")
VlnPlot(Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Seurat[["sample"]] <- "F"
Seurat[["tech"]] <- "Sucrose"


Seurat <- subset(Seurat, subset = nFeature_RNA > 200)#took out mtDNA

Seurat <- SCTransform(Seurat)
#automatically sets default to SCT
Seurat <- RunPCA(Seurat)

Seurat <- RunUMAP(Seurat, dims = 1:50)

sweep.res.list <- paramSweep_v3(Seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]
# Sample(name: pK)|A(10x_s1:0.19)|B(10x_s2:0.13)|C(S2_s1:0.07)|D(S2_s2:0.1)|E(Suc_s1:0.02)|F(Suc_s2:0.02)
###############################################################################################################################################
# Use https://satijalab.org/costpercell/ to estimate multiplet rate and adjust round(0.075*nrow(TP323Seurat@meta.data)); make sure the multiplex rate is 1 not 8 to calculate
saveRDS(Seurat, file = "20221102_Suc_s2_postSX.rds")

# Sample(name: mult_rate, cells)|A(10x_s1: 18.91, 23502)|B(10x_s2:24.42, 30351)|C(S2_s1:12.26,15241)|D(S2_s1:15.33,19053)|E(Suc_s1:12.06,14985)|F(Suc_s2:4.96,6162)
pK = 0.02
mult_rate=.0496
homotypic.prop <- modelHomotypic(Seurat@meta.data$ClusteringResults)
nExp_poi <- round(mult_rate*nrow(Seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset; what we did with costpercell
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Seurat <- doubletFinder_v3(Seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#plot single vs doublet QC
summary(factor(Seurat$DF.classifications_0.25_0.02_306))
Seurat$DF=Seurat$DF.classifications_0.25_0.02_306
# Sample(name: doublet%, final singlet count)|A(10x_s1:18.9%, 19058)|B(10x_s2:24.42)|F(Suc_s2:)
VlnPlot(Seurat, features = c("nFeature_SCT","nCount_SCT","percent.mt"),group.by = "DF")
DimPlot(Seurat,reduction = "umap",group.by = "DF")

Seurat2 <- subset(Seurat, cells=rownames(Seurat@meta.data)[which(Seurat@meta.data$DF == "Singlet")])
#already messed up, saved singlet object only. would have to go back through from postSX and redo DF to get doublets from each..
saveRDS(Seurat2,file="20221103_Suc_s2_postSXandDF.rds")

#################Merge,FilterCluster
A=readRDS("20221031_10x_s1_postSXandDF.rds")
B=readRDS("20221102_10x_s2_postSXandDF.rds")
C=readRDS("20221103_S2_s1_postSXandDF.rds")
D=readRDS("20221103_S2_s2_postSXandDF.rds")
E=readRDS("20221103_Suc_s1_postSXandDF.rds")
F=readRDS("20221103_Suc_s2_postSXandDF.rds")
all.list <- list(A=A,B=B,C=C,D=D,E=E,F=F)
#need to double check that "sct" is being used appropriately.. check F:/SERCA1_6v10
rm(A)
rm(B)
rm(C)
rm(D)
rm(E)
rm(F)

#90535 nuclei in all samples;integrated vs SCT
features = SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = features)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", anchor.features = features)
all.combined.sct = IntegrateData(anchorset = all.anchors, normalization.method = "SCT")
#can remove everything but all.combined.sct
all.combined.sct = RunPCA(all.combined.sct,verbose = F)
#LINE513-535run as job_20221104.sh
all.combined.sct=readRDS("20221104_all.combined.sct.rds")
ElbowPlot(all.combined.sct, ndims = 50)
all.combined.sct = RunUMAP(all.combined.sct, reduction = "pca", dims = 1:46, verbose = F)#selected 46
all.combined.sct = FindNeighbors(all.combined.sct, reduction = "pca", dims = 1:46)
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
Idents(all.combined.sct)<-"integrated_snn_res.0.4"
mca_result <- scMCA(scdata = all.combined.sct@assays$SCT@counts, numbers_plot = 3)#SCT counts, not integrated since it only did var features

#create way to get information of sc_MCA for a cluster, sort by score, take MODE of top 10% of scores; maybe manually look at length of summary of top 10% if there is conflicting
x=mca_result$scMCA_probility
saveRDS(x,file = "20221108_all.combined.sct_mca_result_scMCA_probility.rds")
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
  pdf(paste0("./scMCA_annotation_20221108/","scMCAplot_clust",i,".pdf"))
  print(ggplot(k,aes(x=as.factor(cell_type)))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90)))
  dev.off()
}
#now get cluster markers
Idents(all.combined.sct)<-"integrated_snn_res.0.8"
all.combined.sct.markers <- FindAllMarkers(all.combined.sct,only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.combined.sct.markers,"20221109_all.combined.sct_res0.8_40cluster_markers.csv")



Idents(all.combined.sct)
all.combined.sct<-RenameIdents(all.combined.sct,#scMCA->followed by PanglaoDB refinement
                     `0`="N1",#cajal-retzius or interneuron; N1->N1
                     `1`="O1",#oligo/shwann; O1->O1
                     `2`="A1",#astro/bergmann/RadialGlialCell; A1->A1
                     `3`="U1",#really dirty; neuron/immature neuron; no change (sox5)/interneuron
                     `4`="N2",#interneuron
                     `5`="N3",#interneuron
                     `6`="N4",#neuroblast/interneuron
                     `7`="U2-mt",#precursor/astro (high mt-, high Rp-);N5->U2
                     `8`="N5",#interneuron or cajal
                     `9`="N6-mt",#HIGH mt, CAMK2n1 and Apoe
                     `10`="N7",#neuron
                     `11`="U3",#very unclear; need more research; olig/neuron/RGC
                     `12`="M1-mt",#menigeal,macrophage, high mt, mixed immune? U4 -> M1
                     `13`="N8", #neurons
                     `14`="M2",#microglia (siglech, plxdc2, Syk)
                     `15`="N9", #interneuron
                     `16`="N10", #interneuron
                     `17`="N11",# neuron
                     `18`="U4",# neuron DA: Tenm1 (~20th)
                     `19`="N12",# immature neurons
                     `20`="N13",#neural stem?/neuroblast/neuron
                     `21`="N14", #neuroblast/neuron
                     `22`="U5", #neuron (but other annotations)
                     `23`="N15",# immature neuron (others though)
                     `24`="N16",#neuron/neuroblast
                     `25`="O2-OPC", #OPC
                     `26`="N17",#glutamatergic
                     `27`="N18",#immature neurons
                     `28`="N19",#immatuire neuron/cajal
                     `29`="N20",#neurons/purkinje neruons
                     `30`="N21",#GABAergic
                     `31`="N22-mt",#neuron high mt
                     `32`="N23",#interneuron
                     `33`="U6",#very dirty
                     `34`="U7",#unclear
                     `35`="MG1",#mixed glia;
                     `36`="U8",#unclear
                     `37`="N24",#neuron(other annotation/oligo)
                     `38`="A2",#confirmed
                     `39`="A3"#confirmed

)
x=all.combined.sct.markers[all.combined.sct.markers$cluster=="2",]
head(x)
rownames(x)[1:10]
saveRDS(all.combined.sct, file = "20221110_all.combined.sct_annotated_mt-IN.rds")

DimPlot(all.combined.sct, reduction = "umap", split.by = "tech",label = T)
#need to look at:
#all.combined.sct<-RenameIdents(all.combined.sct,#followed by PanglaoDB refinement
#                     "U1"="M1"
#)
#all.combined.sct<-RenameIdents(all.combined.sct,#followed by PanglaoDB refinement
#                     "05"="O5")

DimPlot(all.combined.sct, reduction = "umap", label = T)
DimPlot(all.combined.sct, reduction = "umap", label = T,split.by = "tech")
#mito figures

#redo after filtering by percent.mt; start new IsolationProcess_afterPercentMT.R