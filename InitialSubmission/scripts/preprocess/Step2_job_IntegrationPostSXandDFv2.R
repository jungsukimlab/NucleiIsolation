setwd("/N/project/kim_lab/Acri_Dom/Isolation_Exp/")
library(Seurat)
#v2's addition is to filter based on percent.mt early. 10% was chosen based on clustering done before initial CMG meeting. Note doublets have already been removed by SX and DF, so no additional QC is necessary
A=readRDS("20221031_10x_s1_postSXandDF.rds")
A=subset(A,percent.mt <= 10)
B=readRDS("20221102_10x_s2_postSXandDF.rds")
B=subset(B,percent.mt <= 10)
C=readRDS("20221103_S2_s1_postSXandDF.rds")
C=subset(C,percent.mt <= 10)
D=readRDS("20221103_S2_s2_postSXandDF.rds")
D=subset(D,percent.mt <= 10)
E=readRDS("20221103_Suc_s1_postSXandDF.rds")
E=subset(E,percent.mt <= 10)
X=readRDS("20221103_Suc_s2_postSXandDF.rds")
X=subset(X,percent.mt <= 10)

all.list <- list(A=A,B=B,C=C,D=D,E=E,X=X)
#need to double check that "sct" is being used appropriately.. check F:/SERCA1_6v10
rm(A)
rm(B)
rm(C)
rm(D)
rm(E)
rm(X)

#nuclei in all samples;integrated vs SCT
features = SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)#in v2 these var features won't include mt- genes. important distinction that affects clustering (DJA-20221112)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = features)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", anchor.features = features)
all.combined.sct = IntegrateData(anchorset = all.anchors, normalization.method = "SCT")
#can remove everything but all.combined.sct
all.combined.sct = RunPCA(all.combined.sct,verbose = F)
saveRDS(all.combined.sct,file = "20221207_all.combined.sct.rds")
