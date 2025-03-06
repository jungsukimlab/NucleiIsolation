set.seed(1221)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

library(Seurat)
library(ggplot2)
library(gridExtra)
#load seurat object
isolation <- readRDS("20231206_snIsolation_collab.rds")

Idents(isolation) <- "tech"
summary(factor(isolation$tech))
iso.10x <- subset(isolation, idents = "10x_kit")
iso.sing <- subset(isolation, idents = "Singulator")
iso.suc <- subset(isolation, idents = "Sucrose")
#Set idents to cell type
Idents(iso.10x) <- "Cluster_Simple"
Idents(iso.sing) <- "Cluster_Simple"
Idents(iso.suc) <- "Cluster_Simple"

#Proportion of Mitochondrial Reads by Cell Type - Kit 
#Mean
a <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Astrocyte"])
b <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Microglia"])
c <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Excitatory Neuron"])
d <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Inhibitory Neuron"])
e <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Oligo"])
f <- mean(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Other"])
#SD
aa <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Astrocyte"])
bb <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Microglia"])
cc <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Excitatory Neuron"])
dd <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Inhibitory Neuron"])
ee <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Oligo"])
ff <- sd(iso.10x$percent.mt[iso.10x$Cluster_Simple=="Other"])
mean.kit <- data.frame(matrix(nrow = 6))
mean.kit$mean <- c(a, b, c, d, e, f)
mean.kit$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
mean.kit$sd <-c(aa, bb, cc, dd, ee, ff)
write.csv(mean.kit, "mt reads kit.csv")

#Proportion of Mitochondrial Reads by Cell Type - Machine
#Mean
g <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Astrocyte"])
h <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Microglia"])
i <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Excitatory Neuron"])
j <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Inhibitory Neuron"])
k <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Oligo"])
l <- mean(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Other"])
#SD
gg <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Astrocyte"])
hh <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Microglia"])
ii <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Excitatory Neuron"])
jj <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Inhibitory Neuron"])
kk <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Oligo"])
ll <- sd(iso.sing$percent.mt[iso.sing$Cluster_Simple=="Other"])
mean.sing <- data.frame(matrix(nrow = 6))
mean.sing$sd <- c(gg,hh,ii,jj,kk,ll)
mean.sing$mean <- c(g,h,i,j,k,l)
mean.sing$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
write.csv(mean.sing, "mt mean machine.csv")

#Proportion of Mitochondrial Reads by Cell Type - Sucrose Centrifugation
#Mean
m <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Astrocyte"])
n <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Microglia"])
o <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Excitatory Neuron"])
p <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Inhibitory Neuron"])
q <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Oligo"])
r <- mean(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Other"])
#SD
mm <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Astrocyte"])
nn <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Microglia"])
oo <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Excitatory Neuron"])
pp <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Inhibitory Neuron"])
qq <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Oligo"])
rr <- sd(iso.suc$percent.mt[iso.suc$Cluster_Simple=="Other"])
mean.suc <- data.frame(matrix(nrow = 6))
mean.suc$mean <- c(m,n,o,p,q,r)
mean.suc$sd <- c(mm,nn,oo,pp,qq,rr)
mean.suc$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
write.csv(mean.suc, "mt mean custom.csv")

#Data was imported into GraphPad for manuscript figures, but can visualize with ggplot2

#Ribosomal Reads
ribo_genes <- rownames(isolation)[grep("^Rp[sl]", rownames(isolation))]
isolation$percent.ribo <- colSums(isolation@assays$RNA@counts[ribo_genes, ])/isolation$nCount_RNA

#Proportion of ribosomal Reads by Cell Type - Kit
#Mean
a.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Astrocyte"])
b.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Microglia"])
c.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Excitatory Neuron"])
d.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Inhibitory Neuron"])
e.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Oligo"])
f.r <- mean(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Other"])
#SD
aa.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Astrocyte"])
bb.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Microglia"])
cc.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Excitatory Neuron"])
dd.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Inhibitory Neuron"])
ee.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Oligo"])
ff.r <- sd(iso.10x$percent.ribo[iso.10x$Cluster_Simple=="Other"])
mean.kit.ribo <- data.frame(matrix(nrow = 6))
mean.kit.ribo$mean <- c(a.r, b.r, c.r, d.r, e.r, f.r)
mean.kit.ribo$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
mean.kit.ribo$sd <-c(aa.r, bb.r, cc.r, dd.r, ee.r, ff.r)
write.csv(mean.kit.ribo, "rb reads kit.csv")

#Proportion of ribosomal Reads by Cell Type - Machine
#Mean
g.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Astrocyte"])
h.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Microglia"])
i.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Excitatory Neuron"])
j.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Inhibitory Neuron"])
k.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Oligo"])
l.r <- mean(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Other"])
#SD
gg.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Astrocyte"])
hh.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Microglia"])
ii.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Excitatory Neuron"])
jj.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Inhibitory Neuron"])
kk.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Oligo"])
ll.r <- sd(iso.sing$percent.ribo[iso.sing$Cluster_Simple=="Other"])
mean.sing.ribo <- data.frame(matrix(nrow=6))
mean.sing.ribo$mean <- c(g.r,h.r,i.r,j.r,k.r,l.r)
mean.sing.ribo$sd <- c(gg.r,hh.r,ii.r,jj.r,kk.r,ll.r)
mean.sing.ribo$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
write.csv(mean.sing.ribo, "rb reads sing.csv")

#Proportion of ribosomal Reads by Cell Type - Sucrose Centrifugation
#Mean
m.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Astrocyte"])
n.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Microglia"])
o.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Excitatory Neuron"])
p.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Inhibitory Neuron"])
q.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Oligo"])
r.r <- mean(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Other"])
#SD
mm.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Astrocyte"])
nn.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Microglia"])
oo.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Excitatory Neuron"])
pp.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Inhibitory Neuron"])
qq.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Oligo"])
rr.r <- sd(iso.suc$percent.ribo[iso.suc$Cluster_Simple=="Other"])
mean.suc.ribo <- data.frame(matrix(nrow=6))
mean.suc.ribo$mean <- c(m.r,n.r,o.r,p.r,q.r,r.r)
mean.suc.ribo$sd <- c(mm.r,nn.r,oo.r,pp.r,qq.r,rr.r)
mean.suc.ribo$celltype <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Ihibitory Neuron", "Oligo","Other")
write.csv(mean.suc.ribo, "rb reads suc.csv")

#Data was imported into GraphPad for manuscript figures, but can visualize with ggplot2
