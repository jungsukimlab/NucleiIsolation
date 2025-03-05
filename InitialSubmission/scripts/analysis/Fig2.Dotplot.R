set.seed(1221)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/N/project/kim_lab/Kersey_Holly/IsolationCollab")

project <- readRDS("20231206_snIsolation_collab.rds")

genes_to_plot = c("Gad1", "Gad2", "Sv2b", "Slc17a7", "Mog", "Mag", "Slc1a3", "Gja1", "Inpp5d", "Cx3cr1", "Pdgfra")
 
exp_mat <- as.matrix(project[["RNA"]]@data[genes_to_plot,])
meta <- project@meta.data %>% 
  select(seurat_clusters, new_annotation)
 
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
meta <- pivot_longer(meta, cols = 3:13, names_to="Gene", values_to="Expression")
 
meta_summary <- meta %>%
group_by(new_annotation, seurat_clusters, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)
 
custom_order_gene <- rev(c("Gja1", 
                       "Slc1a3", 
                       "Cx3cr1", 
                       "Inpp5d", 
                       "Sv2b", 
                       "Slc17a7", 
                       "Gad1", 
                       "Gad2", 
                       "Mog", 
                       "Mag", "Pdgfra"))
 
meta_summary$Gene <- factor(meta_summary$Gene, levels = custom_order_gene)
 
sub_anno_level <- c("A", "M", "eN", "iN", "O", "Ot")
meta_summary$new_annotation <- factor(meta_summary$new_annotation, levels = sub_anno_level)
 
ggplot(meta_summary, aes(x=Gene, y=seurat_clusters)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_radius("% detected", range = c(0,3)) +
  scale_fill_gradient(low = "lightgray", high = "seagreen",
                       guide = guide_colorbar(ticks.colour = "black",  frame.colour = "black"),
                       name = "Average\nexpression") +
  ylab("Cluster") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=6, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title = element_text(size=8)) +
  facet_grid(~new_annotation, scales = "free", space = "free") +
  coord_flip() 
