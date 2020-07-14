#Extract T-cell clusters
#Check cells in 'bad' clusters
T_clusters <- c("CD3+CD8+ T-Cells I", "CD3+CD4+ T-Cells I",
                "CD3+CD8+ T-Cells II", "CD3+CD4+ T-Cells II")

# Update seurat object
all.seur.combined.s3 <- UpdateSeuratObject(all.seur.combined)
T_cells <- WhichCells(object = all.seur.combined.s3, ident = T_clusters)
length(T_cells)

all.seur.combined.T_clusters <- subset(all.seur.combined.s3, cells = T_cells)

#Split based on CD4 and CD8 expression
CD4CD8 <- GetAssayData(all.seur.combined.T_clusters[c("CD4","CD8A"),])
apply(CD4CD8,1,mean) #Means are comparable

#CD4-high cells
CD4.cells <- colnames(CD4CD8[,CD4CD8[1,] > CD4CD8[2,]])
#CD8-high cells
CD8.cells <- colnames(CD4CD8[,CD4CD8[1,] < CD4CD8[2,]])
sum(CD4CD8[1,CD4CD8[1,] == CD4CD8[2,]] != 0)

#Check out the CD4 cells first
all.seur.combined.T_split_clusters.CD4 <- subset(all.seur.combined.T_clusters, cells = CD4.cells)

#Check out the CD8 cells
all.seur.combined.T_split_clusters.CD8 <- subset(all.seur.combined.T_clusters, cells = CD8.cells)

length(WhichCells(all.seur.combined.T_split_clusters.CD8))
length(WhichCells(all.seur.combined.T_split_clusters.CD4))
#Make directory for T-cell_split_CD4_subset results
dir.create("T-cell_split_CD4_subset_results", showWarnings = FALSE)

# t-SNE and Clustering
all.seur.combined.T_split_clusters.CD4 <- RunTSNE(all.seur.combined.T_split_clusters.CD4, reduction = "pca.aligned", dims = 1:15, 
                                        do.fast = T)
all.seur.combined.T_split_clusters.CD4 <- FindNeighbors(all.seur.combined.T_split_clusters.CD4, reduction = "pca.aligned", dims = 1:15)
all.seur.combined.T_split_clusters.CD4 <- FindClusters(all.seur.combined.T_split_clusters.CD4, resolution = 1.2)

# Swap labels to align with v2 analysis
temp.idents <- as.character(Idents(all.seur.combined.T_split_clusters.CD4))
temp.idents[which(temp.idents == 4)] <- "N3"
temp.idents[which(temp.idents == 3)] <- "4"
temp.idents[which(temp.idents == "N3")] <- "3"
temp.idents <- as.factor(temp.idents)
Idents(all.seur.combined.T_split_clusters.CD4) <- temp.idents

TSNEPlot(all.seur.combined.T_split_clusters.CD4, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD4_subset_results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD4, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                   axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD4_subset_results/tSNE patient.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD4, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                               text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                               axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD4_subset_results/tSNE sex.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD4, pt.size = 2, label.size = 10, group.by = "Plate") + theme(panel.background = element_blank(),
                                                                                                 text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                 axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD4_subset_results/tSNE plate.pdf", width = 15, height = 15)

#Define marker genes per cluster
all.all.seur.combined.T_split_clusters.CD4.markers <- FindAllMarkers(object = all.seur.combined.T_split_clusters.CD4, only.pos = TRUE, min.pct = 0.25, 
                                                           thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.T_split_clusters.CD4.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.T_split_clusters.CD4.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#plot the top9 markers per cluster
for(i in 0:4){
  FeaturePlot(object = all.seur.combined.T_split_clusters.CD4,
              features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])) ,
              cols = c("grey", "blue"), 
              reduction = "tsne", pt.size = 2
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("T-cell_split_CD4_subset_results/cluster",i, " markers.pdf", sep = ""), width = 15, height = 15)
  VlnPlot(object = all.seur.combined.T_split_clusters.CD4,
          features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("T-cell_split_CD4_subset_results/cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}


#Plot top markers in a heatmap
top10 <- all.all.seur.combined.T_split_clusters.CD4.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = all.seur.combined.T_split_clusters.CD4, genes.use = top10$gene, slim.col.label = TRUE, remove.key = T) + theme(panel.background = element_blank(),
                                                                                                                        text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                        axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD4_subset_results/clusters.pdf", width = 15, height =15)

#Different view of markers
markers.to.plot <- c("CD4", "CD3E", "CD8A",
                     "CD28","TIGIT","GNLY",
                     "SELL", "IL7R")

DotPlot(all.seur.combined.T_split_clusters.CD4, features = rev(markers.to.plot), 
                      cols = c("dodgerblue", "coral", "mediumorchid4"), dot.scale = 8) + theme(panel.background = element_blank(),
                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                   axis.text        = element_text(size = 14, face = "plain")
                      )
ggsave("T-cell_split_CD4_subset_results/patient dot.r07.pdf", width = 10, height = 20)


#Write var genes per cluster to files
for (i in 0:4){
  x <- all.all.seur.combined.T_split_clusters.CD4.markers[which(all.all.seur.combined.T_split_clusters.CD4.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("T-cell_split_CD4_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:4){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.T_split_clusters.CD4.markers[which(all.all.seur.combined.T_split_clusters.CD4.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "T-cell_split_CD4_subset_results", full_GSEA = F)
}

cluster.averages <- AverageExpression(all.seur.combined.T_split_clusters.CD4)

#-------------------------------------------------
#Make heatmap of a selection of genes
sel.genes <- c("TCF7", "SELL", "LEF1", "CCR7", "IL7R",
               "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4", 
               "IL2", "GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7", "IL32", "IL10",
               "CD28", "TNFRSF14", "ICOS", "TNFRSF9", "CD40LG",
               "EOMES", "HOPX", "TBX21", "ZEB2", "ZNF683", "HIF1A", "ID2", "TOX", "ID3", "RORA", "RORC", "GATA3", "PRDM1", "RUNX3", "ZBTB7B",
               "IL2RA", "FOXP3", "IKZF2",
               "CCL5", "CCL4", "CCR5", "CCR6", "CXCR1", "CD69"
)
length(sel.genes)

#Keep those in expressed in CD4 clusters
sel.genes <- sel.genes[sel.genes %in% row.names(cluster.averages)]
length(sel.genes)
sel.genes[!sel.genes %in% row.names(cluster.averages)]

#Subset the expression data
sel.data <- cluster.averages[sel.genes,]
nrow(sel.data)

#Keep only those genes expressed in more than 0 cells
sel.data[apply(sel.data,1,sum) == 0,] <- 0.0001
nrow(sel.data)

#Make the annotations
T.ann.row <- data.frame(Group=c(rep("Naïve Markers", 5), 
                                rep("Inhibitory receptors", 5), 
                                rep("Cytokine and effector molecules", 10), 
                                rep("Co-stimulatory molecules", 5),
                                rep("Transcription factors", 15),
                                rep("Treg markers", 3),
                                rep("Chemokine receptors and ligands", 6)),
                        row.names = row.names(sel.data)
)
T.ann.row$Group <- as.factor(T.ann.row$Group)
levels(T.ann.row$Group)
T.ann.col <- data.frame(Cluster=c("CD4.0","CD4.1","CD4.2","CD4.3","CD4.4"), row.names = colnames(sel.data))

T.ann.colors <- list(
  Group = c("Naïve Markers"="darkturquoise",
            "Inhibitory receptors"="dodgerblue2",
            "Cytokine and effector molecules"="coral1",
            "Co-stimulatory molecules"="darkgoldenrod1",
            "Transcription factors"="blueviolet",
            "Treg markers"="aquamarine4",
            "Chemokine receptors and ligands"="darkkhaki"),
  Cluster = c("CD4.0"="darkslategray1",
              "CD4.1"="deepskyblue",
              "CD4.2"="dodgerblue",
              "CD4.3"="dodgerblue4",
              "CD4.4"="darkblue")
)

#Draw the plot
pheatmap(sel.data, 
         scale = "row", 
         cluster_cols = F, 
         cluster_rows = F,
         show_colnames = F,
         annotation_row = T.ann.row,
         annotation_col = T.ann.col,
         annotation_colors = T.ann.colors,
         gaps_row = c(5, 10, 19, 24, 39, 42),
         cellheight = 10,
         cellwidth = 10,
         filename = "T-cell_split_CD4_subset_results/select_genes.heatmap.pdf"
)

sel.data.4 <- sel.data

#Add T cells from cluster M.4
sel.data.M <- MT_cells[sel.genes,,drop = F]
colnames(sel.data.M) <- "M.4_CD4"

#Merge and clean up
sel.data.M4 <- merge(sel.data,sel.data.M, by=0, sort = F, all.x = T)
row.names(sel.data.M4) <- sel.data.M4$Row.names
sel.data.M4$Row.names <- NULL
colnames(sel.data.M4) <- c("CD4.0", "CD4.1", "CD4.2", "CD4.3", "CD4.4", "CD4.M.4")

#Make the annotations
T.ann.row <- data.frame(Group=c(rep("Naïve Markers", 5), 
                                rep("Inhibitory receptors", 5), 
                                rep("Cytokine and effector molecules", 10), 
                                rep("Co-stimulatory molecules", 5),
                                rep("Transcription factors", 15),
                                rep("Treg markers", 3),
                                rep("Chemokine receptors and ligands", 6)),
                        row.names = row.names(sel.data)
)
T.ann.row$Group <- as.factor(T.ann.row$Group)
levels(T.ann.row$Group)
T.ann.col <- data.frame(Cluster=c("CD4.0","CD4.1","CD4.2","CD4.3","CD4.4", "CD4.M.4"), row.names = colnames(sel.data.M4))

T.ann.colors <- list(
  Group = c("Naïve Markers"="darkturquoise",
            "Inhibitory receptors"="dodgerblue2",
            "Cytokine and effector molecules"="coral1",
            "Co-stimulatory molecules"="darkgoldenrod1",
            "Transcription factors"="blueviolet",
            "Treg markers"="aquamarine4",
            "Chemokine receptors and ligands"="darkkhaki"),
  Cluster = c("CD4.0"="darkslategray1",
              "CD4.1"="deepskyblue",
              "CD4.2"="dodgerblue",
              "CD4.3"="dodgerblue3",
              "CD4.4"="dodgerblue4",
              "CD4.M.4"="purple")
)

#Draw the plot
pheatmap(sel.data.M4, 
         scale = "row", 
         cluster_cols = F, 
         cluster_rows = F,
         show_colnames = F,
         annotation_row = T.ann.row,
         annotation_col = T.ann.col,
         annotation_colors = T.ann.colors,
         gaps_row = c(5, 10, 19, 24, 39, 42),
         cellheight = 10,
         cellwidth = 10,
         filename = "T-cell_split_CD4_subset_results/select_genes_including_M.4.heatmap.pdf"
)

#------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
#Now check out the CD8 cells
#Make directory for T-cell_split_CD8_subset results
dir.create("T-cell_split_CD8_subset_results", showWarnings = FALSE)

# t-SNE and Clustering
all.seur.combined.T_split_clusters.CD8 <- RunTSNE(all.seur.combined.T_split_clusters.CD8, reduction = "pca.aligned", dims = 1:15, 
                                                  do.fast = T)
all.seur.combined.T_split_clusters.CD8 <- FindNeighbors(all.seur.combined.T_split_clusters.CD8, reduction = "pca.aligned", dims = 1:15)
all.seur.combined.T_split_clusters.CD8 <- FindClusters(all.seur.combined.T_split_clusters.CD8, resolution = 0.4)

TSNEPlot(all.seur.combined.T_split_clusters.CD8, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                                     text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                     axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD8_subset_results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD8, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                             text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                             axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD8_subset_results/tSNE patient.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD8, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                                         text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                         axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD8_subset_results/tSNE sex.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.T_split_clusters.CD8, pt.size = 2, label.size = 10, group.by = "Plate") + theme(panel.background = element_blank(),
                                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD8_subset_results/tSNE plate.pdf", width = 15, height = 15)

#Define marker genes per cluster
all.all.seur.combined.T_split_clusters.CD8.markers <- FindAllMarkers(object = all.seur.combined.T_split_clusters.CD8, only.pos = TRUE, min.pct = 0.25, 
                                                                     thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.T_split_clusters.CD8.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.T_split_clusters.CD8.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#plot the top9 markers per cluster
for(i in 0:4){
  FeaturePlot(object = all.seur.combined.T_split_clusters.CD8,
              features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])) ,
              cols.use = c("grey", "blue"), 
              reduction.use = "tsne", pt.size = 2
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("T-cell_split_CD8_subset_results/cluster",i, " markers.pdf", sep = ""), width = 15, height = 15)
  VlnPlot(object = all.seur.combined.T_split_clusters.CD8,
          features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("T-cell_split_CD8_subset_results/cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}


#Plot top markers in a heatmap
top10 <- all.all.seur.combined.T_split_clusters.CD8.markers %>% group_by(cluster) %>% top_n(10, -p_val_adj)
DoHeatmap(object = all.seur.combined.T_split_clusters.CD8, features = top10$gene) + theme(panel.background = element_blank(),
                                                                                        text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                        axis.text        = element_text(size = 14, face = "plain")
)
ggsave("T-cell_split_CD8_subset_results/clusters.pdf", width = 15, height =15)

#Different view of markers
markers.to.plot <- c("CD4", "CD3E", "CD8A",
                     "CD28","TIGIT","GNLY",
                     "SELL", "IL7R")

sdp <- SplitDotPlotGG(all.seur.combined.T_split_clusters.CD8, genes.plot = rev(markers.to.plot), 
                      cols.use = c("dodgerblue", "coral", "mediumorchid4"), 
                      x.lab.rot = T, plot.legend = T, dot.scale = 8, 
                      do.return = T, grouping.var = "Sex") + theme(panel.background = element_blank(),
                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                   axis.text        = element_text(size = 14, face = "plain")
                      )
ggsave("T-cell_split_CD8_subset_results/patient dot.pdf", width = 10, height = 20)


#Write var genes per cluster to files
for (i in 0:2){
  x <- all.all.seur.combined.T_split_clusters.CD8.markers[which(all.all.seur.combined.T_split_clusters.CD8.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("T-cell_split_CD8_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:2){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.T_split_clusters.CD8.markers[which(all.all.seur.combined.T_split_clusters.CD8.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "T-cell_split_CD8_subset_results", full_GSEA = F)
}

#-------------------------------------------------
#Make heatmap of a selection of genes
cluster.averages <- AverageExpression(all.seur.combined.T_split_clusters.CD8)
sel.genes <- c("TCF7", "SELL", "LEF1", "CCR7", "IL7R",
               "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4", 
               "IL2", "GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7", "IL32", "IL10",
               "CD28", "TNFRSF14", "ICOS", "TNFRSF9", "CD40LG",
               "EOMES", "HOPX", "TBX21", "ZEB2", "ZNF683", "HIF1A", "ID2", "TOX", "ID3", "RORA", "RORC", "GATA3", "PRDM1", "RUNX3", "ZBTB7B",
               "IL2RA", "FOXP3", "IKZF2",
               "CCL5", "CCL4", "CCR5", "CCR6", "CXCR1", "CD69"
)
length(sel.genes)

#Keep those in expressed in CD8 clusters
sel.genes <- sel.genes[sel.genes %in% row.names(cluster.averages)]
length(sel.genes)
sel.genes[!sel.genes %in% row.names(cluster.averages)]

#Subset the expression data
sel.data <- cluster.averages[sel.genes,]
nrow(sel.data)

#Keep only those genes expressed in more than 0 cells
sel.data[apply(sel.data,1,sum) == 0,] <- 0.0001
nrow(sel.data)

#Make the annotations
T.ann.row <- data.frame(Group=c(rep("Naïve Markers", 5), 
                                rep("Inhibitory receptors", 5), 
                                rep("Cytokine and effector molecules", 10), 
                                rep("Co-stimulatory molecules", 5),
                                rep("Transcription factors", 15),
                                rep("Treg markers", 3),
                                rep("Chemokine receptors and ligands", 6)),
                        row.names = row.names(sel.data)
)
T.ann.row$Group <- as.factor(T.ann.row$Group)
levels(T.ann.row$Group)
T.ann.col <- data.frame(Cluster=c("CD8.0","CD8.1","CD8.2","CD8.3"), row.names = colnames(sel.data))

T.ann.colors <- list(
  Group = c("Naïve Markers"="darkturquoise",
            "Inhibitory receptors"="dodgerblue2",
            "Cytokine and effector molecules"="coral1",
            "Co-stimulatory molecules"="darkgoldenrod1",
            "Transcription factors"="blueviolet",
            "Treg markers"="aquamarine4",
            "Chemokine receptors and ligands"="darkkhaki"),
  Cluster = c("CD8.0"="orangered",
              "CD8.1"="firebrick1",
              "CD8.2"="firebrick3",
              "CD8.3"="firebrick4"
              )
)

#Draw the plot
pheatmap(sel.data, 
         scale = "row", 
         cluster_cols = F, 
         cluster_rows = F,
         show_colnames = F,
         annotation_row = T.ann.row,
         annotation_col = T.ann.col,
         annotation_colors = T.ann.colors,
         gaps_row = c(5, 10, 19, 24, 39, 42),
         cellheight = 10,
         cellwidth = 10,
         filename = "T-cell_split_CD8_subset_results/select_genes.heatmap.pdf"
)


#Match with CD4 expression
colnames(sel.data.4) <- c("CD4.0","CD4.1","CD4.2", "CD4.3","CD4.4")
colnames(sel.data) <- c("CD8.0", "CD8.1","CD8.2","CD8.3")
dim(sel.data)
dim(sel.data.4)

#Merge and clean up
sel.data.84 <- merge(sel.data,sel.data.4, by=0, sort = F, all.x = T)
row.names(sel.data.84) <- sel.data.84$Row.names
sel.data.84$Row.names <- NULL

#Make the annotations
T.ann.row <- data.frame(Group=c(rep("Naïve Markers", 5), 
                                rep("Inhibitory receptors", 5), 
                                rep("Cytokine and effector molecules", 10), 
                                rep("Co-stimulatory molecules", 5),
                                rep("Transcription factors", 15),
                                rep("Treg markers", 3),
                                rep("Chemokine receptors and ligands", 6)),
                        row.names = row.names(sel.data)
)
T.ann.row$Group <- as.factor(T.ann.row$Group)
levels(T.ann.row$Group)
T.ann.col <- data.frame(Cluster=c("CD8.0","CD8.1","CD8.2","CD8.3","CD4.0","CD4.1","CD4.2","CD4.3","CD4.4"), row.names = colnames(sel.data.84))

T.ann.colors <- list(
  Group = c("Naïve Markers"="darkturquoise",
            "Inhibitory receptors"="dodgerblue2",
            "Cytokine and effector molecules"="coral1",
            "Co-stimulatory molecules"="darkgoldenrod1",
            "Transcription factors"="blueviolet",
            "Treg markers"="aquamarine4",
            "Chemokine receptors and ligands"="darkkhaki"),
  Cluster = c("CD8.0"="orangered",
              "CD8.1"="firebrick1",
              "CD8.2"="firebrick3",
              "CD8.3"="firebrick4",
              "CD4.0"="darkslategray1",
              "CD4.1"="deepskyblue",
              "CD4.2"="dodgerblue",
              "CD4.3"="dodgerblue3",
              "CD4.4"="dodgerblue4")
)

#Draw the plot
pheatmap(sel.data.84, 
         scale = "row", 
         cluster_cols = F, 
         cluster_rows = F,
         show_colnames = F,
         annotation_col = T.ann.col,
         annotation_row = T.ann.row,
         annotation_colors = T.ann.colors,
         gaps_row = c(5, 10, 20, 25, 40, 43),
         gaps_col = 4,
         cellheight = 10,
         cellwidth = 15,
         filename = "results/CD4_CD8_select_genes.heatmap.pdf"
)



#-----------------------------------------------------------------------------------------------
#Let's get the top 1000 genes for all clusters
dirname <- "T-cell_split_CD8_subset_results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined.T_split_clusters.CD8)

#take the top 1000 from each dataframe and put them in a list for easier handling
df.list <- list()
for(i in 1:ncol(cluster.averages)){
  df <- cluster.averages[order(cluster.averages[,i], decreasing = T),][1:1000,]
  df.list[[i]] <- df
}
names(df.list) <- paste("CD8", colnames(cluster.averages), sep = ".")

#Iterate over the list and:
for (i in 1:length(df.list)){
  #Write a table with the top 1000 genes and their expression
  write.better.table(df.list[[i]], file = paste(dirname, "/top1000_genes_from_", names(df.list[i]), ".expression.txt", sep = ""))
}

#-----------------------------------------------------------------------------------------------
#Let's get the top 1000 genes for all clusters
dirname <- "T-cell_split_CD4_subset_results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined.T_split_clusters.CD4)

#take the top 1000 from each dataframe and put them in a list for easier handling
df.list <- list()
for(i in 1:ncol(cluster.averages)){
  df <- cluster.averages[order(cluster.averages[,i], decreasing = T),][1:1000,]
  df.list[[i]] <- df
}
names(df.list) <- paste("CD4", colnames(cluster.averages), sep = ".")

#Iterate over the list and:
for (i in 1:length(df.list)){
  #Write a table with the top 1000 genes and their expression
  write.better.table(df.list[[i]], file = paste(dirname, "/top1000_genes_from_", names(df.list[i]), ".expression.txt", sep = ""))
}



#---------------------------------------------------------------------------------------------------
# Dotplot CD4 markers

CD4.markers <- c("TBX21",
                  "GATA3",
                  "RORA",
                  "FOXP3",
                  "CD27",
                  "CD28",
                  "RUNX3",
                  "CD69",
                  "SELL",
                  "IL7R",
                  "LEF1",
                  "GZMB",
                  "PRF1",
                  "GZMA",
                  "GZMK")

DotPlot(seurat3.all.seur.combined.T_split_clusters.CD4, features = CD4.markers, cols = c("dodgerblue","firebrick1")) +
  theme(axis.text.x = element_text(angle = 90),
        aspect.ratio = 3/1)
ggsave("T-cell_split_CD4_subset_results/CD4_final_markers.dotplot.pdf")
