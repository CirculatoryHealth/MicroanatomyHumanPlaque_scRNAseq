#Extract endo cell clusters
#Check cells in 'bad' clusters
E_clusters <- c("CD34+ Endothelial Cells", "CD34+ Endothelial Cells II")
E_cells <- WhichCells(object = v3.all.seur.combined, ident = E_clusters)

v3.all.seur.combined.E_clusters <- subset(v3.all.seur.combined, cells = E_cells)

#Make directory for Endothelial_cell_subset_results
dir.create("Endothelial_cell_subset_results", showWarnings = FALSE)

# t-SNE and Clustering
v3.all.seur.combined.E_clusters <- RunTSNE(v3.all.seur.combined.E_clusters, reduction = "pca.aligned", dims = 1:15, 
                                             do.fast = T)
v3.all.seur.combined.E_clusters <- FindNeighbors(v3.all.seur.combined.E_clusters, reduction = "pca.aligned", dims = 1:15)
v3.all.seur.combined.E_clusters <- FindClusters(v3.all.seur.combined.E_clusters, resolution = 0.5)

TSNEPlot(v3.all.seur.combined.E_clusters, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Endothelial_cell_subset_results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.E_clusters, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                   axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Endothelial_cell_subset_results/tSNE patient.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.E_clusters, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                                      text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                      axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Endothelial_cell_subset_results/tSNE sex.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.E_clusters, pt.size = 2, label.size = 10, group.by = "Plate") + theme(panel.background = element_blank(),
                                                                                                 text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                 axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Endothelial_cell_subset_results/tSNE plate.pdf", width = 15, height = 15)

#Define marker genes per cluster
all.all.seur.combined.E_clusters.markers <- FindAllMarkers(object = v3.all.seur.combined.E_clusters, only.pos = TRUE, min.pct = 0.25, 
                                                           thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.E_clusters.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.E_clusters.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#And extract the gene names
c0.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==0),"gene"]))
c1.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==1),"gene"]))
c2.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==2),"gene"]))
c3.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==3),"gene"]))

#plot the top9 markers per cluster
VlnPlot(object = v3.all.seur.combined.E_clusters,
            features = c0.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
)
ggsave("Endothelial_cell_subset_results/cluster0 markers.pdf", width = 15, height = 15)

VlnPlot(object = v3.all.seur.combined.E_clusters,
            features = c1.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
)
ggsave("Endothelial_cell_subset_results/cluster1 markers.pdf", width = 15, height = 15)

VlnPlot(object = v3.all.seur.combined.E_clusters,
            features = c2.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
)
ggsave("Endothelial_cell_subset_results/cluster2 markers.pdf", width = 15, height = 15)

VlnPlot(object = v3.all.seur.combined.E_clusters,
            features = c3.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
)
ggsave("Endothelial_cell_subset_results/cluster3 markers.pdf", width = 15, height = 15)

#Plot top markers in a heatmap
top10 <- all.all.seur.combined.E_clusters.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = v3.all.seur.combined.E_clusters, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + theme(panel.background = element_blank(),
                                                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Endothelial_cell_subset_results/clusters.pdf", width = 15, height =15)

#Replot the top 10 markers with pheatmap
#Subset the expression data
endo.EMT.data <- as.matrix(all.seur.combined.E_clusters@data)[top10$gene,]

#Keep only those genes expressed in more than 0 cells
endo.EMT.data <- endo.EMT.data[apply(endo.EMT.data,1,sum) != 0,]

#Order on cluster
endo.EMT.data <- endo.EMT.data[, names(sort(all.seur.combined.E_clusters@ident))]

#Draw the plot
pheatmap(endo.EMT.data, 
         scale = "row", 
         cluster_cols = F, 
         cluster_rows = T,
         show_colnames = F, 
         gaps_col = endo.cluster.lengths,
         color = colorRampPalette(c("blue","goldenrod","red"))(ncol(endo.EMT.data)),
         breaks = seq(-2,2, length.out = ncol(endo.EMT.data)), 
         cellheight = 10,
         filename = "Endothelial_cell_subset_results/top10_markers.heatmap.pdf"
)

#Write var genes per cluster to files
for (i in 0:3){
  x <- all.all.seur.combined.E_clusters.markers[which(all.all.seur.combined.E_clusters.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Endothelial_cell_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:3){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.E_clusters.markers[which(all.all.seur.combined.E_clusters.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "Endothelial_cell_subset_results", full_GSEA = F)
}

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#Now let's look at some genes/clusters!
#First let's get the top 1000 genes for all clusters
dirname <- "Endothelial_cell_subset_results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined.E_clusters)

#take the top 1000 from each dataframe and put them in a list for easier handling
df.list <- list()
for(i in 1:ncol(cluster.averages)){
  df <- cluster.averages[order(cluster.averages[,i], decreasing = T),][1:1000,]
  df.list[[i]] <- df
}
names(df.list) <- colnames(cluster.averages)

#Iterate over the list and:
for (i in 1:length(df.list)){
  #Plot kmaens clusters and z scores
  pheatmap(df.list[[i]][,1:ncol(cluster.averages)], scale = "row", show_rownames = F, filename = paste(dirname, "/top1000_genes_from_", names(df.list[i]), ".zscores.pdf", sep = ""))
  pheatmap(df.list[[i]][,1:ncol(cluster.averages)], scale = "row", show_rownames = T, kmeans_k = 6, filename = paste(dirname, "/top1000_genes_from_", names(df.list[i]), ".kmeans.pdf", sep = ""))
  
  #Write a table with the top 1000 genes, their expression, and kmeans cluster
  d <- pheatmap(df.list[[i]],scale = "row", kmeans_k = 6, silent = T)
  d <- data.frame(kmeans.cluster = d$kmeans$cluster)
  df.list[[i]] <- cbind(df.list[[i]],d)
  write.better.table(df.list[[i]], file = paste(dirname, "/top1000_genes_from_", names(df.list[i]), ".expressions_and_kmeans_cluster.txt", sep = ""))
  
  #Get pathways per clusters
  x <- data.frame(gene = rownames(df.list[[i]]), avg_logFC = 0)
  get_ontology(x, name = paste("Cluster", names(df.list[i]), "Pathways for top 1000 genes", sep = " "), return.data = F, outdir = dirname, full_GSEA = F)
  
  #Get pathways for all kmeans clusters (for each cluster)
  numclus <- sort(unique(df.list[[i]]$kmeans.cluster))
  subdirname <- paste(dirname, "/top1000_genes_from_", names(df.list[i]), "_pathways_from_kmeans_clusters", sep = "")
  dir.create(subdirname, showWarnings = FALSE)
  
  for (j in numclus){
    x <- data.frame(gene = rownames(df.list[[i]][which(df.list[[i]]$kmeans.cluster == j),]), avg_logFC = 0)
    if(nrow(x)>10){
      get_ontology(x, name = paste("Cluster", names(df.list[i]), "Pathways_for_kmeans_cluster", j, sep = "_"), return.data = F, outdir = subdirname, full_GSEA = F)
    }
  }
}


#==================================================================================================
# Look at some IPA results
# Load data for upstream regulators (UR)
E.0123.UR <- read.table(file = "IPA/E0123 UR.txt", header = T, row.names = 1)
dim(E.0123.UR)
colnames(E.0123.UR) <- c("E.3", "E.2", "E.1", "E.0")
E.0123.UR[grep("IFNG", row.names(E.0123.UR)),]

# Order
E.0.UR <- E.0123.UR[order(E.0123.UR$E.0, decreasing = T),]
E.1.UR <- E.0123.UR[order(E.0123.UR$E.1, decreasing = T),]
E.2.UR <- E.0123.UR[order(E.0123.UR$E.2, decreasing = T),]
E.3.UR <- E.0123.UR[order(E.0123.UR$E.3, decreasing = T),]

# Get the top 100 each and compile a unique list
E.UR.top100 <- E.0123.UR[unique(c(row.names(E.0.UR[1:100,]),
                                      row.names(E.1.UR[1:100,]),
                                      row.names(E.2.UR[1:100,]),
                                      row.names(E.3.UR[1:100,])
)
)
,]

# Extract k-means clusters
set.seed(1)
pheatmap(E.UR.top100,
         scale        = "row", 
         cellwidth    = 30, 
         cluster_rows = F, 
         cluster_cols = F, kmeans_k = 4
)
set.seed(1)
E.0.UR.k <- E.UR.top100[pheatmap(E.UR.top100,
                                     scale        = "row", 
                                     cellwidth    = 30, 
                                     cluster_rows = F, 
                                     cluster_cols = F,
                                     kmeans_k     = 4
)$kmeans$cluster == 1,]
set.seed(1)
E.1.UR.k <- E.UR.top100[pheatmap(E.UR.top100,
                                     scale        = "row", 
                                     cellwidth    = 30, 
                                     cluster_rows = F, 
                                     cluster_cols = F,
                                     kmeans_k     = 4
)$kmeans$cluster == 3,]
set.seed(1)
E.2.UR.k <- E.UR.top100[pheatmap(E.UR.top100,
                                     scale        = "row", 
                                     cellwidth    = 30, 
                                     cluster_rows = F, 
                                     cluster_cols = F,
                                     kmeans_k     = 4
)$kmeans$cluster == 2,]
set.seed(1)
E.3.UR.k <- E.UR.top100[pheatmap(E.UR.top100,
                                     scale        = "row", 
                                     cellwidth    = 30, 
                                     cluster_rows = F, 
                                     cluster_cols = F,
                                     kmeans_k     = 4
)$kmeans$cluster == 4,]

# Order
E.0.UR.k <- E.0.UR.k[order(E.0.UR.k$E.0, decreasing = T),]
E.1.UR.k <- E.1.UR.k[order(E.1.UR.k$E.1, decreasing = T),]
E.2.UR.k <- E.2.UR.k[order(E.2.UR.k$E.2, decreasing = T),]
E.3.UR.k <- E.3.UR.k[order(E.3.UR.k$E.3, decreasing = T),]

# Get the top 10 each and compile a unique list
E.UR.k.top10 <- E.0123.UR[unique(c(row.names(E.0.UR.k[1:10,]),
                                       row.names(E.1.UR.k[1:10,]),
                                       row.names(E.2.UR.k[1:10,]),
                                       row.names(E.3.UR.k[1:10,])
))
,]

E.UR.k.top10 <- E.UR.k.top10[!is.na(E.UR.k.top10[,1]),]
dim(E.UR.k.top10)
write.table(as.data.frame(row.names(E.UR.k.top10))[1:10,, drop = F], file = "IPA/E.0.UR.txt", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(row.names(E.UR.k.top10))[11:20,, drop = F], file = "IPA/E.1.UR.txt", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(row.names(E.UR.k.top10))[21:30,, drop = F], file = "IPA/E.2.UR.txt", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(row.names(E.UR.k.top10))[31:40,, drop = F], file = "IPA/E.3.UR.txt", quote = F, row.names = F, col.names = F)

# Plot
pheatmap(E.UR.k.top10,
         scale        = "row", 
         cellwidth    = 30, 
         cluster_rows = F, 
         cluster_cols = F,
         gaps_row     = c(10,20,30),
         filename     = "IPA/E.01234.UR.top10.heatmap.pdf"
)

# Check expression of the URs in the expression of the subsets
#Update the v2 seaurat object
all.seur.combined.E_clusters.V3 <- UpdateSeuratObject(all.seur.combined.E_clusters)

# E.0 URs
# Get the genes
gene.list <- row.names(E.0.UR.k[1:10,])

# Fix gene names
gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)
gene.list  <- gene.list[gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)]
long.gene.list <- gene.list

# Plot
VlnPlot(all.seur.combined.E_clusters.V3, 
        features = gene.list, 
        ncol     = 3
)
ggsave("IPA/E.0.UR.violin.pdf")

# E.1 URs
# Get the genes
gene.list <- row.names(E.1.UR.k[1:10,])

# Fix gene names
gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)
gene.list  <- gene.list[gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)]
long.gene.list <- c(long.gene.list, gene.list)

# Plot
VlnPlot(all.seur.combined.E_clusters.V3, 
        features = gene.list, 
        ncol     = 3
)
ggsave("IPA/E.1.UR.violin.pdf")

# E.2 URs
# Get the genes
gene.list <- row.names(E.2.UR.k[1:10,])

# Fix gene names
gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)
gene.list  <- gene.list[gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)]
long.gene.list <- c(long.gene.list, gene.list)

# Plot
VlnPlot(all.seur.combined.E_clusters.V3, 
        features = gene.list, 
        ncol     = 3
)
ggsave("IPA/E.2.UR.violin.pdf")

# E.3 URs
# Get the genes
gene.list <- row.names(E.3.UR.k[1:10,])

# Fix gene names
gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)
gene.list  <- gene.list[gene.list %in% row.names(all.seur.combined.E_clusters.V3$RNA)]
long.gene.list <- c(long.gene.list, gene.list)

# Plot
VlnPlot(all.seur.combined.E_clusters.V3, 
        features = gene.list, 
        ncol     = 3
)
ggsave("IPA/E.3.UR.violin.pdf")


dev.off()
DoHeatmap(all.seur.combined.E_clusters.V3, features = long.gene.list)


# Load data for canonical pathways (CP)
E.0123.CP <- read.table(file = "IPA/E0123 CP.txt", header = T, row.names = 1)
dim(E.0123.CP)
colnames(E.0123.CP) <- c("E.3", "E.2", "E.1", "E.0")

# Order
E.0.CP <- E.0123.CP[order(E.0123.CP$E.0, decreasing = T),]
E.1.CP <- E.0123.CP[order(E.0123.CP$E.1, decreasing = T),]
E.2.CP <- E.0123.CP[order(E.0123.CP$E.2, decreasing = T),]
E.3.CP <- E.0123.CP[order(E.0123.CP$E.3, decreasing = T),]

# Get the top 100 each and compile a unique list
E.CP.top100 <- E.0123.CP[unique(c(row.names(E.0.CP[1:100,]),
                                  row.names(E.1.CP[1:100,]),
                                  row.names(E.2.CP[1:100,]),
                                  row.names(E.3.CP[1:100,])
)
)
,]

# Extract k-means clusters
set.seed(1)
pheatmap(E.CP.top100,
         scale        = "row", 
         cellwidth    = 30, 
         cluster_rows = F, 
         cluster_cols = F, kmeans_k = 4
)
set.seed(1)
E.0.CP.k <- E.CP.top100[pheatmap(E.CP.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 4
)$kmeans$cluster == 1,]
set.seed(1)
E.1.CP.k <- E.CP.top100[pheatmap(E.CP.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 4
)$kmeans$cluster == 3,]
set.seed(1)
E.2.CP.k <- E.CP.top100[pheatmap(E.CP.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 4
)$kmeans$cluster == 4,]
set.seed(1)
E.3.CP.k <- E.CP.top100[pheatmap(E.CP.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 4
)$kmeans$cluster == 2,]

# Order
E.0.CP.k <- E.0.CP.k[order(E.0.CP.k$E.0, decreasing = T),]
E.1.CP.k <- E.1.CP.k[order(E.1.CP.k$E.1, decreasing = T),]
E.2.CP.k <- E.2.CP.k[order(E.2.CP.k$E.2, decreasing = T),]
E.3.CP.k <- E.3.CP.k[order(E.3.CP.k$E.3, decreasing = T),]

# Get the top 10 each and compile a unique list
E.CP.k.top10 <- E.0123.CP[unique(c(row.names(E.0.CP.k[1:10,]),
                                   row.names(E.1.CP.k[1:10,]),
                                   row.names(E.2.CP.k[1:10,]),
                                   row.names(E.3.CP.k[1:10,])
))
,]

E.CP.k.top10 <- E.CP.k.top10[!is.na(E.CP.k.top10[,1]),]

# Plot
pheatmap(E.CP.k.top10,
         scale        = "row", 
         cellwidth    = 30, 
         labels_row   = gsub("_", " ", row.names(E.CP.k.top10)),
         cluster_rows = F, 
         cluster_cols = F,
         gaps_row     = c(10,20,30),
         filename     = "IPA/E.01234.CP.top10.heatmap.pdf"
)
