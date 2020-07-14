#Extract T-cell clusters
#Check cells in 'bad' clusters
SMC_clusters <- c("MYH11+ Smooth Muscle Cells")
SMC_cells <- WhichCells(object = v3.all.seur.combined, ident = SMC_clusters)

v3.all.seur.combined.SMC_clusters <- subset(v3.all.seur.combined, cells = SMC_cells)

#Make directory for SMC_subset_results
dir.create("SMC_subset_results", showWarnings = FALSE)

# t-SNE and Clustering
v3.all.seur.combined.SMC_clusters <- RunTSNE(v3.all.seur.combined.SMC_clusters, reduction = "pca.aligned", dims = 1:15, 
                                                  do.fast = T)
v3.all.seur.combined.SMC_clusters <- FindNeighbors(v3.all.seur.combined.SMC_clusters, reduction = "pca.aligned", dims = 1:15)
v3.all.seur.combined.SMC_clusters <- FindClusters(v3.all.seur.combined.SMC_clusters, resolution = 0.3)

TSNEPlot(v3.all.seur.combined.SMC_clusters, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                                     text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                     axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.SMC_clusters, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                     text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                     axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/tSNE patient.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.SMC_clusters, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                                        text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                        axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/tSNE sex.pdf", width = 15, height = 15)
TSNEPlot(v3.all.seur.combined.SMC_clusters, pt.size = 2, label.size = 10, group.by = "Plate") + theme(panel.background = element_blank(),
                                                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                   axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/tSNE Plate.pdf", width = 15, height = 15)

#Define marker genes per cluster
all.all.seur.combined.SMC_clusters.markers <- FindAllMarkers(object = v3.all.seur.combined.SMC_clusters, only.pos = TRUE, min.pct = 0.25, 
                                                           thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.SMC_clusters.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.SMC_clusters.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#And extract the gene names
c0.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==0),"gene"]))
c1.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==1),"gene"]))

#plot the top9 markers per cluster
VlnPlot(object = v3.all.seur.combined.SMC_clusters,
            features = c0.sep.markers
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/cluster0 markers violin.pdf", width = 15, height = 15)

VlnPlot(object = v3.all.seur.combined.SMC_clusters,
            features = c1.sep.markers
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/cluster1 markers violin.pdf", width = 15, height = 15)



VlnPlot(object = v3.all.seur.combined.SMC_clusters,
        features = c("ACTA2", "MYH11", "COL1A1","COL1A2"), nCol = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/some markers violin.pdf", width = 15, height = 15)

#Plot top markers in a heatmap
top10 <- all.all.seur.combined.SMC_clusters.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = v3.all.seur.combined.SMC_clusters, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + theme(panel.background = element_blank(),
                                                                                                                             text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                             axis.text        = element_text(size = 14, face = "plain")
)
ggsave("SMC_subset_results/clusters.pdf", width = 15, height =15)


#Write var genes per cluster to files
for (i in 0:1){
  x <- all.all.seur.combined.SMC_clusters.markers[which(all.all.seur.combined.SMC_clusters.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("SMC_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:1){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.SMC_clusters.markers[which(all.all.seur.combined.SMC_clusters.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "SMC_subset_results",full_GSEA = F)
}

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#Now let's look at some genes/clusters!
#First let's get the top 1000 genes for all clusters
dirname <- "SMC_subset_results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined.SMC_clusters)

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

