#Extract myeloid clusters
#Check cells in 'bad' clusters
M_clusters <- c("CD14+CD68+ Macrophages I", "CD14+CD68+ Macrophages II","CD14+CD68+ Macrophages III", "CD14+CD68+ Macrophages IV")
M_cells <- WhichCells(object = all.seur.combined, ident = M_clusters)

all.seur.combined.M_clusters <- SubsetData(object = all.seur.combined, ident.use = M_clusters)

#Make directory for Myeloid_cell_subset_results
dir.create("Myeloid_cell_subset_results", showWarnings = FALSE)

# t-SNE and Clustering
all.seur.combined.M_clusters <- RunTSNE(all.seur.combined.M_clusters, reduction.use = "pca.aligned", dims.use = 1:15, 
                                        do.fast = T)
all.seur.combined.M_clusters <- FindClusters(all.seur.combined.M_clusters, reduction.type = "pca.aligned", 
                                             resolution = 1.2, dims.use = 1:15, force.recalc = T)

TSNEPlot(all.seur.combined.M_clusters, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.M_clusters, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                   text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                   axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/tSNE patient.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.M_clusters, pt.size = 2, label.size = 10, group.by = "Plate") + theme(panel.background = element_blank(),
                                                                                                      text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                      axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/tSNE plate.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined.M_clusters, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                                 text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                 axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/tSNE sex.pdf", width = 15, height = 15)

# Update to Seurat v3
v3.all.seur.combined.M_clusters <- UpdateSeuratObject(all.seur.combined.M_clusters)

# Update ident labels
Mye.idents <- paste0("Mye.", v3.all.seur.combined.M_clusters@active.ident)
Idents(v3.all.seur.combined.M_clusters) <- Mye.idents

#Update ident order
Idents(v3.all.seur.combined.M_clusters) <- factor(x =Idents(v3.all.seur.combined.M_clusters), levels = c("Mye.0", "Mye.1", "Mye.2", "Mye.3", "Mye.4"))

#Define marker genes per cluster
all.all.seur.combined.M_clusters.markers <- FindAllMarkers(object = all.seur.combined.M_clusters, only.pos = TRUE, min.pct = 0.25, 
                                                           thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.M_clusters.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.M_clusters.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#And extract the gene names
c0.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==0),"gene"]))
c1.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==1),"gene"]))
c2.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==2),"gene"]))
c3.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==3),"gene"]))
c4.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==4),"gene"]))

#plot the top9 markers per cluster
FeaturePlot(object = all.seur.combined.M_clusters,
            features.plot = c0.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/cluster0 markers.pdf", width = 15, height = 15)

FeaturePlot(object = all.seur.combined.M_clusters,
            features.plot = c1.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/cluster1 markers.pdf", width = 15, height = 15)

FeaturePlot(object = all.seur.combined.M_clusters,
            features.plot = c2.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/cluster2 markers.pdf", width = 15, height = 15)

FeaturePlot(object = all.seur.combined.M_clusters,
            features.plot = c3.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/cluster3 markers.pdf", width = 15, height = 15)

FeaturePlot(object = all.seur.combined.M_clusters,
            features.plot = c4.sep.markers ,
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size = 2
) + theme(panel.background = element_blank(),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
         axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/cluster4 markers.pdf", width = 15, height = 15)


#Plot top markers in a heatmap
top10 <- all.all.seur.combined.M_clusters.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = all.seur.combined.M_clusters, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + theme(panel.background = element_blank(),
                                                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/clusters.pdf", width = 15, height =15)

#Write var genes per cluster to files
for (i in 0:4){
  x <- all.all.seur.combined.M_clusters.markers[which(all.all.seur.combined.M_clusters.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Myeloid_cell_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:4){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.M_clusters.markers[which(all.all.seur.combined.M_clusters.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "Myeloid_cell_subset_results", full_GSEA = F)
}

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#Now let's look at some genes/clusters!
#First let's get the top 1000 genes for all clusters
dirname <- "Myeloid_cell_subset_results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined.M_clusters)

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


#Extract "T-cells" cluster 4
all.seur.combined.MT_cluster <- SubsetData(all.seur.combined.M_clusters, ident.use = 4)
all.seur.combined.MT_cluster@data["CD3E",]
all.seur.combined.MT_T.only <- SubsetData(all.seur.combined.M_clusters, cells.use = colnames(all.seur.combined.MT_cluster@data)[!(all.seur.combined.MT_cluster@data["CD3E",] == 0 & all.seur.combined.MT_cluster@data["CD4",] == 0)])
MT_cells <- AverageExpression(all.seur.combined.MT_T.only)[,'4', drop = F]
MT_cells["CD3E",]
MT_cells["CD4",]

#Filter out 'unique' genes from the top 1000 lists. I.e., only keep those that are *not* in the other sub clusters
#Extract the top 1000 dfs
My0.top1000 <- df.list$`0`
My1.top1000 <- df.list$`1`
My2.top1000 <- df.list$`2`

#Keep only those genes that are unique to each
My0.top1000.unique <- My0.top1000[which(!row.names(My0.top1000)[1:1000] %in% c(row.names(My1.top1000)[1:1000], row.names(My2.top1000)[1:1000])), 1:3]
My1.top1000.unique <- My1.top1000[which(!row.names(My1.top1000)[1:1000] %in% c(row.names(My0.top1000)[1:1000], row.names(My2.top1000)[1:1000])), 1:3]
My2.top1000.unique <- My2.top1000[which(!row.names(My2.top1000)[1:1000] %in% c(row.names(My0.top1000)[1:1000], row.names(My1.top1000)[1:1000])), 1:3]

#Write the files
write.better.table(My0.top1000.unique, file = "Myeloid_cell_subset_results/My0.top1000.unique.genes.txt")
write.better.table(My1.top1000.unique, file = "Myeloid_cell_subset_results/My1.top1000.unique.genes.txt")
write.better.table(My2.top1000.unique, file = "Myeloid_cell_subset_results/My2.top1000.unique.genes.txt")

#Visualize with RPKM heatmaps
pheatmap(mat = My0.top1000.unique, 
         scale = "row",
         cluster_cols = T,
         show_rownames = F,
         cellwidth = 20,
         cellheight = 1,
         cutree_rows = 3,
         border_color = "black",
         color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(1000),
         main = "My0 uniquely expressed genes from top 1000",
         filename = "Myeloid_cell_subset_results/My0 uniquely expressed genes from top 1000.pdf"
)

pheatmap(mat = My1.top1000.unique, 
         scale = "row",
         cluster_cols = T,
         show_rownames = F,
         cellwidth = 20,
         cellheight = 1,
         cutree_rows = 3,
         border_color = "black",
         color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(1000),
         main = "My1 uniquely expressed genes from top 1000",
         filename = "Myeloid_cell_subset_results/My1 uniquely expressed genes from top 1000.pdf"
)

pheatmap(mat = My2.top1000.unique, 
         scale = "row",
         cluster_cols = T,
         show_rownames = F,
         cellwidth = 20,
         cellheight = 1,
         cutree_rows = 3,
         border_color = "black",
         color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(1000),
         main = "My2 uniquely expressed genes from top 1000",
         filename = "Myeloid_cell_subset_results/My2 uniquely expressed genes from top 1000.pdf"
)

My.top20.unique <- rbind(My0.top1000.unique[1:20,],My1.top1000.unique[1:20,],My2.top1000.unique[1:20,])
pheatmap(mat = My.top20.unique, 
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,
         cellwidth = 20,
         cellheight = 10,
         gaps_row = c(20,40),
         border_color = "black",
         color = colorRampPalette(c("#FF00FF","#000000","#FFFF00"))(1000),
         main = "Top 20 uniquely expressed genes per subset",
         filename = "Myeloid_cell_subset_results/Top 20 uniquely expressed genes per subset from top 1000.pdf"
)


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
# Look at some IPA results
# Load data for upstream regulators (UR)
My.012.UR <- read.table(file = "IPA/M012 UR.txt", header = T, row.names = 1)
dim(My.012.UR)
colnames(My.012.UR) <- c("My.2", "My.1", "My.0")
My.012.UR <- My.012.UR[,c(3,2,1)]
My.012.UR[grep("IFNG", row.names(My.012.UR)),]

# Order
My.M0.UR <- My.012.UR[order(My.012.UR$My.0, decreasing = T),]
My.M1.UR <- My.012.UR[order(My.012.UR$My.1, decreasing = T),]
My.M2.UR <- My.012.UR[order(My.012.UR$My.2, decreasing = T),]

# Get the top 100 each and compile a unique list
My.UR.top100 <- My.012.UR[unique(c(row.names(My.M0.UR[1:100,]),row.names(My.M1.UR[1:100,]),row.names(My.M2.UR[1:100,]))),]

# Extract k-means clusters
set.seed(1)
pheatmap(My.UR.top100,
         scale        = "row", 
         cellwidth    = 30, 
         cluster_rows = F, 
         cluster_cols = F, kmeans_k = 3
         )
set.seed(1)
My.M0.UR.k <- My.UR.top100[pheatmap(My.UR.top100,
                           scale        = "row", 
                           cellwidth    = 30, 
                           cluster_rows = F, 
                           cluster_cols = F,
                           kmeans_k     = 3
)$kmeans$cluster == 2,]
set.seed(1)
My.M1.UR.k <- My.UR.top100[pheatmap(My.UR.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 3
)$kmeans$cluster == 3,]
set.seed(1)
My.M2.UR.k <- My.UR.top100[pheatmap(My.UR.top100,
                                 scale        = "row", 
                                 cellwidth    = 30, 
                                 cluster_rows = F, 
                                 cluster_cols = F,
                                 kmeans_k     = 3
)$kmeans$cluster == 1,]

# Order
My.M0.UR.k <- My.M0.UR.k[order(My.M0.UR.k$My.0, decreasing = T),]
My.M1.UR.k <- My.M1.UR.k[order(My.M1.UR.k$My.1, decreasing = T),]
My.M2.UR.k <- My.M2.UR.k[order(My.M2.UR.k$My.2, decreasing = T),]

# Get the top 10 each and compile a unique list
My.UR.k.top10 <- My.012.UR[unique(c(row.names(My.M0.UR.k[1:10,]),row.names(My.M1.UR.k[1:10,]),row.names(My.M2.UR.k[1:10,]))),]
dim(My.UR.k.top10)
write.table(as.data.frame(row.names(My.UR.k.top10))[1:10,, drop = F], file = "IPA/My.0.UR.txt", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(row.names(My.UR.k.top10))[11:20,, drop = F], file = "IPA/My.1.UR.txt", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(row.names(My.UR.k.top10))[21:30,, drop = F], file = "IPA/My.2.UR.txt", quote = F, row.names = F, col.names = F)


# Plot
pheatmap(My.UR.k.top10,
         scale        = "row", 
         cellwidth    = 30, 
         cluster_rows = F, 
         cluster_cols = F,
         gaps_row     = c(10,20),
         filename     = "IPA/My.012.UR.top10.heatmap.pdf"
)

# Check expression of the URs in the expression of the subsets
#Update the v2 seaurat object
all.seur.combined.M_clusters.V3 <- UpdateSeuratObject(all.seur.combined.M_clusters)

# My.0 URs
# Get the genes
gene.list <- row.names(My.M0.UR.k[1:10,])

# Fix gene names
row.names(My.M0.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)
gene.list  <- gene.list[row.names(My.M0.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)]
long.gene.list <- gene.list

# Plot
VlnPlot(all.seur.combined.M_clusters.V3, 
        features = gene.list, 
        idents   = c(0,1,2),
        ncol     = 3
        )
ggsave("IPA/My.0.UR.violin.pdf")

# My.1 URs
# Get the genes
gene.list <- row.names(My.M1.UR.k[1:10,])

# Fix gene names
row.names(My.M1.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)
gene.list  <- gene.list[row.names(My.M1.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)]
long.gene.list <- c(long.gene.list, gene.list)

# Plot
VlnPlot(all.seur.combined.M_clusters.V3, 
        features = gene.list, 
        idents   = c(0,1,2),
        ncol     = 3
)
ggsave("IPA/My.1.UR.violin.pdf")

# My.2 URs
# Get the genes
gene.list <- row.names(My.M2.UR.k[1:10,])

# Fix gene names
row.names(My.M2.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)
gene.list  <- gene.list[row.names(My.M2.UR.k[1:10,]) %in% row.names(all.seur.combined.M_clusters.V3$RNA)]
long.gene.list <- c(long.gene.list, gene.list)

# Plot
VlnPlot(all.seur.combined.M_clusters.V3, 
        features = gene.list, 
        idents   = c(0,1,2),
        ncol     = 3
)
ggsave("IPA/My.2.UR.violin.pdf")

dev.off()
DoHeatmap(all.seur.combined.M_clusters.V3, features = long.gene.list, cells = WhichCells(all.seur.combined.M_clusters.V3,idents = c(0,1,2)))

#----------------
# Try looking at the abs top 20 per cluster instead of going the K-means route
My.UR.top20 <- My.012.UR[unique(c(row.names(My.M0.UR[1:20,]),row.names(My.M1.UR[1:20,]),row.names(My.M2.UR[1:20,]))),]
pheatmap(My.UR.top20, 
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         cellwidth = 20, 
         cellheight = 12,
         filename = "IPA/My.012.UR.abs_top20.pdf")

# Do the same for the pathways
# Load data for canonical pathways (CP)
My.012.CP <- read.table(file = "IPA/M012 CP.txt", header = T, row.names = 1)
dim(My.012.CP)
colnames(My.012.CP) <- c("My.2", "My.1", "My.0")
My.012.CP <- My.012.CP[,c(3,2,1)]

# Order
My.M0.CP <- My.012.CP[order(My.012.CP$My.0, decreasing = T),]
My.M1.CP <- My.012.CP[order(My.012.CP$My.1, decreasing = T),]
My.M2.CP <- My.012.CP[order(My.012.CP$My.2, decreasing = T),]

# Get the top 20 each and compile a unique list
My.CP.top20 <- My.012.CP[unique(c(row.names(My.M0.CP[1:20,]),row.names(My.M1.CP[1:20,]),row.names(My.M2.CP[1:20,]))),]

pheatmap(My.CP.top20, 
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         cellwidth = 20, 
         cellheight = 12)

# Try with our own pathway data
# Extract pvals
x <- all.all.seur.combined.M_clusters.markers[which(all.all.seur.combined.M_clusters.markers$cluster == 0),-6]
My.0.PW <- get_ontology(x, return.data = T, full_GSEA = F)

x <- all.all.seur.combined.M_clusters.markers[which(all.all.seur.combined.M_clusters.markers$cluster == 1),-6]
My.1.PW <- get_ontology(x, return.data = T, full_GSEA = F)

x <- all.all.seur.combined.M_clusters.markers[which(all.all.seur.combined.M_clusters.markers$cluster == 2),-6]
My.2.PW <- get_ontology(x, return.data = T, full_GSEA = F)

# Merge into one big table
My.CP <- merge(My.0.PW$Pathways, My.1.PW$Pathways, by = "name", all = T)
My.CP <- merge(My.CP, My.2.PW$Pathways, by = "name", all = T)
colnames(My.CP) <- c("name", "My.0", "My.1", "My.2")

# Remove NA's (turn into 1)
My.CP[is.na(My.CP$My.0),"My.0"] <- 1
My.CP[is.na(My.CP$My.1),"My.1"] <- 1
My.CP[is.na(My.CP$My.2),"My.2"] <- 1

# Take the -log10 of the pvals for plotting
My.CP$My.0 <- -log10(My.CP$My.0)
My.CP$My.1 <- -log10(My.CP$My.1)
My.CP$My.2 <- -log10(My.CP$My.2)

# Clean up the table
row.names(My.CP) <- My.CP$name
My.CP$name <- NULL

# Take the top10 per cluster
My.CP.top10 <- My.CP[unique(c(as.character(My.0.PW$Pathways[1:10,"name"]),
                              as.character(My.1.PW$Pathways[1:10,"name"]),
                              as.character(My.2.PW$Pathways[1:10,"name"])
                              )),]

# And plot
pheatmap(My.CP.top10, 
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         cellwidth = 20, 
         cellheight = 12)


# Check myeloid genes in ACTA2 high cells
Idents(v3.all.seur.combined.M_clusters) <- factor(Idents(v3.all.seur.combined.M_clusters), levels = c("Mye.0","Mye.1","Mye.2","Mye.3", "Mye.4" ))

ACTA2 <- GetAssayData(v3.all.seur.combined.M_clusters[c("ACTA2"),])
ACTA2.cells <- colnames(v3.all.seur.combined.M_clusters)[as.vector(ACTA2 > 0)]
v3.all.seur.combined.M_clusters.ACTA2_cells <- subset(v3.all.seur.combined.M_clusters, cells = ACTA2.cells)

VlnPlot(v3.all.seur.combined.M_clusters.ACTA2_cells, cols = c("#2AB34B"), idents = "Mye.2", features = c("SPI1", "CEBPB", "MYOCD", "MRTFA"), ncol = 2)
ggsave("Myeloid_cell_subset_results/ACTA2.pos_Cells.TFs.pdf")

# Check diff with ACTA2 low cells
no.ACTA2.cells <- colnames(v3.all.seur.combined.M_clusters)[as.vector(ACTA2 == 0)]
v3.all.seur.combined.M_clusters.no_ACTA2_cells <- subset(v3.all.seur.combined.M_clusters, cells = no.ACTA2.cells)

ncol(v3.all.seur.combined.M_clusters.no_ACTA2_cells)
ncol(v3.all.seur.combined.M_clusters.ACTA2_cells)

v3.all.seur.combined.M_clusters@meta.data$ACTA2_level <- as.vector(ACTA2 != 0)
VlnPlot(v3.all.seur.combined.M_clusters, cols = c("#F2756D", "#B6A031", "#2AB34B"), split.by = "ACTA2_level", idents = c("Mye.0", "Mye.1", "Mye.2"), features = c("ACTA2", "KLF4"), ncol = 2)

# Check KLF4
KLF4 <- GetAssayData(v3.all.seur.combined.M_clusters[c("KLF4"),])
v3.all.seur.combined.M_clusters@meta.data$KLF4_level <- as.vector(KLF4 != 0)
Idents(v3.all.seur.combined.M_clusters) <- factor(Idents(v3.all.seur.combined.M_clusters), levels = c("Mye.0","Mye.1","Mye.2","Mye.3", "Mye.4" ))
VlnPlot(v3.all.seur.combined.M_clusters, cols = c("#F2756D", "#B6A031", "#2AB34B"), split.by = "KLF4_level", idents = c("Mye.0", "Mye.1", "Mye.2"), features = c("KLF4", "IL1B", "TNF", "ACTA2"), ncol = 2)

# Check Il1B
IL1B <- GetAssayData(v3.all.seur.combined.M_clusters[c("IL1B"),])
v3.all.seur.combined.M_clusters@meta.data$IL1B_level <- as.vector(IL1B != 0)
VlnPlot(v3.all.seur.combined.M_clusters, cols = c("#F2756D", "#B6A031", "#2AB34B"), split.by = "IL1B_level", idents = c("Mye.0"), features = c("CASP1", "CASP4"), ncol = 2)
ggsave("Myeloid_cell_subset_results/IL1B.pos_Cells.CASPs.pdf")

# Check mac markers in ACTA2 high cells
VlnPlot(v3.all.seur.combined.M_clusters.ACTA2_cells, cols = c("#2AB34B"), idents = "Mye.2", features = c("CD68", "LGALS3"), ncol = 2)
ggsave("Myeloid_cell_subset_results/ACTA2.pos_Cells.mac_markers.pdf")

v3.all.seur.combined.M_clusters@meta.data$ACTA2_level <- as.vector(ACTA2 != 0)
VlnPlot(v3.all.seur.combined.M_clusters, cols = c("#2AB34B"), idents = "Mye.2", split.by = "ACTA2_level", features = c("ACTA2","CD68", "LGALS3"), ncol = 3)
ggsave("Myeloid_cell_subset_results/ACTA2.pos_Cells.mac_markers_split.pdf")

