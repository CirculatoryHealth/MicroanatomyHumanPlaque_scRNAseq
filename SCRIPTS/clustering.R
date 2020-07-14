# t-SNE and Clustering
set.seed(1)
all.seur.combined <- RunTSNE(all.seur.combined, reduction.use = "pca.aligned", dims.use = 1:15, 
                             do.fast = T)
all.seur.combined <- FindClusters(all.seur.combined, reduction.type = "pca.aligned", 
                                  resolution = 1.2, dims.use = 1:15, force.recalc = T)

TSNEPlot(all.seur.combined, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                axis.text        = element_text(size = 14, face = "plain")
                                                                                )
ggsave("results/tSNE.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined, pt.size = 2, label.size = 10, group.by = "ID") + theme(panel.background = element_blank(),
                                                                                        text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                        axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/tSNE ID.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                           axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/tSNE Sex.pdf", width = 15, height = 15)
TSNEPlot(all.seur.combined, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                      text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                      axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/tSNE Patient.pdf", width = 15, height = 15)


#Get cluster colors for subsequent plots
t                      <- TSNEPlot(all.seur.combined, do.return = T)
tbuild                 <- ggplot_build(t)
tdata                  <- tbuild$data[[1]] 
tdata                  <- tdata[order(tdata$group), ]
cluster_colours        <- unique(tdata$colour)
#names(cluster_colours) <- new.ident

#Define marker genes per cluster
all.all.seur.combined.markers <- FindAllMarkers(object = all.seur.combined, only.pos = TRUE, min.pct = 0.25, 
                                                thresh.use = 0.25)

#Show top2 marker genes
all.all.seur.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#Save the top9 markers per cluster
sep.markers <- all.all.seur.combined.markers %>% group_by(cluster) %>% top_n(9, avg_logFC)

#plot the top9 markers per cluster
for(i in 0:13){
  FeaturePlot(object = all.seur.combined,
              features.plot = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])) ,
              cols.use = c("grey", "blue"), 
              reduction.use = "tsne", pt.size = 2
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/cluster",i, " markers.pdf", sep = ""), width = 15, height = 15)
  VlnPlot(object = all.seur.combined,
          features.plot = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(family = "Arial", size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}

#Plot top markers in a heatmap
top10 <- all.all.seur.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = all.seur.combined, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE) + theme(panel.background = element_blank(),
                                                                                                                 text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                 axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/clusters.pdf", width = 15, height =15)
DoHeatmap(object = all.seur.combined, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE) + theme(panel.background = element_blank(),
                                                                                                                text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                                                axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/clusters_key.pdf", width = 15, height =15)


#Get cluster cell types via singleR
#See singleR_analyses.R

#Rename clusters and replot tSNE
new.ident <- c( "CD3+CD8+ T-Cells I", "CD3+CD4+ T-Cells I", "Mixed Cells", #0,1,2
                "CD3+CD8+ T-Cells II", "CD3+CD4+ T-Cells II", "CD14+CD68+ Macrophages I", #3,4,5
                "CD14+CD68+ Macrophages II", "CD14+CD68+ Macrophages III", "MYH11+ Smooth Muscle Cells", #6,7,8
                "CD34+ Endothelial Cells", "CD34+ Endothelial Cells II", "CD79A+ B-Cells", #9,10,11
                "CD14+CD68+ Macrophages IV", "KIT+ Mast Cells") #12,13
for (i in 0:13) {
  all.seur.combined <- RenameIdent(object = all.seur.combined, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

TSNEPlot(all.seur.combined, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/tSNE labels.pdf", width = 15, height = 15)


#Compare cluster size per patient
freq_table <- prop.table(x = table(all.seur.combined@ident, all.seur.combined@meta.data[, "Patient"]), margin = 2)
row.names(freq_table) <- new.ident
m <- melt(freq_table)
m$Var2 <- as.factor(m$Var2)
ggplot(m, aes(x = Var2, fill = Var1, y = value)) +
  geom_col(position = position_fill(reverse = T)) +
  coord_flip() +
  theme_light() + 
  ylab("% of cells") +
  xlab("Patient") +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.50,0.75,1.00), labels = c(0,25,50,75,100)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(), aspect.ratio = 1:1
  )
ggsave("results/freq patient.pdf",width = 10, height = 10)

freq_table <- prop.table(x = table(all.seur.combined@ident, all.seur.combined@meta.data[, "Sex"]), margin = 2)
m <- melt(freq_table)
m$Var2 <- as.factor(m$Var2)
ggplot(m, aes(x = Var2, fill = Var1, y = value)) +
  geom_col(position = position_fill(reverse = T)) +
  coord_flip() +
  theme_light() + 
  ylab("% of cells") +
  xlab("Sex") +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.25,0.50,0.75,1.00), labels = c(0,25,50,75,100)) +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(), aspect.ratio = 1:1
  )
ggsave("results/freq sex.pdf",width = 10, height = 10)

freq_table <- prop.table(x = table(all.seur.combined@ident, all.seur.combined@meta.data[, "Plate"]), margin = 2)
pdf("results/freq Plate.pdf")
barplot(height = freq_table, legend.text = T, horiz = T)
dev.off()
row.names(freq_table) <- new.ident
freq_table

#Different view of markers
markers.to.plot <- c("CD3E", "CD4", "CD8A", "CD14", 
                     "CD68", "CD79A", "MYH11", "ACTA2", 
                     "KIT", "CD34", "OLR1","ABCA1", "CD163", 
                     "SELL", "TIGIT", "FOS", "PRF1", 
                     "GZMK", "FOSB","IL7R" ,"TNF", "IL1B","PTPRC")

DotPlot(v3.all.seur.combined, features = rev(markers.to.plot),
                      cols = c("salmon", "cornflowerblue"), dot.scale = 8, split.by = "Sex") + theme(panel.background = element_blank(),
                                                                       text             = element_text(size = 16, face = "bold"),
                                                                       axis.text        = element_text(size = 14, face = "plain")
                      )
ggsave("results/dot.pdf", width = 10, height = 20)

sdp <- Spl(all.seur.combined, genes.plot = rev(markers.to.plot), 
                      cols.use = colorRampPalette(c("Red","Blue"))(18), 
                      x.lab.rot = T, plot.legend = T, dot.scale = 8, 
                      do.return = T, grouping.var = "Patient") + theme(panel.background = element_blank(),
                                                                          text             = element_text(size = 16, face = "bold"),
                                                                          axis.text        = element_text(size = 14, face = "plain")
                      )
ggsave("results/patient dot.pdf", width = 10, height = 20)


#Write var genes per cluster to files
for (i in 0:13){
  x <- all.all.seur.combined.markers[which(all.all.seur.combined.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

#Pathway analyses per cluster
for (i in 0:13){
  print("")
  print(paste("Working on cluster:", i, sep = " "))
  print("")
  x <- all.all.seur.combined.markers[which(all.all.seur.combined.markers$cluster == i),-6]
  get_ontology(x, name = paste("Cluster", i, sep = "_"), return.data = F, outdir = "results", full_GSEA = F)
}




#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#Now let's look at some genes/clusters!
#First let's get the top 1000 genes for all clusters
dirname <- "results/top1000_expressed_genes"
dir.create(dirname, showWarnings = FALSE)

#Get the average expression  matrix
cluster.averages <- AverageExpression(all.seur.combined)

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
    get_ontology(x, name = paste("Cluster", names(df.list[i]), "Pathways_for_kmeans_cluster", j, sep = "_"), return.data = F, outdir = subdirname, full_GSEA = F)
  }
}

#-----------------------------------------------------------------------------------------------
#Compare expression levels in cluster 2 (weird Myeloid/T-cell cluster) with the rest
#Get all genes expressed over 0 (per cluster)
exp.list <- list()
for (i in 1:ncol(cluster.averages)){
  l <- cluster.averages[which(cluster.averages[,i] > 0),i]
  exp.list[[i]] <- l
}
names(exp.list) <- colnames(cluster.averages)

#Check the stats
lapply(exp.list,summary)
lapply(exp.list,length)

#Box plot for expression levels per cluster
median.order <- order(unlist(lapply(exp.list, median)), decreasing = T)
pdf("results/expression_of_genes_per_cluster_cut-off-over-0.pdf")
par(mar=c(15,4,2,1)+0.1)
boxplot(exp.list[median.order], 
        las = 2, 
        outline = F, 
        notch = T, 
        col = cluster_colours[median.order], 
        ylab = "Normalized Expression",
        main = "All genes per cluster expressed > 0"
        )
dev.off()

#Plot bars for number of genes
d <- melt(lapply(exp.list,length))
colnames(d) <- c("Number","Cluster")

ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of genes epressed > 0", breaks = c(10000,12500,15000,17500,20000), limits = c(0000,20000)) +
  coord_cartesian(ylim = c(10000,20000)) +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()
  )
ggsave("results/number_of_genes_expressed_over_0.pdf")

#Get all genes expressed over 1 (per cluster)
exp.list <- list()
for (i in 1:ncol(cluster.averages)){
  l <- cluster.averages[which(cluster.averages[,i] > 1),i]
  exp.list[[i]] <- l
}
names(exp.list) <- colnames(cluster.averages)

#Check the stats
lapply(exp.list,summary)
lapply(exp.list,length)

#Box plot for expression levels per cluster
median.order <- order(unlist(lapply(exp.list, median)), decreasing = T)
pdf("results/expression_of_genes_per_cluster_cut-off-over-1.pdf")
par(mar=c(15,4,2,1)+0.1)
boxplot(exp.list[median.order], 
        las = 2, 
        outline = F, 
        notch = T, 
        col = cluster_colours[median.order], 
        ylab = "Normalized Expression",
        main = "All genes per cluster expressed > 1"
)
dev.off()

#Plot bars for number of genes
d <- melt(lapply(exp.list,length))
colnames(d) <- c("Number","Cluster")

ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of genes epressed > 1", breaks = c(1000,1500,2000,2500), limits = c(0000,2500)) +
  coord_cartesian(ylim = c(1000,2500)) +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()
  )
ggsave("results/number_of_genes_expressed_over_1.pdf")

#Get all genes expressed over 5 (per cluster)
exp.list <- list()
for (i in 1:ncol(cluster.averages)){
  l <- cluster.averages[which(cluster.averages[,i] > 5),i]
  exp.list[[i]] <- l
}
names(exp.list) <- colnames(cluster.averages)

#Check the stats
lapply(exp.list,summary)
lapply(exp.list,length)

#Box plot for expression levels per cluster
median.order <- order(unlist(lapply(exp.list, median)), decreasing = T)
pdf("results/expression_of_genes_per_cluster_cut-off-over-5.pdf")
par(mar=c(15,4,2,1)+0.1)
boxplot(exp.list[median.order], 
        las = 2, 
        outline = F, 
        notch = T, 
        col = cluster_colours[median.order], 
        ylab = "Normalized Expression",
        main = "All genes per cluster expressed > 5"
)
dev.off()

#Plot bars for number of genes
d <- melt(lapply(exp.list,length))
colnames(d) <- c("Number","Cluster")

ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of genes epressed > 5", breaks = c(100,150,200,250,300,350), limits = c(0000,350)) +
  coord_cartesian(ylim = c(100,350)) +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()
  )
ggsave("results/number_of_genes_expressed_over_5.pdf")

#Get all genes expressed over 10 (per cluster)
exp.list <- list()
for (i in 1:ncol(cluster.averages)){
  l <- cluster.averages[which(cluster.averages[,i] > 10),i]
  exp.list[[i]] <- l
}
names(exp.list) <- colnames(cluster.averages)

#Check the stats
lapply(exp.list,summary)
lapply(exp.list,length)

#Box plot for expression levels per cluster
median.order <- order(unlist(lapply(exp.list, median)), decreasing = T)
pdf("results/expression_of_genes_per_cluster_cut-off-over-10.pdf")
par(mar=c(15,4,2,1)+0.1)
boxplot(exp.list[median.order], 
        las = 2, 
        outline = F, 
        notch = T, 
        col = cluster_colours[median.order], 
        ylab = "Normalized Expression",
        main = "All genes per cluster expressed > 10"
)
dev.off()

#Plot bars for number of genes
d <- melt(lapply(exp.list,length))
colnames(d) <- c("Number","Cluster")

ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of genes epressed > 10", breaks = c(25,50,75,100,125,150), limits = c(0000,150)) +
  coord_cartesian(ylim = c(0,150)) +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()
  )
ggsave("results/number_of_genes_expressed_over_10.pdf")


#-----------------------------------------------------------------------------------------------
#Plot number of cells per cluster
#Get the number of cells in each cluster
cell.numbers <- list()
for(i in new.ident){
  cell.numbers[i] <- length(SubsetData(all.seur.combined, ident.use = i)@cell.names)
}

#Get the total number of cells
total.cells <- length(all.seur.combined@cell.names)

#Sanity check
total.cells == Reduce("+", cell.numbers)

#Plot bars for number or percentage of cells
d <- melt(cell.numbers)
colnames(d) <- c("Number","Cluster")
d$Percentage <- (d$Number / total.cells) * 100
d$Cluster_Number <- as.factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

#Sanity check
Reduce("+", d$Percentage)

#Plot absolute values
ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of cells", breaks = c(0,50,100,150,200,250,300,500,750,1000), limits = c(000,1000)) +
  theme_light() +
  ggtitle("Number of cells per cluster") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_cluster.pdf")

#Plot relative values
ggplot(d, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "% of cells", breaks = c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25), limits = c(00,25)) +
  theme_light() +
  ggtitle("Percentage of cells per cluster") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_cluster_relative.pdf")

#Plot absolute values with numbers
ggplot(d, aes(x = Cluster, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", position = "Dodge") +
  geom_text(aes(label=Number), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), name = "Number of cells", breaks = c(0,50,100,150,200,250,300,500,750,1000), limits = c(000,1000)) +
  theme_light() +
  ggtitle("Number of cells per cluster") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_cluster_show_value.pdf")

#Plot relative values with numbers
ggplot(d, aes(x = Cluster, y = Percentage, fill = Cluster_Number)) +
  geom_bar(stat = "identity", position = "Dodge") +
  geom_text(aes(label=round(Percentage, 1)), position=position_dodge(width = 1), vjust = 1.2, colour = "yellow") +
  scale_fill_manual(values = cluster_colours) +
  scale_x_discrete(limits = d$Cluster[order(d$Number, decreasing = T)]) +
  scale_y_continuous(expand = c(0,0), 
                     name = "% of cells", 
                     breaks = c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25), 
                     limits = c(00,25),
                     sec.axis = sec_axis(trans = ~ (. / 100) * total.cells, name = "# of cells", breaks = c(0,50,100,150,200,250,300,500))) +
  theme_light() +
  ggtitle("Percentage and number of cells per cluster") +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 270, hjust = 0),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size=0.1, colour = "grey")
  )
ggsave("results/number_of_cell_per_cluster_relative_show_values.pdf")


#----------------------------------------------------------
#Calculate GWAS enrichment
gwas.deg.patterns <- data.frame(gwas = scan(file = "gwas.length",what = "list"), 
                                deg = scan(file = "deg.length",what = "list"), 
                                pval = "NA", 
                                stringsAsFactors = F
                                )

gwas.deg.patterns$gwas <- as.numeric(gwas.deg.patterns$gwas)
gwas.deg.patterns$deg <- as.numeric(gwas.deg.patterns$deg)
gwas.deg <- 65
tot.pop <- 3244

for (i in 1:nrow(gwas.deg.patterns)){
  gwas.deg.patterns[i,"pval"] <- phyper(gwas.deg.patterns[i, "gwas"] - 1, 
                                        gwas.deg.patterns[i, "deg"],
                                        tot.pop - gwas.deg.patterns[i, "deg"], 
                                        gwas.deg, 
                                        lower.tail = F 
  )
}


for(i in T_clusters){
  cells <- SubsetData(all.seur.combined, ident.use = i)
  cells <- as.matrix(cells@data)
  write.table(transform(cells, Symbol = rownames(cells))[,c(length(colnames(cells))+1,1:length(colnames(cells)))],
              row.names = F, quote = F, sep = "\t", file = paste("results/cluster_",i,"_normalized_counts.txt",sep = ""))
}


#----------------------------------------------------------
#Export count table
write.table(as.matrix(all.seur.combined@data),
            file = "big_table.relative_expression.normalized_and_log_transformed_counts.txt", 
            quote = F, 
            sep = "\t")



all.all.seur.combined.markers


#---------------------------------------------------------
#Density plots of number of genes per cell
m <- melt(all.seur.combined@meta.data[,c(1,10)])
colnames(m) <- c("Cluster", "var", "value")
m$Cluster <- as.factor(m$Cluster)
m$Cluster <- mapvalues(m$Cluster, from = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), to = levels(all.seur.combined@ident))
m$Cluster <- factor(m$Cluster,levels(all.seur.combined@ident))

mu <- ddply(m, "Cluster", summarise, grp.mean=mean(value))

ggplot(m,aes(x = value, fill = Cluster)) +
  geom_density(alpha = 0.4) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Cluster),
             linetype="dashed") +
  xlab("# of genes per cell") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(500, 6000), breaks = c(500,1000,2000,3000,4000,5000,6000)) +
  theme_light() +
  theme(text = element_text(size = 16, family = "Arial"),
        panel.border = element_rect(size = 1, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1), 
        axis.ticks = element_line(size = 0.5, colour = "black"),
        aspect.ratio = 3/2
        )
ggsave("Number of genes per cell density plot.pdf", width = 10, height = 10)



#write out the normalized count table
write.better.table(as.data.frame(as.matrix(all.seur.combined@data)), file = "normalized_expression_18p.txt")

#Save the object
save(all.seur.combined, file = "all.seur.combined.RData")

# Export myeloid and T sub clusters
Mye.idents <- paste0("Mye.", v3.all.seur.combined.M_clusters@active.ident)
names(Mye.idents) <- names(v3.all.seur.combined.M_clusters@active.ident)
Idents(v3.all.seur.combined.M_clusters) <- Mye.idents

CD4.idents <- paste0("CD4.", v3.all.seur.combined.T_split_clusters.CD4@active.ident)
names(CD4.idents) <- names(v3.all.seur.combined.T_split_clusters.CD4@active.ident)
Idents(v3.all.seur.combined.T_split_clusters.CD4) <- CD4.idents

CD8.idents <- paste0("CD8.", all.seur.combined.T_split_clusters.CD8@active.ident)
names(CD8.idents) <- names(all.seur.combined.T_split_clusters.CD8@active.ident)
Idents(all.seur.combined.T_split_clusters.CD8) <- CD8.idents

# Rename CD8 object for constistency
v3.all.seur.combined.T_split_clusters.CD8 <- all.seur.combined.T_split_clusters.CD8

saveRDS(v3.all.seur.combined.M_clusters, file = "v3.all.seur.combined.M_clusters")
saveRDS(v3.all.seur.combined.T_split_clusters.CD4, file = "v3.all.seur.combined.T_split_clusters.CD4")
saveRDS(v3.all.seur.combined.T_split_clusters.CD8, file = "v3.all.seur.combined.T_split_clusters.CD8")

#Add subcluster names and replot tSNE
v3.all.seur.combined <- UpdateSeuratObject(all.seur.combined)

# We already prepended the idents for these sub clusters
#Mye.idents <- paste0("Mye.", v3.all.seur.combined.M_clusters@active.ident)
names(Mye.idents) <- names(v3.all.seur.combined.M_clusters@active.ident)

#CD4.idents <- paste0("CD4.", v3.all.seur.combined.T_split_clusters.CD4@active.ident)
names(CD4.idents) <- names(v3.all.seur.combined.T_split_clusters.CD4@active.ident)

#CD8.idents <- paste0("CD8.", all.seur.combined.T_split_clusters.CD8@active.ident)
names(CD8.idents) <- names(all.seur.combined.T_split_clusters.CD8@active.ident)

E.idents <- paste0("E.", v3.all.seur.combined.E_clusters@active.ident)
names(E.idents) <- names(v3.all.seur.combined.E_clusters@active.ident)

S.idents <- paste0("SMC.", v3.all.seur.combined.SMC_clusters@active.ident)
names(S.idents) <- names(v3.all.seur.combined.SMC_clusters@active.ident)

B.idents <- rep("CD79A+ B-Cells", length(WhichCells(v3.all.seur.combined, idents = "CD79A+ B-Cells")))
names(B.idents) <- WhichCells(v3.all.seur.combined, idents = "CD79A+ B-Cells")

Mast.idents <- rep("KIT+ Mast Cells", length(WhichCells(v3.all.seur.combined, idents = "KIT+ Mast Cells")))
names(Mast.idents) <- WhichCells(v3.all.seur.combined, idents = "KIT+ Mast Cells")

Mixed.idents <- rep("Mixed Cells", length(WhichCells(v3.all.seur.combined, idents = "Mixed Cells")))
names(Mixed.idents) <- WhichCells(v3.all.seur.combined, idents = "Mixed Cells")

sub.idents <- c(Mye.idents, CD4.idents, CD8.idents, E.idents, S.idents, B.idents, Mast.idents, Mixed.idents)
v3.all.seur.combined <- AddMetaData(v3.all.seur.combined, sub.idents, col.name = "sub.idents")

# We threw some T cells out because we couldnt confidently call them CD4 or CD8, let's call them CD3+
CD3.idents <- rep("CD3+ T-cells", length(v3.all.seur.combined@active.ident[is.na(v3.all.seur.combined@meta.data$sub.idents)]))
names(CD3.idents) <- names(v3.all.seur.combined@active.ident[is.na(v3.all.seur.combined@meta.data$sub.idents)])
sub.idents <- c(sub.idents, CD3.idents)
v3.all.seur.combined <- AddMetaData(v3.all.seur.combined, sub.idents, col.name = "sub.idents")

TSNEPlot(v3.all.seur.combined, group.by = "sub.idents",  pt.size = 2, label.size = 6, label=T) + theme(panel.background = element_blank(),
                                                                                    text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                    axis.text        = element_text(size = 14, face = "plain"),
                                                                                    aspect.ratio     = 1 
)
ggsave("results/tSNE subclusters.pdf", width = 10, height = 10)
saveRDS(v3.all.seur.combined, file = "v3.all.seur.combined.RDS")