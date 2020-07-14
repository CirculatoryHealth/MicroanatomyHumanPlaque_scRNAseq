# Extract number of cells per patient
# Update our object to seurat v3
v3.all.seur.combined <- UpdateSeuratObject(all.seur.combined)

# Populate with total cells per patient
num_cells.per_patient <- data.frame()
for(i in unique(v3.all.seur.combined$Patient)){
  num_cells.per_patient[i, "Total"] <- sum(v3.all.seur.combined$Patient == i)
}
# Display
num_cells.per_patient 

# Add cells per cluster
for(i in unique(v3.all.seur.combined$Patient)){
  for(j in unique(Idents(v3.all.seur.combined))){
    num_cells.per_patient[i, j] <- sum(v3.all.seur.combined$Patient == i & Idents(v3.all.seur.combined) == j)
  }
}

# Order on patient number
num_cells.per_patient <- num_cells.per_patient[order(row.names(num_cells.per_patient)),]

# Sanity check
cat("Rows are sane: ", sum((rowSums(num_cells.per_patient) - num_cells.per_patient$Total) == num_cells.per_patient$Total) == length(unique(v3.all.seur.combined$Patient)), "\n")
cat("Columns are sane: ", sum(colSums(num_cells.per_patient[,2:dim(num_cells.per_patient)[2]])[order(colnames(num_cells.per_patient[,2:dim(num_cells.per_patient)[2]]))] == cell.numbers[order(names(cell.numbers))]) == length(unique(Idents(v3.all.seur.combined))), "\n")

# Display
num_cells.per_patient

# Remove outlier patient 4440 (only 7 cells total o_O)
num_cells.per_patient <- num_cells.per_patient[which(row.names(num_cells.per_patient) != "4440"),]

# move 'Totals column to it's own frame
total_cells.per_patient <- num_cells.per_patient[,"Total", drop = F]
num_cells.per_patient$Total <- NULL

# Normalize to number of cells per patient
x <- colnames(num_cells.per_patient)
num_cells.per_patient <- as.data.frame(apply(num_cells.per_patient,2,function(x) x / total_cells.per_patient))
colnames(num_cells.per_patient) <- x

# Add plaque phenotype metadata
plaque.type <- unique(read.table(file= "Project SCS 18 patients meta data.txt", header=T)[,c(1,3)])
plaque.type <- plaque.type[order(plaque.type$AE),]
row.names(plaque.type) <- plaque.type$AE
plaque.type$AE <- NULL
colnames(plaque.type) <- "Type"
plaque.type <- plaque.type[!row.names(plaque.type) == "4440", ,drop =F] #Remove outlier


# Split per plaque type
A.num_cells.per_patient  <- num_cells.per_patient[which(plaque.type$Type == "Atheromatous"),]
F.num_cells.per_patient  <- num_cells.per_patient[which(plaque.type$Type == "Fibrous"),]
FA.num_cells.per_patient <- num_cells.per_patient[which(plaque.type$Type == "Fibro-atheromatous"),]

# Combine counts per type
summed.num_cells.per_patient <- t(data.frame(Atheromatous = apply(A.num_cells.per_patient, 2, sum), 
                                                 Fibrous = apply(F.num_cells.per_patient, 2, sum),
                                                 Fibro.atheromatous = apply(FA.num_cells.per_patient, 2, sum))
                                        )

# Run Fisher Exact tests
fishere.summed.num_cells.per_patient <- list()
for (i in 1:(dim(summed.num_cells.per_patient)[2]-1)){
  sublist <- list()
  for (j in (i+1):dim(summed.num_cells.per_patient)[2]){
    sublist[[colnames(num_cells.per_patient)[j]]] <- fisher.test(num_cells.per_patient[,i], num_cells.per_patient[,j])
  }
  fishere.summed.num_cells.per_patient[[colnames(num_cells.per_patient)[i]]] <- sublist
}

# Extract p values
pvals.fishere.summed.num_cells.per_patient <- data.frame()
for (i in 1:(dim(summed.num_cells.per_patient)[2]-1)){
  for (j in (i+1):dim(summed.num_cells.per_patient)[2]){
    pvals.fishere.summed.num_cells.per_patient <- rbind(pvals.fishere.summed.num_cells.per_patient , 
                                                           data.frame(First   = colnames(num_cells.per_patient)[i],
                                                                      Second  = colnames(num_cells.per_patient)[j],
                                                                      p.value = fishere.summed.num_cells.per_patient[[colnames(num_cells.per_patient)[i]]][[colnames(num_cells.per_patient)[j]]]$p.value
                                                     )
    )
  }
}

#Adjust p value for multiple testing
pvals.fishere.summed.num_cells.per_patient$padj <- p.adjust(pvals.fishere.summed.num_cells.per_patient$p.value, method = "BH")

# Check for significance
pvals.fishere.summed.num_cells.per_patient[which(pvals.fishere.summed.num_cells.per_patient$padj < 0.05),]
pvals.fishere.summed.num_cells.per_patient[which(pvals.fishere.summed.num_cells.per_patient$p.value < 0.05),]

# Try MANOVA instead
m <- manova(as.matrix(num_cells.per_patient) ~ Type, data = plaque.type)
summary(m)
summary.aov(m)

# Build a plot based on percentages to visualize
perc.summed.num_cells.per_patient <- prop.table(summed.num_cells.per_patient, 2)
write.table(perc.summed.num_cells.per_patient, file="percentage_types.txt", quote = F, sep = "\t")

d <- melt(perc.summed.num_cells.per_patient)
colnames(d) <- c("Type", "Cluster", "Percentage")
d$Percentage <- d$Percentage * 100

ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("indianred3", "darkolivegreen4", "goldenrod")) +
  ylab("Percentage of Cells") +
  ggtitle("Main clusters") +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(size = 12,
                                   angle = 45, 
                                   hjust = 1),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(), 
        aspect.ratio = 1/2
  )
ggsave("pheno.cor.pdf")


#===========================================================================
#===========================================================================
#Let's look at some sub-clusters
# CD4, CD8, Endo subsets: not enough cells in a lot of patients!
#Macs reasonably OK so let's check those

# Extract number of cells per patient
# Update our object to seurat v3
v3.all.seur.combined.M_clusters <- UpdateSeuratObject(all.seur.combined.M_clusters)

# Populate with total cells per patient
num_cells.per_patient.M_clusters <- data.frame()
for(i in unique(v3.all.seur.combined.M_clusters$Patient)){
  num_cells.per_patient.M_clusters[i, "Total"] <- sum(v3.all.seur.combined.M_clusters$Patient == i)
}
# Display
num_cells.per_patient.M_clusters

# Add cells per cluster
for(i in unique(v3.all.seur.combined.M_clusters$Patient)){
  for(j in unique(Idents(v3.all.seur.combined.M_clusters))){
    num_cells.per_patient.M_clusters[i, j] <- sum(v3.all.seur.combined.M_clusters$Patient == i & Idents(v3.all.seur.combined.M_clusters) == j)
  }
}

# Order on patient number
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[order(row.names(num_cells.per_patient.M_clusters)),]

# Display
num_cells.per_patient.M_clusters

# Remove outliers patient 4440 and 4477
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[which(row.names(num_cells.per_patient.M_clusters) != "4440"),]
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[which(row.names(num_cells.per_patient.M_clusters) != "4477"),]

# Remove 'T-cell contaminent' subcluster 4
num_cells.per_patient.M_clusters$Total <- num_cells.per_patient.M_clusters$Total - num_cells.per_patient.M_clusters$`4`
num_cells.per_patient.M_clusters$`4` <- NULL

# move 'Totals column to it's own frame
total_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[,"Total", drop = F]
num_cells.per_patient.M_clusters$Total <- NULL

# Normalize to number of cells per patient
x <- colnames(num_cells.per_patient.M_clusters)
num_cells.per_patient.M_clusters <- as.data.frame(apply(num_cells.per_patient.M_clusters,2,function(x) x / total_cells.per_patient.M_clusters))
colnames(num_cells.per_patient.M_clusters) <- x

# Add plaque phenotype metadata
plaque.type.M_clusters <- plaque.type
plaque.type.M_clusters <- plaque.type.M_clusters[!row.names(plaque.type.M_clusters) == "4477",, drop = F]

# Split per plaque type
A.num_cells.per_patient.M_clusters  <- num_cells.per_patient.M_clusters[which(plaque.type.M_clusters$Type == "Atheromatous"),]
F.num_cells.per_patient.M_clusters  <- num_cells.per_patient.M_clusters[which(plaque.type.M_clusters$Type == "Fibrous"),]
FA.num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[which(plaque.type.M_clusters$Type == "Fibro-atheromatous"),]

# Combine counts per type
summed.num_cells.per_patient.M_clusters <- t(data.frame(Atheromatous = apply(A.num_cells.per_patient.M_clusters, 2, sum), 
                                             Fibrous = apply(F.num_cells.per_patient.M_clusters, 2, sum),
                                             Fibro.atheromatous = apply(FA.num_cells.per_patient.M_clusters, 2, sum))
)

# Run Fisher Exact tests
fishere.summed.num_cells.per_patient.M_clusters <- list()
for (i in 1:(dim(summed.num_cells.per_patient.M_clusters)[2]-1)){
  sublist <- list()
  for (j in (i+1):dim(summed.num_cells.per_patient.M_clusters)[2]){
    sublist[[colnames(num_cells.per_patient.M_clusters)[j]]] <- fisher.test(num_cells.per_patient.M_clusters[,i], num_cells.per_patient.M_clusters[,j])
  }
  fishere.summed.num_cells.per_patient.M_clusters[[colnames(num_cells.per_patient.M_clusters)[i]]] <- sublist
}

# Extract p values
pvals.fishere.summed.num_cells.per_patient.M_clusters <- data.frame()
for (i in 1:(dim(summed.num_cells.per_patient.M_clusters)[2]-1)){
  for (j in (i+1):dim(summed.num_cells.per_patient.M_clusters)[2]){
    pvals.fishere.summed.num_cells.per_patient.M_clusters <- rbind(pvals.fishere.summed.num_cells.per_patient.M_clusters , 
                                                      data.frame(First   = colnames(num_cells.per_patient.M_clusters)[i],
                                                                 Second  = colnames(num_cells.per_patient.M_clusters)[j],
                                                                 p.value = fishere.summed.num_cells.per_patient.M_clusters[[colnames(num_cells.per_patient.M_clusters)[i]]][[colnames(num_cells.per_patient.M_clusters)[j]]]$p.value
                                                      )
    )
  }
}

#Adjust p value for multiple testing
pvals.fishere.summed.num_cells.per_patient.M_clusters$padj <- p.adjust(pvals.fishere.summed.num_cells.per_patient.M_clusters$p.value, method = "BH")

# Check for significance
pvals.fishere.summed.num_cells.per_patient.M_clusters[which(pvals.fishere.summed.num_cells.per_patient.M_clusters$padj < 0.05),]
pvals.fishere.summed.num_cells.per_patient.M_clusters[which(pvals.fishere.summed.num_cells.per_patient.M_clusters$p.value < 0.05),]
### Nothing significant ###

# Try MANOVA instead
m.M_clusters <- manova(as.matrix(num_cells.per_patient.M_clusters) ~ Type, data = plaque.type.M_clusters)
summary(m.M_clusters)
summary.aov(m.M_clusters)
### Nope also nothing ###

# Build a plot based on percentages to visualize
perc.summed.num_cells.per_patient.M_clusters <- prop.table(summed.num_cells.per_patient.M_clusters, 2)

d <- melt(perc.summed.num_cells.per_patient.M_clusters)
colnames(d) <- c("Type", "Cluster", "Percentage")
d$Cluster <- paste("My.", d$Cluster, sep = "")
d$Cluster <- as.factor(d$Cluster)
d$Percentage <- d$Percentage * 100

ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("indianred3", "darkolivegreen4", "goldenrod")) +
  ylab("Percentage of Cells") +
  ggtitle("Myeloid subclusters") +
  theme_light() +
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()
  )
ggsave("pheno.cor.M_clusters.pdf")

#=========================================================================================================
# Get pathways per plaque type
# Fetch cells per type
A.M.Cells  <- WhichCells(object     = v3.all.seur.combined.M_clusters, 
                         expression = Patient %in% row.names(plaque.type.M_clusters[which(plaque.type.M_clusters$Type == "Atheromatous"),,drop = F]))
F.M.Cells  <- WhichCells(object     = v3.all.seur.combined.M_clusters, 
                         expression = Patient %in% row.names(plaque.type.M_clusters[which(plaque.type.M_clusters$Type == "Fibrous"),,drop = F]))
AF.M.Cells <- WhichCells(object     = v3.all.seur.combined.M_clusters, 
                         expression = Patient %in% row.names(plaque.type.M_clusters[which(plaque.type.M_clusters$Type == "Fibro-atheromatous"),,drop = F]))

# Subset the seurat object per type
v3.all.seur.combined.M_clusters.A.cells  <- subset(v3.all.seur.combined.M_clusters, cells = A.M.Cells)
v3.all.seur.combined.M_clusters.F.cells  <- subset(v3.all.seur.combined.M_clusters, cells = F.M.Cells)
v3.all.seur.combined.M_clusters.AF.cells <- subset(v3.all.seur.combined.M_clusters, cells = AF.M.Cells)

# Fetch average expression per cluster (except the T-contiminant cluster 4)
A.M.avg  <- AverageExpression(v3.all.seur.combined.M_clusters.A.cells)$RNA[,-5]
F.M.avg  <- AverageExpression(v3.all.seur.combined.M_clusters.F.cells)$RNA[,-5]
AF.M.avg <- AverageExpression(v3.all.seur.combined.M_clusters.AF.cells)$RNA[,-5]

# Remove genes with no expression in any cluster
A.M.avg  <- A.M.avg[rowSums(A.M.avg) != 0,]
F.M.avg  <- F.M.avg[rowSums(F.M.avg) != 0,]
AF.M.avg <- AF.M.avg[rowSums(AF.M.avg) != 0,]

# k-means cluster average expression values to extract population specific genes
set.seed(666)
k.A.M.avg   <- pheatmap(A.M.avg, scale = "row", kmeans_k = 4)
My0.A.M.avg <- A.M.avg[k.A.M.avg$kmeans$cluster == 3,]
My1.A.M.avg <- A.M.avg[k.A.M.avg$kmeans$cluster == 2,]
My2.A.M.avg <- A.M.avg[k.A.M.avg$kmeans$cluster == 1,]
My3.A.M.avg <- A.M.avg[k.A.M.avg$kmeans$cluster == 4,]

set.seed(666)
k.F.M.avg   <- pheatmap(F.M.avg, scale = "row", kmeans_k = 4)
My0.F.M.avg <- F.M.avg[k.F.M.avg$kmeans$cluster == 1,]
My1.F.M.avg <- F.M.avg[k.F.M.avg$kmeans$cluster == 2,]
My2.F.M.avg <- F.M.avg[k.F.M.avg$kmeans$cluster == 4,]
My3.F.M.avg <- F.M.avg[k.F.M.avg$kmeans$cluster == 3,]

set.seed(666)
k.AF.M.avg   <- pheatmap(AF.M.avg, scale = "row", kmeans_k = 4)
My0.AF.M.avg <- AF.M.avg[k.AF.M.avg$kmeans$cluster == 1,]
My1.AF.M.avg <- AF.M.avg[k.AF.M.avg$kmeans$cluster == 3,]
My2.AF.M.avg <- AF.M.avg[k.AF.M.avg$kmeans$cluster == 4,]
My3.AF.M.avg <- AF.M.avg[k.AF.M.avg$kmeans$cluster == 2,]

# Sort on highest expression
My0.A.M.avg <- My0.A.M.avg[order(My0.A.M.avg$`0`, decreasing = T),]
My1.A.M.avg <- My1.A.M.avg[order(My1.A.M.avg$`1`, decreasing = T),]
My2.A.M.avg <- My2.A.M.avg[order(My2.A.M.avg$`2`, decreasing = T),]
My3.A.M.avg <- My3.A.M.avg[order(My3.A.M.avg$`3`, decreasing = T),]

My0.F.M.avg <- My0.F.M.avg[order(My0.F.M.avg$`0`, decreasing = T),]
My1.F.M.avg <- My1.F.M.avg[order(My1.F.M.avg$`1`, decreasing = T),]
My2.F.M.avg <- My2.F.M.avg[order(My2.F.M.avg$`2`, decreasing = T),]
My3.F.M.avg <- My3.F.M.avg[order(My3.F.M.avg$`3`, decreasing = T),]

My0.AF.M.avg <- My0.AF.M.avg[order(My0.AF.M.avg$`0`, decreasing = T),]
My1.AF.M.avg <- My1.AF.M.avg[order(My1.AF.M.avg$`1`, decreasing = T),]
My2.AF.M.avg <- My2.AF.M.avg[order(My2.AF.M.avg$`2`, decreasing = T),]
My3.AF.M.avg <- My3.AF.M.avg[order(My3.AF.M.avg$`3`, decreasing = T),]

# Collect top 25 per cluster for visualization
top25.A.M.avg <- rbind(My0.A.M.avg[1:25,], My1.A.M.avg[1:25,])
top25.A.M.avg <- rbind(top25.A.M.avg, My2.A.M.avg[1:25,])
top25.A.M.avg <- rbind(top25.A.M.avg, My3.A.M.avg[1:25,])

top25.F.M.avg <- rbind(My0.F.M.avg[1:25,], My1.F.M.avg[1:25,])
top25.F.M.avg <- rbind(top25.F.M.avg, My2.F.M.avg[1:25,])
top25.F.M.avg <- rbind(top25.F.M.avg, My3.F.M.avg[1:25,])

top25.AF.M.avg <- rbind(My0.AF.M.avg[1:25,], My1.AF.M.avg[1:25,])
top25.AF.M.avg <- rbind(top25.AF.M.avg, My2.AF.M.avg[1:25,])
top25.AF.M.avg <- rbind(top25.AF.M.avg, My3.AF.M.avg[1:25,])

top25.M.avg <- rbind(top25.A.M.avg, top25.F.M.avg)
top25.M.avg <- rbind(top25.M.avg, top25.AF.M.avg)

# Plot
pheatmap(top25.A.M.avg, 
         cluster_rows = F, 
         cluster_cols = F,
                scale = "row", 
            cellwidth = 50, 
           cellheight = 10,
             gaps_row = c(25,50,75),
             filename = "pheno.cor.Ath.Mye.top25.heatmap.pdf"
         )

pheatmap(top25.F.M.avg, 
         cluster_rows = F, 
         cluster_cols = F,
                scale = "row", 
            cellwidth = 50, 
           cellheight = 10,
             gaps_row = c(25,50,75),
             filename = "pheno.cor.Fib.Mye.top25.heatmap.pdf"
)

pheatmap(top25.AF.M.avg, 
         cluster_rows = F, 
         cluster_cols = F,
                scale = "row", 
            cellwidth = 50, 
           cellheight = 10,
             gaps_row = c(25,50,75),
             filename = "pheno.cor.Ath-Fib.Mye.top25.heatmap.pdf"
)

pheatmap(top25.M.avg,
         show_rownames = F,
         cluster_rows = F, 
         cluster_cols = F,
         scale = "row", 
         cellwidth = 15, 
         cellheight = 1.5,
         border_color = NA
)

# Collect top 500 per cluster for pathway analysis
top500.A.M.avg <- rbind(My0.A.M.avg[1:500,], My1.A.M.avg[1:500,])
top500.A.M.avg <- rbind(top500.A.M.avg, My2.A.M.avg[1:500,])
top500.A.M.avg <- rbind(top500.A.M.avg, My3.A.M.avg[1:500,])

top500.F.M.avg <- rbind(My0.F.M.avg[1:500,], My1.F.M.avg[1:500,])
top500.F.M.avg <- rbind(top500.F.M.avg, My2.F.M.avg[1:500,])
top500.F.M.avg <- rbind(top500.F.M.avg, My3.F.M.avg[1:500,])

top500.AF.M.avg <- rbind(My0.AF.M.avg[1:500,], My1.AF.M.avg[1:500,])
top500.AF.M.avg <- rbind(top500.AF.M.avg, My2.AF.M.avg[1:500,])
top500.AF.M.avg <- rbind(top500.AF.M.avg, My3.AF.M.avg[1:500,])

top500.M.avg <- rbind(top500.A.M.avg, top500.F.M.avg)
top500.M.avg <- rbind(top500.M.avg, top500.AF.M.avg)

# Plot as sanity check
pheatmap(top500.M.avg,
         show_rownames = F,
          cluster_rows = F, 
          cluster_cols = F,
                 scale = "row", 
             cellwidth = 15, 
            cellheight = 1,
          border_color = NA,
              filename = "pheno.cor.Mye.top500.heatmap.pdf"
)

# Run pathway analysis
ont.My0.A.M.avg <- get_ontology_no_direction(thegenes = row.names(My0.A.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye0.Ath.pheno.cor")
ont.My1.A.M.avg <- get_ontology_no_direction(thegenes = row.names(My1.A.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye1.Ath.pheno.cor")
ont.My2.A.M.avg <- get_ontology_no_direction(thegenes = row.names(My2.A.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye2.Ath.pheno.cor")
ont.My3.A.M.avg <- get_ontology_no_direction(thegenes = row.names(My3.A.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye3.Ath.pheno.cor")

ont.My0.F.M.avg <- get_ontology_no_direction(thegenes = row.names(My0.F.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye0.Fib.pheno.cor")
ont.My1.F.M.avg <- get_ontology_no_direction(thegenes = row.names(My1.F.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye1.Fib.pheno.cor")
ont.My2.F.M.avg <- get_ontology_no_direction(thegenes = row.names(My2.F.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye2.Fib.pheno.cor")
ont.My3.F.M.avg <- get_ontology_no_direction(thegenes = row.names(My3.F.M.avg[1:500,]), 
                                             return.data = T,
                                             outdir = "pheno.cor", 
                                             full_GSEA = F,
                                             name = "Mye3.Fib.pheno.cor")

ont.My0.AF.M.avg <- get_ontology_no_direction(thegenes = row.names(My0.AF.M.avg[1:500,]), 
                                              return.data = T,
                                              outdir = "pheno.cor", 
                                              full_GSEA = F,
                                              name = "Mye0.Ath-Fib.pheno.cor")
ont.My1.AF.M.avg <- get_ontology_no_direction(thegenes = row.names(My1.AF.M.avg[1:500,]), 
                                              return.data = T,
                                              outdir = "pheno.cor", 
                                              full_GSEA = F,
                                              name = "Mye1.Ath-Fib.pheno.cor")
ont.My2.AF.M.avg <- get_ontology_no_direction(thegenes = row.names(My2.AF.M.avg[1:500,]), 
                                              return.data = T,
                                              outdir = "pheno.cor", 
                                              full_GSEA = F,
                                              name = "Mye2.Ath-Fib.pheno.cor")
ont.My3.AF.M.avg <- get_ontology_no_direction(thegenes = row.names(My3.AF.M.avg[1:500,]),
                                              return.data = T,
                                              outdir = "pheno.cor", 
                                              full_GSEA = F,
                                              name = "Mye3.Ath-Fib.pheno.cor")

# Compare the 3 phenotypes
# Add -log10(padj) column to cluster on
# Pathways
ont.My0.A.M.avg$Pathways$log10padj <- -log10(ont.My0.A.M.avg$Pathways$padj)
ont.My1.A.M.avg$Pathways$log10padj <- -log10(ont.My1.A.M.avg$Pathways$padj)
ont.My2.A.M.avg$Pathways$log10padj <- -log10(ont.My2.A.M.avg$Pathways$padj)
ont.My3.A.M.avg$Pathways$log10padj <- -log10(ont.My3.A.M.avg$Pathways$padj)

ont.My0.F.M.avg$Pathways$log10padj <- -log10(ont.My0.F.M.avg$Pathways$padj)
ont.My1.F.M.avg$Pathways$log10padj <- -log10(ont.My1.F.M.avg$Pathways$padj)
ont.My2.F.M.avg$Pathways$log10padj <- -log10(ont.My2.F.M.avg$Pathways$padj)
ont.My3.F.M.avg$Pathways$log10padj <- -log10(ont.My3.F.M.avg$Pathways$padj)

ont.My0.AF.M.avg$Pathways$log10padj <- -log10(ont.My0.AF.M.avg$Pathways$padj)
ont.My1.AF.M.avg$Pathways$log10padj <- -log10(ont.My1.AF.M.avg$Pathways$padj)
ont.My2.AF.M.avg$Pathways$log10padj <- -log10(ont.My2.AF.M.avg$Pathways$padj)
ont.My3.AF.M.avg$Pathways$log10padj <- -log10(ont.My3.AF.M.avg$Pathways$padj)

# Go terms
ont.My0.A.M.avg$GO_terms$log10padj <- -log10(ont.My0.A.M.avg$GO_terms$padj)
ont.My1.A.M.avg$GO_terms$log10padj <- -log10(ont.My1.A.M.avg$GO_terms$padj)
ont.My2.A.M.avg$GO_terms$log10padj <- -log10(ont.My2.A.M.avg$GO_terms$padj)
ont.My3.A.M.avg$GO_terms$log10padj <- -log10(ont.My3.A.M.avg$GO_terms$padj)

ont.My0.F.M.avg$GO_terms$log10padj <- -log10(ont.My0.F.M.avg$GO_terms$padj)
ont.My1.F.M.avg$GO_terms$log10padj <- -log10(ont.My1.F.M.avg$GO_terms$padj)
ont.My2.F.M.avg$GO_terms$log10padj <- -log10(ont.My2.F.M.avg$GO_terms$padj)
ont.My3.F.M.avg$GO_terms$log10padj <- -log10(ont.My3.F.M.avg$GO_terms$padj)

ont.My0.AF.M.avg$GO_terms$log10padj <- -log10(ont.My0.AF.M.avg$GO_terms$padj)
ont.My1.AF.M.avg$GO_terms$log10padj <- -log10(ont.My1.AF.M.avg$GO_terms$padj)
ont.My2.AF.M.avg$GO_terms$log10padj <- -log10(ont.My2.AF.M.avg$GO_terms$padj)
ont.My3.AF.M.avg$GO_terms$log10padj <- -log10(ont.My3.AF.M.avg$GO_terms$padj)

# Build common tables
# Pathways
ont.My0.pathways <- merge(ont.My0.A.M.avg$Pathways, ont.My0.F.M.avg$Pathways, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My0.pathways <- merge(ont.My0.pathways, ont.My0.AF.M.avg$Pathways, by = "name")
colnames(ont.My0.pathways)[6:7] <- c("padj.AthFib", "log10padj.AthFib")
row.names(ont.My0.pathways) <- ont.My0.pathways$name
ont.My0.pathways <- ont.My0.pathways[,c(3,5,7)]

ont.My1.pathways <- merge(ont.My1.A.M.avg$Pathways, ont.My1.F.M.avg$Pathways, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My1.pathways <- merge(ont.My1.pathways, ont.My1.AF.M.avg$Pathways, by = "name")
colnames(ont.My1.pathways)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My1.pathways) <- ont.My1.pathways$name
ont.My1.pathways <- ont.My1.pathways[,c(3,5,7)]

ont.My2.pathways <- merge(ont.My2.A.M.avg$Pathways, ont.My2.F.M.avg$Pathways, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My2.pathways <- merge(ont.My2.pathways, ont.My2.AF.M.avg$Pathways, by = "name")
colnames(ont.My2.pathways)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My2.pathways) <- ont.My2.pathways$name
ont.My2.pathways <- ont.My2.pathways[,c(3,5,7)]

ont.My3.pathways <- merge(ont.My3.A.M.avg$Pathways, ont.My3.F.M.avg$Pathways, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My3.pathways <- merge(ont.My3.pathways, ont.My3.AF.M.avg$Pathways, by = "name")
colnames(ont.My3.pathways)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My3.pathways) <- ont.My3.pathways$name
ont.My3.pathways <- ont.My3.pathways[,c(3,5,7)]

# GO terms
ont.My0.GO_terms <- merge(ont.My0.A.M.avg$GO_terms, ont.My0.F.M.avg$GO_terms, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My0.GO_terms <- merge(ont.My0.GO_terms, ont.My0.AF.M.avg$GO_terms, by = "name")
colnames(ont.My0.GO_terms)[6:7] <- c("padj.AthFib", "log10padj.AthFib")
row.names(ont.My0.GO_terms) <- ont.My0.GO_terms$name
ont.My0.GO_terms <- ont.My0.GO_terms[,c(3,5,7)]

ont.My1.GO_terms <- merge(ont.My1.A.M.avg$GO_terms, ont.My1.F.M.avg$GO_terms, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My1.GO_terms <- merge(ont.My1.GO_terms, ont.My1.AF.M.avg$GO_terms, by = "name")
colnames(ont.My1.GO_terms)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My1.GO_terms) <- ont.My1.GO_terms$name
ont.My1.GO_terms <- ont.My1.GO_terms[,c(3,5,7)]

ont.My2.GO_terms <- merge(ont.My2.A.M.avg$GO_terms, ont.My2.F.M.avg$GO_terms, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My2.GO_terms <- merge(ont.My2.GO_terms, ont.My2.AF.M.avg$GO_terms, by = "name")
colnames(ont.My2.GO_terms)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My2.GO_terms) <- ont.My2.GO_terms$name
ont.My2.GO_terms <- ont.My2.GO_terms[,c(3,5,7)]

ont.My3.GO_terms <- merge(ont.My3.A.M.avg$GO_terms, ont.My3.F.M.avg$GO_terms, by = "name", suffixes = c(".Ath", ".Fib"))
ont.My3.GO_terms <- merge(ont.My3.GO_terms, ont.My3.AF.M.avg$GO_terms, by = "name")
colnames(ont.My3.GO_terms)[6:7] <- c("padj.AthFib", "log10padj.AthFib") 
row.names(ont.My3.GO_terms) <- ont.My3.GO_terms$name
ont.My3.GO_terms <- ont.My3.GO_terms[,c(3,5,7)]

# Extract top 10 specific pathways per pheno type
# Cluster mye 0
dev.off() # Reset xquartz device after pheatmap's ggsave call or it wont display (stupid bug)
set.seed(100)
k.ont.My0.pathways <- pheatmap(ont.My0.pathways, scale = "row", kmeans_k = 5)
ont.My0.top.A.pathways <- ont.My0.pathways[k.ont.My0.pathways$kmeans$cluster == 2,]
ont.My0.top.A.pathways <- ont.My0.top.A.pathways[order(ont.My0.top.A.pathways$log10padj.Ath, decreasing = T),][1:10,]

ont.My0.top.F.pathways <- ont.My0.pathways[k.ont.My0.pathways$kmeans$cluster == 3,]
ont.My0.top.F.pathways <- ont.My0.top.F.pathways[order(ont.My0.top.F.pathways$log10padj.Fib, decreasing = T),][1:10,]

ont.My0.top.AF.pathways <- ont.My0.pathways[k.ont.My0.pathways$kmeans$cluster == 5,]
ont.My0.top.AF.pathways <- ont.My0.top.AF.pathways[order(ont.My0.top.AF.pathways$log10padj.AthFib, decreasing = T),][1:10,]

ont.My0.top.pathways <- rbind(ont.My0.top.A.pathways,ont.My0.top.F.pathways)
ont.My0.top.pathways <- rbind(ont.My0.top.pathways,ont.My0.top.AF.pathways)
pheatmap(ont.My0.top.pathways,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My0.top.diff.pathways.heatmap.pdf")

# Cluster mye 1
dev.off()
set.seed(100)
k.ont.My1.pathways <- pheatmap(ont.My1.pathways, scale = "row", kmeans_k = 5)
ont.My1.top.A.pathways <- ont.My1.pathways[k.ont.My1.pathways$kmeans$cluster == 5,]
ont.My1.top.A.pathways <- ont.My1.top.A.pathways[order(ont.My1.top.A.pathways$log10padj.Ath, decreasing = T),][1:10,]

ont.My1.top.F.pathways <- ont.My1.pathways[k.ont.My1.pathways$kmeans$cluster == 3,]
ont.My1.top.F.pathways <- ont.My1.top.F.pathways[order(ont.My1.top.F.pathways$log10padj.Fib, decreasing = T),][1:10,]

ont.My1.top.AF.pathways <- ont.My1.pathways[k.ont.My1.pathways$kmeans$cluster == 2,]
ont.My1.top.AF.pathways <- ont.My1.top.AF.pathways[order(ont.My1.top.AF.pathways$log10padj.AthFib, decreasing = T),][1:10,]

ont.My1.top.pathways <- rbind(ont.My1.top.A.pathways,ont.My1.top.F.pathways)
ont.My1.top.pathways <- rbind(ont.My1.top.pathways,ont.My1.top.AF.pathways)
pheatmap(ont.My1.top.pathways,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My1.top.diff.pathways.heatmap.pdf")

# Cluster mye 2
dev.off()
set.seed(100)
k.ont.My2.pathways <- pheatmap(ont.My2.pathways, scale = "row", kmeans_k = 5)
ont.My2.top.A.pathways <- ont.My2.pathways[k.ont.My2.pathways$kmeans$cluster == 4,]
ont.My2.top.A.pathways <- ont.My2.top.A.pathways[order(ont.My2.top.A.pathways$log10padj.Ath, decreasing = T),][1:10,]

ont.My2.top.F.pathways <- ont.My2.pathways[k.ont.My2.pathways$kmeans$cluster == 5,]
ont.My2.top.F.pathways <- ont.My2.top.F.pathways[order(ont.My2.top.F.pathways$log10padj.Fib, decreasing = T),][1:10,]

ont.My2.top.AF.pathways <- ont.My2.pathways[k.ont.My2.pathways$kmeans$cluster == 3,]
ont.My2.top.AF.pathways <- ont.My2.top.AF.pathways[order(ont.My2.top.AF.pathways$log10padj.AthFib, decreasing = T),][1:10,]

ont.My2.top.pathways <- rbind(ont.My2.top.A.pathways,ont.My2.top.F.pathways)
ont.My2.top.pathways <- rbind(ont.My2.top.pathways,ont.My2.top.AF.pathways)
pheatmap(ont.My2.top.pathways,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My2.top.diff.pathways.heatmap.pdf")

# Cluster mye 3
dev.off()
set.seed(100)
k.ont.My3.pathways <- pheatmap(ont.My3.pathways, scale = "row", kmeans_k = 5)
ont.My3.top.A.pathways <- ont.My3.pathways[k.ont.My3.pathways$kmeans$cluster == 2,]
ont.My3.top.A.pathways <- ont.My3.top.A.pathways[order(ont.My3.top.A.pathways$log10padj.Ath, decreasing = T),][1:10,]

ont.My3.top.F.pathways <- ont.My3.pathways[k.ont.My3.pathways$kmeans$cluster == 1,]
ont.My3.top.F.pathways <- ont.My3.top.F.pathways[order(ont.My3.top.F.pathways$log10padj.Fib, decreasing = T),][1:10,]

ont.My3.top.AF.pathways <- ont.My3.pathways[k.ont.My3.pathways$kmeans$cluster == 5 | k.ont.My3.pathways$kmeans$cluster == 3,]
ont.My3.top.AF.pathways <- ont.My3.top.AF.pathways[order(ont.My3.top.AF.pathways$log10padj.AthFib, decreasing = T),][1:10,]

ont.My3.top.pathways <- rbind(ont.My3.top.A.pathways,ont.My3.top.F.pathways)
ont.My3.top.pathways <- rbind(ont.My3.top.pathways,ont.My3.top.AF.pathways)
pheatmap(ont.My3.top.pathways,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My3.top.diff.pathways.heatmap.pdf")


# Extract top 10 specific GO_terms per pheno type
# Cluster mye 0
dev.off() # Reset xquartz device after ggsave call or it wont show properly (stupid bug)
set.seed(100)
k.ont.My0.GO_terms <- pheatmap(ont.My0.GO_terms, scale = "row", kmeans_k = 5)
ont.My0.top.A.GO_terms <- ont.My0.GO_terms[k.ont.My0.GO_terms$kmeans$cluster == 3,]
ont.My0.top.A.GO_terms <- ont.My0.top.A.GO_terms[order(ont.My0.top.A.GO_terms$log10padj.Ath, decreasing = T),][1:10,]

ont.My0.top.F.GO_terms <- ont.My0.GO_terms[k.ont.My0.GO_terms$kmeans$cluster == 1,]
ont.My0.top.F.GO_terms <- ont.My0.top.F.GO_terms[order(ont.My0.top.F.GO_terms$log10padj.Fib, decreasing = T),][1:10,]

ont.My0.top.AF.GO_terms <- ont.My0.GO_terms[k.ont.My0.GO_terms$kmeans$cluster == 2,]
ont.My0.top.AF.GO_terms <- ont.My0.top.AF.GO_terms[order(ont.My0.top.AF.GO_terms$log10padj.AthFib, decreasing = T),][1:10,]

ont.My0.top.GO_terms <- rbind(ont.My0.top.A.GO_terms,ont.My0.top.F.GO_terms)
ont.My0.top.GO_terms <- rbind(ont.My0.top.GO_terms,ont.My0.top.AF.GO_terms)
pheatmap(ont.My0.top.GO_terms,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My0.top.diff.GO_terms.heatmap.pdf")

# Cluster mye 1
dev.off()
set.seed(100)
k.ont.My1.GO_terms <- pheatmap(ont.My1.GO_terms, scale = "row", kmeans_k = 5)
ont.My1.top.A.GO_terms <- ont.My1.GO_terms[k.ont.My1.GO_terms$kmeans$cluster == 1,]
ont.My1.top.A.GO_terms <- ont.My1.top.A.GO_terms[order(ont.My1.top.A.GO_terms$log10padj.Ath, decreasing = T),][1:10,]

ont.My1.top.F.GO_terms <- ont.My1.GO_terms[k.ont.My1.GO_terms$kmeans$cluster == 3 | k.ont.My1.GO_terms$kmeans$cluster == 4,]
ont.My1.top.F.GO_terms <- ont.My1.top.F.GO_terms[order(ont.My1.top.F.GO_terms$log10padj.Fib, decreasing = T),][1:10,]

ont.My1.top.AF.GO_terms <- ont.My1.GO_terms[k.ont.My1.GO_terms$kmeans$cluster == 2,]
ont.My1.top.AF.GO_terms <- ont.My1.top.AF.GO_terms[order(ont.My1.top.AF.GO_terms$log10padj.AthFib, decreasing = T),][1:10,]

ont.My1.top.GO_terms <- rbind(ont.My1.top.A.GO_terms,ont.My1.top.F.GO_terms)
ont.My1.top.GO_terms <- rbind(ont.My1.top.GO_terms,ont.My1.top.AF.GO_terms)
pheatmap(ont.My1.top.GO_terms,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My1.top.diff.GO_terms.heatmap.pdf")

# Cluster mye 2
dev.off()
set.seed(100)
k.ont.My2.GO_terms <- pheatmap(ont.My2.GO_terms, scale = "row", kmeans_k = 5)
ont.My2.top.A.GO_terms <- ont.My2.GO_terms[k.ont.My2.GO_terms$kmeans$cluster == 4,]
ont.My2.top.A.GO_terms <- ont.My2.top.A.GO_terms[order(ont.My2.top.A.GO_terms$log10padj.Ath, decreasing = T),][1:10,]

ont.My2.top.F.GO_terms <- ont.My2.GO_terms[k.ont.My2.GO_terms$kmeans$cluster == 3,]
ont.My2.top.F.GO_terms <- ont.My2.top.F.GO_terms[order(ont.My2.top.F.GO_terms$log10padj.Fib, decreasing = T),][1:10,]

ont.My2.top.AF.GO_terms <- ont.My2.GO_terms[k.ont.My2.GO_terms$kmeans$cluster == 1,]
ont.My2.top.AF.GO_terms <- ont.My2.top.AF.GO_terms[order(ont.My2.top.AF.GO_terms$log10padj.AthFib, decreasing = T),][1:10,]

ont.My2.top.GO_terms <- rbind(ont.My2.top.A.GO_terms,ont.My2.top.F.GO_terms)
ont.My2.top.GO_terms <- rbind(ont.My2.top.GO_terms,ont.My2.top.AF.GO_terms)
pheatmap(ont.My2.top.GO_terms,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My2.top.diff.GO_terms.heatmap.pdf")

# Cluster mye 3
dev.off()
set.seed(100)
k.ont.My3.GO_terms <- pheatmap(ont.My3.GO_terms, scale = "row", kmeans_k = 5)
ont.My3.top.A.GO_terms <- ont.My3.GO_terms[k.ont.My3.GO_terms$kmeans$cluster == 4,]
ont.My3.top.A.GO_terms <- ont.My3.top.A.GO_terms[order(ont.My3.top.A.GO_terms$log10padj.Ath, decreasing = T),][1:10,]

ont.My3.top.F.GO_terms <- ont.My3.GO_terms[k.ont.My3.GO_terms$kmeans$cluster == 5,]
ont.My3.top.F.GO_terms <- ont.My3.top.F.GO_terms[order(ont.My3.top.F.GO_terms$log10padj.Fib, decreasing = T),][1:10,]

ont.My3.top.AF.GO_terms <- ont.My3.GO_terms[k.ont.My3.GO_terms$kmeans$cluster == 2 | k.ont.My3.GO_terms$kmeans$cluster == 3,]
ont.My3.top.AF.GO_terms <- ont.My3.top.AF.GO_terms[order(ont.My3.top.AF.GO_terms$log10padj.AthFib, decreasing = T),][1:10,]

ont.My3.top.GO_terms <- rbind(ont.My3.top.A.GO_terms,ont.My3.top.F.GO_terms)
ont.My3.top.GO_terms <- rbind(ont.My3.top.GO_terms,ont.My3.top.AF.GO_terms)
pheatmap(ont.My3.top.GO_terms,
         scale = "row", 
         cutree_rows = 3, 
         cellwidth = 30, 
         cellheight = 10, 
         border_color = NA, filename = "pheno.cor/My3.top.diff.GO_terms.heatmap.pdf")

# Try clustering the top 10 per pheno type, regardless of clustering on exlusivity
# Cluster mye 0
# Atherogenous
pheatmap(ont.My0.pathways[order(ont.My0.pathways$log10padj.Ath, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My0.top_Ath.pathways.pdf"
         )

# Fibrous
pheatmap(ont.My0.pathways[order(ont.My0.pathways$log10padj.Fib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My0.top_Fib.pathways.pdf"
)         

# Atherogenous-Fibrous
pheatmap(ont.My0.pathways[order(ont.My0.pathways$log10padj.AthFib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000) , filename = "pheno.cor/My0.top_Ath-Fib.pathways.pdf"
)

# Cluster mye 1
# Atherogenous
pheatmap(ont.My1.pathways[order(ont.My1.pathways$log10padj.Ath, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My1.top_Ath.pathways.pdf"
)

# Fibrous
pheatmap(ont.My1.pathways[order(ont.My1.pathways$log10padj.Fib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My1.top_Fib.pathways.pdf"
)         

# Atherogenous-Fibrous
pheatmap(ont.My1.pathways[order(ont.My1.pathways$log10padj.AthFib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My1.top_Ath-Fib.pathways.pdf"
)

# Cluster mye 2
# Atherogenous
pheatmap(ont.My2.pathways[order(ont.My2.pathways$log10padj.Ath, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My2.top_Ath.pathways.pdf"
)

# Fibrous
pheatmap(ont.My2.pathways[order(ont.My2.pathways$log10padj.Fib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My2.top_Fib.pathways.pdf"
)         

# Atherogenous-Fibrous
pheatmap(ont.My2.pathways[order(ont.My2.pathways$log10padj.AthFib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My2.top_Ath-Fib.pathways.pdf"
)

# Cluster mye 3
# Atherogenous
pheatmap(ont.My3.pathways[order(ont.My3.pathways$log10padj.Ath, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My3.top_Ath.pathways.pdf"
)

# Fibrous
pheatmap(ont.My3.pathways[order(ont.My3.pathways$log10padj.Fib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My3.top_Fib.pathways.pdf"
)         

# Atherogenous-Fibrous
pheatmap(ont.My3.pathways[order(ont.My3.pathways$log10padj.AthFib, decreasing = T),][1:10,],
         cluster_rows = F, 
         cluster_cols = F, cellwidth = 30, cellheight = 15, border_color = NA,
         color = colorRampPalette(c("yellow","red"))(1000), filename = "pheno.cor/My3.top_Ath-Fib.pathways.pdf"
)



