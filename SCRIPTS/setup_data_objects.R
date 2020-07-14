#-----------------------------------------------------------------------------------------
#Setup data

#Make directory for results
dir.create("results", showWarnings = FALSE)

# merge suerat objects from list
all.seur.combined <- seurat.object.list[[1]]
for (i in 2:length(x = seurat.object.list)) {
  # only filter genes with too low expression during last merge
  if (i == length(seurat.object.list)){
    all.seur.combined <- MergeSeurat(object1 = all.seur.combined, object2 = seurat.object.list[[i]], do.normalize = F,
                          min.cells = 3)
  } else {
    all.seur.combined <- MergeSeurat(object1 = all.seur.combined, object2 = seurat.object.list[[i]], do.normalize = F)
  }
}
# remove object list to save memory
#rm(seurat.object.list)

dim(all.seur.combined@data)

## Filter cells with too low and too high gene counts
# pick thresholds
VlnPlot(object = all.seur.combined, features.plot = c("nGene"))

all.seur.combined <- FilterCells(object = all.seur.combined,
                                 subset.names = c("nGene"),
                                 low.thresholds = 500,
                                 high.thresholds = 10000
)
dim(all.seur.combined@data)

# Normalize data
all.seur.combined <- NormalizeData(object = all.seur.combined,
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000
)

# scale data (regress out UMI)
all.seur.combined <- ScaleData(object = all.seur.combined, 
                               vars.to.regress = c("nUMI"))

# Find variable genes 
all.seur.combined <- FindVariableGenes(all.seur.combined, do.plot = F)
head(all.seur.combined@var.genes, n = 20)

#Run PCA
all.seur.combined <-RunPCA(all.seur.combined)

#visualize results of PCA
DimPlot(object = all.seur.combined, reduction.use = "pca", group.by = "ID", 
        pt.size = 0.5, do.return = F)
VlnPlot(object = all.seur.combined, features.plot = "PC1", group.by = "ID", 
        do.return = F)
PrintDim(object = all.seur.combined, reduction.type = "pca", dims.print = 1:2, 
         genes.print = 10)
DimHeatmap(object = all.seur.combined, reduction.type = "pca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)
pdf("results/PCs.pdf", width = 10, height = 10)
DimHeatmap(object = all.seur.combined, reduction.type = "PCA", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)
dev.off()

#Align PCA subspaces
all.seur.combined <- AlignSubspace(all.seur.combined, reduction.type = "pca", grouping.var = "ID", 
                                   dims.align = 1:15)

length(x = all.seur.combined@var.genes)
