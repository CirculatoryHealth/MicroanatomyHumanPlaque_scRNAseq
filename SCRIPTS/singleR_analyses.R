#Load the package
library(SingleR)

#Make object
#Lotte's fix
singler <- CreateSinglerObject(counts = all.seur.combined@raw.data[,rownames(all.seur.combined@meta.data)],       # somehow this eliminates one issue
                               annot = all.seur.combined@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "Single Cell Plaques",
                               min.genes = 0,                          # if starting from seurat, set min.genes to 0
                               technology = "CEL-seq2",
                               species = "Human",
                               citation = "",
                               ref.list = list(),
                               normalize.gene.length = FALSE,
                               variable.genes = "de",
                               fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
                               do.signatures = TRUE,
                               clusters = all.seur.combined@ident,
                               do.main.types = TRUE,
                               reduce.file.size = TRUE,
                               numCores = 4
)

#So let's try constructing the object manually
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                       all.seur.combined@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                            all.seur.combined@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                       all.seur.combined@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                            all.seur.combined@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined@data)
singler$meta.data$project.name = "Single Cell Plaques"
singler$meta.data$orig.ident = all.seur.combined@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- all.seur.combined

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 80)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = F, top.n = 80)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T,top.n = 50)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = F,top.n = 50)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 50)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters),order.by.clusters = F, top.n = 50)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 50)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = F, top.n = 50)


#Now look only at the CD4 T-cells
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.CD4_clusters@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined.CD4_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.CD4_clusters@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined.CD4_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.CD4_clusters@data, 
                                                     hpca$data, 
                                                     hpca$types, 
                                                     clusters = all.seur.combined.CD4_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.CD4_clusters@data, 
                                                          hpca$data, 
                                                          hpca$main_types, 
                                                          clusters = all.seur.combined.CD4_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.CD4_clusters@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined.CD4_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.CD4_clusters@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined.CD4_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.CD4_clusters@data, 
                                                     blueprint_encode$data, 
                                                     blueprint_encode$types, 
                                                     clusters = all.seur.combined.CD4_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.CD4_clusters@data, 
                                                          blueprint_encode$data, 
                                                          blueprint_encode$main_types, 
                                                          clusters = all.seur.combined.CD4_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined.CD4_clusters@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined.CD4_clusters@data)
singler$meta.data$project.name = "Single Cell Plaques CD4 subset"
singler$meta.data$orig.ident = all.seur.combined.CD4_clusters@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- all.seur.combined.CD4_clusters

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined.CD4_clusters@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined.CD4_clusters@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)


#Now look only at the myeloid cells
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.M_clusters@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined.M_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.M_clusters@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined.M_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.M_clusters@data, 
                                                     hpca$data, 
                                                     hpca$types, 
                                                     clusters = all.seur.combined.M_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.M_clusters@data, 
                                                          hpca$data, 
                                                          hpca$main_types, 
                                                          clusters = all.seur.combined.M_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.M_clusters@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined.M_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.M_clusters@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined.M_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.M_clusters@data, 
                                                     blueprint_encode$data, 
                                                     blueprint_encode$types, 
                                                     clusters = all.seur.combined.M_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.M_clusters@data, 
                                                          blueprint_encode$data, 
                                                          blueprint_encode$main_types, 
                                                          clusters = all.seur.combined.M_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined.M_clusters@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined.M_clusters@data)
singler$meta.data$project.name = "Single Cell Plaques Myeloid subset"
singler$meta.data$orig.ident = all.seur.combined.M_clusters@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- all.seur.combined.M_clusters

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined.M_clusters@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined.M_clusters@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)

SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)

#Now look only at the T-cells
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.T_clusters@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined.T_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.T_clusters@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined.T_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.T_clusters@data, 
                                                     hpca$data, 
                                                     hpca$types, 
                                                     clusters = all.seur.combined.T_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.T_clusters@data, 
                                                          hpca$data, 
                                                          hpca$main_types, 
                                                          clusters = all.seur.combined.T_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.T_clusters@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined.T_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.T_clusters@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined.T_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.T_clusters@data, 
                                                     blueprint_encode$data, 
                                                     blueprint_encode$types, 
                                                     clusters = all.seur.combined.T_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.T_clusters@data, 
                                                          blueprint_encode$data, 
                                                          blueprint_encode$main_types, 
                                                          clusters = all.seur.combined.T_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined.T_clusters@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined.T_clusters@data)
singler$meta.data$project.name = "Single Cell Plaques CD4 subset"
singler$meta.data$orig.ident = all.seur.combined.T_clusters@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- all.seur.combined.T_clusters

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined.T_clusters@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined.T_clusters@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)


#Now look only at cluster2
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.cluster_2@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = all.seur.combined.cluster_2@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.cluster_2@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = all.seur.combined.cluster_2@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.cluster_2@data, 
                                                     hpca$data, 
                                                     hpca$types, 
                                                     clusters = all.seur.combined.cluster_2@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.cluster_2@data, 
                                                          hpca$data, 
                                                          hpca$main_types, 
                                                          clusters = all.seur.combined.cluster_2@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.cluster_2@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = all.seur.combined.cluster_2@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.cluster_2@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = all.seur.combined.cluster_2@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.cluster_2@data, 
                                                     blueprint_encode$data, 
                                                     blueprint_encode$types, 
                                                     clusters = all.seur.combined.cluster_2@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.cluster_2@data, 
                                                          blueprint_encode$data, 
                                                          blueprint_encode$main_types, 
                                                          clusters = all.seur.combined.cluster_2@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined.cluster_2@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined.cluster_2@data)
singler$meta.data$project.name = "Single Cell Plaques cluster2 subset"
singler$meta.data$orig.ident = all.seur.combined.cluster_2@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- all.seur.combined.cluster_2

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined.cluster_2@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined.cluster_2@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)


#Now look only at cluster2 without T cells
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       noCD3_in_cluster2@data, 
                                                       hpca$data, 
                                                       hpca$types, 
                                                       clusters = noCD3_in_cluster2@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            noCD3_in_cluster2@data, 
                                                            hpca$data, 
                                                            hpca$main_types, 
                                                            clusters = noCD3_in_cluster2@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     noCD3_in_cluster2@data, 
                                                     hpca$data, 
                                                     hpca$types, 
                                                     clusters = noCD3_in_cluster2@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          noCD3_in_cluster2@data, 
                                                          hpca$data, 
                                                          hpca$main_types, 
                                                          clusters = noCD3_in_cluster2@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)),
                       list(SingleR.clusters = SingleR(method = "cluster",
                                                       noCD3_in_cluster2@data, 
                                                       blueprint_encode$data, 
                                                       blueprint_encode$types, 
                                                       clusters = noCD3_in_cluster2@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            noCD3_in_cluster2@data, 
                                                            blueprint_encode$data, 
                                                            blueprint_encode$main_types, 
                                                            clusters = noCD3_in_cluster2@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     noCD3_in_cluster2@data, 
                                                     blueprint_encode$data, 
                                                     blueprint_encode$types, 
                                                     clusters = noCD3_in_cluster2@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          noCD3_in_cluster2@data, 
                                                          blueprint_encode$data, 
                                                          blueprint_encode$main_types, 
                                                          clusters = noCD3_in_cluster2@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)))

#Manual object: Add signatures
singler$signatures = calculateSingScores(noCD3_in_cluster2@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(noCD3_in_cluster2@data)
singler$meta.data$project.name = "Single Cell Plaques noCD3 in cluster2 subset"
singler$meta.data$orig.ident = noCD3_in_cluster2@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData = "HPCA"
singler$singler[[2]]$about$organism = "Human"
singler$singler[[2]]$about$Technology = "Cel-seq2"
singler$singler[[2]]$about$Citation = ""
singler$singler[[2]]$about$RefData <- "Blueprint_Encode"
singler$seurat <- noCD3_in_cluster2

#Add tSNE coordinates
singler$meta.data$xy = noCD3_in_cluster2@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = noCD3_in_cluster2@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#Unfortunatly, object thusly created does not work well in the f***ing browser tool
#So let's just try some plots right here
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[2]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T, top.n = 10)



#Load macrophage subtype reference set
load(file = "~/Desktop/AMC/MfTypes.singleR_referenceSet.RData")

#And look only at the myeloid cells
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       all.seur.combined.M_clusters@data, 
                                                       MfTypes.ref$data, 
                                                       MfTypes.ref$types, 
                                                       clusters = all.seur.combined.M_clusters@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            all.seur.combined.M_clusters@data, 
                                                            MfTypes.ref$data, 
                                                            MfTypes.ref$main_types, 
                                                            clusters = all.seur.combined.M_clusters@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     all.seur.combined.M_clusters@data, 
                                                     MfTypes.ref$data, 
                                                     MfTypes.ref$types, 
                                                     clusters = all.seur.combined.M_clusters@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          all.seur.combined.M_clusters@data, 
                                                          MfTypes.ref$data, 
                                                          MfTypes.ref$main_types, 
                                                          clusters = all.seur.combined.M_clusters@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)
                            )
                       )

#Manual object: Add signatures
singler$signatures = calculateSingScores(all.seur.combined.M_clusters@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(all.seur.combined.M_clusters@data)
singler$meta.data$project.name = "Single Cell Plaques Myeloid subset Mf refdata"
singler$meta.data$orig.ident = all.seur.combined.M_clusters@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData <- "Mf Types"
singler$seurat <- all.seur.combined.M_clusters

#Add tSNE coordinates
singler$meta.data$xy = all.seur.combined.M_clusters@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = all.seur.combined.M_clusters@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#And plot
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T,top.n = 50)

#Let's also check cluster 2 (no CD3) with this ref set
singler$singler = list(list(SingleR.clusters = SingleR(method = "cluster",
                                                       noCD3_in_cluster2@data, 
                                                       MfTypes.ref$data, 
                                                       MfTypes.ref$types, 
                                                       clusters = noCD3_in_cluster2@ident,
                                                       genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                       fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                       do.pvals = T, numCores = 4),
                            SingleR.clusters.main = SingleR(method = "cluster",
                                                            noCD3_in_cluster2@data, 
                                                            MfTypes.ref$data, 
                                                            MfTypes.ref$main_types, 
                                                            clusters = noCD3_in_cluster2@ident,
                                                            genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                            fine.tune = T, fine.tune.thres = 0.05, sd.thres = 1,
                                                            do.pvals = T, numCores = 4),
                            SingleR.single = SingleR(method = "single",
                                                     noCD3_in_cluster2@data, 
                                                     MfTypes.ref$data, 
                                                     MfTypes.ref$types, 
                                                     clusters = noCD3_in_cluster2@ident,
                                                     genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                     fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                     do.pvals = T, numCores = 4),
                            SingleR.single.main = SingleR(method = "single",
                                                          noCD3_in_cluster2@data, 
                                                          MfTypes.ref$data, 
                                                          MfTypes.ref$main_types, 
                                                          clusters = noCD3_in_cluster2@ident,
                                                          genes = "de", quantile.use = 0.8, p.threshold = 0.05,
                                                          fine.tune = F, fine.tune.thres = 0.05, sd.thres = 1,
                                                          do.pvals = T, numCores = 4)
)
)

#Manual object: Add signatures
singler$signatures = calculateSingScores(noCD3_in_cluster2@data, species = "Human")
#MAnual object: Add some more info
singler$other = SingleR.CreateKangAnnotations(noCD3_in_cluster2@data)
singler$meta.data$project.name = "Single Cell Plaques cluster 2 Mf refdata"
singler$meta.data$orig.ident = noCD3_in_cluster2@meta.data$orig.ident
singler$singler[[1]]$about$organism = "Human"
singler$singler[[1]]$about$Technology = "Cel-seq2"
singler$singler[[1]]$about$Citation = ""
singler$singler[[1]]$about$RefData <- "Mf Types"
singler$seurat <- noCD3_in_cluster2

#Add tSNE coordinates
singler$meta.data$xy = noCD3_in_cluster2@dr$tsne@cell.embeddings # the tSNE coordinates
#Add cluster identities
singler$meta.data$clusters = noCD3_in_cluster2@ident


#Convert to new object type for export to browser and save
singler.new <- convertSingleR2Browser(singler)
saveRDS(singler.new, file = paste0(singler.new@project.name,'.rds'))

#And plot
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters, clusters = levels(singler$meta.data$clusters), order.by.clusters = T, top.n = 10)
SingleR::SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$clusters),order.by.clusters = T,top.n = 50)


