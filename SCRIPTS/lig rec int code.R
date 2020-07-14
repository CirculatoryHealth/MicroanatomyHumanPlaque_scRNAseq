# Receptor ligand interaction figure building
dir.create("lig.rec.int", showWarnings = F)

# Get interaction score files
lig.rec.score.files.dir <- "Interaction score files"
lig.rec.score.files <- list.files(lig.rec.score.files.dir)

# Fetch interaction scores and stuff 'em in a list
lig.rec.scores <- list()
for (thefile in lig.rec.score.files){
  cat(paste("Importing ", thefile, "...\n", sep = ""))
  lig.rec.scores[[sub('\\.txt$', '', thefile)]] <- read.table(paste(lig.rec.score.files.dir, thefile, sep = "/"), 
                                                              header = T, 
                                                              row.names = 1)
  cat("\t",dim(lig.rec.scores[[sub('\\.txt$', '', thefile)]])[1],"interactions\n")
}
length(lig.rec.scores)

# Add cell type information
for(i in names(lig.rec.scores)){
  lig.rec.scores[[i]][,"ligand.cells"]   <- strsplit(i, "_")[[1]][1]
  lig.rec.scores[[i]][,"receptor.cells"] <- strsplit(i, "_")[[1]][2]
}

# Combine everything in one big data frame
lri.big.table <- rbind(lig.rec.scores[[1]], lig.rec.scores[[2]])
for (i in 3:length(lig.rec.scores)){
  lri.big.table <- rbind(lri.big.table, lig.rec.scores[[i]])
}

# Reset row names because dups will be marked by rbind
row.names(lri.big.table) <- NULL

# Clean up unused columns
lri.big.table$ligands.logical <- NULL
lri.big.table$receptors.logical <- NULL

dim(lri.big.table)
head(lri.big.table)

# Filter out 'mixed' and 'Mye.4' clusters as they are both of dubious quality.
lri.big.table <- lri.big.table[grep("Mixed|My.4", lri.big.table$ligand.cells, invert = T),]
lri.big.table <- lri.big.table[grep("Mixed|My.4", lri.big.table$receptor.cells, invert = T),]
dim(lri.big.table)

# Save HLA pairs to a separate tabel for later
hla.big.table <- lri.big.table[grep("HLA", lri.big.table$ligands, invert = F),]
hla.big.table <- rbind(hla.big.table, lri.big.table[grep("HLA", lri.big.table$receptors, invert = F),])
dim(hla.big.table)

# Remove HLA pairs because they are too generic (actually that's fine now with cellphonedb doing the filtering for us already, so let's keep them in)
# lri.big.table <- lri.big.table[grep("HLA", lri.big.table$ligands, invert = T),]
# lri.big.table <- lri.big.table[grep("HLA", lri.big.table$receptors, invert = T),]
dim(lri.big.table)

# Filter out scores over 1 and sort
score_over_1.lri.big.table <- subset(lri.big.table, interaction.score > 1)
score_over_1.lri.big.table <- score_over_1.lri.big.table[order(score_over_1.lri.big.table$interaction.score, decreasing = T),]
dim(score_over_1.lri.big.table)

# Subset out the sub clusters
sub_clusters_lig.so1.lri.big.table <- list()
for (i in unique(score_over_1.lri.big.table$ligand.cells)){
  sub_clusters_lig.so1.lri.big.table[[i]] <- score_over_1.lri.big.table[grep(i, score_over_1.lri.big.table$ligand.cells),]
}

sub_clusters_rec.so1.lri.big.table <- list()
for (i in unique(score_over_1.lri.big.table$receptor.cells)){
  sub_clusters_rec.so1.lri.big.table[[i]] <- score_over_1.lri.big.table[grep(i, score_over_1.lri.big.table$receptor.cells),]
}

#Define colors for the subclusters
subclus.colors <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#E7298A", # CD79 B
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", # Endo
                    "#35978F", # Mast
                    "#000000", # Mixed
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B", # Myeloid
                    "#8856A7", "#810F7C") # SMC
names(subclus.colors) <- sort(unique(score_over_1.lri.big.table$receptor.cells))

subclus.colors.generic <- subclus.colors
subclus.colors.generic[grep("^CD4", names(subclus.colors.generic))] <- "#FED976"
subclus.colors.generic[grep("^CD8", names(subclus.colors.generic))] <- "#CB181D"
subclus.colors.generic[grep("^CD79", names(subclus.colors.generic))] <- "#FFCCFF"
subclus.colors.generic[grep("^E", names(subclus.colors.generic))] <- "#4292C6"
subclus.colors.generic[grep("^My", names(subclus.colors.generic))] <- "#74C476"
subclus.colors.generic[grep("^KIT", names(subclus.colors.generic))] <- "#CCCCCC"
subclus.colors.generic[grep("^S", names(subclus.colors.generic))] <- "#8856A7"

dir.create("lig.rec.int/ligands_per_population/", showWarnings = F)


# Make overall plots for ALL lig and recs
d <- melt(score_over_1.lri.big.table, 
          id.vars = c("ligands", "receptors", "ligand.cells"), 
          measure.vars = "interaction.score"
)

ggplot(d, aes(x = ligands, y= value, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/3
  )
ggsave(paste("lig.rec.int/all ligands in ligand cells.dots.pdf", sep = ""))

ggplot(d, aes(x = receptors, y= value, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors paired with ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/3
  )
ggsave(paste("lig.rec.int/all receptors paired with ligand cells.dots.pdf", sep = ""))

d <- melt(score_over_1.lri.big.table, 
          id.vars = c("ligands", "receptors", "receptor.cells"), 
          measure.vars = "interaction.score"
)

ggplot(d, aes(x = ligands, y= value, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands paired with receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/3
  )
ggsave(paste("lig.rec.int/all ligands paired with receptor cells.dots.pdf", sep = ""))

ggplot(d, aes(x = receptors, y= value, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/3
  )
ggsave(paste("lig.rec.int/all receptors in receptor cells.dots.pdf", sep = ""))


#Make sub cluster plots
for( i in unique(score_over_1.lri.big.table$ligand.cells)){
  d <- melt(sub_clusters_lig.so1.lri.big.table[[i]], 
            id.vars = c("ligands", "receptors", "receptor.cells"), 
            measure.vars = "interaction.score"
  )
  
  ggplot(d, aes(x = ligands, y= value, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste(i, "expressed ligands")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 2/1.5
       )
  ggsave(paste("lig.rec.int/ligands_per_population/", i, " expressed ligands.dots.pdf", sep = ""))
  
  ggplot(d, aes(x = receptors, y= value, col = receptor.cells)) +
    geom_point(size = 2) +
    scale_fill_manual(values = subclus.colors, aesthetics = "col") +
    ggtitle(paste(i, "linked receptors")) +
    ylab("Interaction score") +
    theme_light() +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 2/1.5
    )
  ggsave(paste("lig.rec.int/ligands_per_population/", i, " linked receptors.dots.pdf", sep = ""))
}


#=============================================================================
# Filter out ligands only expressed in 1 major cell type
CD4.lri.ligs <- unique(score_over_1.lri.big.table[grep("^CD4", score_over_1.lri.big.table$ligand.cells),]$ligands)
CD8.lri.ligs <- unique(score_over_1.lri.big.table[grep("^CD8", score_over_1.lri.big.table$ligand.cells),]$ligands)
T.lri.ligs   <- unique(unlist(list(CD4.lri.ligs, CD8.lri.ligs)))
Mye.lri.ligs <- unique(score_over_1.lri.big.table[grep("^My", score_over_1.lri.big.table$ligand.cells),]$ligands)
E.lri.ligs   <- unique(score_over_1.lri.big.table[grep("^E", score_over_1.lri.big.table$ligand.cells),]$ligands)
SMC.lri.ligs <- unique(score_over_1.lri.big.table[grep("^S", score_over_1.lri.big.table$ligand.cells),]$ligands)
Ma.lri.ligs  <- unique(score_over_1.lri.big.table[grep("^KIT", score_over_1.lri.big.table$ligand.cells),]$ligands)
Mi.lri.ligs  <- unique(score_over_1.lri.big.table[grep("^Mix", score_over_1.lri.big.table$ligand.cells),]$ligands)
B.lri.ligs   <- unique(score_over_1.lri.big.table[grep("^CD79", score_over_1.lri.big.table$ligand.cells),]$ligands)

CD8.uniq.lri.ligs <- CD8.lri.ligs[!CD8.lri.ligs %in% unique(unlist(list(Mye.lri.ligs, E.lri.ligs,   Ma.lri.ligs,  B.lri.ligs)))]
Mye.uniq.lri.ligs <- Mye.lri.ligs[!Mye.lri.ligs %in% unique(unlist(list(CD8.lri.ligs, E.lri.ligs,   Ma.lri.ligs,  B.lri.ligs)))]
E.uniq.lri.ligs   <- E.lri.ligs[    !E.lri.ligs %in% unique(unlist(list(Mye.lri.ligs, CD8.lri.ligs, Ma.lri.ligs,  B.lri.ligs)))]
Ma.uniq.lri.ligs  <- Ma.lri.ligs[  !Ma.lri.ligs %in% unique(unlist(list(Mye.lri.ligs, E.lri.ligs,   CD8.lri.ligs, B.lri.ligs)))]
B.uniq.lri.ligs   <- B.lri.ligs[    !B.lri.ligs %in% unique(unlist(list(Mye.lri.ligs, E.lri.ligs,   Ma.lri.ligs,  CD8.lri.ligs)))]

unique.lri.ligs <- list(CD8.uniq.lri.ligs, Mye.uniq.lri.ligs, E.uniq.lri.ligs, Ma.uniq.lri.ligs, B.uniq.lri.ligs)
names(unique.lri.ligs) = c("CD8","My","E","Ma","B")
dir.create("lig.rec.int/unique_ligands_per_celltype", showWarnings = F)

# Plot per celtype
for (i in names(unique.lri.ligs)){
  for( j in unique(score_over_1.lri.big.table$ligand.cells)[grep(i, unique(score_over_1.lri.big.table$ligand.cells))]){
    d <- melt(sub_clusters_lig.so1.lri.big.table[[j]], 
              id.vars = c("ligands", "receptors", "receptor.cells"), 
              measure.vars = "interaction.score"
    )
    d <- d[d$ligands %in% unique.lri.ligs[[i]],]
    
    ggplot(d, aes(x = ligands, y= value, col = receptor.cells)) +
      geom_point(size = 2) +
      scale_fill_manual(values = subclus.colors, aesthetics = "col") +
      ggtitle(paste(j, "expressed ligands")) +
      ylab("Interaction score") +
      theme_light() +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2/1.5
      )
    ggsave(paste("lig.rec.int/unique_ligands_per_celltype/", j, " expressed ", i, " unique ligands.dots.pdf", sep = ""))
    
    ggplot(d, aes(x = receptors, y= value, col = receptor.cells)) +
      geom_point(size = 2) +
      scale_fill_manual(values = subclus.colors, aesthetics = "col") +
      ggtitle(paste(j, "linked receptors")) +
      ylab("Interaction score") +
      theme_light() +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2/1.5
      )
    ggsave(paste("lig.rec.int/unique_ligands_per_celltype/", j, " linked ", i, " unique  receptors.dots.pdf", sep = ""))
  }
}


#=============================================================================
# Vice versa, filter out ligands whose recptors are only expressed in 1 major cell type
CD4.lri.receptors <- unique(score_over_1.lri.big.table[grep("^CD4", score_over_1.lri.big.table$ligand.cells),]$receptors)
CD8.lri.receptors <- unique(score_over_1.lri.big.table[grep("^CD8", score_over_1.lri.big.table$ligand.cells),]$receptors)
Mye.lri.receptors <- unique(score_over_1.lri.big.table[grep("^My", score_over_1.lri.big.table$ligand.cells),]$receptors)
E.lri.receptors   <- unique(score_over_1.lri.big.table[grep("^E", score_over_1.lri.big.table$ligand.cells),]$receptors)
SMC.lri.receptors <- unique(score_over_1.lri.big.table[grep("^S", score_over_1.lri.big.table$ligand.cells),]$receptors)
Ma.lri.receptors  <- unique(score_over_1.lri.big.table[grep("^KIT", score_over_1.lri.big.table$ligand.cells),]$receptors)
Mi.lri.receptors  <- unique(score_over_1.lri.big.table[grep("^Mix", score_over_1.lri.big.table$ligand.cells),]$receptors)
B.lri.receptors   <- unique(score_over_1.lri.big.table[grep("^CD79", score_over_1.lri.big.table$ligand.cells),]$receptors)

CD8.uniq.lri.receptors <- CD8.lri.receptors[!CD8.lri.receptors %in% unique(unlist(list(CD4.lri.receptors,
                                                                                       Mye.lri.receptors,
                                                                                       SMC.lri.receptors,
                                                                                       Ma.lri.receptors)))]
Mye.uniq.lri.receptors <- Mye.lri.receptors[!Mye.lri.receptors %in% unique(unlist(list(CD4.lri.receptors,
                                                                                       CD8.lri.receptors,
                                                                                       SMC.lri.receptors,
                                                                                       Ma.lri.receptors)))]
SMC.uniq.lri.receptors <- SMC.lri.receptors[!SMC.lri.receptors %in% unique(unlist(list(CD4.lri.receptors,
                                                                                       Mye.lri.receptors,
                                                                                       CD8.lri.receptors,   
                                                                                       Ma.lri.receptors)))]
Ma.uniq.lri.receptors <- Ma.lri.receptors[!Ma.lri.receptors %in% unique(unlist(list(CD4.lri.receptors,
                                                                                       Mye.lri.receptors,
                                                                                       CD8.lri.receptors,   
                                                                                       SMC.lri.receptors)))]

unique.lri.receptors <- list(CD8.uniq.lri.receptors, Mye.uniq.lri.receptors, SMC.uniq.lri.receptors, Ma.uniq.lri.receptors)
names(unique.lri.receptors) = c("CD8","My","SMC","Ma")
dir.create("lig.rec.int/unique_receptors_per_celltype", showWarnings = F)

# Plot per celtype
for (i in names(unique.lri.receptors)){
  for( j in unique(score_over_1.lri.big.table$receptor.cells)[grep(i, unique(score_over_1.lri.big.table$receptor.cells))]){
    d <- melt(sub_clusters_lig.so1.lri.big.table[[j]], 
              id.vars = c("ligands", "receptors", "receptor.cells","ligand.cells"), 
              measure.vars = "interaction.score"
    )
    d <- d[d$receptors %in% unique.lri.receptors[[i]],]
    
    ggplot(d, aes(x = ligands, y= value, col = receptor.cells)) +
      geom_point(size = 2) +
      scale_fill_manual(values = subclus.colors, aesthetics = "col") +
      ggtitle(paste(j, "linked ligands")) +
      ylab("Interaction score") +
      theme_light() +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2/1.5
      )
    ggsave(paste("lig.rec.int/unique_receptors_per_celltype/", j, " linked ligands ", i, " unique receptors.dots.pdf", sep = ""))
    
    ggplot(d, aes(x = receptors, y= value, col = receptor.cells)) +
      geom_point(size = 2) +
      scale_fill_manual(values = subclus.colors, aesthetics = "col") +
      ggtitle(paste(j, "expressed receptors")) +
      ylab("Interaction score") +
      theme_light() +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 2/1.5
      )
    ggsave(paste("lig.rec.int/unique_receptors_per_celltype/", j, " expressed ", i, " unique receptors.dots.pdf", sep = ""))
  }
}


# Merge ligands and receptors in one plot and print per subcluster
# Plot per ligand cell
for( i in unique(score_over_1.lri.big.table$ligand.cells)){
  d <- melt(sub_clusters_lig.so1.lri.big.table[[i]], 
            id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
            measure.vars = "interaction.score"
  )
  d$variable <- NULL
  colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")
  
  d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
  d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
  d <- d[order(d$interaction.score, decreasing = T),]
  
  ggplot(d, aes(x = ligands, y= receptors, col = receptor.cells, size = interaction.score)) +
    geom_point() +
    scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
    ggtitle(paste(i, "expressed ligands")) +
    theme_light() +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste("lig.rec.int/ligands_per_population/Population ", i, " ligands.dots.pdf", sep = ""))
}

# PLot per receptor cell
for( i in unique(score_over_1.lri.big.table$receptor.cells)){
  d <- melt(sub_clusters_rec.so1.lri.big.table[[i]], 
            id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
            measure.vars = "interaction.score"
  )
  d$variable <- NULL
  colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")
  
  d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
  d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
  d <- d[order(d$interaction.score, decreasing = T),]
  
  ggplot(d, aes(x = ligands, y= receptors, col = ligand.cells, size = interaction.score)) +
    geom_point() +
    scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
    ggtitle(paste(i, "expressed receptors")) +
    theme_light() +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 16),
          aspect.ratio = 1)
  ggsave(paste("lig.rec.int/ligands_per_population/Population ", i, " receptors.dots.pdf", sep = ""))
}


# direct myeloid comparison
# Receptor based
d <- rbind(melt(sub_clusters_rec.so1.lri.big.table[["My.0"]], 
          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
          measure.vars = "interaction.score"), melt(sub_clusters_rec.so1.lri.big.table[["My.1"]], 
                                                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                   measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_rec.so1.lri.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_rec.so1.lri.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

library(clipr)
write_clip(as.character(unique(d[d$receptor.cells == "My.0","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.1","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.2","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.3","receptors"])))

ggplot(d, aes(x = ligands, y= receptors, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid expressed receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations receptors.dots.plus_HLA.pdf", sep = ""))

# Ligand based
d <- rbind(melt(sub_clusters_lig.so1.lri.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_lig.so1.lri.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_lig.so1.lri.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_lig.so1.lri.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

library(clipr)
write_clip(as.character(unique(d[d$ligand.cells == "My.0","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.1","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.2","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.3","ligands"])))

ggplot(d, aes(x = ligands, y= receptors, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid expressed ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations ligands.dots.plus_HLA.pdf", sep = ""))

# Only plot interactions unique to one subtype
# Receptor based
d <- rbind(melt(sub_clusters_rec.so1.lri.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_rec.so1.lri.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_rec.so1.lri.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_rec.so1.lri.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

d$interaction <- paste(d$ligands, d$receptors, sep = ".")
d <- d[!duplicated(d$interaction),]

ggplot(d, aes(x = ligands, y= receptors, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors uniquely expressed")) +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations unique receptors.dots.plus_HLA.pdf", sep = ""),width = 20, height =  20)

# Ligand based
d <- rbind(melt(sub_clusters_lig.so1.lri.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_lig.so1.lri.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_lig.so1.lri.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_lig.so1.lri.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

d$interaction <- paste(d$ligands, d$receptors, sep = ".")
d <- d[!duplicated(d$interaction),]

ggplot(d, aes(x = ligands, y= receptors, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands uniquely expressed")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations unique ligands.dots.plus_HLA.pdf", sep = ""), width = 20, height = 20)



#==========================================
# Plot the HLA interactions in the myeloid fraction

# Filter out scores over 1 and sort
score_over_1.hla.big.table <- subset(hla.big.table, interaction.score > 1)
score_over_1.hla.big.table <- score_over_1.hla.big.table[order(score_over_1.hla.big.table$interaction.score, decreasing = T),]
dim(score_over_1.hla.big.table)

dim(hla.big.table[grep("HLA", hla.big.table$receptors, invert = F),])

# Subset out the sub clusters
sub_clusters_lig.so1.hla.big.table <- list()
for (i in unique(score_over_1.hla.big.table$ligand.cells)){
  sub_clusters_lig.so1.hla.big.table[[i]] <- score_over_1.hla.big.table[grep(i, score_over_1.hla.big.table$ligand.cells),]
}

sub_clusters_rec.so1.hla.big.table <- list()
for (i in unique(score_over_1.hla.big.table$receptor.cells)){
  sub_clusters_rec.so1.hla.big.table[[i]] <- score_over_1.hla.big.table[grep(i, score_over_1.hla.big.table$receptor.cells),]
}

# direct myeloid comparison
# Receptor based
d <- rbind(melt(sub_clusters_rec.so1.hla.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_rec.so1.hla.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_rec.so1.hla.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_rec.so1.hla.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

library(clipr)
write_clip(as.character(unique(d[d$receptor.cells == "My.0","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.1","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.2","receptors"])))
write_clip(as.character(unique(d[d$receptor.cells == "My.3","receptors"])))

ggplot(d, aes(x = ligands, y= receptors, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid expressed receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations HLA receptors.dots.pdf", sep = ""))

# Ligand based
d <- rbind(melt(sub_clusters_lig.so1.hla.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_lig.so1.hla.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_lig.so1.hla.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_lig.so1.hla.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

library(clipr)
write_clip(as.character(unique(d[d$ligand.cells == "My.0","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.1","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.2","ligands"])))
write_clip(as.character(unique(d[d$ligand.cells == "My.3","ligands"])))

ggplot(d, aes(x = ligands, y= receptors, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid expressed ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations HLA ligands.dots.pdf", sep = ""))

# Only plot interactions unique to one subtype
# Receptor based
d <- rbind(melt(sub_clusters_rec.so1.hla.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_rec.so1.hla.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_rec.so1.hla.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_rec.so1.hla.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

d$interaction <- paste(d$ligands, d$receptors, sep = ".")
d <- d[!duplicated(d$interaction),]

ggplot(d, aes(x = ligands, y= receptors, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors uniquely expressed")) +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations unique HLA receptors.dots.pdf", sep = ""))

# Ligand based
d <- rbind(melt(sub_clusters_lig.so1.hla.big.table[["My.0"]], 
                id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                measure.vars = "interaction.score"), melt(sub_clusters_lig.so1.hla.big.table[["My.1"]], 
                                                          id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                                                          measure.vars = "interaction.score")
)
d <- rbind(d, melt(sub_clusters_lig.so1.hla.big.table[["My.2"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))
d <- rbind(d, melt(sub_clusters_lig.so1.hla.big.table[["My.3"]], 
                   id.vars = c("ligands", "receptors", "ligand.cells", "receptor.cells"), 
                   measure.vars = "interaction.score"))

d$variable <- NULL
colnames(d) <- c("ligands", "receptors", "ligand.cells", "receptor.cells", "interaction.score")

d$ligands <- factor(d$ligands, levels = unique(sort(d$ligands)))
d$receptors <- factor(d$receptors, levels = unique(sort(d$receptors)))
d <- d[order(d$interaction.score, decreasing = T),]

d$interaction <- paste(d$ligands, d$receptors, sep = ".")
d <- d[!duplicated(d$interaction),]

ggplot(d, aes(x = ligands, y= receptors, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands uniquely expressed")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave(paste("lig.rec.int/ligands_per_population/Myeloid populations unique HLA ligands.dots.pdf", sep = ""))
