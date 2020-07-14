#===========================================================================================================================
#===========================================================================================================================
# Export for cellphonedb
# Get sub cluster idents
sub.idents <- paste0("E.", Idents(v3.all.seur.combined.E_clusters))
sub.idents <- c(sub.idents, paste0("My.", Idents(v3.all.seur.combined.M_clusters)))
sub.idents <- c(sub.idents, paste0("CD8.", Idents(all.seur.combined.T_split_clusters.CD8)))
sub.idents <- c(sub.idents, paste0("CD4.", Idents(all.seur.combined.T_split_clusters.CD4)))
sub.idents <- c(sub.idents, paste0("SMC.", Idents(v3.all.seur.combined.SMC_clusters)))
sub.idents <- c(sub.idents, as.character(Idents(v3.all.seur.combined)[grep("KIT", Idents(v3.all.seur.combined))]))
sub.idents <- c(sub.idents, as.character(Idents(v3.all.seur.combined)[grep("CD79", Idents(v3.all.seur.combined))]))
unique(sub.idents)
length(sub.idents)

# Get cell names
names.idents <- names(Idents(v3.all.seur.combined.E_clusters))
names.idents <- c(names.idents, names(Idents(v3.all.seur.combined.M_clusters)))
names.idents <- c(names.idents, names(Idents(all.seur.combined.T_split_clusters.CD8)))
names.idents <- c(names.idents, names(Idents(all.seur.combined.T_split_clusters.CD4)))
names.idents <- c(names.idents, names(Idents(v3.all.seur.combined.SMC_clusters)))
names.idents <- c(names.idents, names(Idents(v3.all.seur.combined)[grep("KIT",  Idents(v3.all.seur.combined))]))
names.idents <- c(names.idents, names(Idents(v3.all.seur.combined)[grep("CD79", Idents(v3.all.seur.combined))]))
length(names.idents)

# Construct metadata table
cellphone.meta <- data.frame(Cell = names.idents, cell_type = sub.idents)

# Get counts
cellphone.counts <- all.seur.combined@raw.data

# Normalize
cellphone.counts <- apply(cellphone.counts, 2, function(x) (x/sum(x))*10000)

# Remove unwanted clusters (mixed and My.4) (by simply taking the cells we have in our metadata frame)
cellphone.counts <- cellphone.counts[,colnames(cellphone.counts) %in% cellphone.meta$Cell]
sum(!cellphone.meta$Cell %in% colnames(cellphone.counts))
sum(!colnames(cellphone.counts) %in% cellphone.meta$Cell)

# Write the tables to disk for use with cellphonedb (python)
write.table(cellphone.counts, file = "lig.rec.int/cellphonedb/counts.txt", quote = F, sep = "\t")
write.table(cellphone.meta,   file = "lig.rec.int/cellphonedb/meta.txt",   quote = F, sep = "\t", row.names = F)



#===========================================================================================================================
#===========================================================================================================================
dir.create("lig.rec.int/cellphonedb_results", showWarnings = F)

# Import results from cellphonedb
# ### CAUTION: As a prep step, I changed the column names in 'significant_means.txt' to be compliant with R. i.e. I removed | and + chars.
cpdb.means <- read.table(file = "lig.rec.int/cellphonedb/outs/hum_plaq_new/significant_means.txt", row.names = 1, header = T, sep = "\t")
cpdb.means[is.na(cpdb.means)] <- 0

# Remove integrin interactions
cpdb.means <- cpdb.means[cpdb.means$is_integrin == "False",]

# Save all HLA interactions for later
hla.cpdb.means <-cpdb.means[grep("HLA",cpdb.means$interacting_pair),]

# Keep only top ranked interactions
summary(cpdb.means$rank)
cpdb.means <- cpdb.means[cpdb.means$rank <mean(cpdb.means$rank),]

#===========================================================================================================================
# Plot overall interactions
# Shape the data for plotting
s.cpdb.means <- cpdb.means

# Prep ligand and receptor name columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(s.cpdb.means$interacting_pair), "_")){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
s.cpdb.means$ligand <- l
s.cpdb.means$receptor <- r

# Melt the 'frame using ligand and receptor info and the interaction scores
m.s.cpdb.means <- melt(s.cpdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:length(colnames(s.cpdb.means))), variable.name = "cells", value.name = "interaction.score")

# Subset pairs with significant interactions (>0)
exp.m.s.cpdb.means <- m.s.cpdb.means[m.s.cpdb.means$interaction.score > 0,]
dim(exp.m.s.cpdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.s.cpdb.means$cells), split = "_", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.s.cpdb.means$ligand.cells <- l
exp.m.s.cpdb.means$receptor.cells <- r

# Clean up
exp.m.s.cpdb.means$interaction.score <- as.numeric(exp.m.s.cpdb.means$interaction.score)
exp.m.s.cpdb.means <- exp.m.s.cpdb.means[!is.na(exp.m.s.cpdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
exp.m.s.cpdb.means <- exp.m.s.cpdb.means[grep("My.4|Mixed", exp.m.s.cpdb.means$receptor.cells, invert = T),]
exp.m.s.cpdb.means <- exp.m.s.cpdb.means[grep("My.4|Mixed", exp.m.s.cpdb.means$ligand.cells, invert = T),]

# Fix factors
exp.m.s.cpdb.means$receptor       <- droplevels(as.factor(exp.m.s.cpdb.means$receptor))
exp.m.s.cpdb.means$ligand         <- droplevels(as.factor(exp.m.s.cpdb.means$ligand))
exp.m.s.cpdb.means$receptor.cells <- droplevels(as.factor(exp.m.s.cpdb.means$receptor.cells))
exp.m.s.cpdb.means$ligand.cells   <- droplevels(as.factor(exp.m.s.cpdb.means$ligand.cells))
length(levels(exp.m.s.cpdb.means$receptor.cells))

#Define colors for the subclusters
subclus.colors <- c("#E7298A", # CD79 B
                    "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", "#3366FF", "#0000FF", # Endo
                    "#35978F", # Mast
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", # Myeloid
                    "#8856A7", "#810F7C", "#663399") # SMC
names(subclus.colors) <- levels(exp.m.s.cpdb.means$ligand.cells)

# Subset a bit for clarity
summary(exp.m.s.cpdb.means$interaction.score)
sub.exp.m.s.cpdb.means <- exp.m.s.cpdb.means[exp.m.s.cpdb.means$interaction.score > summary(exp.m.s.cpdb.means$interaction.score)[5],]
sub.exp.m.s.cpdb.means <- exp.m.s.cpdb.means[exp.m.s.cpdb.means$interaction.score < 20,]

# plot
ggplot(sub.exp.m.s.cpdb.means, aes(x = ligand, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("lig.rec.int/cellphonedb_results/All.ligands_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.means, aes(x = receptor, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("lig.rec.int/cellphonedb_results/All.receptors_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.means, aes(x = ligand, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("lig.rec.int/cellphonedb_results/All.ligands_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.means, aes(x = receptor, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("lig.rec.int/cellphonedb_results/All.receptors_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

#===========================================================================================================================
# Subset Myeloid ligand interactions
Mye.lig.cpdb.means <- cpdb.means[,c(1:11, grep("My.[0-9].", colnames(cpdb.means)))]

# Subset significant interaction pairs
sig.Mye.lig.cpdb.means <- Mye.lig.cpdb.means[apply(Mye.lig.cpdb.means[,12:length(colnames(Mye.lig.cpdb.means))], 1, sum) != 0,]
dim(sig.Mye.lig.cpdb.means)

# Shape the data for plotting
# Prep ligand and receptor name columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(sig.Mye.lig.cpdb.means$interacting_pair), "_")){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
sig.Mye.lig.cpdb.means$ligand <- l
sig.Mye.lig.cpdb.means$receptor <- r

# Melt the 'frame using ligand and receptor info and the interaction scores
m.sig.Mye.lig.cpdb.means <- melt(sig.Mye.lig.cpdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:length(colnames(sig.Mye.lig.cpdb.means))), variable.name = "cells", value.name = "interaction.score")

# Subset pairs with significant interactions (>0)
exp.m.sig.Mye.lig.cpdb.means <- m.sig.Mye.lig.cpdb.means[m.sig.Mye.lig.cpdb.means$interaction.score > 0,]
dim(exp.m.sig.Mye.lig.cpdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.sig.Mye.lig.cpdb.means$cells), split = "_", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.sig.Mye.lig.cpdb.means$ligand.cells <- l
exp.m.sig.Mye.lig.cpdb.means$receptor.cells <- r

# Clean up
exp.m.sig.Mye.lig.cpdb.means$interaction.score <- as.numeric(exp.m.sig.Mye.lig.cpdb.means$interaction.score)
exp.m.sig.Mye.lig.cpdb.means <- exp.m.sig.Mye.lig.cpdb.means[!is.na(exp.m.sig.Mye.lig.cpdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
exp.m.sig.Mye.lig.cpdb.means <- exp.m.sig.Mye.lig.cpdb.means[grep("My.4|Mixed", exp.m.sig.Mye.lig.cpdb.means$receptor.cells, invert = T),]
exp.m.sig.Mye.lig.cpdb.means <- exp.m.sig.Mye.lig.cpdb.means[grep("My.4", exp.m.sig.Mye.lig.cpdb.means$ligand.cells, invert = T),]

# Fix factors
exp.m.sig.Mye.lig.cpdb.means$receptor       <- droplevels(as.factor(exp.m.sig.Mye.lig.cpdb.means$receptor))
exp.m.sig.Mye.lig.cpdb.means$ligand         <- droplevels(as.factor(exp.m.sig.Mye.lig.cpdb.means$ligand))
exp.m.sig.Mye.lig.cpdb.means$receptor.cells <- droplevels(as.factor(exp.m.sig.Mye.lig.cpdb.means$receptor.cells))
exp.m.sig.Mye.lig.cpdb.means$ligand.cells   <- droplevels(as.factor(exp.m.sig.Mye.lig.cpdb.means$ligand.cells))
length(levels(exp.m.sig.Mye.lig.cpdb.means$receptor.cells))

#Define colors for the subclusters
subclus.colors <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#E7298A", # CD79 B
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", "#3366FF", "#0000FF", # Endo
                    "#35978F", # Mast
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", # Myeloid
                    "#8856A7", "#810F7C", "#663399") # SMC
names(subclus.colors) <- levels(exp.m.sig.Mye.lig.cpdb.means$receptor.cells)

subclus.colors.generic <- subclus.colors
subclus.colors.generic[grep("^CD4", names(subclus.colors.generic))] <- "#FED976"
subclus.colors.generic[grep("^CD8", names(subclus.colors.generic))] <- "#CB181D"
subclus.colors.generic[grep("^B", names(subclus.colors.generic))] <- "#FFCCFF"
subclus.colors.generic[grep("^E", names(subclus.colors.generic))] <- "#4292C6"
subclus.colors.generic[grep("^My", names(subclus.colors.generic))] <- "#74C476"
subclus.colors.generic[grep("^Ma", names(subclus.colors.generic))] <- "#CCCCCC"
subclus.colors.generic[grep("^S", names(subclus.colors.generic))] <- "#8856A7"

# Check dimensions
length(unique(exp.m.sig.Mye.lig.cpdb.means$ligand))
length(unique(exp.m.sig.Mye.lig.cpdb.means$receptor))

# And plot
ggplot(exp.m.sig.Mye.lig.cpdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/Myeloid.ligands.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.sig.Mye.lig.cpdb.means$interaction.score)
sub.exp.m.sig.Mye.lig.cpdb.means <- exp.m.sig.Mye.lig.cpdb.means[exp.m.sig.Mye.lig.cpdb.means$interaction.score > summary(exp.m.sig.Mye.lig.cpdb.means$interaction.score)[5],]

summary(summary(exp.m.sig.Mye.lig.cpdb.means$ligand))
sub.exp.m.sig.Mye.lig.cpdb.means <- sub.exp.m.sig.Mye.lig.cpdb.means[sub.exp.m.sig.Mye.lig.cpdb.means$ligand %in% names(summary(sub.exp.m.sig.Mye.lig.cpdb.means$ligand))[summary(sub.exp.m.sig.Mye.lig.cpdb.means$ligand) < 15],]

# Check dimensions
length(unique(sub.exp.m.sig.Mye.lig.cpdb.means$ligand))
length(unique(sub.exp.m.sig.Mye.lig.cpdb.means$receptor))

# And plot
ggplot(sub.exp.m.sig.Mye.lig.cpdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 14),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/Myeloid.ligands.subset.dots.pdf", width = 10, height = 10)


#===========================================================================================================================
# Subset Myeloid receptor interactions
Mye.rec.cpdb.means <- cpdb.means[,c(1:11, grep(".My.[0-9]", colnames(cpdb.means)))]

# Subset significant interaction pairs
sig.Mye.rec.cpdb.means <- Mye.rec.cpdb.means[apply(Mye.rec.cpdb.means[,12:length(colnames(Mye.rec.cpdb.means))], 1, sum) != 0,]
dim(sig.Mye.rec.cpdb.means)

# Shape the data for plotting
# Prep ligand and receptor name columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(sig.Mye.rec.cpdb.means$interacting_pair), "_")){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
sig.Mye.rec.cpdb.means$ligand <- l
sig.Mye.rec.cpdb.means$receptor <- r

# Melt the 'frame using ligand and receptor info and the interaction scores
m.sig.Mye.rec.cpdb.means <- melt(sig.Mye.rec.cpdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:length(colnames(sig.Mye.rec.cpdb.means))), variable.name = "cells", value.name = "interaction.score")

# Subset pairs with significant interactions (>0)
exp.m.sig.Mye.rec.cpdb.means <- m.sig.Mye.rec.cpdb.means[m.sig.Mye.rec.cpdb.means$interaction.score > 0,]
dim(exp.m.sig.Mye.rec.cpdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.sig.Mye.rec.cpdb.means$cells), split = "_", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.sig.Mye.rec.cpdb.means$ligand.cells <- l
exp.m.sig.Mye.rec.cpdb.means$receptor.cells <- r

# Clean up
exp.m.sig.Mye.rec.cpdb.means$interaction.score <- as.numeric(exp.m.sig.Mye.rec.cpdb.means$interaction.score)
exp.m.sig.Mye.rec.cpdb.means <- exp.m.sig.Mye.rec.cpdb.means[!is.na(exp.m.sig.Mye.rec.cpdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
exp.m.sig.Mye.rec.cpdb.means <- exp.m.sig.Mye.rec.cpdb.means[grep("My.4|Mixed", exp.m.sig.Mye.rec.cpdb.means$ligand.cells, invert = T),]
exp.m.sig.Mye.rec.cpdb.means <- exp.m.sig.Mye.rec.cpdb.means[grep("My.4", exp.m.sig.Mye.rec.cpdb.means$receptor.cells, invert = T),]

# Fix factors
exp.m.sig.Mye.rec.cpdb.means$receptor       <- droplevels(as.factor(exp.m.sig.Mye.rec.cpdb.means$receptor))
exp.m.sig.Mye.rec.cpdb.means$ligand         <- droplevels(as.factor(exp.m.sig.Mye.rec.cpdb.means$ligand))
exp.m.sig.Mye.rec.cpdb.means$receptor.cells <- droplevels(as.factor(exp.m.sig.Mye.rec.cpdb.means$receptor.cells))
exp.m.sig.Mye.rec.cpdb.means$ligand.cells   <- droplevels(as.factor(exp.m.sig.Mye.rec.cpdb.means$ligand.cells))
length(levels(exp.m.sig.Mye.rec.cpdb.means$ligand.cells))

#Define colors for the subclusters
subclus.colors <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#E7298A", # CD79 B
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", "#3366FF", "#0000FF", # Endo
                    "#35978F", # Mast
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", # Myeloid
                    "#8856A7", "#810F7C", "#663399") # SMC
names(subclus.colors) <- levels(exp.m.sig.Mye.rec.cpdb.means$ligand.cells)

subclus.colors.generic <- subclus.colors
subclus.colors.generic[grep("^CD4", names(subclus.colors.generic))] <- "#FED976"
subclus.colors.generic[grep("^CD8", names(subclus.colors.generic))] <- "#CB181D"
subclus.colors.generic[grep("^B", names(subclus.colors.generic))] <- "#FFCCFF"
subclus.colors.generic[grep("^E", names(subclus.colors.generic))] <- "#4292C6"
subclus.colors.generic[grep("^My", names(subclus.colors.generic))] <- "#74C476"
subclus.colors.generic[grep("^Ma", names(subclus.colors.generic))] <- "#CCCCCC"
subclus.colors.generic[grep("^S", names(subclus.colors.generic))] <- "#8856A7"

# Check dimensions
length(unique(exp.m.sig.Mye.rec.cpdb.means$ligand))
length(unique(exp.m.sig.Mye.rec.cpdb.means$receptor))

# And plot
ggplot(exp.m.sig.Mye.rec.cpdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/Myeloid.receptors.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.sig.Mye.rec.cpdb.means$interaction.score)
sub.exp.m.sig.Mye.rec.cpdb.means <- exp.m.sig.Mye.rec.cpdb.means[exp.m.sig.Mye.rec.cpdb.means$interaction.score > summary(exp.m.sig.Mye.rec.cpdb.means$interaction.score)[5],]

summary(summary(exp.m.sig.Mye.rec.cpdb.means$ligand))
sub.exp.m.sig.Mye.rec.cpdb.means <- sub.exp.m.sig.Mye.rec.cpdb.means[sub.exp.m.sig.Mye.rec.cpdb.means$ligand %in% names(summary(sub.exp.m.sig.Mye.rec.cpdb.means$ligand))[summary(sub.exp.m.sig.Mye.rec.cpdb.means$ligand) < 20],]

# Check dimensions
length(unique(sub.exp.m.sig.Mye.rec.cpdb.means$ligand))
length(unique(sub.exp.m.sig.Mye.rec.cpdb.means$receptor))

# And plot
ggplot(sub.exp.m.sig.Mye.rec.cpdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/Myeloid.receptors.subset.dots.pdf", width = 10, height = 10)




#=============================================================================================================
# Investigate random interactions

sub.exp.m.s.cpdb.means[grep("FAM3C", sub.exp.m.s.cpdb.means$receptor),][sub.exp.m.s.cpdb.means[grep("FAM3C", sub.exp.m.s.cpdb.means$receptor),"interaction.score"] > 15,]


#=============================================================================================================
# PLot HLA interactions in and with myeloid cells
#===========================================================================================================================
# Subset Myeloid receptor interactions
Mye.rec.hla.cpdb.means <- hla.cpdb.means[,c(1:11, grep(".My.[0-9]", colnames(hla.cpdb.means)))]
dim(Mye.rec.hla.cpdb.means)

# Shape the data for plotting
# Prep ligand and receptor name columns
l <- vector()
r <- vector()

for (i in strsplit(as.character(Mye.rec.hla.cpdb.means$interacting_pair), "_")){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
Mye.rec.hla.cpdb.means$ligand <- l
Mye.rec.hla.cpdb.means$receptor <- r

# Melt the 'frame using ligand and receptor info and the interaction scores
m.Mye.rec.hla.cpdb.means <- melt(Mye.rec.hla.cpdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:(length(colnames(Mye.rec.hla.cpdb.means))-2)), variable.name = "cells", value.name = "interaction.score")

# Subset pairs with significant interactions (>0)
exp.m.Mye.rec.hla.cpdb.means <- m.Mye.rec.hla.cpdb.means[m.Mye.rec.hla.cpdb.means$interaction.score > 0,]
dim(exp.m.Mye.rec.hla.cpdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.Mye.rec.hla.cpdb.means$cells), split = "_", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.Mye.rec.hla.cpdb.means$ligand.cells <- l
exp.m.Mye.rec.hla.cpdb.means$receptor.cells <- r

# Clean up
exp.m.Mye.rec.hla.cpdb.means$interaction.score <- as.numeric(exp.m.Mye.rec.hla.cpdb.means$interaction.score)
exp.m.Mye.rec.hla.cpdb.means <- exp.m.Mye.rec.hla.cpdb.means[!is.na(exp.m.Mye.rec.hla.cpdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
exp.m.Mye.rec.hla.cpdb.means <- exp.m.Mye.rec.hla.cpdb.means[grep("My.4|Mixed", exp.m.Mye.rec.hla.cpdb.means$ligand.cells, invert = T),]
exp.m.Mye.rec.hla.cpdb.means <- exp.m.Mye.rec.hla.cpdb.means[grep("My.4", exp.m.Mye.rec.hla.cpdb.means$receptor.cells, invert = T),]

# Fix factors
exp.m.Mye.rec.hla.cpdb.means$receptor       <- droplevels(as.factor(exp.m.Mye.rec.hla.cpdb.means$receptor))
exp.m.Mye.rec.hla.cpdb.means$ligand         <- droplevels(as.factor(exp.m.Mye.rec.hla.cpdb.means$ligand))
exp.m.Mye.rec.hla.cpdb.means$receptor.cells <- droplevels(as.factor(exp.m.Mye.rec.hla.cpdb.means$receptor.cells))
exp.m.Mye.rec.hla.cpdb.means$ligand.cells   <- droplevels(as.factor(exp.m.Mye.rec.hla.cpdb.means$ligand.cells))
length(levels(exp.m.Mye.rec.hla.cpdb.means$ligand.cells))

#Define colors for the subclusters
subclus.colors <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#E7298A", # CD79 B
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", "#3366FF", "#0000FF", # Endo
                    "#35978F", # Mast
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", # Myeloid
                    "#8856A7", "#810F7C", "#663399") # SMC
names(subclus.colors) <- levels(exp.m.Mye.rec.hla.cpdb.means$ligand.cells)

subclus.colors.generic <- subclus.colors
subclus.colors.generic[grep("^CD4", names(subclus.colors.generic))] <- "#FED976"
subclus.colors.generic[grep("^CD8", names(subclus.colors.generic))] <- "#CB181D"
subclus.colors.generic[grep("^B", names(subclus.colors.generic))] <- "#FFCCFF"
subclus.colors.generic[grep("^E", names(subclus.colors.generic))] <- "#4292C6"
subclus.colors.generic[grep("^My", names(subclus.colors.generic))] <- "#74C476"
subclus.colors.generic[grep("^Ma", names(subclus.colors.generic))] <- "#CCCCCC"
subclus.colors.generic[grep("^S", names(subclus.colors.generic))] <- "#8856A7"

# Check dimensions
length(unique(exp.m.Mye.rec.hla.cpdb.means$ligand))
length(unique(exp.m.Mye.rec.hla.cpdb.means$receptor))

# And plot
ggplot(exp.m.Mye.rec.hla.cpdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/HLA-only.Myeloid.receptors.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.Mye.rec.hla.cpdb.means$interaction.score)
sub.exp.m.Mye.rec.hla.cpdb.means <- exp.m.Mye.rec.hla.cpdb.means[exp.m.Mye.rec.hla.cpdb.means$interaction.score > summary(exp.m.Mye.rec.hla.cpdb.means$interaction.score)[5],]

summary(summary(exp.m.Mye.rec.hla.cpdb.means$ligand))
sub.exp.m.Mye.rec.hla.cpdb.means <- sub.exp.m.Mye.rec.hla.cpdb.means[sub.exp.m.Mye.rec.hla.cpdb.means$ligand %in% names(summary(sub.exp.m.Mye.rec.hla.cpdb.means$ligand))[summary(sub.exp.m.Mye.rec.hla.cpdb.means$ligand) < 20],]

# Check dimensions
length(unique(sub.exp.m.Mye.rec.hla.cpdb.means$ligand))
length(unique(sub.exp.m.Mye.rec.hla.cpdb.means$receptor))

# And plot
ggplot(sub.exp.m.Mye.rec.hla.cpdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid receptors")) +
  theme_light() +
  ylab("Receptors (in myeloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/HLA-only.Myeloid.receptors.subset.dots.pdf", width = 10, height = 10)




# Subset Myeloid ligand interactions
Mye.lig.hla.cpdb.means <- hla.cpdb.means[,c(1:11, grep("My.[0-9].", colnames(hla.cpdb.means)))]
dim(Mye.lig.hla.cpdb.means)

# Shape the data for plotting
# Prep ligand and receptor name columns
l <- vector()
r <- vector()

for (i in strsplit(as.character(Mye.lig.hla.cpdb.means$interacting_pair), "_")){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
Mye.lig.hla.cpdb.means$ligand <- l
Mye.lig.hla.cpdb.means$receptor <- r

# Melt the 'frame using ligand and receptor info and the interaction scores
m.Mye.lig.hla.cpdb.means <- melt(Mye.lig.hla.cpdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:(length(colnames(Mye.lig.hla.cpdb.means))-2)), variable.name = "cells", value.name = "interaction.score")

# Subset pairs with significant interactions (>0)
exp.m.Mye.lig.hla.cpdb.means <- m.Mye.lig.hla.cpdb.means[m.Mye.lig.hla.cpdb.means$interaction.score > 0,]
dim(exp.m.Mye.lig.hla.cpdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.Mye.lig.hla.cpdb.means$cells), split = "_", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.Mye.lig.hla.cpdb.means$ligand.cells <- l
exp.m.Mye.lig.hla.cpdb.means$receptor.cells <- r

# Clean up
exp.m.Mye.lig.hla.cpdb.means$interaction.score <- as.numeric(exp.m.Mye.lig.hla.cpdb.means$interaction.score)
exp.m.Mye.lig.hla.cpdb.means <- exp.m.Mye.lig.hla.cpdb.means[!is.na(exp.m.Mye.lig.hla.cpdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
exp.m.Mye.lig.hla.cpdb.means <- exp.m.Mye.lig.hla.cpdb.means[grep("My.4|Mixed", exp.m.Mye.lig.hla.cpdb.means$ligand.cells, invert = T),]
exp.m.Mye.lig.hla.cpdb.means <- exp.m.Mye.lig.hla.cpdb.means[grep("My.4", exp.m.Mye.lig.hla.cpdb.means$receptor.cells, invert = T),]

# Fix factors
exp.m.Mye.lig.hla.cpdb.means$receptor       <- droplevels(as.factor(exp.m.Mye.lig.hla.cpdb.means$receptor))
exp.m.Mye.lig.hla.cpdb.means$ligand         <- droplevels(as.factor(exp.m.Mye.lig.hla.cpdb.means$ligand))
exp.m.Mye.lig.hla.cpdb.means$receptor.cells <- droplevels(as.factor(exp.m.Mye.lig.hla.cpdb.means$receptor.cells))
exp.m.Mye.lig.hla.cpdb.means$ligand.cells   <- droplevels(as.factor(exp.m.Mye.lig.hla.cpdb.means$ligand.cells))
length(levels(exp.m.Mye.lig.hla.cpdb.means$ligand.cells))

#Define colors for the subclusters
subclus.colors <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", # CD4 T
                    "#E7298A", # CD79 B
                    "#EF3B2C", "#CB181D", "#A50F15", "#67000D", # CD8 T
                    "#4292C6", "#2171B5", "#08519C", "#08306B", "#3366FF", "#0000FF", # Endo
                    "#35978F", # Mast
                    "#74C476", "#41AB5D", "#238B45", "#006D2C", # Myeloid
                    "#8856A7", "#810F7C", "#663399") # SMC
names(subclus.colors) <- levels(exp.m.Mye.lig.hla.cpdb.means$receptor.cells)

subclus.colors.generic <- subclus.colors
subclus.colors.generic[grep("^CD4", names(subclus.colors.generic))] <- "#FED976"
subclus.colors.generic[grep("^CD8", names(subclus.colors.generic))] <- "#CB181D"
subclus.colors.generic[grep("^B", names(subclus.colors.generic))] <- "#FFCCFF"
subclus.colors.generic[grep("^E", names(subclus.colors.generic))] <- "#4292C6"
subclus.colors.generic[grep("^My", names(subclus.colors.generic))] <- "#74C476"
subclus.colors.generic[grep("^Ma", names(subclus.colors.generic))] <- "#CCCCCC"
subclus.colors.generic[grep("^S", names(subclus.colors.generic))] <- "#8856A7"

# Check dimensions
length(unique(exp.m.Mye.lig.hla.cpdb.means$ligand))
length(unique(exp.m.Mye.lig.hla.cpdb.means$receptor))

# And plot
ggplot(exp.m.Mye.lig.hla.cpdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/HLA-only.Myeloid.ligands.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.Mye.lig.hla.cpdb.means$interaction.score)
sub.exp.m.Mye.lig.hla.cpdb.means <- exp.m.Mye.lig.hla.cpdb.means[exp.m.Mye.lig.hla.cpdb.means$interaction.score > summary(exp.m.Mye.lig.hla.cpdb.means$interaction.score)[5],]

summary(summary(exp.m.Mye.lig.hla.cpdb.means$ligand))
sub.exp.m.Mye.lig.hla.cpdb.means <- sub.exp.m.Mye.lig.hla.cpdb.means[sub.exp.m.Mye.lig.hla.cpdb.means$ligand %in% names(summary(sub.exp.m.Mye.lig.hla.cpdb.means$ligand))[summary(sub.exp.m.Mye.lig.hla.cpdb.means$ligand) < 20],]

# Check dimensions
length(unique(sub.exp.m.Mye.lig.hla.cpdb.means$ligand))
length(unique(sub.exp.m.Mye.lig.hla.cpdb.means$receptor))

# And plot
ggplot(sub.exp.m.Mye.lig.hla.cpdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors.generic, aesthetics = "col") +
  ggtitle(paste("Myeloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in myeloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("lig.rec.int/cellphonedb_results/HLA-only.Myeloid.ligands.subset.dots.pdf", width = 10, height = 10)
