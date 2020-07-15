### 2019-5-16

# Project SCS 18 patients Experiment Netherlands

# single cell vs Bulk again

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Netherlands")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(plyr)          # v0.7.8
library(SingleR)        # v0.2.2
library(org.Hs.eg.db)

# object
load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined
seuset <- UpdateSeuratObject(all.seur.combined)
rm(all.seur.combined)

# data 
single.vs.whole <- read.csv("pseudobulk_vs_bulk.txt", header = T, row.names = 1)
head(single.vs.whole)

# set seed
set.seed(69)

#-----------------------------------------------------------------------------------------
### get all marker genes
## get marker genes/top differentially expressed
# note that marker genes are only positively DEG
#load and save
all.marker.files <- list.files(pattern="*genes.txt")
all.marker.files
marker.genes.list <- list()

for (i in seq_along(all.marker.files)) {
  current.markers <- read.delim(all.marker.files[i], stringsAsFactors = F)
  current.markers <- current.markers[1:10,"gene"]
  marker.genes.list[[i]] <- current.markers
}

names(marker.genes.list) <- paste("Cluster", 0:13, sep = " ")
View(marker.genes.list)

# marker.genes <- as.data.frame(marker.genes.list)
# marker.genes
# melted.marker.genes <- reshape2::melt(t(marker.genes))
# melted.marker.genes

#-----------------------------------------------------------------------------------------
### create count matrix for bulk
# load files
whole.tissue <- read.delim("DH_combined_counts.names.txt", header = TRUE, stringsAsFactors = F)
dim(whole.tissue)
head(whole.tissue)

# remove unwanted genes
unwanted.genes <- c("^ERCC", "^MT-", "^RPL", "^RPS", "^KCNQ1OT1$", "^UGDH-AS1$", "^MALAT1$", "^EEF1A$")
whole.tissue <- whole.tissue[grep(paste(unwanted.genes, collapse = "|"), whole.tissue$X ,invert=TRUE),]
dim(whole.tissue)

## Focus on Bulk 1
# order genes by highest expression of Bulk1
whole.tissue <- whole.tissue[order(whole.tissue$RNASvDL2207.sam.counts, decreasing = T),]
# then remove duplicates (will remove the second entry of each duplicate which is lowest expressed)
whole.tissue <- whole.tissue[!duplicated(whole.tissue$X),]
head(whole.tissue)
dim(whole.tissue)

## update gene symbols
# Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
# Retrieve the mapping
gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = whole.tissue$X, columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
dim(gene.info)

# alias is columX
# remove duplicate aliases
gene.info <- gene.info[!duplicated(gene.info$ALIAS),]
dim(gene.info)

# add symbols to table
whole.tissue$SYMBOL <- gene.info$SYMBOL

# Remove non-mappings (old LOCs and stuff that are not availble anymore)
whole.tissue <- whole.tissue[!is.na(whole.tissue$SYMBOL),]

# remove duplicated symbols
whole.tissue <- whole.tissue[!duplicated(whole.tissue$SYMBOL),]

# remove excess info
whole.tissue <- whole.tissue[,c(2,3,5)]
names(whole.tissue) <- c("Bulk", "Bulk2", "SYMBOL")
head(whole.tissue)
dim(whole.tissue)

#-----------------------------------------------------------------------------------------
### create count matrix for pseudobulk
# load files
all.files <- list.files(pattern="*.tsv")
single.cell <- read.delim(all.files[1])
for (df in c(2:length(all.files))) {
  next.df <- read.delim(all.files[df])
  single.cell <- full_join(single.cell, next.df, by = "GENEID")
}

# replace NA with 0
single.cell[is.na(single.cell)] <- 0
# fix GENEID (get rid of chromosome info)
single.cell$GENEID <- sub("__.*", "" , single.cell$GENEID)
dim(single.cell)

# remove unwanted genes
unwanted.genes <- c("^ERCC", "^MT-", "^RPL", "^RPS", "^KCNQ1OT1$", "^UGDH-AS1$", "^MALAT1$", "^EEF1A$")
single.cell <- single.cell[grep(paste(unwanted.genes, collapse = "|"), single.cell$GENEID ,invert=TRUE),]
dim(single.cell)

# create pseudobulk
single.cell$Pseudobulk <- round(rowSums(single.cell[,-which(names(single.cell) == "GENEID")]))
single.cell <- single.cell[,c("Pseudobulk", "GENEID")]
head(single.cell)

# order genes by highest expression
single.cell <- single.cell[order(single.cell$Pseudobulk, decreasing = T),]
head(single.cell)

## update gene symbols
# Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
# Retrieve the mapping
gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = single.cell$GENEID, columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
dim(gene.info)

# alias is columX
# remove duplicate aliases
gene.info <- gene.info[!duplicated(gene.info$ALIAS),]
dim(gene.info)

# add symbols to table
single.cell$SYMBOL <- gene.info$SYMBOL

# Remove non-mappings (old LOCs and stuff that are not availble anymore)
single.cell <- single.cell[!is.na(single.cell$SYMBOL),]

# remove duplicated symbols
single.cell <- single.cell[!duplicated(single.cell$SYMBOL),]

# select columns
single.cell <- single.cell[,c("Pseudobulk","SYMBOL")]
head(single.cell)

## sanity check
sum(melted.marker.genes$value %in% single.cell$SYMBOL)

#-----------------------------------------------------------------------------------------
### combine single and whole tissue data
single.vs.whole <- full_join(single.cell, whole.tissue, by = "SYMBOL")
head(single.vs.whole)

# replace NA with 0
single.vs.whole[is.na(single.vs.whole)] <- 0

# make it pretty
rownames(single.vs.whole) <- single.vs.whole$SYMBOL
single.vs.whole <- single.vs.whole[,c("Pseudobulk","Bulk", "Bulk2")]
head(single.vs.whole)
dim(single.vs.whole)

#write.csv(single.vs.whole, "pseudobulk_vs_bulk.txt")

#-----------------------------------------------------------------------------------------
# define plotting data
plot.pseudo <- log2(single.vs.whole[,c(1,2)] + 0.1)

# get colours
colour.palette <- scales::hue_pal()(14)

# input
ggplot.input <- unlist(marker.genes.list, recursive=FALSE)

## plot all markers together
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[ggplot.input,], aes(plot.pseudo[ggplot.input,1], plot.pseudo[ggplot.input,2]), 
             colour = "red", size = 4) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("Top 10 DEG of each cell cluster", sep = " ")) 
#geom_text(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2], label = all.top.10))
ggsave(paste(result.folder, "/", Sys.Date(), " Pseudobulk vs. Bulk all marker genes in red.pdf", sep = ""))

top.10s <- marker.genes
q <- ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) 

q

lapply(marker.genes.list, function(x){
  q <- q +
      geom_point(data = plot.pseudo[x,], aes(plot.pseudo[x,1], plot.pseudo[x,2]),
                 colour = "black", size = 5, shape = 21, fill = "red")
  return(q)
})




## plot all markers seperate
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2]), colour = "black", size = 1) +
  # CD3 T-Cells I
  geom_point(data = plot.pseudo[top.10s$Cluster.0,], aes(plot.pseudo[top.10s$Cluster.0,1], plot.pseudo[top.10s$Cluster.0,2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[1]) + 
  # CD3+ T-Cells II
  geom_point(data = plot.pseudo[top.10s$Cluster.1,], aes(plot.pseudo[top.10s$Cluster.1,1], plot.pseudo[top.10s$Cluster.1,2]),
           colour = "black", size = 5, shape = 21, fill = colour.palette[2]) +
  # Mixed cells
  geom_point(data = plot.pseudo[top.10s[,3],], aes(plot.pseudo[top.10s[,3],1], plot.pseudo[top.10s[,3],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[3] ) +
  # CD3+ T-Cells III
  geom_point(data = plot.pseudo[top.10s[,4],], aes(plot.pseudo[top.10s[,4],1], plot.pseudo[top.10s[,4],2]),
             colour = "black", shape = 21, fill = colour.palette[4], size = 5) +
  # CD3+ T-cells IV
  geom_point(data = plot.pseudo[top.10s[,5],], aes(plot.pseudo[top.10s[,5],1], plot.pseudo[top.10s[,5],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[5]) +
  # CD68+ MYeloid Cells I
  geom_point(data = plot.pseudo[top.10s[,6],], aes(plot.pseudo[top.10s[,6],1], plot.pseudo[top.10s[,6],2]),
             colour = "black", shape = 21, fill = colour.palette[6], size = 5) +
  # CD14+CD68+ Macrophages II
  geom_point(data = plot.pseudo[top.10s[,7],], aes(plot.pseudo[top.10s[,7],1], plot.pseudo[top.10s[,7],2]),
             colour = "black", shape = 21, fill = colour.palette[7], size = 5) +
  # CCD68+ MYeloid Cells III
  geom_point(data = plot.pseudo[top.10s[,8],], aes(plot.pseudo[top.10s[,8],1], plot.pseudo[top.10s[,8],2]),
             colour = "black", shape = 21, fill = colour.palette[8], size = 5) +
  # ACTA2+ Smooth muscle cells
  geom_point(data = plot.pseudo[top.10s[,9],], aes(plot.pseudo[top.10s[,9],1], plot.pseudo[top.10s[,9],2]),
             colour = "black", shape = 21, fill = colour.palette[9], size = 5) +
  # CD34+ endothelial Cells I
  geom_point(data = plot.pseudo[top.10s[,10],], aes(plot.pseudo[top.10s[,10],1], plot.pseudo[top.10s[,10],2]),
             colour = "black", shape = 21, fill = colour.palette[10], size = 5) +
  # CD34+ endothelial Cells II
  geom_point(data = plot.pseudo[top.10s[,11],], aes(plot.pseudo[top.10s[,11],1], plot.pseudo[top.10s[,11],2]),
             colour = "black", shape = 21, fill = colour.palette[11], size = 5) +
  # CD79A+ B cells
  geom_point(data = plot.pseudo[top.10s[,12],], aes(plot.pseudo[top.10s[,12],1], plot.pseudo[top.10s[,12],2]),
             colour = "black", shape = 21, fill = colour.palette[12], size = 5) +
  # CD68+ MYeloid Cells IV
  geom_point(data = plot.pseudo[top.10s[,13],], aes(plot.pseudo[top.10s[,13],1], plot.pseudo[top.10s[,13],2]),
             colour = "black", shape = 21, fill = colour.palette[13], size = 5) +
  # KIT+ Mast Cells
  geom_point(data = plot.pseudo[top.10s[,14],], aes(plot.pseudo[top.10s[,14],1], plot.pseudo[top.10s[,14],2]),
             colour = "black", shape = 21, fill = colour.palette[14], size = 5) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2]))


ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2]), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[top.10s[,2],], aes(plot.pseudo[top.10s[,2],1], plot.pseudo[top.10s[,2],2]),
            colour = "black", size = 5, shape = 21, fill = colour.palette[2])


#-----------------------------------------------------------------------------------------
head(plot.pseudo)

cor.test(plot.pseudo[,1], plot.pseudo[,2], method = "pearson")

cor.test(plot.pseudo[melted.marker.genes$value,1], plot.pseudo[melted.marker.genes$value,2], method = "pearson")
