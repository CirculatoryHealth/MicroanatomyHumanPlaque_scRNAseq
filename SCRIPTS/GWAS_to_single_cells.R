### 2019-5-16

# Project SCS 18 patients Experiment Kenya

# DEG and gene pattern analysis

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Kenya")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.8
library(org.Hs.eg.db)   # v3.7.0
library(SingleR)        # v0.2.2

# object
seuset <-  load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined
average.expression <- readRDS("seuset_average_expression_DEG.RDS")
Kmeans.object <- readRDS("seuset_Kmeans_object_average_expression_DEG.RDS")

# set seed
set.seed(15)

# check tSNE
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "check")

#-----------------------------------------------------------------------------------------
### Differential expression
## pairwise comparison of differential expression data between clusters

# fetch identity class names
fetch.identity <- FetchData(object = seuset,
                            vars.all = "ident")

# convert factor to string/characters
identity.classes <- as.character((unique(fetch.identity$ident)))
# cluster number
cluster.number <- length(identity.classes)

## loop over different clusters and save differential exression data
# create dataframe
differential.expression.clusters <- data.frame()

# prepare all possible combinations
choices <- choose(cluster.number, 2)        # identify the number of possible pairs
combinations <- combn(cluster.number, 2)    # create those pairs

# pairwise comparison for differential expression
for (choice in 1:choices) {
  # select the cluster numbers to run
  cluster1 <- combinations[,choice][1] # the first index of the combination
  cluster2 <- combinations[,choice][2] # the second index of the combination
  
  # input function must be character
  cluster1 <- identity.classes[cluster1]
  cluster2 <- identity.classes[cluster2]
  
  # run formula
  pairwise.comparison <- FindMarkers(object = seuset,
                                     ident.1 = cluster1 ,
                                     ident.2 = cluster2,
                                     logfc.threshold = 0.6, #default 0.25
                                     min.pct = 0.1,
                                     test.use = "wilcox") # default wilcox
  
  # add gene names for later
  pairwise.comparison$gene <- rownames(pairwise.comparison)
  # add clusters
  pairwise.comparison$cluster1 <- cluster1
  pairwise.comparison$cluster2 <- cluster2
  # select only values with p adjusted of certain value
  pairwise.comparison <- pairwise.comparison[(pairwise.comparison$p_val_adj <= .05), ]
  # bind to dataframe
  differential.expression.clusters <- rbind(differential.expression.clusters, pairwise.comparison)
  
}

dim(differential.expression.clusters)

## differential expression findAllMarkers
differential.expression.all <- FindAllMarkers(object = seuset,
                                              logfc.threshold = 0.6,
                                              min.pct = 0.1,
                                              only.pos = FALSE,
                                              test.use = "wilcox")

dim(differential.expression.all)
differential.expression.all <- differential.expression.all[(differential.expression.all$p_val_adj <= 0.05),]
dim(differential.expression.all)

## combine differential expression based on gene names, not rownames
differential.expression.data <- union(differential.expression.all$gene, differential.expression.clusters$gene)
length(differential.expression.data)
length(unique(differential.expression.data))

## save the table with results
# write.table(differential.expression.clusters, "DEG_one_vs_one.txt")
# write.table(differential.expression.all, "DEG_one_vs_all.txt")
# write.table(differential.expression.data, "All_DEG_(gene_names).txt")

#-----------------------------------------------------------------------------------------
### K means
## K means with average expression
  
# calculate average expression with the genes from differential.expression.data that are present in the object
average.expression <-  AverageExpression(object = seuset,
                                         genes.use = differential.expression.data,
                                         return.seurat = TRUE,
                                         show.progress = T)

# save object
#saveRDS(average.expression, "seuset_average_expression_DEG.RDS")

## K means
geneclusters <- 15
seed <- 67

## calculate Kmeans
# option 1 (seurat heatmap)
Kmeans.object <- DoKMeans(average.expression,
                          genes.use = differential.expression.data,
                          k.genes = geneclusters,
                          do.plot = T,
                          k.seed = seed) #67
ggsave(paste(result.folder, "/", Sys.Date(), " Kmeans ", geneclusters, " seed = ", seed, " .pdf", sep = ""), height = 8, limitsize = FALSE)

# save object
#saveRDS(average.expression, "seuset_Kmeans_object_average_expression_DEG.RDS")

#-----------------------------------------------------------------------------------------
### prepare reference data
reference.data <- read.delim2("reference_REFSEQ_UCSC_2015.txt", header = F, stringsAsFactors = F)
dim(reference.data)
# collapse duplicates
reference.data <- unique(reference.data)
dim(reference.data)
reference.data <- reference.data$V1

## remove unwanted genes
unwanted.genes <- c("^ERCC", "^MT-","^RPL", "^RPS","^UGDH-AS1$","^PGM2P2$", "^LOC100131257$", "^KCNQ1OT1$", 
"^MALAT1$", "^PGM5P2$", "^MAB21L3$","^EEF1A1$")

reference.data <- reference.data[grep(paste(unwanted.genes, collapse = "|"),reference.data,invert=TRUE)]
length(reference.data)

## update gene names
# update gene symbols to newest version
library(org.Hs.eg.db)

update.symbols <- function(my.counts.table = my.counts.table){
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping
  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = my.counts.table, columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- unique(gene.info)
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- gene.info$SYMBOL
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  #Overwrite old symbols in row names and clean up
  row.names(my.counts.table) <- my.counts.table$SYMBOL
  my.counts.table$SYMBOL <- NULL
  
  return(my.counts.table)
}

reference.data <- update.symbols(reference.data)
length(reference.data)


#-----------------------------------------------------------------------------------------
### FUMA
## Which proportion of the FUMA genes is expressed in each k means pattern?

# differentially expressed genes (load from object)
differential.expression.data <- rownames(Kmeans.object@raw.data)

# extract which gene belongs to which K means gene cluster:
genes.distributed.across.kmeans <- Kmeans.object@kmeans@gene.kmeans.obj$cluster
# save
write.table(genes.distributed.across.kmeans, "All_DEG_distributed_across_Kmeans_12062019.txt")

# prepare files
all.files <- list.files(pattern="*.txt")
all.files <- all.files[-1]
all.files
FUMA.names <- c("ArteryBlood.withMAGMA", "ArteryBlood", "UKBB_CAD", "UKBB_CAD.withMAGMA")

# FUMA info
FUMA.variables <- c("total.FUMA.genes", "possible.FUMA.genes", "FUMA.genes.after.magma.filtering", "FUMA.differentially.expressed")
FUMA.info <- data.frame(matrix(nrow = length(FUMA.names), ncol = length(FUMA.variables)))
colnames(FUMA.info) <- FUMA.variables
rownames(FUMA.info) <- FUMA.names

# to write
differentially.expressed.per.FUMA.list <- list()
FUMA.files <- list()

for (df in c(1:length(all.files))) {
  next.df <- read.table((all.files[df]), header = TRUE)
  
  # check if MAGMA file by checking if there is a column named P
  if ("P" %in% colnames(next.df) == T) { #if true
    
    # take only the rows that are present in the reference dataset
    intersecting.genes <- next.df[is.element(next.df$symbol, reference.data),]
    
    # filtering restrictions
    filtered.genes <- intersecting.genes[(intersecting.genes$P <= 0.05), ]
    
  } else { # if false
    
    # take only the rows that are present in the dataset
    intersecting.genes <- next.df[is.element(next.df$symbol, reference.data),]
    
    # filtering restrictions
    filtered.genes <- intersecting.genes[(intersecting.genes$eqtlMapSNPs > 0 | intersecting.genes$ciMap == "Yes"), ]
  } 
  
  # save information
  FUMA.files[[df]] <- filtered.genes
  
  # assingn new name
  new.name <- FUMA.names[df]
  # write.csv(filtered.genes, file = paste(result.folder, "/", new.name, ".csv", sep = ""))
  assign(new.name, unique(as.character(filtered.genes$symbol)))
  
  ## Complete FUMA.info
  # total FUMA genes
  FUMA.info$total.FUMA.genes[df] <- dim(next.df)[1]
  
  # FUMA genes in dataset
  left.after.intersecting <- dim(intersecting.genes)[1]
  FUMA.info$possible.FUMA.genes[df] <- left.after.intersecting
  
  # FUMA genes after filtering
  left.after.filtering <- dim(filtered.genes)[1]
  FUMA.info$FUMA.genes.after.magma.filtering[df] <- left.after.filtering
  
  # FUMA differentially expressed (filtered genes only)
  differentially.expressed.genes <- intersect(filtered.genes$symbol, differential.expression.data)
  FUMA.differentially.expressed <- length(differentially.expressed.genes)
  
  FUMA.info$FUMA.differentially.expressed[df] <- FUMA.differentially.expressed
  
  # #write all differentially expressed genes to file
  write.table(intersecting.genes, file = paste(result.folder, "/", new.name, " Possible FUMA genes.txt", sep = ""),
              sep = "\t")
  write.table(filtered.genes, file = paste(result.folder, "/", new.name, " Filtered with MAGMA.txt", sep = ""),
              sep = "\t")
  
  # save all differentailly expressed genes
  differentially.expressed.per.FUMA.list[[df]] <- differentially.expressed.genes
  
}

# check
#View(differentially.expressed.per.FUMA.list)
FUMA.info
write.table(FUMA.info, file = paste(result.folder, "/", new.name, " FUMA_info.txt", sep = ""),
            sep = "\t")

## Individual expression of differential genes per FUMA entry per genepattern
FUMA.genecluster.counts <- data.frame(row.names = FUMA.names)

# for each FUMA entry
for (listing in 1:length(differentially.expressed.per.FUMA.list)) {
  # select the list of genes that are differentially expressed
  my.list <- differentially.expressed.per.FUMA.list[[listing]]
  # create subset specific for my.list (as matrix is needed)
  current.df <- as.matrix(genes.distributed.across.kmeans[my.list])
  
  write.table(current.df, file = paste(result.folder, "/", FUMA.names[listing], " all DEG with genepattern.txt", sep = ""))
  
  for (cluster in 1:geneclusters) {
    # calculate the amount of genes corresponding to a gene cluster
    current.length <- length(current.df[current.df[,1] == cluster, ])
    # add this to data frame
    FUMA.genecluster.counts[listing,cluster] <- current.length
  }
}

# check
FUMA.genecluster.counts
rowSums(FUMA.genecluster.counts)
FUMA.info$FUMA.differentially.expressed

# write.table(FUMA.info, file = paste(result.folder, "/", new.name, " FUMA_genecluster_counts.txt", sep = ""),
#             sep = "\t")

#-----------------------------------------------------------------------------------------
### CAD expression vs random genes seperately
# select random data (which is the same as reference data)
random.data <- reference.data

# set parameters
sample.size <- FUMA.info[1,3]
nr.of.random.samples <- 75000
sample.data <- matrix(nrow = sample.size, ncol = nr.of.random.samples)


# generate random samples from total dataset
for (number in c(1:nr.of.random.samples) ) {
  current.genes <- sample(x = random.data, size = sample.size)
  sample.data[,number] <- current.genes
}


## Get (differentially) expressed
random.data.variables <- c("total.random.genes", "random.genes.after.magma.filtering", "fraction.random.differentially.expressed")
random.data.info <- data.frame(matrix(nrow = dim(sample.data)[2], ncol = length(random.data.variables)))
colnames(random.data.info)  <- random.data.variables

differentially.expressed.per.random.dataset <- list()

for (list in c(1:nr.of.random.samples)) {
  next.list <- sample.data[,list]
  
  # take only the rows that are present in the reference dataset
  intersecting.genes <- intersect(next.list, reference.data)
  
  ## Complete random.data.info
  # total random genes
  random.data.info$total.random.genes[list] <- length(next.list)
  # FUMA genes after filtering
  left.after.filtering <- length(intersecting.genes)
  random.data.info$random.genes.after.magma.filtering[list] <- left.after.filtering
  # fraction fuma differentially expressed
  differentially.expressed.genes <- intersect(intersecting.genes, differential.expression.data)
  random.differentially.expressed <- length(differentially.expressed.genes)
  #FUMA.info$fraction.FUMA.differentially.expressed[df] <- FUMA.differentially.expressed
  random.data.info$fraction.random.differentially.expressed[list] <- random.differentially.expressed/random.data.info$total.random.genes[list]
  
  # save the differentially expressed genes
  differentially.expressed.per.random.dataset[[list]] <- differentially.expressed.genes
  if(list %% 5000 == 0){
    print(list)
  }
  
}

# check
head(random.data.info)

## Individual expression of differential genes per random entry
random.genes.genecluster.counts <- data.frame(matrix(ncol = geneclusters, nrow = nr.of.random.samples))

# for each random entry
for (list in 1:length(differentially.expressed.per.random.dataset)) {
  # select the list of genes that are differentially expressed
  my.list <- differentially.expressed.per.random.dataset[[list]]
  # create subset specific for my.list (as matrix is needed)
  current.df <- as.matrix(genes.distributed.across.kmeans[my.list])
  
  for (cluster in 1:geneclusters) {
    # calculate the amount of genes corresponding to a gene cluster
    current.length <- length(current.df[current.df[,1] == cluster, ])
    # add this to data frame
    random.genes.genecluster.counts[list, cluster] <- current.length
  }
  if(list %% 5000 == 0){
    print(list)
  }
}

head(rowSums(random.genes.genecluster.counts))
saveRDS(random.genes.genecluster.counts, file = "random.genes.genecluster.counts")

### enrichment
## get fold enrichment per cluster
CAD <- FUMA.genecluster.counts[1,]
random.sample.mean <- apply(random.genes.genecluster.counts, 2, mean)
enrichment <- CAD

# only positive enrichment
for (i in c(1:length(CAD))) {
  result <- CAD[[i]] / random.sample.mean[[i]]
  enrichment[i] <- result
}

enrichment

p.values <- c()

# calculate p.value
for (i in 1:geneclusters) {
  current.value <- CAD[[i]]
  count = 0
  # check if random value is same or higher than current value
  for (j in 1:nr.of.random.samples) {
    if (random.genes.genecluster.counts[j,i] >= current.value) {
      count = count + 1
    }
  } 
  # divide by nr.of.random.samples to get p value
  p.value = count/nr.of.random.samples
  # add to vector
  p.values <- c(p.values, p.value)
}

p.values

## get fold enrichment for total DEG
DEG.enrichment <- c()
FUMA.total.DEG <- FUMA.info[1,4]
random.total.DEG <- rowSums(random.genes.genecluster.counts)

# enrichment
for (i in 1:nr.of.random.samples) {
  result <- FUMA.total.DEG / random.total.DEG[[i]]
  DEG.enrichment[i] <- result
}

# total DEG p.values
total.DEG.p.value <- sum(random.total.DEG >= FUMA.total.DEG)
total.DEG.p.value <- total.DEG.p.value / nr.of.random.samples
total.DEG.p.value

# save values
# write.table(rbind(enrichment, p.values), file = paste(result.folder, "/", new.name, " Enrichment results_1.txt", sep = ""), sep = "\t")
# write.table(total.DEG.p.value, file = paste(result.folder, "/", new.name, " Enrichment results_2.txt", sep = ""), sep = "\t")


#-----------------------------------------------------------------------------------------
### plot fraction of whole and absolute counts with colour based on enrichment

FUMA.genecluster.counts.ind <- FUMA.genecluster.counts[1,]

FUMA.genecluster.counts.ind$sums <- rowSums(FUMA.genecluster.counts.ind)
FUMA.genecluster.counts.ind <- FUMA.genecluster.counts.ind[,1:geneclusters,]/FUMA.genecluster.counts.ind$sums
FUMA.genecluster.counts.ind[is.na(FUMA.genecluster.counts.ind)] <- 0

# reverse enrichment for the plot
enrichment <- rev(enrichment)

# colours for barplot
test.cols <- c()
for (i in 1:length(enrichment)){
  j <- enrichment[1,i]
  if (j >= 2) {
    test.cols <- c(test.cols, "black")
  } else if (j <= 1) {
    test.cols <- c(test.cols, "white")
  } else if (1 < j && j < 1.5) {
    test.cols <- c(test.cols, "lightgrey")
  } else if ( j >= 1.5 && j < 2) {
    test.cols <- c(test.cols, "darkgrey")
  }
}

# barplot to match genepatterns with heatmap
pdf(paste(result.folder, "/", Sys.Date(), " Distribution differentially expressed genes", FUMA.names[1], ".pdf", sep = ""))
brplt <- barplot(height = t(as.matrix(rev(FUMA.genecluster.counts.ind))),
                 main = "Distribution of differentially expressed genes across gene clusters fraction",
                 ylab = "Fraction",
                 las = 2,
                 width = 1,
                 col = test.cols)
dev.off()
brplt

# absolute counts
pdf(paste(result.folder, "/", Sys.Date(), " Distribution differentially expressed genes absolute counts", FUMA.names[1], ".pdf", sep = ""))
brplt <- barplot(height = t(as.matrix(rev(FUMA.genecluster.counts[1,]))),
                 main = "Distribution of differentially expressed genes across gene clusters",
                 ylab = "Counts",
                 las = 2,
                 col = test.cols)
dev.off()

#-----------------------------------------------------------------------------------------
### make boxplot for each Kmeans gene pattern instead of the big heatmap

## the heatmap
# sort the clusters by number
genes.distributed.across.kmeans.sorted <- as.matrix(sort(Kmeans.object@kmeans@gene.kmeans.obj$cluster))

DoHeatmap(object = Kmeans.object,
          genes.use = rownames(genes.distributed.across.kmeans.sorted),
          slim.col.label = TRUE,
          group.label.rot = TRUE,
          cex.row = 3,
          title = "Average expression differentially expressed genes",
          do.plot = TRUE)
ggsave(paste(result.folder, "/", Sys.Date(), " Average expression differentially expressed genes ordererd by genecluster", ".pdf", sep = ""), height = 11, width = 6, limitsize = FALSE)


genes.distributed.across.kmeans <- as.matrix(Kmeans.object@kmeans@gene.kmeans.obj$cluster)
# extract celltypes from meta data (starts at 0)
seurat.meta.data <- cbind(seuset@meta.data, fetch.identity)
numeric.celltype <- sort(as.integer(unique(seurat.meta.data$res.1.2)))

## subset all genes per celltype
# to store seurat data per celltype
scale.data.per.celltype <- list()

for (celltype in numeric.celltype) {
  # extract all samples corresponding to celltype
  meta.data <- seuset@meta.data[seuset@meta.data$res.1.2 == celltype,]
  samples <- rownames(meta.data)
  
  # use these samples to subset the right columns from dataset
  df <- seuset@scale.data[,samples]
  scale.data.per.celltype[[celltype + 1]] <- df
}

## Plot gene expression for every genepattern in every celltype
genes.per.genepattern <- list()

for (genepattern in c(1:geneclusters)) {
  # subset data
  data <- as.matrix(genes.distributed.across.kmeans[genes.distributed.across.kmeans[,1] ==  genepattern,])
  genes <- rownames(data)
  genes.per.genepattern[[genepattern]] <- genes
  
  # create list that resets for each genepattern
  current.list <- list()
  current.names <- c()
  
  # dataframe that resets for each genepattern as input for violin plot
  current.vln.input <- data.frame(matrix(nrow = length(genes)))
  
  # subset each celltype with the genes from the genepattern and plot this
  for (df in c(1:length(scale.data.per.celltype))) {
    # get dataset
    current.df <- scale.data.per.celltype[[df]]
    # select only genes from gene pattern
    current.df <- current.df[genes,]
    # average expression
    current.df <- apply(current.df, 1, mean)
    # add to list
    current.list[[df]] <- current.df
    # for violin plot
    current.vln.input[,df] <- current.df
    
    # get name of celltype
    current.name <- as.character(seurat.meta.data[seurat.meta.data$res.1.2 == (df-1), "ident"])
    current.names <- c(current.names, current.name[1])
  }
  
  # prepare for violinplot
  colnames(current.vln.input) <- current.names
  current.vln.input <- reshape2::melt(current.vln.input, id.vars = NULL)
  
  # plot and save it
  pdf(paste(result.folder, "/", Sys.Date(), " Boxplot gene pattern ", genepattern, ".pdf", sep = ""))
  boxplot(current.list, main = paste("pattern #", genepattern, sep = ""), xlab = "celltypes", ylab = "Average expression", names = current.names, las = 2)
  dev.off()
  ggplot(current.vln.input, aes(x = variable, y = value)) + geom_violin() + ggtitle(paste("pattern #", genepattern, sep = "")) +
    xlab("celltypes") + xlab("Average expression") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5, colour = "darkgrey") + 
    geom_boxplot(width=0.08, colour = "red")
  ggsave(paste(result.folder, "/", Sys.Date(), " Violinplot gene pattern ", genepattern, ".pdf", sep = ""))
}

ggplot(current.vln.input, aes(x = variable, y = value)) + geom_violin() + ggtitle("Hoi") + xlab("celltypes") + 
  ylab("Average expressionnnnn") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5, colour = "black") + 
  geom_boxplot(width=0.08, colour = "red")

#-----------------------------------------------------------------------------------------
### Explore Enriched
# select the genes
select.FUMA <- differentially.expressed.per.FUMA.list[[1]]
genes.distributed.across.kmeans <- as.matrix(Kmeans.object@kmeans@gene.kmeans.obj$cluster)
explore.enriched <- as.data.frame(genes.distributed.across.kmeans)
head(explore.enriched)

explore.enriched <- genes.distributed.across.kmeans[select.FUMA,,drop = F]
head(explore.enriched)
explore.enriched <- explore.enriched[order(explore.enriched[,1]),,drop = F]

explore.enriched <- rownames(subset(explore.enriched, explore.enriched[,1] == 3 |
                                      explore.enriched[,1] == 8 |
                                      explore.enriched[,1] == 14))


DoHeatmap(object = Kmeans.object,
          genes.use = explore.enriched,
          group.label.rot = T,
          do.plot = T,
          slim.col.label = T,
          cex.row = 5)
ggsave(paste(result.folder, "/", Sys.Date(), " Average expression of enriched genes", ".pdf", sep = ""))
