#Load in the data
#Setup the matrixes and perform initial filtering
#-----------------------------------------------------------------------------------------

#Load packages
library(Seurat)
library(plyr)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(EGSEA)
library(ggplot2)
library(extrafont)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(SingleR)
library(genefilter)
library(monocle3)
source("functions.R")

#Get Arial font
font_import(pattern = "Arial.ttf", prompt = F, recursive = F)
loadfonts()
loadfonts(device = "postscript")

#Set number of cells / plate.
n = 384

#-----------------------------------------------------------------------------------------
#Set up gene ontology objects for pathway analyses
ont_sets <- readList("cp_and_h.gmt.txt")

#Initialize bioMarts
marts <- listMarts()
mart <- marts[1,1]

human <- useMart(mart, dataset = "hsapiens_gene_ensembl")
mouse <- useMart(mart, dataset = "mmusculus_gene_ensembl")

#-----------------------------------------------------------------------------------------
### prep meta data and files
# load meta data
meta.data <- read.table('Project SCS 18 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
# rename first column
names(meta.data)[1] <- "Patient"
# add ID column
meta.data$ID <- paste(meta.data$Patient, ".P",meta.data$Plate, sep = "")
# transform patient colum from integer to character
meta.data$Patient <- as.character(meta.data$Patient)

sum(duplicated(meta.data$ID))
View(meta.data)

# check files
all.files <- list.files(pattern="*.tsv")
all.files
all.files %in% meta.data$File
length(all.files) == dim(meta.data)[1]

## read in and process
seurat.object.list <- list(c)
unwanted.genes <- c("^ERCC", "^MT-", "^RPL", "^RPS", "^KCNQ1OT1$", "^UGDH-AS1$", "^MALAT1$", "^EEF1A$")

for (df in seq_along(all.files)) {
  current.df <- read.delim(all.files[df])
  
  # fix gene names and rownames
  current.df$GENEID <- sub("__.*", "" ,  current.df$GENEID)
  rownames(current.df) <- current.df$GENEID
  current.df <- current.df[,-1]
  # remove unwanted genes
  current.df <- current.df[grep(paste(unwanted.genes, collapse = "|"),rownames(current.df),invert=TRUE),]
  
  # match AE and plate number to file
  current.meta.data <- meta.data[meta.data$File == all.files[df],]
  colnames(current.df) <- paste(current.meta.data$ID, ".", 1:length(current.df), sep = "")
  
  # update gene names
  current.df <- update.symbols(current.df, n)
  
  ## create object
  current.object <- CreateSeuratObject(current.df, project = current.meta.data$ID)
  # add corresponding meta data
  current.object@meta.data <- cbind(current.object@meta.data, current.meta.data)
  # add to list
  seurat.object.list[[df]] <- current.object
}