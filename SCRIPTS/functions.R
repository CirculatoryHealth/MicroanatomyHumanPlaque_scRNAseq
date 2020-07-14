merge.on.rownames <- function(x, y){
  m <- merge(x, y, by = 0)
  row.names(m) <- m$Row.names
  m$Row.names <- NULL
  return(m)
}

#-----------------------------------------------------------------------------------------

update.symbols <- function(my.counts.table = my.counts.table, n = n){
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping
  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(my.counts.table), columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- gene.info[!d,]
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- gene.info$SYMBOL
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  #Keep only highest expressed SYMBOL if multiple old aliases are now merged
  #Get maximum number of dups per gene
  num_dups <- max(table(my.counts.table$SYMBOL))
  
  #Loop until only one entry per gene left
  while(num_dups > 1){
    #First retrieve the list of duplicated symbols
    o <- order(my.counts.table$SYMBOL)
    my.counts.table <- my.counts.table[o,]
    d <- duplicated(my.counts.table$SYMBOL)
    
    #Then keep only the highest expressed one (this part will fail if we happen upon 2 rows with equal expression, but for now that hasn't happened yet so meh ;)
    h <- rowSums(my.counts.table[d,1:n]) > rowSums(my.counts.table[which(d==T)-1,1:n])
    h <- c(h, rowSums(my.counts.table[which(d==T)-1,1:n]) > rowSums(my.counts.table[d,1:n]))
    h <- h[which(!h)]
    my.counts.table <- my.counts.table[!row.names(my.counts.table) %in% names(h),]
    
    num_dups <- num_dups - 1 #One dup removed
  }
  
  #Overwrite old symbols in row names and clean up
  row.names(my.counts.table) <- my.counts.table$SYMBOL
  my.counts.table$SYMBOL <- NULL
  
  return(my.counts.table)
}


get_ontology <- function(res, name = "cluster", outdir = ".", return.data = F, full_GSEA = T){
  #Set up our data
  #Get background of all expressed genes in our dataset
  universe <- row.names(all.seur.combined@raw.data)
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID
  
  #Extract genes from res
  genes <- data.frame(Symbol = res[,"gene"])
  
  #Return emptyhanded if empty cluster
  if (nrow(genes)==0){
    return("NA")
  }

  #Retrieve entrez IDs
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = as.vector(genes$Symbol), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  colnames(genes) <- c("symbol", "entrezID")
  genes <- genes[,c(2,1)]
  
  #Add logFC column
  genes$logFC <- res[,"avg_logFC"]
  
  #-------------------------------------------------------
  #First run the custom collection and GO to make neat bar graphs
  #CP & H
  #Create annotation table
  gs.annots <- buildCustomIdx(genes$entrezID,
                              species = "human",
                              gsets = ont_sets,
                              label = "CP_and_H",
                              name = "CP_and_H"
  )
  
  #Get pathway enrichments
  gsa.pathways <- egsea.ora(geneIDs = genes$entrezID,
                            title = "x", 
                            universe = universe,
                            gs.annots = gs.annots,
                            symbolsMap = genes,
                            sort.by = "p.adj",
                            report.dir = ".",
                            num.threads = 4,
                            report = F
  )
  
  #GO
  #Create annotation table
  gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezID, 
                              species = "human", 
                              geneSets = "c5"
  )
  
  #Get GO enrichments
  gsa.go <- egsea.ora(geneIDs = genes$entrezID,
                      gs.annots = gs.annots,
                      title = "x",
                      universe = universe,
                      symbolsMap = genes,
                      sort.by = "p.adj",
                      report.dir = ".",
                      num.threads = 4,
                      report = F
  )
  
  #Extract top 20 enriched sets
  top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)),
                             padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
  top.GO       <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), 
                             padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
  
  #Check if any sets are significant, otherwise return emptyhanded
  if (!sum(unlist(lapply(c(top.pathways$padj, top.GO$padj), function(x) any(x<0.1))))){
    cat("No significant sets found! Exiting...\n")
    return("NA")
  }
  
  #Plot!
  plot_GO_from_table(top.pathways, name = paste(outdir, "/", name, ".pathways.pdf", sep = ""))
  plot_GO_from_table(top.GO,       name = paste(outdir, "/", name, ".GO_terms.pdf", sep = ""))
  
  #-------------------------------------------------------
  #Now run the full GSEA suite and make a report dir, if wanted
  if(full_GSEA){
    #Build annotation index
    gs.annots <- buildIdx(entrezIDs = genes$entrezID,
                          species = "human",
                          msigdb.gsets = "all",
                          gsdb.gsets = "all",
                          go.part = T,
                          kegg.updated = T
    )
    
    #Run EGSEA
    egsea.ora(geneIDs = genes$entrezID,
              title = name,
              universe = universe,
              gs.annots = gs.annots,
              logFC = genes$logFC,
              symbolsMap = genes,
              display.top = 20,
              sort.by = "p.adj",
              report.dir = paste(outdir, "/", name, "_EGSEA", sep = ""),
              kegg.dir = paste(outdir, "/", name, "_EGSEA/kegg-dir", sep = ""),
              num.threads = 4,
              interactive = F,
              report = T,
              verbose = F
    )
  }
  #Return the tables if wanted
  if(return.data){
    cat("Returning Data...\n")
    pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                           padj = topSets(gsa.pathways, number = 500, names.only = F)$p.adj)
    goterms  <- data.frame(name = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                           padj = topSets(gsa.go, number = 500, names.only = F)$p.adj)
    returnList <- list(pathways,goterms)
    names(returnList) <- c("Pathways", "GO_terms")
    return(returnList)
  }
} #get_ontology(res)

get_ontology_no_direction <- function(thegenes, name = "cluster", outdir = ".", return.data = F, full_GSEA = T){
  #Set up our data
  #Get background of all expressed genes in our dataset
  universe <- row.names(all.seur.combined@raw.data)
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID
  
  #Extract genes from res
  genes <- data.frame(Symbol = thegenes)
  
  #Return emptyhanded if empty cluster
  if (nrow(genes)==0){
    return("NA")
  }
  
  #Retrieve entrez IDs
  genes <- AnnotationDbi::select(org.Hs.eg.db, keys = as.vector(genes$Symbol), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  colnames(genes) <- c("symbol", "entrezID")
  genes <- genes[,c(2,1)]
  
  #-------------------------------------------------------
  #First run the custom collection and GO to make neat bar graphs
  #CP & H
  #Create annotation table
  gs.annots <- buildCustomIdx(genes$entrezID,
                              species = "human",
                              gsets = ont_sets,
                              label = "CP_and_H",
                              name = "CP_and_H"
  )
  
  #Get pathway enrichments
  gsa.pathways <- egsea.ora(geneIDs = genes$entrezID,
                            title = "x", 
                            universe = universe,
                            gs.annots = gs.annots,
                            symbolsMap = genes,
                            sort.by = "p.adj",
                            report.dir = ".",
                            num.threads = 4,
                            report = F
  )
  
  #GO
  #Create annotation table
  gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezID, 
                              species = "human", 
                              geneSets = "c5"
  )
  
  #Get GO enrichments
  gsa.go <- egsea.ora(geneIDs = genes$entrezID,
                      gs.annots = gs.annots,
                      title = "x",
                      universe = universe,
                      symbolsMap = genes,
                      sort.by = "p.adj",
                      report.dir = ".",
                      num.threads = 4,
                      report = F
  )
  
  #Extract top 20 enriched sets
  top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)),
                             padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
  top.GO       <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), 
                             padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
  
  #Check if any sets are significant, otherwise return emptyhanded
  if (!sum(unlist(lapply(c(top.pathways$padj, top.GO$padj), function(x) any(x<0.1))))){
    return("NA")
  }
  
  #Plot!
  plot_GO_from_table(top.pathways, name = paste(outdir, "/", name, ".pathways.pdf", sep = ""))
  plot_GO_from_table(top.GO,       name = paste(outdir, "/", name, ".GO_terms.pdf", sep = ""))
  
  #-------------------------------------------------------
  #Now run the full GSEA suite and make a report dir, if wanted
  if(full_GSEA){
    #Build annotation index
    gs.annots <- buildIdx(entrezIDs = genes$entrezID,
                          species = "human",
                          msigdb.gsets = "all",
                          gsdb.gsets = "all",
                          go.part = T,
                          kegg.updated = T
    )
    
    #Run EGSEA
    egsea.ora(geneIDs = genes$entrezID,
              title = name,
              universe = universe,
              gs.annots = gs.annots,
              symbolsMap = genes,
              display.top = 20,
              sort.by = "p.adj",
              report.dir = paste(outdir, "/", name, "_EGSEA", sep = ""),
              kegg.dir = paste(outdir, "/", name, "_EGSEA/kegg-dir", sep = ""),
              num.threads = 4,
              interactive = F,
              report = T,
              verbose = F
    )
  }
  #Return the tables if wanted
  if(return.data){
    cat("Returning Data...\n")
    pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 500, names.only = F)),
                           padj = topSets(gsa.pathways, number = 500, names.only = F)$p.adj)
    goterms  <- data.frame(name = row.names(topSets(gsa.go, number = 500, names.only = F)), 
                           padj = topSets(gsa.go, number = 500, names.only = F)$p.adj)
    returnList <- list(pathways,goterms)
    names(returnList) <- c("Pathways", "GO_terms")
    return(returnList)
  }
} #get_ontology_no_direction(thegenes)

#-----------------------------------------------------------------------------------------
#Calculate GO and pathway enrichment
#Define a function to plot the GO graphs
plot_GO_from_table <- function(go_data, name="GO_terms.pdf", sig_cut=0.05){
  #Prep our vars
  go_data[,1]       <-  gsub("_", " ", go_data[,1])
  go_data[,3]       <- -log10(go_data[,2])
  colnames(go_data) <-  c("term","padj","log10_FDR_qvalue")
  sig_cut           <- -log10(sig_cut)
  
  #Calculate how wide the labels will be
  label_width <- 1 + (max(strwidth(go_data$term, units = "inches")) * 2.54) #Converting to cm
  #print(label_width) #old debug line
  
  #Make the plot!
  ggplot(go_data,aes(x = reorder(term, log10_FDR_qvalue)  , y = log10_FDR_qvalue)) +
    geom_bar(stat="identity", fill = "Dodgerblue") + 
    coord_flip() +
    ylab("-log10(FDR q-value)") +
    xlab("") +
    ylim(0,roundCeiling(x = max(go_data$log10_FDR_qvalue))) +
    geom_hline(yintercept = sig_cut, col = "red", linetype = "dashed") +
    theme(panel.background = element_blank(),
          axis.ticks.y     = element_blank(),
          plot.margin      = unit(c(1,1,1,1), "cm"),
          text             = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
    )
    
  #And save it
  ggsave(name,width = (label_width / 2.54) + 7.36)
} #plot_GO_from_table()

#Define helper function for rounding the axis of the plots
roundCeiling <- function(x) {
  if(length(x) != 1) stop("'x' must be of length 1")
  if(x < 5){
    return(5)
  }
  else if(x < 10){
    return(10)
  }
  else if(x < 15){
    return(15)
  }
  else if(x < 100){
    round(x+5,-1)
  }
  else if(x < 1000){
    round(x+50,-2)
  }
  else{
    round(x+500,-3)
  }
} #roundCeiling()

#wrapper for write.table so that it makes a header without blank tab at the start when using row names
write.better.table <- function(data, file = "myfile.txt"){
  write.table(transform(data,Symbol = rownames(data)
                       )[,c(length(colnames(data))+1,
                            1:length(colnames(data))
                           )
                        ],
              file      = file, 
              quote     = F, 
              row.names = F,
              sep       = "\t"
              )
}

#Retrieve all cell names (not) expressing a certain gene
GetCellNames <- function(myObject = all.seur.combined, theGene, treshold = 0, invert = F){
  m        <- as.matrix(myObject@data)
  theRows  <- grep(theGene, 
                   row.names(m)
                   )
  if(invert){
    theCells <- colnames(m[,
                           m[theRows,] <= treshold
                          ]
                      )
  } else {
    theCells <- colnames(m[,
                           m[theRows,] > treshold
                          ]
    )
  }
  return(theCells)
}
