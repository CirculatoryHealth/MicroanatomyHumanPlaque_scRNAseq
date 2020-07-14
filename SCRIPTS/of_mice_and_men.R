#Read top20 markers from the mouse studies
mouse.markers <- read.table(file = "mouse_top20_markers.txt", header = T)

#Convert to human
human.mouse.markers <- list()
for(i in 1:dim(mouse.markers)[2]){
  print("Working on mouse set:")
  print(i)
  human.mouse.markers[colnames(mouse.markers)[i]] <- getLDS(attributes = "mgi_symbol",
                                        filters = "mgi_symbol",
                                         values = mouse.markers[, i],
                                           mart = mouse,
                                    attributesL = "hgnc_symbol",
                                          martL = human
  )[2]
}

#Retrieve top markers from our mac populations
human.markers <- list(My.0 = subset(all.all.seur.combined.M_clusters.markers, cluster == 0 & p_val_adj < 0.1)[,"gene"],
                      My.1 = subset(all.all.seur.combined.M_clusters.markers, cluster == 1 & p_val_adj < 0.1)[,"gene"],
                      My.2 = subset(all.all.seur.combined.M_clusters.markers, cluster == 2 & p_val_adj < 0.1)[,"gene"])

#Compute overlaps
mouse.human.marker.overlap <- t(as.data.frame(lapply(human.markers, function(x)sum(grepl(paste(human.mouse.markers[[1]], collapse = "|"), x)))))
for (i in 2:length(human.mouse.markers)){
  mouse.human.marker.overlap <- cbind(mouse.human.marker.overlap, t(as.data.frame(lapply(human.markers, function(x)sum(grepl(paste(human.mouse.markers[[i]], collapse = "|"), x))))))
}
colnames(mouse.human.marker.overlap) <- names(human.mouse.markers)
mouse.human.marker.overlap

#Compute significance of overlap
total.genes <- as.character(mouse.markers[,1])
for (i in 2:length(human.mouse.markers)){
  total.genes <- c(total.genes, as.character(mouse.markers[,i]))
}
total.genes <- c(total.genes, human.markers[[1]])
total.genes <- c(total.genes, human.markers[[2]])
total.genes <- c(total.genes, human.markers[[3]])
total.genes <- toupper(total.genes)
#total.genes <- unique(total.genes)
total.genes <- length(total.genes)
#total.genes <- nrow(all.seur.combined@raw.data)
total.genes

p.mouse.human.marker.overlap <- data.frame()
for (i in row.names(mouse.human.marker.overlap)){
  for (j in colnames(mouse.human.marker.overlap)){
    cat(i, j, length(mouse.markers[,j]), length(human.markers[[i]]), "\n")
    p.mouse.human.marker.overlap[i,j] <- phyper(mouse.human.marker.overlap[i,j] - 1, 
                                                length(mouse.markers[,j]), 
                                                total.genes - length(mouse.markers[,j]), 
                                                length(human.markers[[i]]),
                                                lower.tail = F
                                               )
  }
}
p.mouse.human.marker.overlap
rowSums(p.mouse.human.marker.overlap < 0.05)
write.table(p.mouse.human.marker.overlap, file = "mice.overlap.all.clusters.txt", quote = F, sep = "\t")

#Export lists
library(clipr)
write_clip(human.markers$My.0)
write_clip(human.markers$My.1)
write_clip(human.markers$My.2)


#Make 'go-like' table from enriched populations
#mice.My0 <- read.table(file="mice.My.0.txt", header =T, sep = "\t")
#plot_GO_from_table(mice.My0, name = "mice.My0.pdf")

mice.My1 <- read.table(file="mice.My.1.unique.txt", header =T, sep = "\t")
plot_GO_from_table(mice.My1, name = "mice.My1.unique.pdf")

mice.My2 <- read.table(file="mice.My.2.unique.txt", header =T, sep = "\t")
plot_GO_from_table(mice.My2, name = "mice.My2.unique.pdf")

#Prep our vars
#mice.My0$hum_clus <- "My 0"
mice.My1$hum_clus <- "My 1"
mice.My2$hum_clus <- "My 2"
#m <- mice.My0
#m <- rbind(m, mice.My1)
m <- mice.My1
m <- rbind(m, mice.My2)
m[,4]       <- -log10(m[,2])

colnames(m) <-  c("clus","pval","hum_clus","log10_pval")
sig_cut           <- -log10(0.05)

#Calculate how wide the labels will be
label_width <- 1 + (max(strwidth(m$clus, units = "inches")) * 2.54) #Converting to cm

#Make the plot!
ggplot(m, aes(x = reorder(clus, log10_pval)  , y = log10_pval, fill = hum_clus)) +
  geom_bar(stat="identity", position = "dodge") + 
  #scale_fill_manual(values = c("#F2756D", "#B6A031", "#2AB34B")) +
  scale_fill_manual(values = c("#B6A031", "#2AB34B")) +
  coord_flip() +
  ylab("-log10(p-value)") +
  xlab("") +
  #ylim(0,6) +
  geom_hline(yintercept = sig_cut, col = "red", linetype = "dashed") +
  theme(panel.background = element_blank(),
        axis.ticks.y     = element_blank(),
        plot.margin      = unit(c(1,1,1,1), "cm"),
        text             = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain")
  )

  #And save it
ggsave("mice.My.unique.pdf",width = (label_width / 2.54) + 7.36)


#===================================================
#Make 'go-like' table from ALL populations
mice.My0 <- read.table(file="Mice.My.0.all.ref.txt", header =T, sep = "\t")
plot_GO_from_table(mice.My0, name = "mice.My0.all.ref.pdf")

mice.My1 <- read.table(file="mice.My.1.all.ref.txt", header =T, sep = "\t")
plot_GO_from_table(mice.My1, name = "mice.My1.all.ref.pdf")

mice.My2 <- read.table(file="mice.My.2.all.ref.txt", header =T, sep = "\t")
plot_GO_from_table(mice.My2, name = "mice.My2.all.ref.pdf")

#Prep our vars
mice.My0$hum_clus <- "My 0"
mice.My1$hum_clus <- "My 1"
mice.My2$hum_clus <- "My 2"
m <- mice.My0
m <- rbind(m, mice.My1)
m <- rbind(m, mice.My2)
m[,4]       <- -log10(m[,2])
colnames(m) <-  c("clus","pval","hum_clus","log10_pval")
sig_cut           <- -log10(0.05)

#Calculate how wide the labels will be
label_width <- 1 + (max(strwidth(m$clus, units = "inches")) * 2.54) #Converting to cm

#Make the plot!
ggplot(m, aes(x = reorder(clus, log10_pval)  , y = log10_pval, fill = hum_clus)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("#F2756D", "#B6A031", "#2AB34B")) +
  coord_flip() +
  ylab("-log10(p-value)") +
  xlab("") +
  #ylim(0,6) +
  geom_hline(yintercept = sig_cut, col = "red", linetype = "dashed") +
  theme(panel.background = element_blank(),
        axis.ticks.y     = element_blank(),
        plot.margin      = unit(c(1,1,1,1), "cm"),
        text             = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain")
  )

#And save it
ggsave("mice.My.all.pdf",width = (label_width / 2.54) + 7.36)

