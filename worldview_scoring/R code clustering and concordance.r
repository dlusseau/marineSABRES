#############################
### clustering and concordance with CT and between scorers
### 11 Jul 2025
#############################
library(clue)  # consensus clustering
library(irr)   #  Fleiss' Kappa
library(cluster)  # For clustering algorithms
library(vegan) # for count data distances
library(stringr) # deal with text
library(fossil) #Rand index
library(aricode) #mutual information

scoring<-read.csv("/SABRES/worldviews ms/scoring.csv",header=T)
scoring$scenario<-sub(" > .*","",scoring$text)

scorers<-unique(scoring$coder)
scorer1 <- subset(scoring,coder==scorers[1])
scorer2 <- subset(scoring,coder==scorers[2])
scorer3 <- subset(scoring,coder==scorers[3])
scorer4 <- subset(scoring,coder==scorers[4])

scorer1 <-scorer1[,2:58]
scorer2 <-scorer2[,2:58]
scorer3 <-scorer3[,2:58]
scorer4 <-scorer4[,2:58]

pcscorer1<-scorer1[,names(which(colSums(scorer1[2:57])!=0))]
pcscorer2<-scorer2[,names(which(colSums(scorer2[2:57])!=0))]
pcscorer3<-scorer3[,names(which(colSums(scorer3[2:57])!=0))]
pcscorer4<-scorer4[,names(which(colSums(scorer4[2:57])!=0))]



# PCA
pca1 <- prcomp(pcscorer1,scale.=T)
pca2 <- prcomp(pcscorer2,scale.=T)
pca3 <- prcomp(pcscorer3,scale.=T)
pca4 <- prcomp(pcscorer4,scale.=T)

#not super let's go for Hellinger distance and hierarchical clustering with Ward (accounts for different variance of counts)


scorer1_hellinger <- decostand(scorer1[,2:57], method = "hellinger")  
row.names(scorer1_hellinger)<-scorer1[,1]
scorer1_dist <- dist(scorer1_hellinger, method = "euclidean")
scorer1_hc <- hclust(scorer1_dist, method = "ward.D2")
plot(scorer1_hc, main = "Hellinger")

scorer2_hellinger <- decostand(scorer2[,2:57], method = "hellinger")  
row.names(scorer2_hellinger)<-scorer2[,1]
scorer2_dist <- dist(scorer2_hellinger, method = "euclidean")
scorer2_hc <- hclust(scorer2_dist, method = "ward.D2")
plot(scorer2_hc, main = "Hellinger")

scorer3_hellinger <- decostand(scorer3[,2:57], method = "hellinger")  
row.names(scorer3_hellinger)<-scorer3[,1]
scorer3_dist <- dist(scorer3_hellinger, method = "euclidean")
scorer3_hc <- hclust(scorer3_dist, method = "ward.D2")
plot(scorer3_hc, main = "Hellinger")

scorer4_hellinger <- decostand(scorer4[,2:57], method = "hellinger")  
row.names(scorer4_hellinger)<-scorer4[,1]
scorer4_dist <- dist(scorer4_hellinger, method = "euclidean")
scorer4_hc <- hclust(scorer4_dist, method = "ward.D2")
plot(scorer4_hc, main = "Hellinger")

#let's cut in 4
clust1 <- cutree(scorer1_hc, k = 4)
clust2 <- cutree(scorer2_hc, k = 4)
clust3 <- cutree(scorer3_hc, k = 4)
clust4 <- cutree(scorer4_hc, k = 4)

# Combine cluster assignments 
cluster_assignments <- matrix(c(clust1, clust2, clust3, clust4), nrow = length(clust1), ncol = 4)

#  Fleiss' Kappa
fleiss_kappa <- kappam.fleiss(cluster_assignments,detail=TRUE)
fleiss_kappa$detail
fleiss_kappa
#average kappa is 0.48 so moderate agreement - good!
#even better for scorer 3 and 4

#  pairwise Adjusted Rand Index

ari_12 <- adj.rand.index(clust1, clust2)
ari_13 <- adj.rand.index(clust1, clust3)
ari_14 <- adj.rand.index(clust1, clust4)
ari_23 <- adj.rand.index(clust2, clust3)
ari_24 <- adj.rand.index(clust2, clust4)
ari_34 <- adj.rand.index(clust3, clust4)

# Average ARI across all pairs
ari_values <- c(ari_12, ari_13, ari_14, ari_23, ari_24, ari_34)
mean_ari <- mean(ari_values)
#not far from Kappa
#scorer 1 and 2 more agreement, scorer 1 and 3 least agreement

# Optional: Compute Normalized Mutual Information (NMI) using 'aricode' package

nmi_12 <- NMI(clust1, clust2)
nmi_13 <- NMI(clust1, clust3)
nmi_14 <- NMI(clust1, clust4)
nmi_23 <- NMI(clust2, clust3)
nmi_24 <- NMI(clust2, clust4)
nmi_34 <- NMI(clust3, clust4)

# Average NMI
nmi_values <- c(nmi_12, nmi_13, nmi_14, nmi_23, nmi_24, nmi_34)
mean_nmi <- mean(nmi_values)

#same mean information
# scorer 1 and 4 the closest here

### Kappa is best to use here


#################################
##### consensus ensemble clustering


hc_ensemble<-list(scorer1_hc,scorer2_hc,scorer3_hc,scorer4_hc)
names(hc_ensemble)<-c("scorer1","scorer2","scorer3","scorer4")

ensemble_ready <- cl_ensemble(list = hc_ensemble)

consensus_cluster<-cl_consensus(ensemble_ready,method="cophenetic") #we have euclidean on hellinger so should be good
plot(consensus_cluster, main = "Consensus Dendrogram",horiz=TRUE)


sum(cl_dissimilarity(ensemble_ready, consensus_cluster, "cophenetic")^2)
#alright
class_ids <- cutree(as.hclust(consensus_cluster), k = 4)
table(class_ids)

library(ggdendro)
library(ggplot2)
library(dendextend)

colours<-c("#000000","#E69F00","#56B4E9","#009E73","#0072B2","#D55E00")

attr(consensus_cluster[[1]],"Labels")<-sub("_Be", "",attr(consensus_cluster[[1]],"Labels"))
attr(consensus_cluster[[1]],"Labels")<-sub(">", "-", attr(consensus_cluster[[1]],"Labels"))


thecluster<- as.dendrogram(as.hclust(consensus_cluster))

order.dendrogram(thecluster)


colourss<-colours[as.numeric(factor(sub("-.*", "", thecluster$labels)))]

labs<-thecluster$labels[thecluster$order[seq(length(thecluster$order),1,-1)]]
labs.col<-colourss[thecluster$order[seq(length(thecluster$order),1,-1)]]
names(labs.col)<-sub("-.*", "", labs)

# lab_html <- sprintf("<span style='color:%s;'>%s</span>", labs.col, labs)
# html_unescape <- function(x) {
  # x <- gsub("&amp;", "&", x, fixed = TRUE)  # do &amp; first, in case others were double-escaped
  # x <- gsub("&lt;",  "<", x, fixed = TRUE)
  # x <- gsub("&gt;",  ">", x, fixed = TRUE)
  # x
# }
# lab_html <- html_unescape(lab_html)
# lab_html <- c(rep("<span style='color:red;'>Test</span>"),28)

thecluster<-color_labels(thecluster, labs.col)


dendro<-ggdendrogram(thecluster, rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Consensus Dendrogram - Fleiss Kappa = 0.47 (4 scorers)", x = "Distance")
  
p_col <-dendro+
scale_y_discrete(labels = lab_html) +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.y.left = element_text(size=10)) 
		
png("/SABRES/worldviews ms/colour_consensus dendogram.png",width=30,height=15,units="cm",res=300)
p_col
 dev.off()






png("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/worldviews ms/consensus dendogram.png",width=30,height=15,units="cm",res=300)
ggdendrogram(as.hclust(consensus_cluster), rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Consensus Dendrogram - Fleiss Kappa = 0.47 (4 scorers)", x = "Distance")
  
dev.off()


png("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/worldviews ms/scorer1 dendogram.png",width=30,height=15,units="cm",res=300)
ggdendrogram(scorer1_hc, rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Dendrogram - scorer 1", x = "Distance")
dev.off()

png("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/worldviews ms/scorer2 dendogram.png",width=30,height=15,units="cm",res=300)
ggdendrogram(scorer2_hc, rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Dendrogram - scorer 2", x = "Distance")
dev.off()

png("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/worldviews ms/scorer3 dendogram.png",width=30,height=15,units="cm",res=300)
ggdendrogram(scorer3_hc, rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Dendrogram - scorer 3", x = "Distance")
dev.off()

png("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/worldviews ms/scorer4 dendogram.png",width=30,height=15,units="cm",res=300)
ggdendrogram(scorer4_hc, rotate = TRUE, theme_dendro = TRUE) +
  labs(title = "Dendrogram - scorer 4", x = "Distance")
dev.off()


