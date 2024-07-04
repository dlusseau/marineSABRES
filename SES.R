library(igraph)

links<-as.data.frame(read_excel("Simplified final map SES_TurismSimplified.xlsx", 
                                sheet = "Connections"))
elements<-as.data.frame(read_excel("Simplified final map SES_TurismSimplified.xlsx", 
                                   sheet = "Elements"))
net<-graph_from_edgelist(as.matrix(links[,1:2]))

plot(net)
V(net)$Description<-elements[match(V(net)$name, elements$Label),]$Description
links$strength2<-NA
links[grepl("Strong", links$strength),]$strength2<-"Strong"
links[grepl("Medium", links$strength),]$strength2<-"Medium"
links[grepl("Weak", links$strength),]$strength2<-"Weak"

E(net)$strength<-links$strength2
E(net)$sign<-links$Label
E(net)$confidence<-sample(c(1,2,3), size=length(as.vector(E(net)$sign)), replace=T)
pal<-paletteer_d("RColorBrewer::Spectral")[c(2,3,5,7,9,10)]
library(ggplot2)
library(ggnetwork)
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend, fill=Description, color=as.factor(sign), label=name, size=as.factor(strength), linetype=as.factor(confidence))) +
  geom_edges(arrow = arrow(length = unit(4, "pt"), type = "closed"), curvature=0.1) +
  geom_nodes(color="black",size = 10, pch=21) +
  geom_text(color="black",size = 3,lineheight=0.1)+
  scale_fill_manual(values=pal)+
  scale_size_manual(values=c(1,2,3))+
  scale_color_manual(values=c( "#E44244","lightblue"))+
  theme_blank()+theme(legend.title = element_blank(), legend.position = "bottom")

igraph::degree(net)

library(LoopAnalyst)
mat<-as.matrix(as_adjacency_matrix(net))
loops<-enumerate.loops(mat)
