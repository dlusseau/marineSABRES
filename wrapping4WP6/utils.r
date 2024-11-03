
#######
## list of utility functions to walk through interventions SABRES analytical pipelines
### Oct 2024
################

make_matrix<-function(from,to,weight) {

elements<-unique(c(unique(from),unique(to)))

SES.mat<-matrix(0,length(elements),length(elements))
colnames(SES.mat)<-elements
row.names(SES.mat)<-elements

for (i in 1:length(from)) {

SES.mat[which(row.names(SES.mat)==from[i]),which(colnames(SES.mat)==to[i])]<-weight[i]

}
return(SES.mat)
}


######################################


SES.laplacian<-function(SES.mat,from=c("rows","cols")) {	
if (from[1]=="rows") {				
SES.Lap<-diag(rowSums(t(SES.mat))) -t(SES.mat) 
} else {
SES.Lap<-diag(rowSums((SES.mat))) -(SES.mat) 
}
lap.e<-eigen(SES.Lap)$value
names(lap.e)<-row.names(SES.mat)

return(lap.e)

}

######################################

plot.SES<-function(SES.mat,folder,filename,title,w=80,h=40,layout=layout_with_fr,label.cex=1.5,vsize=20,eweight=10) {
require(igraph)

SES.net<-graph_from_adjacency_matrix(
  SES.mat,
  mode = "directed",
  weighted = TRUE)

V(SES.net)$color<-"white"
V(SES.net)$label.cex<-label.cex
V(SES.net)$size<-vsize
E(SES.net)$sign<-sign( E(SES.net)$weight)
E(SES.net)$color<-"blue"
E(SES.net)$color[E(SES.net)$sign<0]<-"red"
E(SES.net)$weight<-abs( E(SES.net)$weight)
E(SES.net)$curved <- 0.2
E(SES.net)$width <- E(SES.net)$weight * eweight
E(SES.net)$arrow.size <-E(SES.net)$weight * (eweight/2)


l <- layout(SES.net)

png(paste0(folder,filename,".png"),width=w,height=h,units="cm",res=200)
plot(SES.net,layout=l,ylab = title)
dev.off()

return(SES.net)

}

#####################################################

boolean_file_creation<-function(SES.mat,folder,filename) {

SES.bin<-sign(SES.mat)
colnames(SES.bin)<-gsub("\\s*\\([^\\)]+\\)","",colnames(SES.bin))
row.names(SES.bin)<-gsub("\\s*\\([^\\)]+\\)","",row.names(SES.bin))
colnames(SES.bin)<-gsub("[^[:alnum:]]","",colnames(SES.bin))
row.names(SES.bin)<-gsub("[^[:alnum:]]","",row.names(SES.bin))

colnames(SES.bin)<-gsub(" ","_",colnames(SES.bin))
row.names(SES.bin)<-gsub(" ","_",row.names(SES.bin))


boolean.df<-data.frame(targets=factor(colnames(SES.bin)),factors=NA) # to is the columns

for (i in 1:ncol(SES.bin)) {
poss<-names(which(SES.bin[,i]==1))
negs<-names(which(SES.bin[,i]==-1))
if (length(negs)>0) {
negs<-paste0("!",negs)
}
all<-c(poss,negs)

boolean.df$factors[i]<-paste(all,collapse="|")
}
write.csv(boolean.df,file=paste0(folder,filename,".csv"),row.names = F,quote=FALSE)
}


##############################################################


boolean_analyses<-function(boolean.net,folder,filename) {
require(BoolNet)
require(igraph)
require(stringr)
 
states<-getAttractors(Tuscany_boolean,returnTable=TRUE)

#############################################
## graphing begins
state.map<-plotStateGraph(states, layout = layout.fruchterman.reingold,plotIt = FALSE)
vertex_attr(state.map,index=1)
table(V(state.map)$color)
V(state.map)$degree<-degree(state.map)
V(state.map)$attractor<-"no"
V(state.map)$attractor[V(state.map)$frame.color=="#000000FF"]<-"yes"
V(state.map)$color2<-str_sub(V(state.map)$color, end=-3)
E(state.map)$width<-1
E(state.map)$lty<-1
E(state.map)$color2<-str_sub(E(state.map)$color, end=-3)
core.state<-subgraph(state.map,which(V(state.map)$degree>1)) # pruning the leaves
#plot(core.state,layout=layout.fruchterman.reingold,vertex.size=2,edge.arrow.size=0.5,edge.width=1)
write_graph(
  core.state,
  file=paste0(folder,filename,".graphml"),
  format = "graphml")

## graphing ends
##############################################  
n_states<-length(states$stateInfo$table)
n_attractors<-length(states$attractors)

attractors<-list(n_attractors)
basin_attraction<-array(0,n_attractors)

for (i in 1:n_attractors) {
attractors[[i]]<-as.data.frame(getAttractorSequence(states,i))
basin_attraction[i]<-states$attractors[[i]]$basinSize
}
 
return(list(states=n_states,N_attractors=n_attractors,basins=basin_attraction,attractors=attractors))
 
}

######################################################


#### plot exisintg weighted dynamics


SES.simulate<-function(SES.mat,iter,folder,filename,title) {
require(reshape2)
require(grid)
require(gridExtra)
require(Polychrome)
require(ggplot2)
require(ggfortify)
require(ggpubr)


SES.sim<-matrix(NA,nrow(SES.mat),iter)
SES.sim[,1]<-runif(nrow(SES.mat),0,1)

for ( i in 2:iter) {
SES.sim[,i]<-t(SES.mat)%*%matrix((SES.sim[,i-1]), ncol = 1)  #column are from/jacobian
}
row.names(SES.sim)<-row.names(SES.mat)
sim.melt<-melt(SES.sim)

fillcol<-glasbey.colors(nrow(SES.mat))
names(fillcol)<-unique(sim.melt$Var1)

p1a<-ggplot(sim.melt,aes(x=log10(Var2),y=(value),colour=(Var1)))+
 geom_path(size=2)+
 labs(colour="Element", y="Progress",x="log(time)")+
 scale_colour_manual(values=fillcol)+
 theme_minimal(base_size=18)+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
 
pca<-prcomp(t(SES.sim),scale=FALSE)

### autoplot looks like
p1b<-ggplot2::autoplot(pca, label = FALSE, loadings.colour="black",colour="blue",loadings.label = TRUE,loadings.label.size = 6,loadings.label.colour = "black")+
geom_path(aes(x=PC1,y=PC2),colour="blue",
		arrow = arrow(type = "closed", angle = 30, length = unit(5, "mm"))
		)+
		labs(x="PC1",y="PC2")+
theme_minimal()

text1 <- textGrob(title, just="centre",rot=90,gp = gpar(fontsize = 20))


png(paste0(folder,filename,".png"),width=60,height=35,units="cm",res=200)
ggarrange(text1,p1a,p1b,nrow=1,ncol=3,widths=c(.1,1,1))
dev.off()

return(SES.sim)

}

###########################################
## partcipation ratio

jacobian<-function(SES) {
SES.j<-t(SES)
return(SES.j)
}
lefteigen<-function(mat) {
left<-eigen(t(mat))$vector
return(left)
}
righteigen<-function(mat) {
right<-eigen((mat))$vector
return(right)
}


participation_ratio<-function(SES,folder,filename,title) { 
require(Polychrome)
require(ggplot2)
#Duan et al. 2024
#careful again colSums is the equivalent of sum(, axis=0) in python NOT rowsums
# and careful again * in python is elementwise multiplication not matrix multiplication proper
### here I implement the code in DUan's github (distance.py line 166) however I don't think I reconcile it with the equation in the methods (Eq 5)!
## actually line 166 is really messed up, / x^-1 ?
##back to the equation instead
SES.j<-jacobian(SES)
LJ<-lefteigen(SES.j)
RJ<-righteigen(SES.j)
PR=rowSums(RJ*LJ)^2/rowSums((RJ*LJ)^2)

PR.df<-data.frame(group=title,components=colnames(SES),PR=PR)
fillcol<-glasbey.colors(nrow(SES))
names(fillcol)<-colnames(SES)
plot1<-ggplot(PR.df,aes(y=Re(PR),x=components,fill=components))+
geom_bar(stat="identity")+
scale_fill_manual(values=fillcol)+
theme_minimal(base_size=20)+
theme(legend.position="none")+
xlab(title)+
ylab("Participation ratio")+
theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(paste0(folder,filename,".png"),plot=plot1,device="png",width=60,height=30,units="cm",dpi=200,bg="white")
#dev.off()

return(PR.df)


}


#####################################################################
#### simulate addition (more to come on simulations overall


simulate.measure<-function(mat,measure,affected,indicators,lower,upper) { #generalise so we can have a mixture of + & - values for lowerS and upperS
measure.ef<-runif(length(affected),lower,upper)
newrow<-matrix(rep(0,ncol(mat)),nrow=1)
newrow[match(affected,colnames(mat))]<-measure.ef
mat<-rbind(mat,newrow)
measure.af<-runif(length(indicators),lower,upper)
newcol<-matrix(rep(0,nrow(mat)),ncol=1)
newcol[match(indicators,row.names(mat))]<-measure.af
mat<-cbind(mat,newcol)

row.names(mat)[nrow(mat)]<-colnames(mat)[ncol(mat)]<-measure

return(mat)
}

