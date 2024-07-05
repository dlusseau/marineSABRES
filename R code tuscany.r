#let's start with tuscany 
tuscany<-read.csv("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/tuscany_kumu_v27.csv",header=T)

net.size<-length(unique(c(tuscany$From,tuscany$To)))
nodes<-unique(c(tuscany$From,tuscany$To))

weights<-unique(tuscany$strength)

tuscany$weight<-1
tuscany$weight[tuscany$strength==weights[3]]<-.5
tuscany$weight[tuscany$strength==weights[4]]<-.25
tuscany$weight[tuscany$strength==weights[7]]<-.25 #typo in original data
tuscany$weight[tuscany$strength==weights[6]]<-.25 #NAs

tuscany$weight[tuscany$strength==weights[5]]<-(-.25)
tuscany$weight[tuscany$strength==weights[1]]<-(-.5)
tuscany$weight[tuscany$strength==weights[2]]<-(-1)
tuscany$weight[tuscany$sign=="-"&tuscany$strength==""]<-(-.25)

tuscany$edge<-1

tuscany$edge[tuscany$sign=="-"]<-(-1)



CommunityMatrix <- matrix(0, net.size, net.size, dimnames = list(nodes,nodes))

#i is To j is From X[i,j]

for (r in 1:nrow(tuscany)) {

CommunityMatrix[which(rownames(CommunityMatrix)==tuscany$To[r]),which(colnames(CommunityMatrix)==tuscany$From[r])]<-tuscany$edge[r]
}

CommunityMatrix<-t(CommunityMatrix) #row produce for columns


A<-matrix(rep(0.00001,nrow(CommunityMatrix)), ncol = 1)

sims<-matrix(0,ncol=nrow(CommunityMatrix),nrow=1000)
colnames(sims)<-colnames(CommunityMatrix)
sims[1,]<-A

for (i in 1:999) {
sims[i+1,]<-CommunityMatrix%*%matrix((sims[i,]), ncol = 1)
}

sim.sign<-sign(sims)
library(reshape2)
sim.melt<-melt(sims)

library(ggplot2)


ggplot(sim.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_path()




CM.weight <- matrix(0, net.size, net.size, dimnames = list(nodes,nodes))

for (r in 1:nrow(tuscany)) {

CM.weight[which(rownames(CM.weight)==tuscany$To[r]),which(colnames(CM.weight)==tuscany$From[r])]<-tuscany$weight[r]
}

CM.weight<-t(CM.weight)


A<-matrix(rep(1,nrow(CommunityMatrix)), ncol = 1)
simw<-matrix(0,ncol=nrow(CommunityMatrix),nrow=1000)
colnames(simw)<-colnames(CM.weight)
simw[1,]<-A

for (i in 1:999) {
simw[i+1,]<-CM.weight%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_path()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
hlabs(colour="SES components")+
theme_minimal()



sim.sub.original<-subset(simw.melt,Var2=="P.oceanica meadows (km2)" | Var2=="Area of P. oceanica habitat disturbed or lost" | Var2=="Recreation (number of tourist presences)")


ggplot(sim.sub.original,aes(x=(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line(size=1.5)+
xlab("time steps")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()+
xlim(950,1000)




pca<-prcomp(simw,scale=FALSE)
biplot(pca)


ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2)))+
geom_point()+
geom_path()+
theme_minimal()





library(ggrgl)


p<-ggplot(as.data.frame(pca$x))+
geom_path_3d(aes(x=(PC1),y=(PC2),z=(PC3)))

devoutrgl::rgldev(fov = 30, view_angle = -30)
p
invisible(dev.off())


try<-svd(simw)

CM.beh<-eigen(CommunityMatrix)
CMw.beh<-eigen(CM.weight)



eigen.df<-data.frame(real=Re(CMw.beh$values),imaginary=Im(CMw.beh$values),abs.dist=sqrt(Re(CMw.beh$values)^2+Im(CMw.beh$values)^2),label=paste0("eigenvalue ",1:27))

windows()
ggplot(eigen.df,aes(x=real,y=imaginary))+
geom_point()+
geom_text_repel(aes(label = label), show.legend = FALSE, size = 4,max.overlaps = Inf)+
theme_minimal()+
theme(text = element_text(size = 20))

library(ggnet)
library(igraph)
ComMatnet<-graph_from_adjacency_matrix(
  (CM.weight),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$weights<-abs( E(ComMatnet)$weight)*2
E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"

ggnet2(ComMatnet,label=TRUE,label.size=4,arrow.size=15,arrow.gap=0.02,edge.color="color",edge.size="weights")

##intervention
#correcting an error: in order to land fish, it has to be fished?

CM.intervention<-CM.weight

CM.intervention["Artisanal fishing activity","Fish landed for human consumption"]<-1
#living systems make more of themselves
CM.intervention["Food web dynamics","Food web dynamics"]<-1
CM.intervention["Biomass of reef fishes","Biomass of reef fishes"]<-1
CM.intervention["Area of P. oceanica habitat disturbed or lost","P.oceanica meadows (km2)"]<-(-1)

#loss of space increases rent
CM.intervention["Extension of natural habitats","Rent cost for housing"]<-.25

#carbon sequestration decreases net co2 emission
CM.intervention["Carbon sequestration","CO2 emission"]<-(-.25)

#######################
ComMatnet<-graph_from_adjacency_matrix(
  (CM.intervention),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$weights<-abs( E(ComMatnet)$weight)*2
E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"

ggnet2(ComMatnet,label=TRUE,label.size=4,arrow.size=15,arrow.gap=0.02,edge.color="color",edge.size="weights")

#######################
#A<-matrix(rep(10^100,nrow(CM.intervention)), ncol = 1)

A<-matrix(runif(nrow(CM.intervention),-1,1), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.intervention),nrow=1000)
colnames(simw)<-colnames(CM.intervention)
simw[1,]<-A

for (i in 1:999) {
simw[i+1,]<-CM.intervention%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()


pca<-prcomp(simw,scale=TRUE)

biplot(pca,col=c("red","black"),cex=c(.2,.75),xlim=c(-.2,.2),ylim=c(-.2,.2))


ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2)))+
geom_point()+
geom_path()+
theme_minimal()





ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2)))+
geom_point()+
geom_path()+
theme_minimal()+
ylim(.07996,.07997)+
xlim(.28602,.28604)


library(ggrgl)


p<-ggplot(as.data.frame(pca$x))+
geom_path_3d(aes(x=(PC1),y=(PC2),z=(PC3)))

devoutrgl::rgldev(fov = 30, view_angle = -30)
p
invisible(dev.off())


sim.sub<-subset(simw.melt,Var2=="P.oceanica meadows (km2)" | Var2=="Area of P. oceanica habitat disturbed or lost" | Var2=="Recreation (number of tourist presences)")


ggplot(sim.sub,aes(x=(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line(size=1.5)+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()+
xlim(950,1000)


sim.sub2<-subset(simw.melt,Var2=="P.oceanica status" | Var2=="Potential anchorages" | Var2=="Turbidity (mg/l)" | Var2=="Artisanal fishing activity")


ggplot(sim.sub2,aes(x=(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line(size=1.5)+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()+
xlim(950,1000)


CMw.complete<-eigen(CM.intervention)


eigen.df<-data.frame(real=Re(CMw.complete$values),imaginary=Im(CMw.complete$values),abs.dist=sqrt(Re(CMw.complete$values)^2+Im(CMw.complete$values)^2),label=paste0("eigenvalue ",1:27))

windows()
ggplot(eigen.df,aes(x=real,y=imaginary))+
geom_point()+
geom_text_repel(aes(label = label), show.legend = FALSE, size = 4,max.overlaps = Inf)+
theme_minimal()+
theme(text = element_text(size = 20))

eigen.df.1<-data.frame(real=Re(CMw.complete$vectors[,2]),imaginary=Im(CMw.complete$vectors[,2]),label=colnames(CM.intervention))


ggplot(eigen.df.1,aes(x=real,y=imaginary))+
geom_point()+
geom_text_repel(aes(label = label), show.legend = FALSE, size = 3,max.overlaps = Inf)+
theme_minimal()+
theme(text = element_text(size = 20))


##################################################################################################################
##################################################################################################################
### restrict vessel movement in MPA
CM.MPA<-CM.weight
CM.MPA["Recreation (number of tourist presences)","Vessel movements in the MPA"]<-0
CM.MPA["Resident population on the islands","Vessel movements in the MPA"]<-0

A<-matrix(runif(nrow(CM.MPA),-1,1), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.MPA),nrow=1000)
colnames(simw)<-colnames(CM.MPA)
simw[1,]<-A

for (i in 1:999) {
simw[i+1,]<-CM.MPA%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()


pca<-prcomp(simw,scale=FALSE)

biplot(pca,col=c("red","black"),cex=c(.2,.75),xlim=c(-.2,.2),ylim=c(-.2,.2))


ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2)))+
geom_point()+
#geom_text(hjust=0, vjust=0)+
geom_path()+
theme_minimal()


CMw.complete<-eigen(CM.MPA)


eigen.df<-data.frame(real=Re(CMw.complete$values),imaginary=Im(CMw.complete$values),abs.dist=sqrt(Re(CMw.complete$values)^2+Im(CMw.complete$values)^2),label=paste0("eigenvalue ",1:27))

windows()
ggplot(eigen.df,aes(x=real,y=imaginary))+
geom_point()+
geom_text_repel(aes(label = label), show.legend = FALSE, size = 4,max.overlaps = Inf)+
theme_minimal()+
theme(text = element_text(size = 20))


sim.sub<-subset(simw.melt,Var2=="P.oceanica meadows (km2)" | Var2=="Area of P. oceanica habitat disturbed or lost" | Var2=="Recreation (number of tourist presences)")


ggplot(sim.sub,aes(x=(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line(size=1.5)+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()+
theme(legend.position="bottom")+
xlim(950,1000)



##################################
## removal approach

CM.MPA<-CM.weight
CM.MPA<-CM.MPA[-(which(row.names(CM.MPA)=="Vessel movements in the MPA")),]
CM.MPA<-CM.MPA[,-(which(colnames(CM.MPA)=="Vessel movements in the MPA"))]

CM.MPA<-CM.MPA[-(which(row.names(CM.MPA)=="Licenses for diving access to  MPA")),]
CM.MPA<-CM.MPA[,-(which(colnames(CM.MPA)=="Licenses for diving access to  MPA"))]

CM.MPA<-CM.MPA[-(which(row.names(CM.MPA)=="Potential anchorages")),]
CM.MPA<-CM.MPA[,-(which(colnames(CM.MPA)=="Potential anchorages"))]


A<-matrix(runif(nrow(CM.MPA),-1000,1000), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.MPA),nrow=1000)
colnames(simw)<-colnames(CM.MPA)
simw[1,]<-A

for (i in 1:999) {
simw[i+1,]<-CM.MPA%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()


pca<-prcomp(simw,scale=TRUE)

biplot(pca,col=c("red","black"),cex=c(.2,.75),xlim=c(-.2,.2),ylim=c(-.2,.2))


ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2),label=1:nrow(pca$x)))+
geom_point()+
geom_text(hjust=0, vjust=0)+
geom_path()+
theme_minimal()+
xlim(-.5,.5)+
ylim(-.5,.5)


sim.sub<-subset(simw.melt,Var2=="P.oceanica meadows (km2)" | Var2=="Area of P. oceanica habitat disturbed or lost" | Var2=="Recreation (number of tourist presences)")


ggplot(sim.sub,aes(x=(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line(size=1.5)+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()+
theme(legend.position="bottom")+
xlim(900,950)


loops<-enumerate.loops(CM.MPA)

max(unlist(lapply(loops,length)))

F1<-Fn(1,CM.MPA)
F2<-Fn(2,CM.MPA)
F3<-Fn(3,CM.MPA)
F4<-Fn(4,CM.MPA)
F5<-Fn(5,CM.MPA)
F6<-Fn(6,CM.MPA)
F7<-Fn(7,CM.MPA)

###############################################################
###############################################################
###############################################################
#loop analysis

library(LoopAnalyst)
library(pracma)
library(RConics)


Fn<-function(n,CM) {(-1)*charpoly(CM)[[n+1]]}	

F1<-Fn(1,CM.weight)
F2<-Fn(2,CM.weight)
F3<-Fn(3,CM.weight)
F4<-Fn(4,CM.weight)
F5<-Fn(5,CM.weight)
F6<-Fn(6,CM.weight)
F7<-Fn(7,CM.weight)
F10<-Fn(10,CM.weight)


CR.adj<-make.adjoint(CM.weight)				#table of predictions
CR.wp<-weighted.predictions(CM.weight)			#weighted predictions
CEM<-make.cem(CM.weight)

loops<-enumerate.loops(CM.weight)



F1<-Fn(1,CM.intervention)
F2<-Fn(2,CM.intervention)
F3<-Fn(3,CM.intervention)
F4<-Fn(4,CM.intervention)
F5<-Fn(5,CM.intervention)
F6<-Fn(6,CM.intervention)
F7<-Fn(7,CM.intervention)
F10<-Fn(10,CM.intervention)


CR.adj<-make.adjoint(sign(CM.intervention))				#table of predictions
CR.wp<-weighted.predictions(CM.intervention)			#weighted predictions

CEM<-make.cem(CM.intervention)

loops<-enumerate.loops(CM.intervention)


####################################################################
##### shrink to activities - ecosystem services
label<-rownames(CM.weight)

activities<-label[c(24,3,21,18,22)]
ES<-label[c(9,6,10,12,27)]

retain<-c(which(label%in%activities),which(label%in%ES))

CM.SES<-CM.weight[retain,retain]


for (i in 1:nrow(CM.SES)) {
 for (j in 1:ncol(CM.SES)) {
#	if (i!=j) { #avoid self loops grabbed from original matrix if they were there given the subsetting above
	removed<-retain[-c(i,j)]
	
temp<-CM.weight[-removed,-removed] #now we removed all paths which go through an activity or ecosystem service

temp.graph<-graph_from_adjacency_matrix(
  (temp),
  mode = "directed",
  weighted = TRUE)

paths<-all_simple_paths(temp.graph,from=label[retain[i]],to=label[retain[j]], mode="out")

element<-0
if (length(paths)>0) {
all.element<-array(0,length(paths))
for (p in 1:length(paths)) {

edges = rep(names(paths[[p]]), each=2)[-1]

edges = edges[-length(edges)]
all.element[p]<-prod(E(temp.graph)$weight[get.edge.ids(temp.graph, edges)])

}
element<-mean(all.element)
}

CM.SES[i,j]<-element
#} #if loop
} #j
} #i

ComMatnet<-graph_from_adjacency_matrix(
  (CM.SES),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$weights<-abs( E(ComMatnet)$weight)*4
E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"

ggnet2(ComMatnet,label=TRUE,label.size=4,arrow.size=15,arrow.gap=0.02,edge.color="color",edge.size="weights")

CM.SESself<-CM.SES
CM.SESself[6,6]<-1
CM.SESself[8,8]<-1
CM.SESself[10,5]<-(-1) #carbon sequestration

A<-matrix(runif(nrow(CM.SESself),-10^10,10^10), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.SESself),nrow=100000)
colnames(simw)<-colnames(CM.SESself)
simw[1,]<-A

for (i in 1:99999) {
simw[i+1,]<-CM.SESself%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()


pca<-prcomp(simw,scale=TRUE)

biplot(pca,col=c("red","black"),cex=c(.01,.75),xlim=c(-.05,.06),ylim=c(-.05,.05))


ggplot(as.data.frame(pca$x),aes(x=(PC1),y=(PC2),label=1:nrow(pca$x)))+
geom_point()+
#geom_text(hjust=0, vjust=0)+
geom_path()+
theme_minimal()+
xlim(-.125,-.1)+
ylim(-1.2,-1.1)

loops<-enumerate.loops(CM.SESself)

max(unlist(lapply(loops,length)))

Fn<-function(n,CM) {(-1)*charpoly(CM)[[n+1]]}	

F0<-Fn(0,CM.SESself)
F1<-Fn(1,CM.SESself)
F2<-Fn(2,CM.SESself)
F3<-Fn(3,CM.SESself)
F4<-Fn(4,CM.SESself)
F5<-Fn(5,CM.SESself)
F6<-Fn(6,CM.SESself)
F7<-Fn(7,CM.SESself)
F8<-Fn(8,CM.SESself)
F9<-Fn(9,CM.SESself)
F10<-Fn(10,CM.SESself)#SES with 10 components has 10 max characteristic polynomials 


F0<-Fn(0,sign(CM.SESself))
F1<-Fn(1,sign(CM.SESself))
F2<-Fn(2,sign(CM.SESself))
F3<-Fn(3,sign(CM.SESself))
F4<-Fn(4,sign(CM.SESself))
F5<-Fn(5,sign(CM.SESself))
F6<-Fn(6,sign(CM.SESself))
F7<-Fn(7,sign(CM.SESself))
F8<-Fn(8,sign(CM.SESself))
F9<-Fn(9,sign(CM.SESself))
F10<-Fn(10,sign(CM.SESself))#SES with 10 components has 10 max characteristic polynomials 

hurwitzII<-matrix(c(
c(F1,F0,0,0,0,0,0,0,0,0),
c(F3,F2,F1,F0,0,0,0,0,0,0),
c(F5,F4,F3,F2,F1,F0,0,0,0,0),
c(F7,F6,F5,F4,F3,F2,F1,F0,0,0),
c(F9,F8,F7,F6,F5,F4,F3,F2,F1,F0),
c(0,F10,F9,F8,F7,F6,F5,F4,F3,F2),
c(0,0,0,F10,F9,F8,F7,F6,F5,F4),
c(0,0,0,0,0,F10,F9,F8,F7,F6),
c(0,0,0,0,0,0,0,F10,F9,F8),
c(0,0,0,0,0,0,0,0,0,F10)
),nrow=10,ncol=10)

det(hurwitzII)>0

CR.adj<-make.adjoint(-sign(CM.SESself))				#table of predictions
CR.wp<-weighted.predictions(CM.SESself)			#weighted predictions


############################################################
############################################################
### moving to Boolean approach advocated by Kristensen (MEE 2019)


###recalculating SES net


####################################################################
##### shrink to activities - ecosystem services
activities<-label[c(24,3,21,18,22)]
ES<-label[c(9,6,10,12,27)]

retain<-c(which(label%in%activities),which(label%in%ES))

CM.SES<-CM.weight[retain,retain]


for (i in 1:nrow(CM.SES)) {
 for (j in 1:ncol(CM.SES)) {
#	if (i!=j) { #avoid self loops grabbed from original matrix if they were there given the subsetting above
	removed<-retain[-c(i,j)]
	
temp<-CM.weight[-removed,-removed] #now we removed all paths which go through an activity or ecosystem service
#temp<-CM.weight
temp.graph<-graph_from_adjacency_matrix(
  (temp),
  mode = "directed",
  weighted = TRUE)

paths<-all_simple_paths(temp.graph,from=label[retain[i]],to=label[retain[j]], mode="out")

element<-0
if (length(paths)>0) {
all.element<-array(0,length(paths))
for (p in 1:length(paths)) {

edges = rep(names(paths[[p]]), each=2)[-1]

edges = edges[-length(edges)]
all.element[p]<-prod(E(temp.graph)$weight[get.edge.ids(temp.graph, edges)])

}
element<-mean(all.element)
}

CM.SES[i,j]<-element
#} #if loop
} #j
} #i

CM.SESself<-CM.SES
CM.SESself[6,6]<-1
CM.SESself[8,8]<-1
CM.SESself[10,5]<-(-1) #carbon sequestration
CM.SESself[3,3]<-(-1)
CM.SESself[2,2]<-(-1)
CM.SESself[1,1]<-(-1)
CM.SESself[1,5]<-(1)

cSESself<-sign(t(CM.SESself))  #careful here row col swapped
diag(cSESself)[diag(cSESself)==0]<-(-1) ### temp fix to deal with singular matrix

det(cSESself)

c.adj<-RConics::adjoint(-cSESself)
absSESself<-abs(cSESself)

library(LoopAnalyst)
T<-make.T(cSESself)

# T<-cSESself

# for (i in 1:nrow(T)) {
	# for (j in 1:ncol(T)) {
		# T[i,j]<-permanent(absSESself[-i,-j]) 
	# }
# }

W<-abs(c.adj)/T

row.names(c.adj)<-row.names(CM.SESself)
colnames(c.adj)<-row.names(CM.SESself)

adj.melt<-melt(c.adj[8,])
adj.melt$Component=rownames(adj.melt)

WPosidonia<-melt(W[8,])
WPosidonia$Component=rownames(WPosidonia)

adj.melt$weight<-WPosidonia$value
ggplot(adj.melt,aes(x=Component,y=value,fill=weight))+
geom_bar(stat="identity")+
ylab("Posidonia response prediction")+
xlab("SES Component")+
theme_minimal()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))

#######################################################################
#### Boolean approach
## exactly the way I approach it
## let's try R boolean capacities

library(Matrix)

mat<-Matrix(sample(c(TRUE,FALSE),25,TRUE),5,5)

A<-Matrix(sample(c(TRUE,FALSE),5,TRUE),5,1)
#nope let's try Matrix instead

matn<-as(mat,"nMatrix")
An<-as(A,"nMatrix")
matn%&%An #boolean matrix product


###################################################
### ok looking good let's try on our SES
SES<-as(matrix(0,10,10),"nMatrix")
SES<-Matrix(CM.SESself>0)
SES[CM.SESself>0]<-TRUE
SES[CM.SESself==0]<-TRUE

SES@x[which(CM.SESself[-which(CM.SESself<0)]==0)]<-FALSE



##SES is boolean now

A<-Matrix(matrix(sample(c(TRUE,FALSE),10,TRUE),10,1))


SES%&%A

simw<-Matrix(TRUE,ncol=nrow(SES),nrow=1000)
colnames(simw)<-colnames(SES)

simw[1,]<-A

for (i in 1:999) {
As<-Matrix(matrix(simw[i,], ncol = 1))
simw[i+1,]<-(SES)%&%As

}


library(reshape2)
simw.melt<-melt(as.matrix(simw))

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=Var2,colour=value))+
geom_point()+
# xlab("time steps (log10)")+
# ylab("production of each component (log10 transformed)")+
# hlabs(colour="SES components")+
theme_minimal()

##################################################
#### repeat with original SES

CM.weight 
SES<-as(matrix(0,27,27),"nMatrix")
SES<-Matrix(CM.weight>0)
SES[CM.weight>0]<-TRUE
SES[CM.weight==0]<-TRUE

SES@x[which(CM.weight[-which(CM.weight<0)]==0)]<-FALSE

A<-Matrix(matrix(sample(c(TRUE,FALSE),27,TRUE),27,1))


simw<-Matrix(TRUE,ncol=nrow(SES),nrow=1000)
colnames(simw)<-colnames(SES)

simw[1,]<-A

for (i in 1:999) {
As<-Matrix(matrix(simw[i,], ncol = 1))
simw[i+1,]<-t(SES)%&%As

}


library(reshape2)
simw.melt<-melt(as.matrix(simw))

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=Var2,colour=value))+
geom_point()+
 xlab("time steps (log10)")+
 ylab("production of each component (boolean)")+
# hlabs(colour="SES components")+
theme_minimal()




#### we can actually complete the Boolean work outline in Kristensen in R
### the espresso algorithm is available in LogicOpt/versions/1.0.0/topics/logicopt
### defunct!!

#using C code instead and call it from R
## ok easiest is to run python code via reticulate

###########################################
### env built on workstation (older python version there functools needs >=3.9
library(reticulate)
myenvs=conda_list()

envname=myenvs$name[2] #py39
use_condaenv(envname, required = TRUE)
#Sys.setenv(RETICULATE_PYTHON = "C:/ProgramData/Anaconda3/envs/py39")
#Sys.setenv(RETICULATE_PYTHON = "C:/Users/davlu/AppData/Local/miniconda3")
#use_condaenv('py39')

#use espresso in pyeda
py_config()
conda_install(envname=envname,packages="pyeda",pip=TRUE)
conda_install(envname=envname,packages="networkx",pip=TRUE)
conda_install(envname=envname,packages="matplotlib",pip=TRUE)
conda_install(envname=envname,packages="ipython",pip=TRUE)
conda_install(envname=envname,packages="bidict",pip=TRUE)
source_python("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/Kristensen/functools.py")
### the excellent stuff from Kristensen
source_python("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/Kristensen/qualmod.py")
source_python("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/Kristensen/findpcu.py")


####DONE!
## next the tutorial is the best template

CM<-sign(CM.weight)
np<-import("numpy")
nplinalg<-import("numpy.linalg")

nplinalg$inv(CM)

# good [well bad but consistent] singular matrix

nx<-import("networkx")
plt<-import(" matplotlib.pyplot")

library(igraph)
ComMatnet<-graph_from_adjacency_matrix(
  (CM),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"

neg_edges<- attr((E(ComMatnet)[E(ComMatnet)$weight<0]),"vnames")
pos_edges<- attr((E(ComMatnet)[E(ComMatnet)$weight>0]),"vnames")

neg_list<-list()
for (i in 1:length(neg_edges)) {
temp<-unlist(strsplit(neg_edges[i],"[|]"))
neg_list[[i]]<-temp[2]
names(neg_list[[i]])<-temp[1]
}


pos_list<-list()
for (i in 1:length(pos_edges)) {
temp<-unlist(strsplit(pos_edges[i],"[|]"))
pos_list[[i]]<-temp[2]
names(pos_list[[i]])<-temp[1]
}

web = initialise_foodweb(pos_list, neg_list)

#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
# actually the downstream boolean algebra can be done natively in R

####################################################################
#####################################################################
library(BoolNet)

#use igf as template to reconstruct the boolean network in the right format
#https://cran.r-project.org/web/packages/BoolNet/vignettes/BoolNet_package_vignette.pdf

tuscany_boolean<-loadNetwork("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/tuscany_boolean_ses.csv")
#simulateSymbolicModel
tuscany_boolean_alt<-loadNetwork("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/tuscany_boolean_ses_alt.csv")
statesa<-getAttractors(tuscany_boolean_alt)

states<-getAttractors(tuscany_boolean)

# print(states, activeOnly=TRUE)
# Attractor 1 is a simple attractor consisting of 1 state(s) and has a basin of 3656704 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, BATHING, DISCHARGE, EMPLOYMENT, FISH_BIOMASS, FISH_LANDED, FISHING, FOOD_WEB, MPA_ACCESS, MPA_VESSELS, NATURAL_HABITAT, POSIDONIA, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, RESIDENTS, TURBIDITY

# Attractor 2 is a simple attractor consisting of 1 state(s) and has a basin of 5376 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, BATHING, DISCHARGE, FISH_LANDED, FISHING, MPA_ACCESS, MPA_VESSELS, NATURAL_HABITAT, POSIDONIA, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, RESIDENTS, TURBIDITY

# Attractor 3 is a simple attractor consisting of 1 state(s) and has a basin of 16128 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, BATHING, DISCHARGE, EMPLOYMENT, FISH_BIOMASS, FISH_LANDED, FISHING, MPA_ACCESS, MPA_VESSELS, NATURAL_HABITAT, POSIDONIA, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, RESIDENTS, TURBIDITY

# Attractor 4 is a simple attractor consisting of 1 state(s) and has a basin of 508928 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, DISCHARGE, EMPLOYMENT, FISH_BIOMASS, FISH_LANDED, FISHING, FOOD_WEB, MPA_ACCESS, MPA_VESSELS, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, RESIDENTS, TURBIDITY

# Attractor 5 is a simple attractor consisting of 1 state(s) and has a basin of 1792 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, DISCHARGE, FISHING, MPA_VESSELS, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, TURBIDITY

# Attractor 6 is a simple attractor consisting of 1 state(s) and has a basin of 5376 state(s).
# Active genes in the attractor state(s):
# State 1: ALIEN, ANCHORAGES, DISCHARGE, EMPLOYMENT, FISH_BIOMASS, FISH_LANDED, FISHING, MPA_ACCESS, MPA_VESSELS, POSIDONIA_LOST, RECREATION, REMOVED_ANIMALS, REMOVED_BIOMASS, RESIDENTS, TURBIDITY


library(igraph)
 
 state.map<-plotStateGraph(states, layout = layout.fruchterman.reingold,plotIt = FALSE)
vertex_attr(state.map,index=1)
table(V(state.map)$color)


V(state.map)$degree<-degree(state.map)
V(state.map)$attractor<-"no"
V(state.map)$attractor[V(state.map)$frame.color=="#000000FF"]<-"yes"

library(stringr)

V(state.map)$color2<-str_sub(V(state.map)$color, end=-3)

E(state.map)$width<-1
E(state.map)$lty<-1

E(state.map)$color2<-str_sub(E(state.map)$color, end=-3)

core.state<-subgraph(state.map,which(V(state.map)$degree>1))

plot(core.state,layout=layout.fruchterman.reingold, vertex.label=NA,vertex.size=0.2,edge.arrow.size=0.05,edge.width=.7)


write_graph(
  state.map,
  file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/tuscany_state.graphml",
  format = "graphml"
)


trans.tab<-getTransitionTable(states)
#plotStateGraph(states)
print(getBasinOfAttraction(states, 1))
 

A1<-plotAttractors(statesa, subset=1)
state.def<-data.frame(components=rownames(A1$"1"),state=1,value=as.numeric(A1$"1"))

for (i in 2:6) {

A<-plotAttractors(statesa, subset=i)
temp<-data.frame(components=rownames(A$"1"),state=i,value=as.numeric(A$"1"))
state.def<-rbind(state.def,temp)

}

state.def$state<-factor(state.def$state)


ggplot(state.def, aes(x = state, y = components)) + 
  geom_tile(aes(fill=value),colour="black",show.legend =FALSE) + 
  scale_fill_gradient(low="red2", high="green4") +
  labs(x="attractor", y="SES components", title="") +
  theme_bw() 
  


path<-getPathToAttractor(tuscany_boolean,sample(c(0,1),22,replace=TRUE))

#plotSequence(sequence=path,drawLegend=FALSE)
plotSequence(network=tuscany_boolean,startState=sample(c(0,1),22,replace=TRUE),drawLegend=FALSE)

 table.seq<-plotSequence(network=tuscany_boolean,startState=sample(c(0,1),22,replace=TRUE),drawLegend=FALSE)

library(ggplot2)
library(reshape2)
tab.melt<-melt(table.seq)
tab.melt$Var2<-tab.melt$Var2-1
tab.melt$Var2<-factor(tab.melt$Var2)

ggplot(tab.melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value),colour="black",show.legend =FALSE) + 
  scale_fill_gradient(low="red2", high="green4") +
  labs(x="time steps", y="SES components", title="") +
  theme_bw() 
  
   sim <- markovSimulation(tuscany_boolean)

### what happens if stop mpa access and plant posidonia and stop anchorage
ko<-fixGenes(tuscany_boolean_alt,c("ANCHORAGES"),0)


   state.ko<-getAttractors(ko)
   
   #2 states

A1<-plotAttractors(state.ko, subset=1)
kstate.def<-data.frame(components=rownames(A1$"1"),state=1,value=as.numeric(A1$"1"))

for (i in 2:2) {

A<-plotAttractors(state.ko, subset=i)
temp<-data.frame(components=rownames(A$"1"),state=i,value=as.numeric(A$"1"))
kstate.def<-rbind(kstate.def,temp)

}

kstate.def$state<-factor(kstate.def$state)


ggplot(kstate.def, aes(x = state, y = components)) + 
  geom_tile(aes(fill=value),colour="black",show.legend =FALSE) + 
  scale_fill_gradient(low="red2", high="green4") +
  labs(x="attractor", y="SES components", title="") +
  theme_bw() 

   
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
  
  
#############################################################
#############################################################
#############################################################
#### optimisation using departure from normality as a target for minimisation
#### some of the non-normal matrix maths from Asllani

SESnorm<-norm(CM.SES,type="F") #Frobenius norm
sum(abs(eigen(CM.SES)$value)^2)
SES.spectrum<-eigen(CM.SES)$value

SES.hermitian<-(CM.SES+t(CM.SES))/2
omega<-max(Re(eigen(SES.hermitian)$value))
alpha<-max(Re(eigen(CM.SES)$value))
omega-alpha



weightnorm<-norm(CM.weight,type="F") #Frobenius norm
dF_std<-sqrt(weightnorm^2-sum(abs(eigen(CM.weight)$value)^2))/weightnorm

weight.spectrum<-eigen(CM.weight)$value

weight.hermitian<-(CM.weight+t(CM.weight))/2
omega<-max(Re(eigen(weight.hermitian)$value))
alpha<-max(Re(eigen(CM.weight)$value))
omega-alpha
assymetry<-abs(sum(CM.weight[lower.tri(CM.weight)])-sum(CM.weight[upper.tri(CM.weight)]))/(sum(CM.weight[lower.tri(CM.weight)])+sum(CM.weight[upper.tri(CM.weight)]))

#CM.weight[4,21]

dF<-function(x) {
Xnorm<-norm(x,type="F")
dF<-sqrt(Xnorm^2-sum(abs(eigen(x)$value)^2))/Xnorm

return(dF)
}


D.anchorage<-c(-1000,-100,-10,-5,-2,seq(-1,1,.1),2,5,10,100,1000)
pseudostability<-array(NA,length(D))

temp<-CM.weight
for (i in 1:length(D.anchorage)) {
temp[4,21]<-D.anchorage[i]
pseudostability[i]<-dF(temp)
}

plot(pseudostability~(D.anchorage),xlab="strength anchorage->tourists",ylab="stability departure",ylim=c(.75,.8),xlim=c(-2,2))

############################
############################
############################3
############################3
#### optimisation sign respected
dF.optim<-function(x) {
temp<-CM.weight
temp[temp!=0]<-x
temp<-sign(CM.weight)*abs(temp) # lacks efficiency
dF<-dF(temp)
return(dF)

}

try<-optim(CM.weight[CM.weight!=0],dF.optim,lower=0.000000001,upper=1,method="L-BFGS-B")

temp<-CM.weight
temp[temp!=0]<-try$par
CM.optim<-sign(CM.weight)*abs(temp)

CM.optim


library(ggnet)
library(igraph)
ComMatnet<-graph_from_adjacency_matrix(
  (CM.optim),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$weights<-abs( E(ComMatnet)$weight)*2
E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"
E(ComMatnet)$real.weight<-E(ComMatnet)$weight
E(ComMatnet)$weight<-abs( E(ComMatnet)$weight)*5

plot(ComMatnet,layout=layout_with_lgl)

ggnet2(ComMatnet,label=TRUE,label.size=4,arrow.size=15,arrow.gap=0.02,edge.color="color",edge.size="weights")




A<-matrix(runif(nrow(CM.optim),-1,1), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.optim),nrow=10000)
colnames(simw)<-colnames(CM.optim)
simw[1,]<-A

for (i in 1:9999) {
simw[i+1,]<-CM.optim%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()




###########################
###########################
###########################
############################
#### optimisation sign not respected
dF.optim<-function(x) {
temp<-CM.weight
temp[temp!=0]<-x
#temp<-sign(CM.weight)*abs(temp) # lacks efficiency
dF<-dF(temp)
return(dF)
}

try<-optim(CM.weight[CM.weight!=0],dF.optim,lower=-10,upper=10,method="L-BFGS-B")

temp<-CM.weight
temp[temp!=0]<-try$par
CM.optim<-(temp)

CM.optim


library(ggnet)
library(igraph)
ComMatnet<-graph_from_adjacency_matrix(
  (CM.optim),
  mode = "directed",
  weighted = TRUE)

E(ComMatnet)$weights<-abs( E(ComMatnet)$weight)*2
E(ComMatnet)$sign<-sign( E(ComMatnet)$weight)
E(ComMatnet)$color<-"blue"
E(ComMatnet)$color[E(ComMatnet)$sign<0]<-"red"

ggnet2(ComMatnet,label=TRUE,label.size=4,arrow.size=15,arrow.gap=0.02,edge.color="color",edge.size="weights")




A<-matrix(runif(nrow(CM.optim),-1,1), ncol = 1)

simw<-matrix(0,ncol=nrow(CM.optim),nrow=10000)
colnames(simw)<-colnames(CM.optim)
simw[1,]<-A

for (i in 1:9999) {
simw[i+1,]<-CM.optim%*%matrix((simw[i,]), ncol = 1)
}

simw.sign<-sign(simw)
library(reshape2)
simw.melt<-melt(simw)

library(ggplot2)

ggplot(simw.melt,aes(x=log10(Var1),y=sign(value)*log10(abs(value)),colour=Var2))+
geom_line()+
xlab("time steps (log10)")+
ylab("production of each component (log10 transformed)")+
labs(colour="SES components")+
theme_minimal()
