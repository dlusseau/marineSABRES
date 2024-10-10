library(igraph)
library(data.table)


tuscany<-fread("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD.csv")
tuscany_element<-fread("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD_element.csv")
tuscany_element<-tuscany_element[,1:2]

unique(tuscany$strength)

tuscany$weight<-1
tuscany$weight[tuscany$strength=="Medium Positive"]<-.5
tuscany$weight[tuscany$strength=="Medium Negative"]<-(-.5)
tuscany$weight[tuscany$strength=="Strong Negative"]<-(-1)
tuscany$sign<-sign(tuscany$weight)


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

tuscany.SES<-make_matrix(from=tuscany$From,to=tuscany$To,weight=tuscany$weight)

tuscany.laplacian<-SES.laplacian(tuscany.SES,from="rows")


###################################################
#### graph

plot.SES<-function(SES.mat,folder,filename,title,w=40,h=40,layout=layout_with_fr,label.cex=1.5,vsize=25,eweight=10) {
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


minC <- rep(-Inf, vcount(SES.net))
maxC <- rep(Inf, vcount(SES.net))
minC[1] <- maxC[1] <- 0
co <- layout(SES.net, minx=minC, maxx=maxC,
                                  miny=minC, maxy=maxC)

png(paste0(folder,filename,".png"),width=w,height=h,units="cm",res=200)
plot(SES.net,layout=co,rescale=FALSE,xlim=range(co[,1]), ylim=range(co[,2]),ylab = title)
dev.off()

return(SES.net)

}

tuscany.graph<-plot.SES(tuscany.SES,
				folder="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/",
				filename="tuscany_cld_graph",
				title="Tuscany SES (CLD)",
				label.cex=1.1,
				vsize=30
				)
				
				
#########################################################
##### boolean

### create file in right format


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


folder<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename<-"Tuscany_boolean"

boolean_file_creation(tuscany.SES,folder,filename)


library(BoolNet)

Tuscany_boolean<-loadNetwork(paste0(folder,filename,".csv"))

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
core.state<-subgraph(state.map,which(V(state.map)$degree>1))
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


folder<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename<-"Tuscany_boolean2"

tuscany_bool_ana<-boolean_analyses(Tuscany_boolean,folder,filename)

attractor1<-data.frame(name=colnames(tuscany.SES),boolean=as.numeric(tuscany_bool_ana[[4]][[1]]))

#################################
## is boolean all that we have? then we do the matrix projection on the boolean matrix. let's make sure to 
## to use boolean algebra (see first code file)



###############################################################################
#### quantitative projection

tuscany.SES


iter=1000

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


#png(paste0(folder,filename,".png"),width=60,height=35,units="cm",res=200)
plot1<-ggarrange(text1,p1a,p1b,nrow=1,ncol=3,widths=c(.1,1,1))
ggsave(paste0(folder,filename,".png"),plot=plot1,device="png",width=60,height=35,units="cm",dpi=200,bg="white")
#dev.off()

return(SES.sim)

}

folder<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"

tuscany.sim<-SES.simulate(tuscany.SES,1000,folder,filename="tuscany_weighted_projection",title="Tuscany SES (CLD)")

sim.melt<-melt(tuscany.sim)

fillcol<-glasbey.colors(nrow(tuscany.sim))
names(fillcol)<-unique(sim.melt$Var1)

ggplot(subset(sim.melt,Var2>50 & Var2<70),aes(x=(Var2),y=(value),colour=(Var1)))+
 geom_path(size=2)+
 labs(colour="Element", y="Progress",x="(time)")+
 scale_colour_manual(values=fillcol)+
 theme_minimal(base_size=18)+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
 theme(legend.position="none")
 
 
####################################################################################################
##############################
### current SES participation ratio
tuscany.jacobian<-t(tuscany.SES)

lefteigen<-function(mat) {
left<-eigen(t(mat))$vector
return(left)
}

righteigen<-function(mat) {
right<-eigen((mat))$vector
return(right)
}

jacobian<-function(SES) {
SES.j<-t(SES)
return(SES.j)
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

folder="C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"

filename="Tuscany_participation_ratio"
title="Tuscany (CLD)"
tuscany.PR<-participation_ratio(tuscany.SES,folder=folder,filename=filename,title=title) #Politics has the broader reach


###############################################################################################

targets<-tuscany_element$Label[which(tuscany_element$Description=="Good and Benefit" |tuscany_element$Description=="Ecosystem Service")]
## for some reason elements are not the same as the graphs

targets<-targets[targets%in%row.names(tuscany.SES)]

#### what weight is important?

simulate.mat<-function(mat) {
mat.sim<-1*(mat!=0)*sign(mat)*matrix(runif(prod(dim(mat)),0,1),nrow(mat),ncol(mat))
return(mat.sim)
}


time.simulate<-function(mat,iter,starting.value) {
#mat is the original orientation of the matrix
SES.sim<-matrix(NA,nrow(mat),iter)
SES.sim[,1]<-starting.value

for ( i in 2:iter) {

SES.sim[,i]<-t(mat)%*%matrix((SES.sim[,i-1]), ncol = 1)

}

return(SES.sim)
}

greed<-1000000
iter=500
mat<-tuscany.SES
starting.value<-runif(nrow(mat),0,1)

SES.obs<-time.simulate(mat,iter,starting.value)
state.obs<-SES.obs[,iter]
names(state.obs)<-colnames(tuscany.SES)

state.sim<-array(NA,dim=c(length(state.obs),greed))
rownames(state.sim)<-colnames(tuscany.SES)

mat.m<-melt(mat)
mat.sim.df<-array(NA,dim=c(nrow(mat.m),greed+2))
mat.sim.df[,1]<-mat.m$Var1
mat.sim.df[,2]<-mat.m$Var2 #can't deal with character and numbers in array, that's ok we keep it to factor levels, array is faster than df

tic<-Sys.time()

for (i in 1:greed) {

mat.sim<-simulate.mat(mat)

SES<-time.simulate(mat.sim,iter,starting.value)

mat.sim.df[,i+2]<-melt(mat.sim)$value
#state.sim[,i]<-SES[,iter]-matrix(state.obs,ncol=1)
state.sim[,i]<-apply(sign(SES[,(iter-100):iter]),1,prod) # here we capture behaviour over last 101 steps, importantly an odd number
}
toc<-Sys.time()-tic

#save!

save(state.sim,mat.sim.df,mat.m,mat,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_greedy_simulation.Rdata")
load("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_greedy_simulation.Rdata")

tol=0.000001
state.sim.bin<-state.sim
state.sim.bin[abs(state.sim.bin)<tol]<-0
state.sim.bin<-sign(state.sim.bin)


desirable.outcomes<-which(colSums(state.sim.bin[match(targets,row.names(state.sim.bin)),])==length(targets))
mat.works<-mat.sim.df[,c(1,2,(desirable.outcomes+2))]

##################################################################
sim.ana<-mat.sim.df[which(mat.m$value!=0),]
sim.ana.df<-t(sim.ana[,3:(greed+2)])
colnames(sim.ana.df)<-apply(sim.ana[,1:2],1,function(x) paste0(row.names(mat)[x[1]]," to ",row.names(mat)[x[2]]))
sim.ana.df<-as.data.frame(sim.ana.df)

sim.ana.df$outcomes<-0
sim.ana.df$outcomes[desirable.outcomes]<-1


library(randomForest)
library(randomForestExplainer)
library(reprtree)
library(cluster)

sim.ana.df$outcomes<-factor(sim.ana.df$outcomes)
colnames(sim.ana.df)<-gsub("\\s*\\([^\\)]+\\)","",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub(" ","_",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub("\\.","",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub("\\-","",colnames(sim.ana.df))


# forest.tuscany <- ranger(outcomes ~ ., data = sim.ana.df, ntree=500,importance="impurity")
# save(forest.tuscany,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_forest_ranger.Rdata")

forest.tuscany <- randomForest(outcomes ~ ., data = sim.ana.df, ntree=500,localImp=TRUE,proximity=FALSE)
save(forest.tuscany,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_forest.Rdata")

###########################

gc()

importance_frame.tuscany <- measure_importance(forest.tuscany)
save(importance_frame.tuscany,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_importance.Rdata")

pa<-plot_multi_way_importance(importance_frame.tuscany, x_measure = "gini_decrease", y_measure = "accuracy_decrease", size_measure = "p_value", main="Tuscany (CLD)")+labs(x="mean decrease of GINI coefficient", y="mean decrease in accuracy")

folder="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename="Tuscany_randomforest_variable_importance"
ggsave(paste0(folder,filename,".png"),plot=pa,device="png",width=40,height=40,units="cm",dpi=200,bg="white")


########################################
########################################
########################################
########################################

###
##### introduce measures

library(igraph)
library(data.table)


tuscany<-fread("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD.csv")
tuscany_element<-fread("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD_element.csv")
tuscany_element<-tuscany_element[,1:2]

unique(tuscany$strength)

tuscany$weight<-1
tuscany$weight[tuscany$strength=="Medium Positive"]<-.5
tuscany$weight[tuscany$strength=="Medium Negative"]<-(-.5)
tuscany$weight[tuscany$strength=="Strong Negative"]<-(-1)
tuscany$sign<-sign(tuscany$weight)



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


tuscany.SES<-make_matrix(from=tuscany$From,to=tuscany$To,weight=tuscany$weight)


targets<-tuscany_element$Label[which(tuscany_element$Description=="Good and Benefit" |tuscany_element$Description=="Ecosystem Service")]
## for some reason elements are not the same as the graphs

targets<-targets[targets%in%row.names(tuscany.SES)]

affected<-row.names(tuscany.SES)[2]  #choose by end

indicators<-tuscany_element$Label[which(tuscany_element$Description=="Good and Benefit")]
indicators<-indicators[indicators%in%row.names(tuscany.SES)]

######################
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

tuscany.alt<-simulate.measure(tuscany.SES,"tourism throttling",affected,indicators,lower=-1,upper=0)

############################
### load functions

simulate.mat<-function(mat) {
mat.sim<-1*(mat!=0)*sign(mat)*matrix(runif(prod(dim(mat)),0,1),nrow(mat),ncol(mat))
return(mat.sim)
}


time.simulate<-function(mat,iter,starting.value) {
#mat is the original orientation of the matrix
SES.sim<-matrix(NA,nrow(mat),iter)
SES.sim[,1]<-starting.value

for ( i in 2:iter) {

SES.sim[,i]<-t(mat)%*%matrix((SES.sim[,i-1]), ncol = 1)

}

return(SES.sim)
}


###################################
### now we do this may time and ask when we get the desired outcomes


##### this is going to be a bit different, we are going to try to be memory efficient and only retain the added row and col
library(reshape2)
greed<-1000000
iter=500
mat<-tuscany.SES
starting.value<-runif(nrow(mat),0,1)


targets<-tuscany_element$Label[which(tuscany_element$Description=="Good and Benefit" |tuscany_element$Description=="Ecosystem Service")]
## for some reason elements are not the same as the graphs

targets<-targets[targets%in%row.names(tuscany.SES)]

affected<-row.names(tuscany.SES)[2]  #choose by end

indicators<-tuscany_element$Label[which(tuscany_element$Description=="Good and Benefit")]
indicators<-indicators[indicators%in%row.names(tuscany.SES)]
measure="tourism throttling"
newnames<-c(row.names(mat),measure)

focus.mat.df<-data.frame(from=factor(c(rep(measure,length(newnames)),newnames)),to=factor(c(newnames,rep(measure,length(newnames)))))


state.sim<-array(NA,dim=c(nrow(mat)+1,greed))
rownames(state.sim)<-c(row.names(mat),measure)

mat.sim.df<-array(NA,dim=c(nrow(focus.mat.df),greed+2))
mat.sim.df[,1]<-focus.mat.df$from
mat.sim.df[,2]<-focus.mat.df$to #can't deal with character and numbers in array, that's ok we keep it to factor levels, array is faster than df


starting.value<-c(starting.value,0) #we start without measures
tic<-Sys.time()

for (i in 1:greed) {

mat.sim<-simulate.measure(mat,measure,affected,indicators,lower=-1,upper=0)


SES<-time.simulate(mat.sim,iter,starting.value)

mat.sim.df[,i+2]<-c(mat.sim[(nrow(mat)+1):nrow(mat.sim),],mat.sim[,(ncol(mat)+1):ncol(mat.sim)])  #like that we account for multiple measures

state.sim[,i]<-apply(sign(SES[,(iter-100):iter]),1,prod) # here we capture behaviour over last 101 steps, importantly an odd number

}
toc<-Sys.time()-tic
toc

#save!
save(state.sim,mat.sim.df,mat.m,mat,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_measure_application_greedy_simulation.Rdata")

load("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_measure_application_greedy_simulation.Rdata")

tol=0.000001
state.sim.bin<-state.sim
state.sim.bin[abs(state.sim.bin)<tol]<-0
state.sim.bin<-sign(state.sim.bin)


desirable.outcomes<-which(colSums(state.sim.bin[match(targets,row.names(state.sim.bin)),])==length(targets))
mat.works<-mat.sim.df[,c(1,2,(desirable.outcomes+2))]

##################################################################

sim.ana.df<-t(mat.sim.df[,3:(greed+2)])
colnames(sim.ana.df)<-apply(mat.sim.df[,1:2],1,function(x) paste0(levels(focus.mat.df$from)[x[1]]," to ",levels(focus.mat.df$to)[x[2]]))
sim.ana.df<-as.data.frame(sim.ana.df)

sim.ana.df$outcomes<-0
sim.ana.df$outcomes[desirable.outcomes]<-1

sim.ana.df$outcomes<-factor(sim.ana.df$outcomes)
select<-colSums(sim.ana.df[,1:42])
which(select!=0)

sim.ana.df<-sim.ana.df[,c(which(select!=0),43)]

library(GGally)
success<-subset(sim.ana.df,outcomes==1)
colnames(success)[1]<-"tourismthrottlingTOrecreation"
colnames(success)[4]<-"waterqualityTOtourismthrottling"


pc1<-ggplot(success,
aes(x=tourismthrottlingTOrecreation,
y=waterqualityTOtourismthrottling))+
geom_density_2d_filled(fill='#1f78b4',aes(alpha=(..level..)))+
#stat_density_2d(aes(alpha=(..level..)),bins=4,fill='#b2df8a',geom="polygon",colour="white")+
labs(title=paste(""),alpha="density")+
#ylim(0,1)+
#xlim(0,1)+
geom_abline(intercept = 0, slope = 1,color = "black", linewidth=1.5)+
theme_minimal(base_size=16)+
theme(legend.position="none")

#
ggpairs(subset(sim.ana.df,outcomes==1),                 
        columns = 1:5,
		lower="blank"
        )   


ggpairs(sim.ana.df,                 
        columns = 1:5,      
        aes(color = outcomes, 
            alpha = 0.5))   
			
			
library(randomForest)
library(randomForestExplainer)
library(reprtree)
library(cluster)

colnames(sim.ana.df)<-gsub("\\s*\\([^\\)]+\\)","",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub(" ","_",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub("\\.","",colnames(sim.ana.df))
colnames(sim.ana.df)<-gsub("\\-","",colnames(sim.ana.df))


# forest.tuscany <- ranger(outcomes ~ ., data = sim.ana.df, ntree=500,importance="impurity")
# save(forest.tuscany,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_forest_ranger.Rdata")

forest.tuscany_measure <- randomForest(outcomes ~ ., data = sim.ana.df, ntree=500,localImp=TRUE,proximity=FALSE)
save(forest.tuscany_measure,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_forest.Rdata")


importance_frame.tuscany_measure <- measure_importance(forest.tuscany_measure)
save(importance_frame.tuscany_measure,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_measure_importance.Rdata")

pa<-plot_multi_way_importance(importance_frame.tuscany_measure, x_measure = "gini_decrease", y_measure = "accuracy_decrease", size_measure = "p_value", main="Tuscany (CLD)")+labs(x="mean decrease of GINI coefficient", y="mean decrease in accuracy")

folder="C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename="Tuscany_measure_randomforest_variable_importance"
ggsave(paste0(folder,filename,".png"),plot=pa,device="png",width=40,height=40,units="cm",dpi=200,bg="white")

load("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_measure_importance.Rdata")






