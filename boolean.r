####### first, simple simulation
library(Matrix)
library(reshape2)
library(ggplot2)

SES<-as(matrix(0,10,10),"nMatrix")
SES<-Matrix(CM.SESself>0)
SES[CM.SESself>0]<-TRUE
SES[CM.SESself==0]<-TRUE

SES@x[which(CM.SESself[-which(CM.SESself<0)]==0)]<-FALSE

##SES is boolean now

A<-Matrix(matrix(sample(c(TRUE,FALSE),10,TRUE),10,1))

simw<-Matrix(TRUE,ncol=nrow(SES),nrow=1000)
colnames(simw)<-colnames(SES)

simw[1,]<-A

for (i in 1:999) {
As<-Matrix(matrix(simw[i,], ncol = 1))
simw[i+1,]<-(SES)%&%As
}

simw.melt<-melt(as.matrix(simw))

ggplot(simw.melt,aes(x=log10(Var1),y=Var2,colour=value))+
geom_point()+
theme_minimal()

#########################################################################
#########################################################################
#### Boolean network appraisal
library(BoolNet)
library(igraph)
library(ggplot2)
library(reshape2)

### check out this very good vignette
#https://cran.r-project.org/web/packages/BoolNet/vignettes/BoolNet_package_vignette.pdf

tuscany_boolean<-loadNetwork("/tuscany_boolean_ses.csv") #comma separated txt file 

states<-getAttractors(tuscany_boolean)

 print(states, activeOnly=TRUE)
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
  file="/tuscany_state.graphml",
  format = "graphml"
)


trans.tab<-getTransitionTable(states)
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

plotSequence(network=tuscany_boolean,startState=sample(c(0,1),22,replace=TRUE),drawLegend=FALSE)

 table.seq<-plotSequence(network=tuscany_boolean,startState=sample(c(0,1),22,replace=TRUE),drawLegend=FALSE)

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
