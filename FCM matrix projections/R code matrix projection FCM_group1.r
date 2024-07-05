##############################################################
#### R code #marine SABRES T5.3 FCM map projection############
##############################################################
####   3 July 2024

# Load packages ###############################

# library(ggplot2)
# library(igraph)
# library(ggnet)
library(BoolNet)
library(stringr)
# library(reshape2)

source('Funcs.R')

# Load data ###############################

# !!!! Select right folder

# folder <- "C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/PESTEL analyses/"
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/PESTEL analyses/"

FCM1 <- read.csv(paste0(folder,"FCM_network_group1.csv"), header = T)

# Quantitative analysis ##################################################

## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements
elements <- c(unique(FCM1$ELEMENT))

# Extract starting condition of each PESTLE element 
starting.value <- FCM1$ELEMENT.VALUE[!duplicated(FCM1$ELEMENT)]
names(starting.value) <- FCM1$ELEMENT[!duplicated(FCM1$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM1)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM1$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM1$INFLUENCE[i])] <- FCM1$INFLUENCE.VALUE[i]
}

## Run simulation ###############################

iter <- 1000
FCM1.sim <- matrix(NA, length(elements), iter) # Matrix with NAs for each iteration
FCM1.sim[, 1] <- starting.value # First values are the input from stakeholders

for (i in 2:1000) {
  # Each iteration the PESTLE matrix is multiplied with the previous outcome
  FCM1.sim[, i] <- PESTLE.mat %*% matrix((FCM1.sim[, i - 1]), ncol = 1)
}

## Figures ###############################

plot.network(PESTLE.mat)
plot.time.prog(FCM1.sim)
plot.PCA(FCM1.sim)

## Stability of graph Laplacian a la Brownski #################################

PESTLE.Lap <- PESTLE.mat - diag(rowSums(PESTLE.mat)) ## following Bronski & Deville 2014 SIAM Appl Math (signed graphs) contrary to the usual L=D-A
## and it works as sumRows(sdglowL)= array(0)

# need to get it the right way around the diag sum is in-strength

PESTLE.Lap <- diag(rowSums(t(PESTLE.mat))) - t(PESTLE.mat)
eigen(PESTLE.Lap)

# The dominant (largest) eigenvalue indicates the dominant state
# The corresponding eigenvector indicates the rank of the PESTLE factors

# Boolean analysis ##################################################

## Transform to binary dataset ##################################################

# Create a binary matrix rather than a weighted one
PESTLE.bin <- sign(PESTLE.mat)
boolean.df <- data.frame(targets = factor(colnames(PESTLE.bin)), factors = NA)

# Create csv with 'or' (|) an negations (!)
for (i in 1:ncol(PESTLE.bin)) {
  poss <- names(which(PESTLE.bin[, i] == 1))
  negs <- names(which(PESTLE.bin[, i] == -1))
  if (length(negs) > 0) {
    negs <- paste0("!", negs)
  }
  all <- c(poss, negs)

  boolean.df$factors[i] <- paste(all, collapse = "|")
}

filename <- "PESTLE_bool_1"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean <- loadNetwork(paste0(folder, filename, ".csv"))
states <- getAttractors(pestle_boolean)

state.map <- plotStateGraph(states, layout = layout.fruchterman.reingold, plotIt = FALSE)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder,"pestle1_boolean.graphml"),
  format = "graphml"
)

# Simple graph
plot.state.map(state.map)

# Print states
trans.tab <- getTransitionTable(states)
plotStateGraph(states)
print(getBasinOfAttraction(states, 1))

##################################################
### back to the weighted matrix
### systematic pressed disturbance on systematic
### is the system behaving as we expect?

#######
#######
####### this is the code to take to do the matrix projection simulation
#######
### in the way the matrix multiplication takes place here
### we need to make sure that we deal with an assortment
### where columns contributes to rows (arrow direction)
### column is the starting point of the arrow
### and row is end point of the arrow

matrix.dist(starting.value, PESTLE.mat, 
            dist.element = 'P', dist.magnitude = 0.5, 
            iter = 10000)
  

######################################################


##### pressed disturbance

iter=10000
FCM1.sim<-matrix(NA,length(elements),iter)
FCM1.sim[,1]<-starting.value

for ( i in 2:10000) {
  
  
  disturbed<-FCM1.sim[,i-1]
  disturbed[6]<-disturbed[6]+.5
  disturbed[3]<-disturbed[3]+.5
  
  FCM1.sim[,i]<-t(PESTLE.mat)%*%matrix((disturbed), ncol = 1)
  #FCM1.sim[,i]<-t(PESTLE.mat)%*%matrix((FCM1.sim[,i-1]), ncol = 1)
  
}


row.names(FCM1.sim)<-elements
sim.melt<-melt(FCM1.sim)

ggplot(sim.melt,aes(x=log10(Var2),y=((value)),colour=factor(Var1)))+
  geom_path()



ggplot(subset(sim.melt,Var2>900&Var2<990),aes(x=(Var2),y=((value)),colour=factor(Var1)))+
  geom_path()


#####################
###sensitivity
sen.iter<-10000
FCM1.sen.init<-data.frame(P=rep(NA,sen.iter),E=rep(NA,sen.iter),S=rep(NA,sen.iter),T=rep(NA,sen.iter),L=rep(NA,sen.iter),EN=rep(NA,sen.iter))

for (j in 1:sen.iter) {
  
  iter=1000
  FCM1.sim<-matrix(NA,length(elements),iter)
  FCM1.sim[,1]<-runif(6,0,1)
  
  for ( i in 2:iter) {
    
    
    disturbed<-FCM1.sim[,i-1]
    disturbed[6]<-disturbed[6]+.5
    
    FCM1.sim[,i]<-t(PESTLE.mat)%*%matrix((disturbed), ncol = 1)
    #FCM1.sim[,i]<-t(PESTLE.mat)%*%matrix((FCM1.sim[,i-1]), ncol = 1)
    
  }
  
  
  FCM1.sen.init[j,]<-FCM1.sim[,iter]
  
}

######


