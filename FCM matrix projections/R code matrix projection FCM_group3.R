##############################################################
#### R code #marine SABRES T5.3 FCM map projection############
#### Group 3                                      ############
##############################################################
####   4 July 2024

# This group is divided into two options:
# a. excluding the 'connection, but no values'
# b. including those factors only in Boolean analysis, and there assume they are 1

# Load packages ###############################

# library(ggplot2)
# library(igraph)
# library(ggnet)
library(BoolNet)
library(stringr)
# library(reshape2)

set.seed(123)

source(here::here('FCM matrix projections', 'Funcs.R'))

# Load data ###############################

# !!!! Select right folder

# folder <- "C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/PESTEL analyses/"
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/PESTEL analyses/"

FCM3 <- read.csv(paste0(folder,"FCM_network_group3.csv"),
                 header = T, sep = ";", dec = ","
)

FCM3a <- FCM3[which(FCM3$INFLUENCE.VALUE != 0),]
FCM3b <- FCM3
FCM3b$INFLUENCE.VALUE[which(FCM3b$INFLUENCE.VALUE == 0)] <- 0.001 # Only for Boolean analysis

############################################################################
## OPTION A
############################################################################


group <- 'group3a'
# Create directory where results are saved
dir.create(paste0('./FCM matrix projections/res',group))


## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements for FCM3a
elements <- c(unique(FCM3$ELEMENT))

# Extract starting condition of each PESTLE element 
starting.value <- FCM3$ELEMENT.VALUE[!duplicated(FCM3$ELEMENT)]
names(starting.value) <- FCM3$ELEMENT[!duplicated(FCM3$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM3a)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM3a$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM3a$INFLUENCE[i])] <- FCM3a$INFLUENCE.VALUE[i]
}

# a. Boolean analysis ##################################################

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

filename <- "PESTLE_bool_3a"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean3a <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle3a <- getAttractors(pestle_boolean3a)

# Simple graph
plot.state.map(states = states.pestle3a,group = group)


state.map <- plotStateGraph(states.pestle3a, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder,"pestle3a_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))

# a. Quantitative analysis ##################################################

## Run simulation ###############################

FCM3a.sim <- simu(starting.values = starting.value, matrix.elems= PESTLE.mat)

## Figures ###############################

p.net.gr3a <- plot.network(PESTLE.mat, group)
p.time.gr3a <- plot.time.prog(sim.output = FCM3a.sim, group)

# PCA
pca.FCM3a <- prcomp(t(FCM3a.sim), scale = FALSE)

p.pca.gr3a <- plot.PCA(pca = pca.FCM3a, group)

p.pca.gr3a[[1]]
p.pca.gr3a[[2]]

## Stability of graph Laplacian a la Brownski #################################

PESTLE.Lap <- t(PESTLE.mat) - diag(rowSums(t(PESTLE.mat))) ## following Bronski & Deville 2014 SIAM Appl Math (signed graphs) contrary to the usual L=D-A
## and it works as sumRows(sdglowL)= array(0)

# need to get it the right way around the diag sum is in-strength

# PESTLE.Lap <- diag(rowSums(t(PESTLE.mat))) - t(PESTLE.mat) # This is another notation of the laplacian, where the sign is the opposite of the previous Laplacian. This is however the 'formal' formulation
eigen(PESTLE.Lap)

# We don't have any antagonistic PESTLE elements (the eigenvalues are all negative or zero)
# Therefore everything moves into the same direction; any disturbance will 
# likely have very little effect on the system.

# The dominant (largest) eigenvalue indicates the dominant state
# The corresponding eigenvector indicates the rank of the PESTLE factors

## Resilience #################################

FCM3a.sim.resilience <- resilience.detracting.node.exp(FCM.sim = FCM3a.sim, logged = FALSE)

## Sensitivity initial conditions ###############################

initial.cond.sens.group3a <- initial.cond.sens(matrix.elems = PESTLE.mat,
                                              original.res = FCM3a.sim.resilience,
                                              log.trans = FALSE)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group3a <- resilience.sens(starting.values = starting.value, 
                                          matrix.elems = PESTLE.mat,
                                          original.res = FCM3a.sim.resilience,
                                          log.trans = FALSE)


############################################################################
## OPTION B
############################################################################


group <- 'group3b'
# Create directory where results are saved
dir.create(paste0('./FCM matrix projections/res',group))

## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements for FCM3b
elements <- c(unique(FCM3b$ELEMENT))

# Extract starting condition of each PESTLE element 
starting.value <- FCM3b$ELEMENT.VALUE[!duplicated(FCM3b$ELEMENT)]
names(starting.value) <- FCM3b$ELEMENT[!duplicated(FCM3b$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM3b)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM3b$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM3b$INFLUENCE[i])] <- FCM3b$INFLUENCE.VALUE[i]
}

# b. Boolean analysis ##################################################

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

filename <- "PESTLE_bool_3b"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean3b <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle3b <- getAttractors(pestle_boolean3b)

# Simple graph
plot.state.map(states = states.pestle3b,group = group)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder,"pestle3b_boolean.graphml"),
  format = "graphml"
)

# Print states
trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
print(getBasinOfAttraction(states, 1))

# Quantitative analysis ##################################################

## Run simulation ###############################

iter <- 1000
FCM3b.sim <- matrix(NA, length(elements), iter) # Matrix with NAs for each iteration
FCM3b.sim[, 1] <- starting.value # First values are the input from stakeholders

for (i in 2:1000) {
  # Each iteration the PESTLE matrix is multiplied with the previous outcome
  FCM3b.sim[, i] <- PESTLE.mat %*% matrix((FCM3b.sim[, i - 1]), ncol = 1)
}

## Figures ###############################

plot.network(PESTLE.mat, group)
plot.time.prog(sim.output = FCM3b.sim, group)
plot.PCA(FCM3b.sim, group)

## Stability of graph Laplacian a la Brownski #################################

PESTLE.Lap <- t(PESTLE.mat) - diag(rowSums(t(PESTLE.mat))) ## following Bronski & Deville 2014 SIAM Appl Math (signed graphs) contrary to the usual L=D-A
## and it works as sumRows(sdglowL)= array(0)

# need to get it the right way around the diag sum is in-strength

# PESTLE.Lap <- diag(rowSums(t(PESTLE.mat))) - t(PESTLE.mat) # This is another notation of the laplacian, where the sign is the opposite of the previous Laplacian. This is however the 'formal' formulation
eigen(PESTLE.Lap)

# We don't have any antagonistic PESTLE elements (the eigenvalues are all negative or zero)
# Therefore everything moves into the same direction; any disturbance will 
# likely have very little effect on the system.

# The dominant (largest) eigenvalue indicates the dominant state
# The corresponding eigenvector indicates the rank of the PESTLE factors


