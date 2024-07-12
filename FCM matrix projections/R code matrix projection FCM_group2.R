##############################################################
#### R code #marine SABRES T5.3 FCM map projection############
#### Group 2                                      ############
##############################################################
####   4 July 2024

# Load packages ###############################

# library(ggplot2)
# library(igraph)
# library(ggnet)
library(BoolNet)
library(stringr)
# library(reshape2)

set.seed(123)

source(here::here("FCM matrix projections", "Funcs.R"))

group <- "group2"

# Load data ###############################

# !!!! Select right folder

# folder <- "C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/PESTEL analyses/"
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/PESTEL analyses/"


FCM2 <- read.csv(paste0(folder, "FCM_network_group2.csv"),
  header = T, sep = ";", dec = ","
)
# Create directory where results are saved
dir.create(paste0("./FCM matrix projections/res", group))


## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements
elements <- c(unique(FCM2$ELEMENT))

# Extract starting condition of each PESTLE element
starting.value <- FCM2$ELEMENT.VALUE[!duplicated(FCM2$ELEMENT)]
names(starting.value) <- FCM2$ELEMENT[!duplicated(FCM2$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM2)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM2$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM2$INFLUENCE[i])] <- FCM2$INFLUENCE.VALUE[i]
}

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

filename <- "PESTLE_bool_2"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean2 <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle2 <- getAttractors(pestle_boolean2)

# Simple graph
plot.state.map(states = states.pestle2, group = group)

state.map <- plotStateGraph(states.pestle2, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle2_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))


# Quantitative analysis ##################################################
## Run simulation ###############################

FCM2.sim <- simu(starting.values = starting.value, matrix.elems = PESTLE.mat)

## Figures ###############################

p.net.gr2 <- plot.network(PESTLE.mat, group)
p.time.gr2 <- plot.time.prog(sim.output = FCM2.sim, group, xlog = T, ylog = T)
p.time.gr2.subs <- plot.time.prog(sim.output = FCM2.sim, group, xlog = F, ylog = T, xlims = c(0, 200), file.name = "progress.subs")
p.time.gr2.subs

# PCA
FCM2.sim.trans <- t(FCM2.sim)
colnames(FCM2.sim.trans) <- elements
pca.FCM2 <- prcomp(FCM2.sim.trans, scale = FALSE)

p.pca.gr2 <- plot.PCA(pca = pca.FCM2, group)


p.pca.gr2[[1]]
p.pca.gr2[[2]]

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

# TODO: make sure that logged = TRUE here

FCM2.sim.resilience <- resilience.detracting.node.exp(FCM.sim = FCM2.sim, logged = TRUE)

## Sensitivity initial conditions ###############################

initial.cond.sens.group2 <- initial.cond.sens(
  matrix.elems = PESTLE.mat,
  original.res = FCM2.sim.resilience,
  log.trans = TRUE
)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group2 <- resilience.sens(
  starting.values = starting.value,
  matrix.elems = PESTLE.mat,
  original.res = FCM2.sim.resilience,
  log.trans = TRUE
)
