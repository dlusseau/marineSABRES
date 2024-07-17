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

source(here::here("FCM matrix projections", "Funcs.R"))

# Load data ###############################

# !!!! Select right folder

# folder <- "C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/PESTEL analyses/"
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/PESTEL analyses/"

FCM3 <- read.csv(paste0(folder, "FCM_network_group3.csv"),
  header = T, sep = ";", dec = ","
)

FCM3a <- FCM3[which(FCM3$INFLUENCE.VALUE != 0), ]
FCM3b <- FCM3
FCM3b$INFLUENCE.VALUE[which(FCM3b$INFLUENCE.VALUE == 0)] <- 0.001 # Very small positive link
FCM3c <- FCM3
FCM3c$INFLUENCE.VALUE[which(FCM3c$INFLUENCE.VALUE == 0)] <- 0.2 # Positive link
FCM3d <- FCM3
FCM3d$INFLUENCE.VALUE[which(FCM3d$INFLUENCE.VALUE == 0)] <- 0.8 # Positive link


############################################################################
## OPTION A
############################################################################


group <- "group3a"
# Create directory where results are saved
dir.create(paste0("./FCM matrix projections/res", group))


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
plot.state.map(states = states.pestle3a, group = group)


state.map <- plotStateGraph(states.pestle3a, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle3a_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))

# a. Quantitative analysis ##################################################

## Run simulation ###############################

FCM3a.sim <- simu(starting.values = starting.value, matrix.elems = PESTLE.mat)

## Figures ###############################

p.net.gr3a <- plot.network(PESTLE.mat, group)
p.time.gr3a <- plot.time.prog(sim.output = FCM3a.sim, group, xlog = TRUE)
p.time.gr3a.subs <- plot.time.prog(sim.output = FCM3a.sim, group, xlog = F, xlims = c(0, 40))

# PCA
FCM3a.sim.trans <- t(FCM3a.sim)
colnames(FCM3a.sim.trans) <- elements
pca.FCM3a <- prcomp(FCM3a.sim.trans, scale = FALSE)

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

initial.cond.sens.group3a <- initial.cond.sens(
  matrix.elems = PESTLE.mat,
  original.res = FCM3a.sim.resilience,
  log.trans = FALSE
)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group3a <- resilience.sens(
  starting.values = starting.value,
  matrix.elems = PESTLE.mat,
  original.res = FCM3a.sim.resilience,
  log.trans = FALSE
)


############################################################################
## OPTION B
############################################################################


group <- "group3b"
# Create directory where results are saved
dir.create(paste0("./FCM matrix projections/res", group))

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
state.map <- plotStateGraph(states.pestle3b, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Simple graph
plot.state.map(states = states.pestle3b, group = group)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle3b_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))

# b. Quantitative analysis ##################################################

## Run simulation ###############################

FCM3b.sim <- simu(starting.values = starting.value, matrix.elems = PESTLE.mat)

## Figures ###############################

p.net.gr3b <- plot.network(PESTLE.mat, group)
p.time.gr3b <- plot.time.prog(sim.output = FCM3b.sim, group, xlog = TRUE)
p.time.gr3b.subs <- plot.time.prog(sim.output = FCM3b.sim, group, xlog = F, xlims = c(0, 40))

# PCA
FCM3b.sim.trans <- t(FCM3b.sim)
colnames(FCM3b.sim.trans) <- elements
pca.FCM3b <- prcomp(FCM3b.sim.trans, scale = FALSE)

p.pca.gr3b <- plot.PCA(pca = pca.FCM3b, group)

p.pca.gr3b[[1]]
p.pca.gr3b[[2]]

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

FCM3b.sim.resilience <- resilience.detracting.node.exp(FCM.sim = FCM3b.sim, logged = FALSE)

## Sensitivity initial conditions ###############################

initial.cond.sens.group3b <- initial.cond.sens(
  matrix.elems = PESTLE.mat,
  original.res = FCM3b.sim.resilience,
  log.trans = FALSE
)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group3b <- resilience.sens(
  starting.values = starting.value,
  matrix.elems = PESTLE.mat,
  original.res = FCM3b.sim.resilience,
  log.trans = FALSE
)

############################################################################
## OPTION C
############################################################################


group <- "group3c"
# Create directory where results are saved
dir.create(paste0("./FCM matrix projections/res", group))

## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements for FCM3c
elements <- c(unique(FCM3c$ELEMENT))

# Extract starting condition of each PESTLE element
starting.value <- FCM3c$ELEMENT.VALUE[!duplicated(FCM3c$ELEMENT)]
names(starting.value) <- FCM3c$ELEMENT[!duplicated(FCM3c$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM3c)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM3c$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM3c$INFLUENCE[i])] <- FCM3c$INFLUENCE.VALUE[i]
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

filename <- "PESTLE_bool_3c"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean3c <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle3c <- getAttractors(pestle_boolean3c)
state.map <- plotStateGraph(states.pestle3c, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Simple graph
plot.state.map(states = states.pestle3c, group = group)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle3c_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))

# b. Quantitative analysis ##################################################

## Run simulation ###############################

FCM3c.sim <- simu(starting.values = starting.value, matrix.elems = PESTLE.mat)

## Figures ###############################

p.net.gr3c <- plot.network(PESTLE.mat, group)
p.time.gr3c <- plot.time.prog(sim.output = FCM3c.sim, group, xlog = TRUE)
p.time.gr3c.subs <- plot.time.prog(sim.output = FCM3c.sim, group, xlog = F, xlims = c(0, 40))

# PCA
FCM3c.sim.trans <- t(FCM3c.sim)
colnames(FCM3c.sim.trans) <- elements
pca.FCM3c <- prcomp(FCM3c.sim.trans, scale = FALSE)

p.pca.gr3c <- plot.PCA(pca = pca.FCM3c, group)

p.pca.gr3c[[1]]
p.pca.gr3c[[2]]

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

FCM3c.sim.resilience <- resilience.detracting.node.exp(FCM.sim = FCM3c.sim, logged = FALSE)

## Sensitivity initial conditions ###############################

initial.cond.sens.group3c <- initial.cond.sens(
  matrix.elems = PESTLE.mat,
  original.res = FCM3c.sim.resilience,
  log.trans = FALSE
)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group3c <- resilience.sens(
  starting.values = starting.value,
  matrix.elems = PESTLE.mat,
  original.res = FCM3c.sim.resilience,
  log.trans = FALSE
)

############################################################################
## OPTION D
############################################################################


group <- "group3d"
# Create directory where results are saved
dir.create(paste0("./FCM matrix projections/res", group))

## Create initial conditions and PESTLE matrix ###############################

# Data originates from input from stakeholders
# "ELEMENT" = Originating node
# "ELEMENT.VALUE" = Present conditions of desirable PESTLE states, ranked on a scale between 0 and 10
# "INFLUENCE" = Influenced node
# "INFLUENCE.VALUE" = Agreed, signed and weighted influence from each PESTLE element onto the others

# Create PESTLE elements for FCM3d
elements <- c(unique(FCM3d$ELEMENT))

# Extract starting condition of each PESTLE element
starting.value <- FCM3d$ELEMENT.VALUE[!duplicated(FCM3d$ELEMENT)]
names(starting.value) <- FCM3d$ELEMENT[!duplicated(FCM3d$ELEMENT)]
starting.value <- t(t(starting.value))

# Create a matrix of zeros with 'elements'
PESTLE.mat <- matrix(0, length(elements), length(elements))
colnames(PESTLE.mat) <- elements
row.names(PESTLE.mat) <- elements

# Fill in the connections beteen PESTLE elements in a matrix
for (i in 1:nrow(FCM3d)) {
  PESTLE.mat[which(row.names(PESTLE.mat) == FCM3d$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM3d$INFLUENCE[i])] <- FCM3d$INFLUENCE.VALUE[i]
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

filename <- "PESTLE_bool_3d"
write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

## Load network and obtain states ##################################################

pestle_boolean3d <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle3d <- getAttractors(pestle_boolean3d)
state.map <- plotStateGraph(states.pestle3d, layout = layout.fruchterman.reingold, plotIt = FALSE)


# Simple graph
plot.state.map(states = states.pestle3d, group = group)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle3d_boolean.graphml"),
  format = "graphml"
)

# Print states
# trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
# print(getBasinOfAttraction(states, 1))

# b. Quantitative analysis ##################################################

## Run simulation ###############################

FCM3d.sim <- simu(starting.values = starting.value, matrix.elems = PESTLE.mat)

## Figures ###############################

p.net.gr3d <- plot.network(PESTLE.mat, group)
p.time.gr3d <- plot.time.prog(sim.output = FCM3d.sim, group, xlog = TRUE)
p.time.gr3d.subs <- plot.time.prog(sim.output = FCM3d.sim, group, xlog = F, xlims = c(0, 40))

# PCA
FCM3d.sim.trans <- t(FCM3d.sim)
colnames(FCM3d.sim.trans) <- elements
pca.FCM3d <- prcomp(FCM3d.sim.trans, scale = FALSE)

p.pca.gr3d <- plot.PCA(pca = pca.FCM3d, group)

p.pca.gr3d[[1]]
p.pca.gr3d[[2]]

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

FCM3d.sim.resilience <- resilience.detracting.node.exp(FCM.sim = FCM3d.sim, logged = FALSE)

## Sensitivity initial conditions ###############################

initial.cond.sens.group3d <- initial.cond.sens(
  matrix.elems = PESTLE.mat,
  original.res = FCM3d.sim.resilience,
  log.trans = FALSE
)

## Sensitivity of PESTLE matrix elements ###############################

# TODO:
# Compare the difference in the resilience between the different elements

resilience.sens.group3d <- resilience.sens(
  starting.values = starting.value,
  matrix.elems = PESTLE.mat,
  original.res = FCM3d.sim.resilience,
  log.trans = FALSE
)



############################################################################
### Compare results from group 3
############################################################################

mat3a <- FCM.mat(FCM3a)
jacobian3a <- t(mat3a$PESTLE.mat)

mat3b <- FCM.mat(FCM3b)
jacobian3b <- t(mat3b$PESTLE.mat)

mat3c <- FCM.mat(FCM3c)
jacobian3c <- t(mat3c$PESTLE.mat)

mat3d <- FCM.mat(FCM3d)
jacobian3d <- t(mat3d$PESTLE.mat)



#### Duan
lefteigen <- function(mat) {
  left <- eigen(t(mat))$vector
  return(left)
}


righteigen <- function(mat) {
  right <- eigen((mat))$vector
  return(right)
}

RJ1 <- righteigen(jacobian3a)
LJ1 <- lefteigen(jacobian3a)
EJ1 <- eigen(jacobian3a)$values

RJ2 <- righteigen(jacobian3b)
LJ2 <- lefteigen(jacobian3b)
EJ2 <- eigen(jacobian3b)$values

RJ3 <- righteigen(jacobian3c)
LJ3 <- lefteigen(jacobian3c)
EJ3 <- eigen(jacobian3c)$values

RJ4 <- righteigen(jacobian3d)
LJ4 <- lefteigen(jacobian3d)
EJ4 <- eigen(jacobian3d)$values

participation_ratio <- function(LJ, RJ) { # Duan et al. 2024
  # careful again colSums is the equivalent of sum(, axis=0) in python NOT rowsums
  # and careful again * in python is elementwise multiplication not matrix multiplication proper
  ### here I implement the code in DUan's github (distance.py line 166) however I don't think I reconcile it with the equation in the methods (Eq 5)!
  ## actually line 166 is really messed up, / x^-1 ?
  ## back to the equation instead

  PR <- rowSums(RJ * LJ)^2 / rowSums((RJ * LJ)^2)
  return(PR)
}

PR1 <- participation_ratio(LJ1, RJ1) # Politics has the broader reach
# Legal is not connected
PR1[6] <- PR1[5]
PR1[5] <- 0
PR2 <- participation_ratio(LJ2, RJ2) # Economy has the broader reach
PR3 <- participation_ratio(LJ3, RJ3) # Economy has the broader reach
PR4 <- participation_ratio(LJ4, RJ4) # Economy and Social have a broad reach

# all reals obviously because of the squaring, need to redeclare them as Real

PR.df <- data.frame(group = rep(c("group3a", "group3b", "group3c", "group3d"), each = 6), components = rep(c("P", "E", "S", "T", "L", "EN"), 4), PR = Re(c(PR1, PR2, PR3, PR4)))
PR.df$components <- factor(PR.df$components, levels = c("P", "E", "S", "T", "L", "EN"))

ggplot(PR.df, aes(x = group, y = PR, fill = components)) +
  geom_col(position = "dodge") +
  theme_minimal()
