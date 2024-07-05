##############################################################
#### R code #marine SABRES T5.3 FCM map projection############
#### Group 2                                      ############
##############################################################
####   4 July 2024

# Load packages ###############################

library(ggplot2)
library(igraph)
# library(ggnet)
library(BoolNet)
library(stringr)
library(reshape2)

# Load data ###############################

# !!!! Select right folder
# folder <- "C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/SABRES/PESTEL analyses/"
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/PESTEL analyses/"

FCM2 <- read.csv(paste0(folder,"FCM_network_group2.csv"),
                 header = T, sep = ";", dec = ","
)

# Quantitative analysis ##################################################

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

## Run simulation ###############################

iter <- 1000
FCM2.sim <- matrix(NA, length(elements), iter) # Matrix with NAs for each iteration
FCM2.sim[, 1] <- starting.value # First values are the input from stakeholders

for (i in 2:1000) {
  # Each iteration the PESTLE matrix is multiplied with the previous outcome
  FCM2.sim[, i] <- t(PESTLE.mat) %*% matrix((FCM2.sim[, i - 1]), ncol = 1)
}

# head(FCM2.sim)

## Explore the outcomes ###############################

sim.melt <- melt(FCM2.sim)
row.names(FCM2.sim) <- elements

# This is now a dataframe with the elements (Var1) and values per iteration (Var2)

## Plot of system dynamics over time ###############################

ggplot(sim.melt, aes(x = log10(Var2), y = sign(value) * log10(abs(value)), colour = factor(Var1))) +
  geom_path()

# Plot of the cyclical fluctuations of all variables at the end of the timeseries
ggplot(subset(sim.melt, Var2 > 950), aes(x = (Var2), y = sign(value) * log10(abs(value)), colour = factor(Var1))) +
  geom_path()

## Dimensions of the basin of attraction (PCA) ###############################

# PCA
pca <- prcomp(t(FCM2.sim), scale = FALSE)

# Plot of PCA with arrows
biplot(pca)

# Plot to see if the 
ggplot(as.data.frame(pca$x), aes(x = (PC1), y = (PC2))) +
  geom_point() +
  geom_path(arrow = arrow()) + #
  theme_minimal()
# This is a 'detracting node', the opposite form group 1 (attracting spiral)

## Plot the network - not working for now#####################

#### some some is going on with ggnet using network and not recognising igraph input

FCM2.net <- graph_from_adjacency_matrix(
  (PESTLE.mat),
  mode = "directed",
  weighted = TRUE
)

E(FCM2.net)$weights <- abs(E(FCM2.net)$weight) * 2
E(FCM2.net)$sign <- sign(E(FCM2.net)$weight)
E(FCM2.net)$color <- "blue"
E(FCM2.net)$color[E(FCM2.net)$sign < 0] <- "red"

ggnet2(FCM2.net, label = TRUE, label.size = 4, arrow.size = 15, arrow.gap = 0.02, edge.color = "color", edge.size = "weights")

## Stability of graph Laplacian a la Brownski #################################
### 

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

pestle_boolean <- loadNetwork(paste0(folder, filename, ".csv"))
states <- getAttractors(pestle_boolean)

state.map <- plotStateGraph(states, layout = layout.fruchterman.reingold, plotIt = FALSE)
vertex_attr(state.map, index = 1)
table(V(state.map)$color)

V(state.map)$degree <- degree(state.map)
V(state.map)$attractor <- "no"
V(state.map)$attractor[V(state.map)$frame.color == "#000000FF"] <- "yes"

V(state.map)$color2 <- str_sub(V(state.map)$color, end = -3)

E(state.map)$width <- 1
E(state.map)$lty <- 1

E(state.map)$color2 <- str_sub(E(state.map)$color, end = -3)

core.state <- subgraph(state.map, which(V(state.map)$degree > 0))

plot(state.map, layout = layout.fruchterman.reingold, vertex.size = 2, edge.arrow.size = 0.5, edge.width = 1)

# Write graph: fial figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder,"pestle2_boolean.graphml"),
  format = "graphml"
)


trans.tab <- getTransitionTable(states)
# plotStateGraph(states)
print(getBasinOfAttraction(states, 1))

