
########################################################################
#### Helper functions marine SABRES T5.3 FCM map projection ############
#### BMJV                                                   ############
########################################################################
####   4 July 2024

## Plot the network #####################

plot.network <- function (origin.mat = PESTLE.mat){
  require(igraph)
  require(GGally)
  
  FCM1.net <- graph_from_adjacency_matrix(
    (origin.mat),
    mode = "directed",
    weighted = TRUE
  )
  
  E(FCM1.net)$weights <- abs(E(FCM1.net)$weight) * 2
  E(FCM1.net)$sign <- sign(E(FCM1.net)$weight)
  E(FCM1.net)$color <- "blue"
  E(FCM1.net)$color[E(FCM1.net)$sign < 0] <- "red"
  
  ggnet2(FCM1.net, label = TRUE, label.size = 4, 
         arrow.size = 15, arrow.gap = 0.02, 
         edge.color = "color", edge.size = "weights")
}

## Plot of system dynamics over time ###############################

plot.time.prog <- function(sim.output = FCM1.sim){
  require(reshape2)
  require(ggplot2)
  
  sim.melt <- melt(sim.output)
  row.names(sim.output) <- elements
  
  # This is now a dataframe with the elements (Var1) and values per iteration (Var2)
  
  ggplot(sim.melt, aes(x = log10(Var2), y = sign(value) * log10(abs(value)), colour = factor(Var1))) +
    geom_path()+
    xlab('log(Time)')+
    ylab('Progress')+
    theme_minimal()
  
}

## Dimensions of the basin of attraction (PCA) ###############################

plot.PCA <- function(sim.output = FCM1.sim){
  
  # PCA
  pca <- prcomp(t(sim.output), scale = FALSE)
  
  # Plot of PCA with arrows
  biplot(pca)
  
  # TODO: Different function, check ggplot for better graphical representation
  
  # Plot of the evolution of the state
  ggplot(as.data.frame(pca$x), aes(x = (PC1), y = (PC2))) +
    geom_point() +
    geom_path(arrow = arrow()) +
    theme_minimal()
  
}

## Simple state map ###############################

plot.state.map <- function(states = state.map){
  
  vertex_attr(states, index = 1)
  table(V(states)$color)
  
  V(states)$degree <- degree(states)
  V(states)$attractor <- "no"
  V(states)$attractor[V(states)$frame.color == "#000000FF"] <- "yes"
  
  V(states)$color2 <- str_sub(V(states)$color, end = -3)
  
  E(states)$width <- 1
  E(states)$lty <- 1
  
  E(states)$color2 <- str_sub(E(states)$color, end = -3)
  
  core.state <- subgraph(states, which(V(states)$degree > 0))
  
  plot(states, layout = layout.fruchterman.reingold, vertex.size = 2, 
       edge.arrow.size = 0.5, edge.width = 1)
  
}

## Matrix disturbances ###############################

# TODO: Make this function applicable to multiple factors and magnitudes per factor

matrix.dist <- function(elements = c("P","E","S","T","L","EN") ,
                        starting.value, PESTLE.mat, 
                        dist.element, dist.magnitude = 0.5, iter = 10000){

  dist.element.index <- which(elements %in% dist.element)
  
  FCM1.sim <- matrix(NA,length(elements),iter)
  FCM1.sim[,1] <- starting.value
  
  for (i in 2:iter) {

    disturbed <- FCM1.sim[,i-1]
    disturbed[dist.element.index] <- disturbed[dist.element.index] + dist.magnitude
    
    FCM1.sim[,i]<-t(PESTLE.mat)%*%matrix((disturbed), ncol = 1)
 
  }
  
  print(paste('Matrix element',dist.element,'is distrubed with a magnitude of',dist.magnitude))
  
  plot.network(PESTLE.mat)
  plot.time.prog(FCM1.sim)
  plot.PCA(FCM1.sim)
  
}


