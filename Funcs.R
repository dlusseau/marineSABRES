
########################################################################
#### Helper functions marine SABRES T5.3 FCM map projection ############
#### BMJV                                                   ############
########################################################################
####   4 July 2024



## Plot the network #####################

PlotNetwork <- function (origin.mat = PESTLE.mat){
  require(igraph)
  
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

PlotTimeProg <- function(sim.output = FCM1.sim){
  require(reshape2)
  require(ggplot2)
  
  sim.melt <- melt(sim.output)
  row.names(sim.output) <- elements
  
  # This is now a dataframe with the elements (Var1) and values per iteration (Var2)
  
  ggplot(sim.melt, aes(x = log10(Var2), y = sign(value) * log10(abs(value)), colour = Var1)) +
    geom_path()+
    xlab('log(Time)')+
    ylab('Progress')+
    theme_minimal()
  
}

## Dimensions of the basin of attraction (PCA) ###############################

PlotPCA <- function(sim.output = FCM1.sim){
  
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
