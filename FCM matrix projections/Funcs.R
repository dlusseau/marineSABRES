
########################################################################
#### Helper functions marine SABRES T5.3 FCM map projection ############
#### BMJV                                                   ############
########################################################################
####   4 July 2024


## Simulation  #####################

simu <- function (iter = 1000, starting.values, matrix.elems,elements = c("P","E","S","T","L","EN")){
  
  FCM.sim <- matrix(NA, length(elements), iter) # Matrix with NAs for each iteration
  FCM.sim[, 1] <- starting.values # First values are the input from stakeholders
  
  for (i in 2:iter) {
    # Each iteration the PESTLE matrix is multiplied with the previous outcome
    # FCM1.sim[, i] <- PESTLE.mat %*% matrix((FCM1.sim[, i - 1]), ncol = 1)
    FCM.sim[,i] <- t(matrix.elems) %*% matrix((FCM.sim[,i-1]), ncol = 1)
    
  }
  
  return(FCM.sim)
}

## Resilience  #####################

resilience.detracting.node.exp<-function(FCM.sim,elements =  c("P","E","S","T","L","EN"),tol=10^-5,logged=TRUE) {
  require(reshape2)
  
  iter<-ncol(FCM.sim)
  if (logged==TRUE) {
    diff.mat<-log10(FCM.sim[,2:iter])-log10(FCM.sim[,1:(iter-1)])
  } else {
    diff.mat<-(FCM.sim[,2:iter])-(FCM.sim[,1:(iter-1)])
  }
  
  stable.rate<-diff.mat[,(iter-1)]
  res.mat<-diff.mat-stable.rate
  res.mat.log<-(abs(res.mat)<tol)
  resilience.elements<-apply(res.mat.log,1,function(x) which(x==TRUE)[1])
  names(resilience.elements)<-elements
  system.resilience<-max(resilience.elements)
  
  return(list(resilience.elements=resilience.elements, resilience = system.resilience,state=(FCM.sim[,iter])))
} 


## Plot the network #####################

plot.network <- function (origin.mat = PESTLE.mat, group, save.plot = TRUE){
  
  require(igraph)
  require(GGally)
  
  FCM1.net <- graph_from_adjacency_matrix(
    (origin.mat),
    mode = "directed",
    weighted = TRUE
  )
  
  E(FCM1.net)$weights <- abs(E(FCM1.net)$weight) * 2
  E(FCM1.net)$sign <- sign(E(FCM1.net)$weight)
  E(FCM1.net)$color <- "dodgerblue"
  E(FCM1.net)$color[E(FCM1.net)$sign < 0] <- "darkorange1"
  
  output.folder <- paste0('./FCM matrix projections/res',group)
  file.name <- 'network'
  tiff.dim <- c(1500, 1000)
  
  if(save.plot) tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                     width = tiff.dim[1], height = tiff.dim[2],
                     units = "px", pointsize = 12,  res = 300,
                     compression=c("lzw"))
  
  print(
  ggnet2(FCM1.net, label = TRUE, label.size = 4, 
         arrow.size = 15, arrow.gap = 0.02, 
         edge.color = "color", edge.size = "weights")
  )
  
  if(save.plot) dev.off()
}

## Plot of system dynamics over time ###############################

plot.time.prog <- function(sim.output, group, save.plot = TRUE, 
                           elements =  c("P","E","S","T","L","EN") ){
  
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer) 
  
  row.names(sim.output) <- elements
  sim.melt <- melt(sim.output)
  
  # This is now a dataframe with the elements (Var1) and values per iteration (Var2)
  
  output.folder <- paste0('./FCM matrix projections/res',group)
  file.name <- 'progress'
  tiff.dim <- c(1500, 1000)
  
  if(save.plot) tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                     width = tiff.dim[1], height = tiff.dim[2],
                     units = "px", pointsize = 12,  res = 300,
                     compression=c("lzw"))
  
    print(
    ggplot(sim.melt, aes(x = log10(Var2), 
                         y = sign(value) * log10(abs(value)), 
                         colour = factor(Var1))) +
      geom_path()+
      scale_colour_brewer(palette = "Set3") +
      xlab('log(Time)')+
      ylab('Progress')+
      labs(colour='Element')+
      theme_minimal()
  )
  
  if(save.plot) dev.off()
  
}

## Dimensions of the basin of attraction (PCA) ###############################

plot.PCA <- function(pca, group, save.plot = TRUE){
  # library("devtools")
  # install_github("kassambara/factoextra")
  
  require(factoextra)
  
  output.folder <- paste0('./FCM matrix projections/res',group)
  file.name <- 'pca'
  tiff.dim <- c(1500, 1000)
  
  if(save.plot) tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                     width = tiff.dim[1], height = tiff.dim[2],
                     units = "px", pointsize = 12,  res = 300,
                     compression=c("lzw"))
  
  # Plot of PCA with arrows
  print(
  factoextra::fviz_pca_biplot(pca, title = NULL, xlab = 'PC1', ylab = 'PC2')
  )  # biplot(pca)
  
  if(save.plot) dev.off()
  
  # TODO: Different function, check ggplot for better graphical representation
  
  output.folder <- paste0('./FCM matrix projections/res',group)
  file.name <- 'state'
  tiff.dim <- c(1500, 1000)
  
  if(save.plot) tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                     width = tiff.dim[1], height = tiff.dim[2],
                     units = "px", pointsize = 12,  res = 300,
                     compression=c("lzw"))
  
  # Plot of the evolution of the state
  print(
    ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) +
    geom_point(colour = "darkslategrey") +
    geom_path(colour = "darkslategrey", #,
              arrow =
                arrow(angle = 15,
                      ends = "last",
                      type = "closed")
              ) + 
    theme_minimal()
  )
  
  if(save.plot) dev.off()
  
}


## Simple state map ###############################

plot.state.map <- function(states, group, save.plot = TRUE){
  
  output.folder <- paste0('./FCM matrix projections/res',group)
  file.name <- 'boolean_states'
  tiff.dim <- c(900, 900)
  
  if(save.plot) tiff(file.path(output.folder, paste0(file.name, ".tiff")),
                     width = tiff.dim[1], height = tiff.dim[2],
                     units = "px", pointsize = 12,  res = 300,
                     compression=c("lzw"))
  
    plotStateGraph(states, layout = layout.fruchterman.reingold, 
                   colorSet = c("darkslategrey","darkorange1"),
                   plotIt = T)
  
  if(save.plot) dev.off()
  
  
  # vertex_attr(states, index = 1)
  # table(V(states)$color)
  # 
  # V(states)$degree <- degree(states)
  # V(states)$attractor <- "no"
  # V(states)$attractor[V(states)$frame.color == "#000000FF"] <- "yes"
  # 
  # V(states)$color2 <- str_sub(V(states)$color, end = -3)
  # 
  # E(states)$width <- 1
  # E(states)$lty <- 1
  # 
  # E(states)$color2 <- str_sub(E(states)$color, end = -3)
  # 
  # core.state <- subgraph(states, which(V(states)$degree > 0))
  # 
  # plot(states, layout = layout.fruchterman.reingold, vertex.size = 2, 
  #      edge.arrow.size = 0.5, edge.width = 1)
  
}

## Matrix disturbances ###############################

# TODO: Make this function applicable to multiple factors and magnitudes per factor

matrix.dist <- function(elements = c("P","E","S","T","L","EN") ,
                        starting.value, PESTLE.mat, group, 
                        dist.element, dist.magnitude = 0.5, iter = 10000){
  
  # Create specific result directory
  dir.create(paste0('./FCM matrix projections/res/',group,'/sens/elem_',dist.element,'_magn_',dist.magnitude))

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


