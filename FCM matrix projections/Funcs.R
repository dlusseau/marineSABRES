########################################################################
#### Helper functions marine SABRES T5.3 FCM map projection ############
#### BMJV                                                   ############
########################################################################
####   4 July 2024

FCM.mat <- function(FCM) {
  elements <- c(unique(FCM$ELEMENT))

  starting.value <- FCM$ELEMENT.VALUE[!duplicated(FCM$ELEMENT)]
  names(starting.value) <- FCM$ELEMENT[!duplicated(FCM$ELEMENT)]
  starting.value <- t(t(starting.value))

  PESTLE.mat <- matrix(0, length(elements), length(elements))
  colnames(PESTLE.mat) <- elements
  row.names(PESTLE.mat) <- elements

  for (i in 1:nrow(FCM)) {
    PESTLE.mat[which(row.names(PESTLE.mat) == FCM$ELEMENT[i]), which(colnames(PESTLE.mat) == FCM$INFLUENCE[i])] <- FCM$INFLUENCE.VALUE[i]
  }
  return(list(PESTLE.mat = PESTLE.mat, starting.value = starting.value))
}

## Simulate matrix  #####################

simulate.mat <- function(mat) {
  mat.sim <- 1 * (mat != 0) * sign(mat) * matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat))
  
  return(mat.sim)
}

## Simulation  #####################

simu <- function(iter = 1000, starting.values, matrix.elems, elements = c("P", "E", "S", "T", "L", "EN")) {
  FCM.sim <- matrix(NA, length(elements), iter) # Matrix with NAs for each iteration
  FCM.sim[, 1] <- starting.values # First values are the input from stakeholders

  for (i in 2:iter) {
    # Each iteration the PESTLE matrix is multiplied with the previous outcome
    # FCM1.sim[, i] <- PESTLE.mat %*% matrix((FCM1.sim[, i - 1]), ncol = 1)
    FCM.sim[, i] <- t(matrix.elems) %*% matrix((FCM.sim[, i - 1]), ncol = 1)
  }

  return(FCM.sim)
}

## Boolean analysis  #####################

## Transform to binary dataset ##################################################
bin.transform <- function(mat, folder, group){
  
  # Create a binary matrix rather than a weighted one
  PESTLE.bin <- sign(mat)
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
  
  filename <- paste0("PESTLE_bool_",group)
  write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)
  
  return(boolean.df)
}

## Load network and obtain states ##################################################

pestle_boolean1 <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle1 <- getAttractors(pestle_boolean1)

# Simple graph
plot.state.map(states = states.pestle1,group = group)

state.map <- plotStateGraph(states.pestle1, layout = layout.fruchterman.reingold, plotIt = FALSE)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder,"pestle1_boolean.graphml"),
  format = "graphml"
)

# Print states
trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
print(getBasinOfAttraction(states, 1))


## Resilience  #####################

resilience.detracting.node.exp <- function(FCM.sim, elements = c("P", "E", "S", "T", "L", "EN"), tol = 10^-5, logged = TRUE) {
  require(reshape2)

  iter <- ncol(FCM.sim)

  if (logged == TRUE) {
    diff.mat <- log10(FCM.sim[, 2:iter]) - log10(FCM.sim[, 1:(iter - 1)])
  } else {
    diff.mat <- (FCM.sim[, 2:iter]) - (FCM.sim[, 1:(iter - 1)])
  }

  stable.rate <- diff.mat[, (iter - 1)]
  res.mat <- diff.mat - stable.rate
  res.mat.log <- (abs(res.mat) < tol)
  resilience.elements <- apply(res.mat.log, 1, function(x) which(x == TRUE)[1])
  names(resilience.elements) <- elements
  system.resilience <- max(resilience.elements)

  return(list(resilience.elements = resilience.elements, resilience = system.resilience, state = (FCM.sim[, iter])))
}


## Plot the network #####################

plot.network <- function(origin.mat = PESTLE.mat, group, save.plot = TRUE) {
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

  p.net <- ggnet2(FCM1.net,
    label = TRUE, label.size = 4,
    arrow.size = 15, arrow.gap = 0.02,
    edge.color = "color", edge.size = "weights"
  )

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "network"
  tiff.dim <- c(1500, 1000)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }
  print(
    p.net
  )
  if (save.plot) dev.off()

  return(p.net)
}

## Plot of system dynamics over time ###############################

plot.time.prog <- function(sim.output, group, save.plot = TRUE,
                           xlims = NULL, ylims = NULL,
                           xlog = FALSE, ylog = FALSE,
                           elements = c("P", "E", "S", "T", "L", "EN")) {
  require(reshape2)
  require(ggplot2)
  require(RColorBrewer)

  row.names(sim.output) <- elements
  sim.melt <- melt(sim.output)

  # This is now a dataframe with the elements (Var1) and values per iteration (Var2)
  p.time <- ggplot(sim.melt, aes(
    x = Var2,
    y = value,
    colour = factor(Var1)
  )) +
    geom_path() +
    scale_colour_brewer(palette = "Set3") +
    xlab("Time") +
    ylab("Progress") +
    labs(colour = "Element") +
    theme_minimal()

  if (xlog) {
    p.time <- p.time + scale_x_continuous(trans = "log10") + xlab("log(Time)")
  }
  if (ylog) {
    p.time <- p.time + scale_y_continuous(trans = "log10") + ylab("log(Progress)")
  }

  if (!is.null(xlims)) {
    p.time <- p.time + xlim(xlims)
  }
  if (!is.null(ylims)) {
    p.time <- p.time + ylim(ylims)
  }

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "progress"
  tiff.dim <- c(1500, 1000)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  print(
    p.time
  )

  if (save.plot) dev.off()

  return(p.time)
}

## Dimensions of the basin of attraction (PCA) ###############################

plot.PCA <- function(pca, group,
                     xlims = NULL, ylims = NULL,
                     xlog = NULL, ylog = NULL,
                     save.plot = TRUE) {
  # library("devtools")
  # install_github("kassambara/factoextra")

  require(factoextra)

  biplot <- factoextra::fviz_pca_biplot(pca,
    title = " ",
    xlab = "PC1",
    ylab = "PC2",
    geom = c("point"), # , "text"
    alpha.ind = 0.8,
    col.ind = "darkgrey",
    col.var = "dodgerblue",
    label = "var",
    ggtheme = theme_minimal()
  ) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 2)) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 2))

  # Manipulation of axes
  if (!is.null(xlog)) {
    biplot <- biplot + scale_x_continuous(breaks = scales::breaks_extended(n = 2), trans = "log10") + xlab("log(Time)")
  }
  if (!is.null(ylog)) {
    biplot <- biplot + scale_y_continuous(breaks = scales::breaks_extended(n = 2), trans = "log10") + ylab("log(Progress)")
  }

  if (!is.null(xlims)) {
    biplot <- biplot + xlim(xlims)
  }
  if (!is.null(ylims)) {
    biplot <- biplot + ylim(ylims)
  }


  # biplot
  #
  # layer_scales(biplot)$y$get_limits() #-3.319612e+181  7.967070e+182
  # layer_scales(biplot)$x$get_limits() #-1.511836e+196  5.516300e+198
  #
  # layer_scales(biplot)$x$get_limits()[1]


  trajectory <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) +
    geom_point(colour = "dodgerblue") +
    geom_path(
      colour = "dodgerblue" # ,
      # arrow =
      #   arrow(angle = 15,
      #         ends = "last",
      #         type = "closed")
    ) +
    theme_minimal()

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "pca"
  tiff.dim <- c(1500, 1000)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  # Plot of PCA with arrows
  print(
    biplot
  ) # biplot(pca)

  if (save.plot) dev.off()

  # TODO: Different function, check ggplot for better graphical representation

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "state"
  tiff.dim <- c(1500, 1000)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  # Plot of the evolution of the state
  print(
    trajectory
  )

  if (save.plot) dev.off()

  return(list(biplot, trajectory))
}


## Simple state map ###############################

plot.state.map <- function(states, group, save.plot = TRUE) {
  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "boolean_states"
  tiff.dim <- c(900, 900)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  plotStateGraph(states,
    layout = layout.fruchterman.reingold,
    colorSet = c("darkslategrey", "darkorange1"),
    plotIt = T
  )

  if (save.plot) dev.off()
}

## Initial condition resilience sensitivity ###############################

initial.cond.sens <- function(matrix.elems, original.res, log.trans,
                              elements = c("P", "E", "S", "T", "L", "EN"), save.plot = TRUE) {
  sens.initial.conditions <- list()
  res.initial.conditions <- list()

  # Include also the PESTLE elements in the df? Keep it simple for now
  stability.initial.conditions.df <- data.frame(matrix(NA, nrow = 1000, ncol = 1))
  colnames(stability.initial.conditions.df) <- "resilience.diff"

  for (m in 1:1000) {
    sens.initial.conditions[[m]] <- simu(starting.values = runif(6, 0, 10), matrix.elems = matrix.elems)
    res.initial.conditions[[m]] <- resilience.detracting.node.exp(FCM.sim = sens.initial.conditions[[m]], logged = log.trans)

    # Save resilience results in dataframe
    stability.initial.conditions.df[m, 1] <- res.initial.conditions[[m]]$resilience
  }

  # Plot first run
  # plot.time.prog(sens.initial.conditions[[1]], diff = F, group, save.plot = T,
  #                xlims = c(0,50),file.name = 'sens_initial_cond_progress')
  # plot.time.prog(sens.initial.conditions[[1]], diff = T, group, save.plot = T,
  #                xlims = c(0,50),file.name = 'sens_initial_cond_progress_diff')

  # Plot distribution of resilience changes
  diff.main.init.conditions <- stability.initial.conditions.df - original.res$resilience

  # Density plot
  ini.cond.resilience <- ggplot(data = diff.main.init.conditions, aes(x = resilience.diff)) +
    geom_density(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.7) +
    theme_minimal() +
    xlab("Resilience (difference compared to baseline)") +
    ylab("Density")

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "initial.condition.resilience"
  tiff.dim <- c(2000, 1500)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  print(ini.cond.resilience)

  if (save.plot) dev.off()

  return(list(ini.cond.resilience, sens.initial.conditions, res.initial.conditions, stability.initial.conditions.df))
}

## System resilience sensitivity ###############################

resilience.sens <- function(starting.values, matrix.elems, original.res, log.trans,
                            elements = c("P", "E", "S", "T", "L", "EN"), save.plot = TRUE) {
  # Create empty list to save results
  sens.res <- list()
  sens.res.elem <- list()

  # vector of all PESTLE matrix elements to disturb
  dist.elems <- matrix.elems[which(matrix.elems != 0)]
  # There are duplicates in this vector -- loop over positions of those elements
  duplicated(dist.elems)

  dist.elems.pos <- which(matrix.elems != 0, arr.ind = T)

  # Create df to save resilience results
  sens.res.resilience.df <- data.frame(matrix(NA, nrow = 1000, ncol = length(dist.elems)))
  colnames(sens.res.resilience.df) <- seq(1:length((dist.elems)))


  for (i in 1:length(dist.elems)) {
    # print(i)
    dist <- dist.elems[i] + rnorm(1000, 0, 0.01)

    for (j in 1:length(dist)) {
      # print(j)

      PESTLE.mat.dist <- matrix.elems

      PESTLE.mat.dist[dist.elems.pos[i, ][1], dist.elems.pos[i, ][2]] <- dist[j]

      # Simulation
      simu.res <- simu(starting.values = starting.values, matrix.elems = PESTLE.mat.dist)
      sens.res.elem[[j]] <- simu.res

      # Resilience
      sens.res.resilience <- resilience.detracting.node.exp(FCM.sim = simu.res, logged = log.trans)

      # Save resilience results in dataframe
      sens.res.resilience.df[j, i] <- sens.res.resilience$resilience

      # TODO:
      # Save resilience results in dataframe PER PESTLE ELEMENT - report in one matrix per element
      # sens.res.resilience.df[j,i] <-  sens.res.resilience$resilience
    }
    # Save all simu outcomes in a list
    sens.res[[i]] <- sens.res.elem
  }

  # Difference between original and new matrix resilience.
  diff.main <- sens.res.resilience.df - original.res$resilience
  diff.main.melted <- melt(diff.main)

  # Include matrix positions
  diff.main.melted$ELEMENT <- dist.elems.pos[diff.main.melted$variable, ][, 1]
  diff.main.melted$INFLUENCE <- dist.elems.pos[diff.main.melted$variable, ][, 2]

  # Rename columns
  colnames(diff.main.melted) <- c("dist.elem", "value", "ELEMENT", "INFLUENCE")

  # Rename ELEMENT and INFLUENCE according to PESTLE elements
  diff.main.melted$ELEMENT <- elements[diff.main.melted$ELEMENT]
  diff.main.melted$INFLUENCE <- elements[diff.main.melted$INFLUENCE]

  res.sens <- ggplot(diff.main.melted, aes(x = value)) +
    geom_density(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.7) +
    theme_minimal() +
    xlab("Resilience (difference compared to baseline)") +
    ylab("Density") +
    facet_grid(factor(ELEMENT, levels = elements) ~ factor(INFLUENCE, levels = elements),
      scales = "free_y"
    )

  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "system.resilience"
  tiff.dim <- c(2000, 1500)

  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  print(res.sens)

  if (save.plot) dev.off()

  return(list(res.sens, sens.res, sens.res.resilience.df))
}

## Greedy approach to state shift  #####################

state.shift <- function(greed, iter, mat, tol = 0.000001, elements = c("P", "E", "S", "T", "L", "EN")){
  
  require(reshape2)
  
  starting.value <- runif(nrow(mat), 0, 1) # Create random starting values between 0 and 1
  
  FCM.obs <- simu(matrix.elems = mat, iter, starting.value) # Simulate outcomes
  
  state.obs <- FCM.obs[, iter]
  names(state.obs) <- elements 
  
  state.sim <- array(NA, dim = c(length(state.obs), greed))
  rownames(state.sim) <- elements
  
  mat.m <- melt(mat)
  mat.sim.df <- array(NA, dim = c(nrow(mat.m), greed + 2))
  mat.sim.df[, 1] <- mat.m$Var1
  mat.sim.df[, 2] <- mat.m$Var2
  
  tic <- Sys.time()
  
  for (i in 1:greed) {
    mat.sim <- simulate.mat(mat)
    
    FCM <- simu(matrix.elems = mat.sim, iter, starting.value)
    
    mat.sim.df[, i + 2] <- melt(mat.sim)$value
    state.sim[, i] <- FCM[, iter] - matrix(state.obs, ncol = 1)
  }
  toc <- Sys.time() - tic
  
  state.sim.bin <- state.sim
  state.sim.bin[abs(state.sim.bin) < tol] <- 0
  state.sim.bin <- sign(state.sim.bin)
  
  return(list(state.sim = state.sim, state.sim.bin = state.sim.bin, mat.sim.df = as.data.frame(mat.sim.df)))
}
