########################################################################
####                                                        ############
#### Helper functions marineSABRES T5.3 network analyses   ############
####                                                        ############
########################################################################
####
####   davlu
####   bmjv
####   Start; October 2024
####   list of utility functions to walk through interventions SABRES analytical pipelines
####   See Demo.R for a demo workflow
####
########################################################################


# Function to create a matrix from 'from', 'to', and 'weight' vectors
make_matrix <- function(from, to, weight) {
  # Combine unique elements from 'from' and 'to' to create a list of all elements
  elements <- unique(c(unique(from), unique(to)))

  # Initialize an empty square matrix with rows and columns named by the elements
  SES.mat <- matrix(0, length(elements), length(elements))
  colnames(SES.mat) <- elements
  rownames(SES.mat) <- elements

  # Populate the SES matrix with weights based on 'from' and 'to' connections
  for (i in 1:length(from)) {
    SES.mat[rownames(SES.mat) == from[i], colnames(SES.mat) == to[i]] <- weight[i]
  }
  return(SES.mat) # Return the completed SES matrix
}

# Function to plot the SES network
plot.SES <- function(SES.mat, folder, filename, title, w = 80, h = 40, layout = layout_with_fr, label.cex = 1.5, vsize = 20, eweight = 10) {
  # Create a directed graph from the SES matrix
  SES.graph <- graph_from_adjacency_matrix(SES.mat, mode = "directed", weighted = TRUE)
  
  # Set vertex properties
  V(SES.graph)$color <- "white" # Set vertex color
  V(SES.graph)$label.cex <- label.cex # Set label size
  V(SES.graph)$size <- vsize # Set vertex size
  
  # Set edge properties
  E(SES.graph)$sign <- sign(E(SES.graph)$weight) # Determine the sign of edge weights
  E(SES.graph)$color <- "blue" # Default edge color
  E(SES.graph)$color[E(SES.graph)$sign < 0] <- "red" # Change color for negative edges
  E(SES.graph)$weight <- abs(E(SES.graph)$weight) # Use absolute values for edge weights
  E(SES.graph)$curved <- 0.2 # Set edge curvature
  E(SES.graph)$width <- E(SES.graph)$weight * eweight # Scale edge width by weight
  E(SES.graph)$arrow.size <- E(SES.graph)$weight * (eweight / 2) # Set arrow size proportional to weight
  
  # Define layout for the plot
  l <- layout(SES.graph)
  
  # Save the plot as a PNG file
  png(paste0(folder, filename, ".png"), width = w, height = h, units = "cm", res = 200)
  plot(SES.graph, layout = l, ylab = title) # Plot the network with defined layout and y-axis label
  dev.off() # Close the PNG device
  
  return(SES.graph) # Return the constructed network
}


# Function to load data and optionally plot the SES network
data.load <- function(df, folder, graph = FALSE, graph.name = NULL, graph.title = NULL) {
  
  # Assign default weight of 1 to all connections
  df$weight <- 1
  
  # Adjust weights based on the strength of the connections
  df$weight[df$strength == "Medium Positive"] <- 0.5
  df$weight[df$strength == "Medium Negative"] <- -0.5
  df$weight[df$strength == "Strong Negative"] <- -1
  
  # Determine the sign of the weights
  df$sign <- sign(df$weight)
  
  # Create an adjacency matrix from the 'from', 'to', and 'weight' columns
  SES.mat <- make_matrix(from = df$from, to = df$to, weight = df$weight)
  
  # If graph plotting is requested, plot the SES network
  if (graph) {
    SES.graph <- plot.SES(
      SES.mat,
      folder = folder,
      filename = graph.name,
      title = graph.title,
      label.cex = 1.1,
      vsize = 30
    )
  }
  
  # Return the adjacency matrix
  return(SES.mat)
}


############ Qualitative (signed) analyses ############ 

# Function to compute the Laplacian of the SES matrix
SES.laplacian <- function(SES.mat, from = c("rows", "cols")) {
  # Check if we are using row sums or column sums for the Laplacian
  if (from[1] == "rows") {
    SES.Lap <- diag(rowSums(t(SES.mat))) - t(SES.mat) # Laplacian using row sums
  } else {
    SES.Lap <- diag(rowSums((SES.mat))) - (SES.mat) # Laplacian using column sums
  }

  # Calculate eigenvalues of the Laplacian matrix
  lap.e <- eigen(SES.Lap)$value
  names(lap.e) <- rownames(SES.mat) # Name the eigenvalues

  return(lap.e) # Return eigenvalues
}

# Function to create a boolean file from the SES matrix
boolean.file.creation <- function(SES.mat, folder, filename) {
  
  # Convert the SES matrix to a binary matrix based on sign
  SES.bin <- sign(SES.mat)

  # Clean column and row names by removing special characters
  colnames(SES.bin) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(SES.bin))
  rownames(SES.bin) <- gsub("\\s*\\([^\\)]+\\)", "", rownames(SES.bin))
  colnames(SES.bin) <- gsub("[^[:alnum:]]", "", colnames(SES.bin))
  rownames(SES.bin) <- gsub("[^[:alnum:]]", "", rownames(SES.bin))
  colnames(SES.bin) <- gsub(" ", "_", colnames(SES.bin))
  rownames(SES.bin) <- gsub(" ", "_", rownames(SES.bin))

  # Initialize a data frame for the boolean factors
  boolean.df <- data.frame(targets = factor(colnames(SES.bin)), factors = NA)

  # Construct boolean factors for each target
  for (i in 1:ncol(SES.bin)) {
    poss <- names(which(SES.bin[, i] == 1)) # Positive influences
    negs <- names(which(SES.bin[, i] == -1)) # Negative influences
    if (length(negs) > 0) {
      negs <- paste0("!", negs) # Prefix negative influences with '!'
    }
    all <- c(poss, negs) # Combine positive and negative influences
    boolean.df$factors[i] <- paste(all, collapse = "|") # Create a combined factor string
  }

  # Write the boolean data frame to a CSV file
  write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = F, quote = FALSE)

  return(boolean.df)
}

# Function to analyze boolean networks
boolean.analyses <- function(boolean.net, folder, filename) {
  # Retrieve attractors of the boolean network
  states <- getAttractors(boolean.net, returnTable = TRUE)

  # Graphing begins
  state.map <- plotStateGraph(states, layout = layout.fruchterman.reingold, plotIt = FALSE)

  # Extract vertex attributes for analysis
  V(state.map)$degree <- degree(state.map) # Calculate the degree of each vertex
  V(state.map)$attractor <- "no" # Initialize attractor status
  V(state.map)$attractor[V(state.map)$frame.color == "#000000FF"] <- "yes" # Mark attractors in black

  # Create a subgraph excluding leaves
  core.state <- subgraph(state.map, which(V(state.map)$degree > 1))

  # Save the core state graph as a GraphML file
  write_graph(core.state, file = paste0(folder, filename, ".graphml"), format = "graphml")
  # Graphing ends

  # Calculate the number of states and attractors
  n_states <- length(states$stateInfo$table)
  n_attractors <- length(states$attractors)

  # Initialize lists for attractors and basin attraction sizes
  attractors <- vector("list", n_attractors)
  basin_attraction <- array(0, n_attractors)

  # Gather attractor sequences and basin sizes
  for (i in 1:n_attractors) {
    attractors[[i]] <- as.data.frame(getAttractorSequence(states, i))
    basin_attraction[i] <- states$attractors[[i]]$basinSize
  }

  # Return a summary list containing states, attractors, and basin sizes
  return(list(states = n_states, N_attractors = n_attractors, basins = basin_attraction, attractors = attractors))
}


# Integrated function to perform qualitative analyses on an SES matrix
qualitative.analyses <- function(SES.mat, folder, filename.boolean.csv, filename.boolean.graph) {
  
  # Calculate the Laplacian of the SES matrix using row sums
  laplacian <- SES.laplacian(SES.mat, from = "rows") 
  
  # Create a boolean file from the SES matrix and save it as a CSV file
  # Note: This step is necessary to use the BoolNet package for boolean network analysis
  boolean.df <- boolean.file.creation(SES.mat, folder, filename.boolean.csv)
  
  # Load the boolean network from the saved CSV file
  boolean.network <- loadNetwork(paste0(folder, filename.boolean.csv, ".csv"))
  
  # Perform boolean analyses on the loaded boolean network
  boolean.results <- boolean.analyses(boolean.network, folder, filename.boolean.graph)
  
  # Return a list containing the Laplacian, the boolean network, and the boolean analysis results
  return(list(laplacian = laplacian, boolean.network = boolean.network, boolean.results = boolean.results))
}

############ Quantitative analyses ############ 

# Function to simulate SES dynamics
SES.simulate <- function(SES.mat, iter, save.fig = FALSE, folder, fig.filename, fig.title) {
  
  require(Polychrome)
  require(ggplot2)
  require(ggfortify)
  require(ggpubr)
  require(grid)
  
  # Initialize a simulation matrix
  SES.sim <- matrix(NA, nrow(SES.mat), iter)
  SES.sim[, 1] <- runif(nrow(SES.mat), 0, 1) # Random initial values for the first iteration

  # Perform the simulation over the specified number of iterations
  for (i in 2:iter) {
    SES.sim[, i] <- t(SES.mat) %*% matrix((SES.sim[, i - 1]), ncol = 1) # Update values based on the SES matrix
  }

  rownames(SES.sim) <- rownames(SES.mat) # Set row names of the simulation matrix
  
  if(save.fig){
  sim.melt <- reshape2::melt(SES.sim) # Melt the simulation matrix for ggplot

  # Define colors for the plot
  fillcol <- glasbey.colors(nrow(SES.mat))
  names(fillcol) <- unique(sim.melt$Var1)

  # Create a line plot of the simulation results
  p1a <- ggplot(sim.melt, aes(x = log10(Var2), y = value, colour = Var1)) +
    geom_path(linewidth = 2) +
    labs(colour = "Element", y = "Progress", x = "log(time)") +
    scale_colour_manual(values = fillcol) +
    theme_minimal(base_size = 18) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  # Perform PCA on the simulation data
  pca <- prcomp(t(SES.sim), scale = FALSE) # Transpose for PCA analysis

  # Create a PCA plot
  p1b <- autoplot(pca, label = FALSE, loadings.colour = "black", colour = "blue", loadings.label = TRUE, loadings.label.size = 6, loadings.label.colour = "black") +
    geom_path(aes(x = PC1, y = PC2), colour = "blue", arrow = arrow(type = "closed", angle = 30, length = unit(5, "mm"))) +
    labs(x = "PC1", y = "PC2") +
    theme_minimal()

  # Create a text grob for the title on the side
  text1 <- textGrob(fig.title, just = "centre", rot = 90, gp = gpar(fontsize = 20))
  
  plot <- ggarrange(text1, p1a, p1b, nrow = 1, ncol = 3, widths = c(.1, 1, 1))

  # Save the combined plots as a PNG file
  png(paste0(folder, fig.filename, ".png"), width = 60, height = 35, units = "cm", res = 200)
  print(plot) # Doesn't work with ggarrange directly into png
  dev.off() # Close the PNG device

  }
  
  return(SES.sim) # Return the simulation results
}

# Function to calculate the Jacobian matrix
jacobian <- function(SES.mat) {
  SES.j <- t(SES.mat) # Transpose the SES matrix to compute Jacobian
  return(SES.j) # Return the Jacobian matrix
}

# Function to calculate the left eigenvectors of a matrix
lefteigen <- function(SES.mat) {
  left <- eigen(t(SES.mat))$vector # Compute left eigenvectors using the transpose
  return(left) # Return left eigenvectors
}

# Function to calculate the right eigenvectors of a matrix
righteigen <- function(SES.mat) {
  right <- eigen(SES.mat)$vector # Compute right eigenvectors
  return(right) # Return right eigenvectors
}

# Function to calculate the participation ratio
participation_ratio <- function(SES.mat, folder, filename, title) {
  # Ensure necessary libraries are loaded
  require(Polychrome)
  require(ggplot2)

  # Create the Jacobian matrix of the SES
  SES.j <- jacobian(SES.mat)

  # Calculate the left and right eigenvectors of the Jacobian matrix
  LJ <- lefteigen(SES.j)
  RJ <- righteigen(SES.j)

  # Calculate the Participation Ratio (PR)
  PR <- rowSums(RJ * LJ)^2 / rowSums((RJ * LJ)^2)

  # Create a data frame for the results
  PR.df <- data.frame(group = title, components = colnames(SES.mat), PR = PR)

  # Generate a color palette for the plot
  fillcol <- glasbey.colors(nrow(SES.mat))
  names(fillcol) <- colnames(SES.mat)

  # Create a bar plot for the participation ratios
  plot1 <- ggplot(PR.df, aes(y = Re(PR), x = components, fill = components)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fillcol) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none") +
    xlab(title) +
    ylab("Participation ratio") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Save the plot as a PNG file
  ggsave(paste0(folder, filename, ".png"), plot = plot1, device = "png", width = 60, height = 30, units = "cm", dpi = 200, bg = "white")

  return(PR.df) # Return the participation ratio data frame
}

# Function to simulate a matrix with random values based on the input matrix
# Note: 'change' should be a number indicating the element, ranging between 1 and prod(dim(mat))
# Note: The analyses using Automated Bayesian Computing (ABC) are still under development

simulate.mat <- function(SES.mat, all = TRUE, change, type = c('uniform', 'ordinal')) {
  
  # Define possible values for ordinal simulation
  draw <- c(0, 0.25, 0.5, 1) # Corresponds to the strong, medium, and weak link definitions in the SES
  
  if (all == TRUE) {
    # Simulate all elements of the matrix
    if (type == 'uniform') {
      # Generate a matrix of random numbers uniformly distributed between 0 and 1
      newmat <- matrix(runif(prod(dim(SES.mat)), 0, 1), nrow(SES.mat), ncol(SES.mat))
    } else if (type == 'ordinal') {
      # Generate a matrix of random values from the defined ordinal set
      newmat <- matrix(sample(draw, prod(dim(SES.mat)), replace = TRUE), nrow(SES.mat), ncol(SES.mat))
    }
    
    # Create the simulated matrix by preserving the structure of the input matrix
    mat.sim <- 1 * (SES.mat != 0) * sign(SES.mat) * newmat
    
  } else {
    # Simulate only specific elements of the matrix
    mat.sim <- SES.mat
    
    if (type == 'uniform') {
      # Generate random numbers uniformly distributed between 0 and 1 for specified elements
      mat.sim[as.matrix(change)] <- sign(SES.mat[as.matrix(change)]) * runif(length(SES.mat[as.matrix(change)]), 0, 1)
    } else if (type == 'ordinal') {
      # Generate random values from the defined ordinal set for specified elements
      mat.sim[as.matrix(change)] <- sign(SES.mat[as.matrix(change)]) * sample(draw, length(SES.mat[as.matrix(change)]), replace = TRUE)
    }
  }
  
  # Return the simulated matrix
  return(mat.sim)
}


# Function to simulate effects on the matrix based on specified measures
simulate.measure <- function(SES.mat, measure, affected, indicators, lower, upper) {
  # Generate random effects for the affected variables
  measure.ef <- runif(length(affected), lower, upper)
  
  # Initialize a new row for the matrix
  newrow <- matrix(rep(0, ncol(SES.mat)), nrow = 1)
  newrow[match(affected, colnames(SES.mat))] <- measure.ef # Assign effects to affected variables
  
  # Add the new row to the original matrix
  SES.mat <- rbind(SES.mat, newrow)
  
  # Generate random effects for the indicators
  measure.af <- runif(length(indicators), lower, upper) 
  
  # Initialize a new column for the matrix
  newcol <- matrix(rep(0, nrow(SES.mat)), ncol = 1)
  newcol[match(indicators, row.names(SES.mat))] <- measure.af # Assign effects to indicators
  
  # Add the new column to the matrix
  SES.mat <- cbind(SES.mat, newcol)
  
  # Name the new row and column appropriately
  row.names(SES.mat)[nrow(SES.mat)] <- colnames(SES.mat)[ncol(SES.mat)] <- measure
  
  return(SES.mat) # Return the modified matrix
}


# Function to observe state shifts in the SES dynamics
state.shift <- function(SES.mat, greed, iter, type = c('uniform', 'ordinal'),
                        all = TRUE, change = NULL,
                        folder, file) {
  
  require(reshape2) # Load reshape2 for data manipulation
  
  newnames <- row.names(SES.mat) # Get the row names of the input matrix
  
  # Simulate the existing matrix
  SES.obs <- SES.simulate(SES.mat = SES.mat, iter, save.fig = FALSE, folder = NULL, fig.filename = NULL, fig.title = NULL)
  state.obs <- SES.obs[, iter]
  names(state.obs) <- newnames
  
  state.sim <- array(NA, dim = c(length(state.obs), greed))
  rownames(state.sim) <- newnames
  
  mat.m <- melt(SES.mat)
  mat.sim.df <- array(NA, dim = c(nrow(mat.m), greed + 2))
  mat.sim.df[, 1] <- mat.m$Var1
  mat.sim.df[, 2] <- mat.m$Var2 
  
  tic <- Sys.time() # Start the timer
  
  # Perform simulations based on the greed parameter
  for (i in 1:greed) {
    
    # Change the weights of all links (uniform or ordinal) or change selected links
    mat.sim <- simulate.mat(SES.mat = SES.mat, type = type, all, change) # Simulate a new matrix
    
    SES <- SES.simulate(SES.mat = mat.sim, iter, save.fig = FALSE, folder = NULL, fig.filename = NULL, fig.title = NULL)
    
    mat.sim.df[, i + 2] <- melt(mat.sim)$value # Store the simulated matrix values
    state.sim[, i] <- apply(sign(SES[, (iter - 100):iter]), 1, prod) # Capture behavior over the last 101 steps
    
    if (i %in% seq(100000, 900000, by = 100000)) { 
      print(i) 
      flush.console()
    }
  }
  
  toc <- Sys.time() - tic # Record the elapsed time
  
  # Save the results
  save(state.sim, mat.sim.df, mat.m, SES.mat, file = paste(folder, file, sep = "/")) 
  
  # Return simulation results
  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, SES.mat = SES.mat)) 
}


# Function to prepare data for random forest analysis using the greedy approach
RF.prep <- function(state.shift.res, tolerance = 0.000001, targets, greed) {
  
  # Extract the state simulation results
  bin <- state.shift.res$state.sim
  
  # Apply a tolerance threshold to the simulation results
  bin[abs(bin) < tolerance] <- 0
  bin <- sign(bin)
  
  # Identify desirable outcomes based on the specified targets
  desirable.outcomes <- which(colSums(bin[match(targets, row.names(bin)), ]) == length(targets))
  
  # Extract the relevant part of the simulation data frame
  sim.ana <- state.shift.res$mat.sim.df[which(state.shift.res$mat.m$value != 0), ]
  
  # Transpose the simulation data frame and set column names
  sim.ana.df <- t(sim.ana[, 3:(greed + 2)])
  colnames(sim.ana.df) <- apply(sim.ana[, 1:2], 1, function(x) paste0(row.names(state.shift.res$SES.mat)[x[1]], " to ", row.names(state.shift.res$SES.mat)[x[2]]))
  
  # Convert the transposed data frame to a standard data frame
  sim.ana.df <- as.data.frame(sim.ana.df)
  
  # Initialize the outcomes column
  sim.ana.df$outcomes <- 0
  
  # Mark the desirable outcomes
  sim.ana.df$outcomes[desirable.outcomes] <- 1
  
  # Return the prepared data frame
  return(sim.ana.df)
}

# Function to perform random forest analysis
random.forest <- function(sim.outcomes, ntree = 500, folder, file.RF) {
  require(randomForest)
  require(randomForestExplainer)
  require(reprtree)
  require(cluster)
  
  # Prepare data for random forest analysis
  sim.outcomes$outcomes <- factor(sim.outcomes$outcomes) # Convert outcomes to factor
  
  # Clean column names for compatibility
  colnames(sim.outcomes) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(sim.outcomes)) # Remove parentheses and their contents
  colnames(sim.outcomes) <- gsub(" ", "_", colnames(sim.outcomes)) # Replace spaces with underscores
  colnames(sim.outcomes) <- gsub("\\.", "", colnames(sim.outcomes)) # Remove periods
  colnames(sim.outcomes) <- gsub("\\-", "", colnames(sim.outcomes)) # Remove hyphens
  
  # Train the random forest model
  forest <- randomForest(outcomes ~ ., data = sim.outcomes, ntree = ntree, localImp = TRUE, proximity = FALSE)
  
  # Save the trained forest model
  save(forest, file = paste(folder, file.RF, sep = "/"))
  
  # Return the trained random forest model
  return(forest)
}

# Function to analyze and visualize random forest results
random.forest.res <- function(input, folder, filename1.RF, filename2.RF) {
  
  require(ggplot2) # Load ggplot2 for plotting
  
  # Measure the importance of each feature in the random forest model
  importance_frame <- measure_importance(input)
  
  # Save the importance frame to a file
  save(importance_frame, file = paste0(folder, filename1.RF))
  
  # Create a multi-way importance plot
  pa <- plot_multi_way_importance(importance_frame, x_measure = "gini_decrease", y_measure = "accuracy_decrease", size_measure = "p_value", main = "Tuscany (CLD)") +
    labs(x = "mean decrease of GINI coefficient", y = "mean decrease in accuracy")
  
  # Save the plot as a PNG file
  ggsave(paste0(folder, filename2.RF, ".png"), plot = pa, device = "png", width = 40, height = 40, units = "cm", dpi = 200, bg = "white")
  
  # Return the importance frame
  return(importance_frame)
}

# Function that combines all quantitative anaylses
quantitative.analyses <- function(SES.mat, greed = 1000, iter = 500, 
                                  indicators,
                                  folder, 
                                  filename.simu.figure, 
                                  title.simu,
                                  filename.PR,
                                  filename.greedy.res, 
                                  file.RF, 
                                  filename1.RF, 
                                  filename2.RF){
  
  sim.res <- SES.simulate(SES.mat = SES.mat, iter = iter, save.fig = T, folder, fig.filename = filename.simu.figure, 
                          fig.title = title.simu)
  
  PR <- participation_ratio(SES.mat = SES.mat, folder = folder, filename = filename.PR , title = title.simu) 
  
  state.shift.res <- state.shift(SES.mat = SES.mat, greed, iter, type = 'uniform', folder = folder, file = filename.greedy.res)
  
  sim.outcomes.res <- RF.prep(state.shift.res = state.shift.res, targets = indicators, greed = greed)
  
  forest <- random.forest(sim.outcomes = sim.outcomes.res, 
                          ntree = 500, folder = folder, 
                          file.RF = file.RF)
  
  gc()
  importance_frame <- random.forest.res(input = forest, folder = folder, 
                                        filename1.RF = filename1.RF, 
                                        filename2.RF = filename2.RF)
  
  return(list(sim.res = sim.res, PR = PR, state.shift.res = state.shift.res, sim.outcomes = sim.outcomes, forest = forest, importance_frame = importance_frame))
}


############ Introduce measures ############ 


# Function to observe state shifts with measures in the SES dynamics
state.shift.measure <- function(measure, affected, SES.mat, greed, iter, indicators, 
                                folder, file.measure) {
  
  require(reshape2) # Load reshape2 for data manipulation
  
  # Combine row names of the matrix with the measure
  newnames <- c(row.names(SES.mat), measure)
  
  # Create a data frame for the new measure
  focus.mat.df <- data.frame(from = factor(c(rep(measure, length(newnames)), newnames)), to = factor(c(newnames, rep(measure, length(newnames)))))
  
  # Initialize an array to store the state simulation results
  state.sim <- array(NA, dim = c(nrow(SES.mat) + 1, greed))
  rownames(state.sim) <- c(row.names(SES.mat), measure)
  
  # Initialize an array to store the simulated matrix values
  mat.sim.df <- array(NA, dim = c(nrow(focus.mat.df), greed + 2))
  mat.sim.df[, 1] <- focus.mat.df$from
  mat.sim.df[, 2] <- focus.mat.df$to # Use factor levels for faster array operations
  
  tic <- Sys.time() # Start the timer
  
  # Perform simulations based on the greed parameter
  for (i in 1:greed) {
    # Simulate the measure's effect on the matrix
    mat.sim <- simulate.measure(SES.mat, measure, affected, indicators, lower = -1, upper = 0)
    
    # Simulate the SES dynamics
    # SES <- time.simulate(mat.sim, iter, starting.value)
    
    SES <- SES.simulate(SES.mat = mat.sim, iter, save.fig = FALSE, folder = NULL, fig.filename = NULL, fig.title = NULL)
    
    # Store the simulated matrix values
    mat.sim.df[, i + 2] <- c(mat.sim[(nrow(mat) + 1):nrow(mat.sim), ], mat.sim[, (ncol(mat) + 1):ncol(mat.sim)]) # Account for multiple measures
    
    # Capture behavior over the last 101 steps
    state.sim[, i] <- apply(sign(SES[, (iter - 100):iter]), 1, prod) # Use an odd number of steps to capture behavior
  }
  
  toc <- Sys.time() - tic # Record the elapsed time
  toc
  
  # Save the results
  save(state.sim, mat.sim.df, mat.m, mat, file = paste(folder, file.measure, sep = "/"))
  
  # Return simulation results
  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, mat = mat.m))
}


###########################################################################################################################





#############################################
# With EasyABC 

state.shift.ABC <- function(x, iter = 500, desirable.outcomes = indicators) {
  
  # mat.sim <- simulate.mat(mat, type = type, all, change, ABC = TRUE, prior = my_prior) # Simulate a new matrix -- this is where the prior comes into play
  
  # Classify the drawn values into weak-strong values
  b <- c(0, 0.25, 0.5, 1)
  
  x1 <- sapply(x, function(x, b) {b[which.min(abs(x-b))]}, b)

  mat.sim <- 1 * (mat != 0) * sign(mat) * matrix(x1, nrow(mat), ncol(mat))
  
    # Initialize a simulation matrix
  SES.sim <- matrix(NA, nrow(mat.sim), iter)
  SES.sim[, 1] <- runif(nrow(mat.sim), 0, 1) # Random initial values for the first iteration
  
  # Perform the simulation over the specified number of iterations
  for (i in 2:500) {
    SES.sim[, i] <- t(mat.sim) %*% matrix((SES.sim[, i - 1]), ncol = 1) # Update values based on the SES matrix
  }
  
  # This is not necessary I think
  rownames(SES.sim) <- rownames(mat) # Set row names of the simulation matrix
  
  # SES <- SES.simulate(SES.mat = mat.sim, iter, save.fig = F, folder = NULL, fig.filename = NULL, fig.title = NULL)
  
  res <- apply(sign(SES.sim[, (iter - 100):iter]), 1, prod) # Capture behavior over the last 101 steps
  
  # Select here the desirable outcomes
  res <- res[desirable.outcomes]
  
  return(res) # Return simulation results
}


##############################################

