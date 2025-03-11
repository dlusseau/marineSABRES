########################################################################
####                                                        ############
#### Helper functions marineSABRES T5.3 network analyses   ############
####                                                        ############
########################################################################
####
####   davlu
####   bmjv
####   October 2024
####   list of utility functions to walk through interventions SABRES analytical pipelines
####   See Demo.R for a demo workflow
####
########################################################################

# TODO:

# Make state.shift funciton so that it can deal with both the 'simple' greedy approacha nd with measures

# Think of what to save in state shift


############ Intro ############

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



data.load <- function(df, folder, graph = FALSE, graph.name = NULL, graph.title = NULL){
  
  df$weight <- 1
  df$weight[df$strength == "Medium Positive"] <- .5
  df$weight[df$strength == "Medium Negative"] <- (-.5)
  df$weight[df$strength == "Strong Negative"] <- (-1)
  df$sign <- sign(df$weight)
  
  SES.mat <- make_matrix(from = df$from, to = df$to, weight = df$weight)
  
  if(graph){
    SES.graph <- plot.SES(SES.mat,
                          folder = folder,
                          filename = graph.name,
                          title = graph.title,
                          label.cex = 1.1,
                          vsize = 30
    )
  }
  
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


qualitative.analyses <- function(SES.mat, folder, filename.boolean.csv, filename.boolean.graph) {
  
  laplacian <- SES.laplacian(SES.mat, from = "rows") 
  # The only way for now is by saving the network as a csv and then re-open it with BoolNet - I went down a rabbit hole to mirror this function with a df as input but no success  
  boolean.df <- boolean.file.creation(SES.mat, folder, filename.boolean.csv)
  boolean.network <- loadNetwork(paste0(folder, filename.boolean.csv, ".csv"))
  
  boolean.results <- boolean.analyses(boolean.network, folder, filename.boolean.graph)
  
  return(list(laplacian = laplacian,boolean.network = boolean.network, boolean.results = boolean.results))
  
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
jacobian <- function(SES) {
  SES.j <- t(SES) # Transpose the SES matrix to compute Jacobian
  return(SES.j) # Return the Jacobian matrix
}

# Function to calculate the left eigenvectors of a matrix
lefteigen <- function(mat) {
  left <- eigen(t(mat))$vector # Compute left eigenvectors using the transpose
  return(left) # Return left eigenvectors
}

# Function to calculate the right eigenvectors of a matrix
righteigen <- function(mat) {
  right <- eigen(mat)$vector # Compute right eigenvectors
  return(right) # Return right eigenvectors
}

# Function to calculate the participation ratio
participation_ratio <- function(SES, folder, filename, title) {
  # Ensure necessary libraries are loaded
  require(Polychrome)
  require(ggplot2)

  # Create the Jacobian matrix of the SES
  SES.j <- jacobian(SES)

  # Calculate the left and right eigenvectors of the Jacobian matrix
  LJ <- lefteigen(SES.j)
  RJ <- righteigen(SES.j)

  # Calculate the Participation Ratio (PR)
  PR <- rowSums(RJ * LJ)^2 / rowSums((RJ * LJ)^2)

  # Create a data frame for the results
  PR.df <- data.frame(group = title, components = colnames(SES), PR = PR)

  # Generate a color palette for the plot
  fillcol <- glasbey.colors(nrow(SES))
  names(fillcol) <- colnames(SES)

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

# Change should be a number indicating the element, ranging between 1 and prod(dim(mat))

simulate.mat <- function(mat, all = TRUE, change, type = c('uniform','ordinal'), ABC = FALSE, prior = NULL) {
  #   # Create a matrix to store simulated values based on the input matrix
  #   # Step-by-step explanation of the operations:
  #   # 1. (mat != 0): Creates a boolean matrix where TRUE indicates non-zero elements in 'mat'.
  #   # 2. 1 * (mat != 0): Converts the boolean matrix to a numeric matrix of 1s and 0s.
  #   # 3. sign(mat): Returns the sign of each element in 'mat' (-1, 1, or 0).
  #   # 4. matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat)): Creates a matrix of random numbers uniformly distributed between 0 and 1.
  #   # 5. The final product yields a matrix where:
  #   #    - Elements corresponding to zeros in 'mat' remain zero.
  #   #    - Non-zero elements are multiplied by random numbers between 0 and 1, preserving their sign.
  #   
  draw <- c(0, 0.25, 0.5, 1) # Corresponds to the strong, medium and weak link definition in the SES
  
  # if(ABC == TRUE){
  #   # For now ABC can only be performed on the full matrix and with a uniform distribution
  #   
  #   mat.sim <- 1 * (mat != 0) * sign(mat) * matrix(prior, nrow(mat), ncol(mat))
  #   
  # } else {
    
    if(all == TRUE){
      
      if(type == 'uniform'){
        newmat <- matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat))
      } else if (type == 'ordinal'){
        newmat <- matrix(sample(draw, prod(dim(mat)), replace = TRUE), nrow(mat), ncol(mat))
      }
      
      mat.sim <- 1 * (mat != 0) * sign(mat) * newmat
      
    } else {
      
      mat.sim <- mat
      
      if(type == 'uniform'){
        mat.sim[as.matrix(change)] <- sign(mat[as.matrix(change)]) * runif(length(mat[as.matrix(change)]), 0, 1)
        
      } else if (type == 'ordinal'){
        
        mat.sim[as.matrix(change)] <- sign(mat[as.matrix(change)]) * sample(draw,length(mat[as.matrix(change)]),replace=T)
        
      }
    }
  # }
  
  return(mat.sim) # Return the simulated matrix
}

# Function to simulate effects on the matrix based on specified measures
simulate.measure <- function(mat, measure, affected, indicators, lower, upper) {
  # Generate random effects for the affected variables
  measure.ef <- runif(length(affected), lower, upper)
  
  # Initialize a new row for the matrix
  newrow <- matrix(rep(0, ncol(mat)), nrow = 1)
  newrow[match(affected, colnames(mat))] <- measure.ef # Assign effects to affected variables
  
  # Add the new row to the original matrix
  mat <- rbind(mat, newrow)
  
  # Generate random effects for the indicators
  measure.af <- runif(length(indicators), lower, upper) 
  
  # Initialize a new column for the matrix
  newcol <- matrix(rep(0, nrow(mat)), ncol = 1)
  newcol[match(indicators, row.names(mat))] <- measure.af # Assign effects to indicators
  
  # Add the new column to the matrix
  mat <- cbind(mat, newcol)
  
  # Name the new row and column appropriately
  row.names(mat)[nrow(mat)] <- colnames(mat)[ncol(mat)] <- measure
  
  return(mat) # Return the modified matrix
}


# Function to observe state shifts in the SES dynamics
state.shift <- function(mat, greed, iter, type = c('uniform','ordinal'),
                        all = TRUE, change = NULL,
                        measure = FALSE, affected = NULL, indicators = NULL, lower = -1, upper = 0,
                        folder, file) {
  
  require(reshape2) # Load reshape2 for data manipulation
  
  newnames <- row.names(mat)
  
  if(measure){
    newnames <- c(newnames, measure) 
    focus.mat.df <- data.frame(from = factor(c(rep(measure, length(newnames)), newnames)), to = factor(c(newnames, rep(measure, length(newnames)))))
    
    state.sim <- array(NA, dim = c(nrow(mat) + 1, greed))
    rownames(state.sim) <- c(row.names(mat), measure)
    
    mat.sim.df <- array(NA, dim = c(nrow(focus.mat.df), greed + 2))
    mat.sim.df[, 1] <- focus.mat.df$from
    mat.sim.df[, 2] <- focus.mat.df$to 
    
  } else {
    
    SES.obs <- SES.simulate(SES.mat = mat, iter, save.fig = F, folder = NULL, fig.filename = NULL, fig.title = NULL)
    state.obs <- SES.obs[, iter]
    names(state.obs) <- newnames
    
    state.sim <- array(NA, dim = c(length(state.obs),greed))
    rownames(state.sim) <- newnames
    
    mat.m <- melt(mat)
    mat.sim.df <- array(NA, dim = c(nrow(mat.m), greed + 2))
    mat.sim.df[,1] <- mat.m$Var1
    mat.sim.df[,2] <- mat.m$Var2 
  }
  
  tic <- Sys.time() # Start the timer
  
  # Perform simulations based on the greed parameter
  for (i in 1:greed) {
    
    # Either simulate a new measure (new link), change the weights of all links (uniform or ordinal) or change selected links
    if(measure){
      mat.sim <- simulate.measure(mat, measure, affected, indicators, lower = -1, upper = 0)
    } else {
      mat.sim <- simulate.mat(mat, type = type, all, change) # Simulate a new matrix
    }
    
    SES <- SES.simulate(SES.mat = mat.sim, iter, save.fig = F, folder = NULL, fig.filename = NULL, fig.title = NULL)
    
    mat.sim.df[, i + 2] <- melt(mat.sim)$value # Store the simulated matrix values
    state.sim[, i] <- apply(sign(SES[, (iter - 100):iter]), 1, prod) # Capture behavior over the last 101 steps
    
    if (i %in% c(100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000)) { 
      print(i) 
      flush.console()
    }
  }
  
  toc <- Sys.time() - tic # Record the elapsed time
  
  # To test
  save(state.sim, mat.sim.df, mat.m, mat, file = paste(folder, file, sep = "/")) # Save the results
  
  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, mat = mat)) # Return simulation results
}

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

# Function when using greedy approach
RF.prep <- function(state.shift.res, tolerance = 0.000001, targets, greed = 100){
  
  bin <- state.shift.res$state.sim
  bin[abs(bin) < tolerance] <- 0
  bin <- sign(bin)
  
  desirable.outcomes <- which(colSums(bin[match(targets, row.names(bin)), ]) == length(targets))
  
  sim.ana <- state.shift.res$mat.sim.df[which(state.shift.res$mat.m$value != 0),]
  sim.ana.df <- t(sim.ana[, 3:(greed + 2)])
  colnames(sim.ana.df) <- apply(sim.ana[, 1:2], 1, function(x) paste0(row.names(state.shift.res$mat)[x[1]], " to ", row.names(state.shift.res$mat)[x[2]]))
  sim.ana.df <- as.data.frame(sim.ana.df)
  
  sim.ana.df$outcomes <- 0
  sim.ana.df$outcomes[desirable.outcomes] <- 1
  
  return(sim.ana.df)
}

# Function to perform random forest analysis
random.forest <- function(sim.outcomes, ntree = 500, folder, file) {
  require(randomForest)
  require(randomForestExplainer)
  require(reprtree)
  require(cluster)

  # Prepare data for random forest analysis
  sim.outcomes$outcomes <- factor(sim.outcomes$outcomes) # Convert outcomes to factor
  # Clean column names for compatibility
  colnames(sim.outcomes) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(sim.outcomes))
  colnames(sim.outcomes) <- gsub(" ", "_", colnames(sim.outcomes))
  colnames(sim.outcomes) <- gsub("\\.", "", colnames(sim.outcomes))
  colnames(sim.outcomes) <- gsub("\\-", "", colnames(sim.outcomes))

  # Train the random forest model
  forest <- randomForest(outcomes ~ ., data = sim.outcomes, ntree = ntree, localImp = TRUE, proximity = FALSE)

  # Save the trained forest model
  save(forest, file = paste(folder, file, sep = "/"))

  return(forest) # Return the trained random forest model
}

# Function to analyze and visualize random forest results
random.forest.res <- function(input, folder, filename1, filename2) {
  
  require(ggplot2)
  
  # Measure the importance of each feature in the random forest model
  importance_frame <- measure_importance(input)

  # Save the importance frame to a file
  save(importance_frame, file = paste0(folder, filename1))

  # Create a multi-way importance plot
  pa <- plot_multi_way_importance(importance_frame, x_measure = "gini_decrease", y_measure = "accuracy_decrease", size_measure = "p_value", main = "Tuscany (CLD)") +
    labs(x = "mean decrease of GINI coefficient", y = "mean decrease in accuracy")

  # Save the plot as a PNG file
  ggsave(paste0(folder, filename2, ".png"), plot = pa, device = "png", width = 40, height = 40, units = "cm", dpi = 200, bg = "white")

  return(importance_frame) # Return the importance frame
}




# From here not clean yet......




#  Use this for the greedy approach again


state.shift.measure <- function(measure, greed, iter, mat, starting.value, indicators) {
  require(reshape2) # Load reshape2 for data manipulation

  newnames <- c(row.names(mat), measure)
  focus.mat.df <- data.frame(from = factor(c(rep(measure, length(newnames)), newnames)), to = factor(c(newnames, rep(measure, length(newnames)))))

  state.sim <- array(NA, dim = c(nrow(mat) + 1, greed))
  rownames(state.sim) <- c(row.names(mat), measure)

  mat.sim.df <- array(NA, dim = c(nrow(focus.mat.df), greed + 2))
  mat.sim.df[, 1] <- focus.mat.df$from
  mat.sim.df[, 2] <- focus.mat.df$to # can't deal with character and numbers in array, that's ok we keep it to factor levels, array is faster than df


  tic <- Sys.time()

  for (i in 1:greed) {
    mat.sim <- simulate.measure(mat, measure, affected, indicators, lower = -1, upper = 0)


    SES <- time.simulate(mat.sim, iter, starting.value)

    mat.sim.df[, i + 2] <- c(mat.sim[(nrow(mat) + 1):nrow(mat.sim), ], mat.sim[, (ncol(mat) + 1):ncol(mat.sim)]) # like that we account for multiple measures

    state.sim[, i] <- apply(sign(SES[, (iter - 100):iter]), 1, prod) # here we capture behaviour over last 101 steps, importantly an odd number
  }

  # Looks at the signs of the last 100 elements, to see if they're all the same

  toc <- Sys.time() - tic
  toc

  # To test
  save(state.sim, mat.sim.df, mat.m, mat, file = paste(folder, file, sep = "/")) # Save the results

  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, mat = mat.m)) # Return simulation results
}


