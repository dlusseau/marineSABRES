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
# Test state shift functions, there are some irregularities where there's still references to Tuscany
# Run with Macaronesia as test
# Make state.shift funciton so that it can deal with both the 'simple' greedy approacha nd with measures


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

# Function to plot the SES network
plot.SES <- function(SES.mat, folder, filename, title, w = 80, h = 40, layout = layout_with_fr, label.cex = 1.5, vsize = 20, eweight = 10) {
  # Create a directed graph from the SES matrix
  SES.net <- graph_from_adjacency_matrix(SES.mat, mode = "directed", weighted = TRUE)

  # Set vertex properties
  V(SES.net)$color <- "white" # Set vertex color
  V(SES.net)$label.cex <- label.cex # Set label size
  V(SES.net)$size <- vsize # Set vertex size

  # Set edge properties
  E(SES.net)$sign <- sign(E(SES.net)$weight) # Determine the sign of edge weights
  E(SES.net)$color <- "blue" # Default edge color
  E(SES.net)$color[E(SES.net)$sign < 0] <- "red" # Change color for negative edges
  E(SES.net)$weight <- abs(E(SES.net)$weight) # Use absolute values for edge weights
  E(SES.net)$curved <- 0.2 # Set edge curvature
  E(SES.net)$width <- E(SES.net)$weight * eweight # Scale edge width by weight
  E(SES.net)$arrow.size <- E(SES.net)$weight * (eweight / 2) # Set arrow size proportional to weight

  # Define layout for the plot
  l <- layout(SES.net)

  # Save the plot as a PNG file
  png(paste0(folder, filename, ".png"), width = w, height = h, units = "cm", res = 200)
  plot(SES.net, layout = l, ylab = title) # Plot the network with defined layout and y-axis label
  dev.off() # Close the PNG device

  return(SES.net) # Return the constructed network
}

# Function to create a boolean file from the SES matrix
boolean_file_creation <- function(SES.mat, folder, filename) {
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
}

# Function to analyze boolean networks
boolean_analyses <- function(boolean.net, folder, filename) {
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

# Function to simulate SES dynamics
SES.simulate <- function(SES.mat, iter, folder, filename, title) {
  # Initialize a simulation matrix
  SES.sim <- matrix(NA, nrow(SES.mat), iter)
  SES.sim[, 1] <- runif(nrow(SES.mat), 0, 1) # Random initial values for the first iteration

  # Perform the simulation over the specified number of iterations
  for (i in 2:iter) {
    SES.sim[, i] <- t(SES.mat) %*% matrix((SES.sim[, i - 1]), ncol = 1) # Update values based on the SES matrix
  }

  rownames(SES.sim) <- rownames(SES.mat) # Set row names of the simulation matrix
  sim.melt <- melt(SES.sim) # Melt the simulation matrix for ggplot

  # Define colors for the plot
  fillcol <- glasbey.colors(nrow(SES.mat))
  names(fillcol) <- unique(sim.melt$Var1)

  # Create a line plot of the simulation results
  p1a <- ggplot(sim.melt, aes(x = log10(Var2), y = value, colour = Var1)) +
    geom_path(size = 2) +
    labs(colour = "Element", y = "Progress", x = "log(time)") +
    scale_colour_manual(values = fillcol) +
    theme_minimal(base_size = 18) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  # Perform PCA on the simulation data
  pca <- prcomp(t(SES.sim), scale = FALSE) # Transpose for PCA analysis

  # Create a PCA plot
  p1b <- ggplot2::autoplot(pca, label = FALSE, loadings.colour = "black", colour = "blue", loadings.label = TRUE, loadings.label.size = 6, loadings.label.colour = "black") +
    geom_path(aes(x = PC1, y = PC2), colour = "blue", arrow = arrow(type = "closed", angle = 30, length = unit(5, "mm"))) +
    labs(x = "PC1", y = "PC2") +
    theme_minimal()

  # Create a text grob for the title
  text1 <- textGrob(title, just = "centre", rot = 90, gp = gpar(fontsize = 20))

  # Save the combined plots as a PNG file
  png(paste0(folder, filename, ".png"), width = 60, height = 35, units = "cm", res = 200)
  ggarrange(text1, p1a, p1b, nrow = 1, ncol = 3, widths = c(.1, 1, 1))
  dev.off() # Close the PNG device

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
simulate.mat <- function(mat) {
  # Create a matrix to store simulated values based on the input matrix
  # Step-by-step explanation of the operations:
  # 1. (mat != 0): Creates a boolean matrix where TRUE indicates non-zero elements in 'mat'.
  # 2. 1 * (mat != 0): Converts the boolean matrix to a numeric matrix of 1s and 0s.
  # 3. sign(mat): Returns the sign of each element in 'mat' (-1, 1, or 0).
  # 4. matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat)): Creates a matrix of random numbers uniformly distributed between 0 and 1.
  # 5. The final product yields a matrix where:
  #    - Elements corresponding to zeros in 'mat' remain zero.
  #    - Non-zero elements are multiplied by random numbers between 0 and 1, preserving their sign.

  mat.sim <- 1 * (mat != 0) * sign(mat) * matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat))

  return(mat.sim) # Return the simulated matrix
}

# Function to simulate the dynamics of the SES over a specified number of iterations
time.simulate <- function(mat, iter, starting.value) {
  # Initialize a matrix to store the simulated values with specified iterations
  SES.sim <- matrix(NA, nrow(mat), iter)
  SES.sim[, 1] <- starting.value # Set the starting values

  # Loop through the number of iterations to compute the SES dynamics
  for (i in 2:iter) {
    # Update values based on the SES matrix
    SES.sim[, i] <- t(mat) %*% matrix((SES.sim[, i - 1]), ncol = 1)
  }

  return(SES.sim) # Return the simulated SES matrix
}

# Function to observe state shifts in the SES dynamics
state.shift <- function(greed, iter, mat, starting.value) {
  require(reshape2) # Load reshape2 for data manipulation

  # Simulate SES over specified iterations
  SES.obs <- time.simulate(mat, iter, starting.value)
  state.obs <- SES.obs[, iter] # Capture the final state
  names(state.obs) <- colnames(tuscany.SES) # Name the state observations

  # Initialize an array to store state simulations
  state.sim <- array(NA, dim = c(length(state.obs), greed))
  rownames(state.sim) <- colnames(tuscany.SES) # Set row names for states

  # Melt the original matrix for easier manipulation
  mat.m <- melt(mat)
  mat.sim.df <- array(NA, dim = c(nrow(mat.m), greed + 2))
  mat.sim.df[, 1] <- mat.m$Var1 # Store the first variable
  mat.sim.df[, 2] <- mat.m$Var2 # Store the second variable

  tic <- Sys.time() # Start the timer

  # Perform simulations based on the greed parameter
  for (i in 1:greed) {
    mat.sim <- simulate.mat(mat) # Simulate a new matrix

    SES <- time.simulate(mat.sim, iter, starting.value) # Simulate SES dynamics

    mat.sim.df[, i + 2] <- melt(mat.sim)$value # Store the simulated matrix values
    state.sim[, i] <- apply(sign(SES[, (iter - 100):iter]), 1, prod) # Capture behavior over the last 101 steps
  }

  toc <- Sys.time() - tic # Record the elapsed time

  # To test
  save(state.sim, mat.sim.df, mat.m, mat, file = paste(folder, file, sep = "/")) # Save the results

  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, mat = mat.m)) # Return simulation results
}

# Function to perform random forest analysis
random.forest <- function(data, ntree = 500, folder, file) {
  require(randomForest)
  require(randomForestExplainer)
  require(reprtree)
  require(cluster)

  # Prepare data for random forest analysis
  sim.ana.df$outcomes <- factor(sim.ana.df$outcomes) # Convert outcomes to factor
  # Clean column names for compatibility
  colnames(sim.ana.df) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(sim.ana.df))
  colnames(sim.ana.df) <- gsub(" ", "_", colnames(sim.ana.df))
  colnames(sim.ana.df) <- gsub("\\.", "", colnames(sim.ana.df))
  colnames(sim.ana.df) <- gsub("\\-", "", colnames(sim.ana.df))

  # Train the random forest model
  forest <- randomForest(outcomes ~ ., data = sim.ana.df, ntree = ntree, localImp = TRUE, proximity = FALSE)

  # Save the trained forest model
  save(forest, file = paste(folder, file, sep = "/"))

  return(forest) # Return the trained random forest model
}

# Function to analyze and visualize random forest results
random.forest.res <- function(input, folder, filename1, filename2) {
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


# From here not clean yet......




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
  toc <- Sys.time() - tic
  toc

  # To test
  save(state.sim, mat.sim.df, mat.m, mat, file = paste(folder, file, sep = "/")) # Save the results

  return(list(state.sim = state.sim, mat.sim.df = mat.sim.df, mat.m = mat.m, mat = mat.m)) # Return simulation results
}











##### Berthe's old code


FCM.mat <- function(FCM) {
  # Extract unique elements from the FCM data frame
  elements <- c(unique(FCM$ELEMENT))

  # Get the initial values for each unique element
  # Use the first occurrence of each element (non-duplicated)
  starting.value <- FCM$ELEMENT.VALUE[!duplicated(FCM$ELEMENT)]
  # Assign element names to these starting values
  names(starting.value) <- FCM$ELEMENT[!duplicated(FCM$ELEMENT)]
  # Transform the starting values to a column vector
  starting.value <- t(t(starting.value))

  # Initialize a square matrix with dimensions equal to the number of unique elements
  # Fill the matrix with zeros
  PESTLE.mat <- matrix(0, length(elements), length(elements))
  # Set the column names of the matrix to the unique elements
  colnames(PESTLE.mat) <- elements
  # Set the row names of the matrix to the unique elements
  row.names(PESTLE.mat) <- elements

  # Fill the matrix with influence values
  # Iterate over each row of the FCM data frame
  for (i in 1:nrow(FCM)) {
    # Find the row index corresponding to the current ELEMENT
    row_idx <- which(row.names(PESTLE.mat) == FCM$ELEMENT[i])
    # Find the column index corresponding to the current INFLUENCE
    col_idx <- which(colnames(PESTLE.mat) == FCM$INFLUENCE[i])
    # Assign the influence value to the corresponding cell in the matrix
    PESTLE.mat[row_idx, col_idx] <- FCM$INFLUENCE.VALUE[i]
  }

  # Return a list containing the influence matrix and the starting values vector
  return(list(PESTLE.mat = PESTLE.mat, starting.value = starting.value))
}


## Simulate matrix  #####################

simulate.mat <- function(mat) {
  # Create a matrix to store the simulated values
  # Step-by-step explanation of the matrix operations:
  # 1. (mat != 0): This creates a boolean matrix where each element is TRUE if the corresponding element in 'mat' is not zero, and FALSE otherwise.
  # 2. 1 * (mat != 0): This converts the boolean matrix to a numeric matrix where TRUE is converted to 1 and FALSE is converted to 0.
  # 3. sign(mat): This returns a matrix of the same dimensions as 'mat' where each element is the sign of the corresponding element in 'mat' (-1 for negative, 1 for positive, 0 for zero).
  # 4. matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat)): This creates a matrix of random numbers uniformly distributed between 0 and 1, with the same dimensions as 'mat'.
  # 5. The product of these matrices results in a matrix where:
  #    - Elements corresponding to zero elements in 'mat' are zero.
  #    - Non-zero elements are random numbers between 0 and 1, multiplied by the sign of the corresponding element in 'mat'.
  mat.sim <- 1 * (mat != 0) * sign(mat) * matrix(runif(prod(dim(mat)), 0, 1), nrow(mat), ncol(mat))

  # Return the simulated matrix
  return(mat.sim)
}

## Simulation  #####################

simu <- function(iter = 1000, starting.values, matrix.elems, elements = c("P", "E", "S", "T", "L", "EN")) {
  # Initialize a matrix to store the simulation results
  # Rows correspond to elements (e.g., P, E, S, T, L, EN)
  # Columns correspond to iterations
  FCM.sim <- matrix(NA, length(elements), iter) # Matrix filled with NAs initially

  # Set the initial values for the first iteration
  FCM.sim[, 1] <- starting.values # The first column is filled with the starting values

  # Loop over the number of iterations starting from the second iteration
  for (i in 2:iter) {
    # Each iteration multiplies the transpose of the PESTLE matrix (matrix.elems) with the previous outcome
    # FCM.sim[, i] <- PESTLE.mat %*% matrix((FCM.sim[, i - 1]), ncol = 1)
    FCM.sim[, i] <- t(matrix.elems) %*% matrix((FCM.sim[, i - 1]), ncol = 1)
  }

  # Return the simulation results matrix
  return(FCM.sim)
}

## Boolean analysis  #####################

## Transform to binary dataset ##################################################

bin.transform <- function(mat, folder, group) {
  # Create a binary matrix where elements are 1, -1, or 0
  # 1 indicates positive influence, -1 indicates negative influence, 0 indicates no influence
  PESTLE.bin <- sign(mat)

  # Initialize a data frame to store the boolean expressions for each target element
  boolean.df <- data.frame(targets = factor(colnames(PESTLE.bin)), factors = NA)

  # Loop through each column of the binary matrix
  for (i in 1:ncol(PESTLE.bin)) {
    # Get the names of elements with positive influence on the current target
    poss <- names(which(PESTLE.bin[, i] == 1))
    # Get the names of elements with negative influence on the current target
    negs <- names(which(PESTLE.bin[, i] == -1))

    # If there are negative influences, prepend '!' to each name to indicate negation
    if (length(negs) > 0) {
      negs <- paste0("!", negs)
    }

    # Combine positive and negative influences into a single vector
    all <- c(poss, negs)

    # Concatenate all influences with '|' to create a boolean expression
    boolean.df$factors[i] <- paste(all, collapse = "|")
  }

  # Create the filename for the CSV file
  filename <- paste0("PESTLE_bool_", group)
  # Write the data frame to a CSV file in the specified folder
  write.csv(boolean.df, file = paste0(folder, filename, ".csv"), row.names = FALSE, quote = FALSE)

  # Return the data frame with boolean expressions
  return(boolean.df)
}


## Load network and obtain states ##################################################

pestle_boolean1 <- loadNetwork(paste0(folder, filename, ".csv"))
states.pestle1 <- getAttractors(pestle_boolean1)

# Simple graph
plot.state.map(states = states.pestle1, group = group)

state.map <- plotStateGraph(states.pestle1, layout = layout.fruchterman.reingold, plotIt = FALSE)

# Write graph: final figure will be made in Gephi
write_graph(
  state.map,
  file = paste0(folder, "pestle1_boolean.graphml"),
  format = "graphml"
)

# Print states
trans.tab <- getTransitionTable(states)
# plot.state.graph(states)
print(getBasinOfAttraction(states, 1))


## Resilience  #####################

resilience.detracting.node.exp <- function(FCM.sim, elements = c("P", "E", "S", "T", "L", "EN"), tol = 10^-5, logged = TRUE) {
  # Load the 'reshape2' package, which is required for the function
  require(reshape2)

  # Determine the number of iterations (columns) in the simulation matrix
  iter <- ncol(FCM.sim)

  # Calculate the difference matrix based on whether logging is applied or not
  if (logged == TRUE) {
    # If logging is true, compute the logarithmic difference between consecutive iterations
    diff.mat <- log10(FCM.sim[, 2:iter]) - log10(FCM.sim[, 1:(iter - 1)])
  } else {
    # Otherwise, compute the simple difference between consecutive iterations
    diff.mat <- (FCM.sim[, 2:iter]) - (FCM.sim[, 1:(iter - 1)])
  }

  # Calculate the stable rate using the second-to-last iteration's difference
  stable.rate <- diff.mat[, (iter - 1)]
  # Compute the resilience matrix by subtracting the stable rate from the difference matrix
  res.mat <- diff.mat - stable.rate
  # Create a logical matrix indicating where the absolute values of res.mat are less than the tolerance
  res.mat.log <- (abs(res.mat) < tol)
  # Determine the first iteration where each element meets the tolerance criteria
  resilience.elements <- apply(res.mat.log, 1, function(x) which(x == TRUE)[1])
  # Assign names to the resilience elements based on the input elements
  names(resilience.elements) <- elements
  # Determine the system resilience as the maximum iteration where any element becomes stable
  system.resilience <- max(resilience.elements)

  # Return a list containing the resilience elements, the system resilience, and the state at the last iteration
  return(list(resilience.elements = resilience.elements, resilience = system.resilience, state = (FCM.sim[, iter])))
}


## Plot the network #####################

plot.network <- function(origin.mat = PESTLE.mat, group, save.plot = TRUE) {
  # Load the necessary libraries
  require(igraph)
  require(GGally)

  # Create a directed graph object from the adjacency matrix 'origin.mat'
  FCM1.net <- graph_from_adjacency_matrix(
    (origin.mat),
    mode = "directed",
    weighted = TRUE
  )

  # Set the weights of the edges in the graph to be twice the absolute value of the original weights
  E(FCM1.net)$weights <- abs(E(FCM1.net)$weight) * 2
  # Determine the sign of each edge (positive or negative)
  E(FCM1.net)$sign <- sign(E(FCM1.net)$weight)
  # Set the default color of edges to 'dodgerblue'
  E(FCM1.net)$color <- "dodgerblue"
  # Change the color of negative edges to 'darkorange1'
  E(FCM1.net)$color[E(FCM1.net)$sign < 0] <- "darkorange1"

  # Generate a network plot using ggnet2
  p.net <- ggnet2(FCM1.net,
    label = TRUE, label.size = 4,
    arrow.size = 15, arrow.gap = 0.02,
    edge.color = "color", edge.size = "weights"
  )

  # Define the output folder and file name for the plot
  output.folder <- paste0("./FCM matrix projections/res", group)
  file.name <- "network"
  tiff.dim <- c(1500, 1000) # Dimensions for the TIFF file

  # Save the plot as a TIFF file if save.plot is TRUE
  if (save.plot) {
    tiff(file.path(output.folder, paste0(file.name, ".tiff")),
      width = tiff.dim[1], height = tiff.dim[2],
      units = "px", pointsize = 12, res = 300,
      compression = c("lzw")
    )
  }

  # Print the network plot
  print(p.net)

  # Close the TIFF device if the plot was saved
  if (save.plot) dev.off()

  # Return the generated plot
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

state.shift <- function(greed, iter, mat, tol = 0.000001, elements = c("P", "E", "S", "T", "L", "EN")) {
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
