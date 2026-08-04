###########################################
##########################################
#### analysis of multi-nation consensus
##########################################

library(BoolNet)
source("R/utils_final.R")
source("R/helpers.R")
source("R/boolean_file_creation_no_write.R")
cld3<-read.csv("arctic/CLD_3nations.csv")


##Boolean approach


cld3.mat <- make_matrix(cld3$From, cld3$To, cld3$weight)


colnames(cld3.mat) <- gsub(" ", "_", colnames(cld3.mat))
rownames(cld3.mat) <- gsub(" ", "_", rownames(cld3.mat))


folder <- "arctic/"
filename <- "cld3nation_boolean"
filegraph <- "cld3nation_statespace"


boolean_file_creation(cld3.mat, folder, filename)

cld3.bool <- loadNetwork(paste0(folder, filename, ".csv"))

#cld3.bool<-boolean_file_creation_no_write(cld3.mat)

states <- getAttractors(cld3.bool , method = "sat.exhaustive", returnTable = TRUE)

# large graphs use "sat.exhaustive"


attractors.df <- as.data.frame(matrix(NA, 1, 1 + length(states$stateInfo$genes), ))
colnames(attractors.df) <- c("attractors", states$stateInfo$genes)


n_attractors <- length(states$attractors)
m <- 1
for (i in 1:n_attractors) {
  A1 <- states$attractors[[i]]$involvedStates
  if (length(A1) == 1) {
    attractors.df[m, 1] <- paste0("attractor ", i)
    attractors.df[m, 2:ncol(attractors.df)] <- dec2bin(A1, length(states$stateInfo$genes))
    
    m <- m + 1
  } else {
    for (j in 1:length(A1)) {
      A1b <- states$attractors[[i]]$involvedStates[j]
      attractors.df[m, 1] <- paste0("attractor ", i)
      attractors.df[m, 2:ncol(attractors.df)] <- dec2bin(A1b, length(states$stateInfo$genes))
      
      m <- m + 1
    }
  }
  print(i)
  flush.console()
}
attractors.df

write.csv(attractors.df, file = "arctic/cld3nation_attractors_bool.csv")




