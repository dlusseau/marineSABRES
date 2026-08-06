###########################################
##########################################
#### analysis of multi-nation consensus
##########################################

library(BoolNet)

library(data.table)

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


components<- states$stateInfo$genes
initialstates<-sapply(states$attractors, function(x) is.null(x$initialStates))

bin_attractors <- lapply(seq_along(states$attractors), function(A) {
  getattractorbin(
    states = states,
    A = A,
    initialstate = initialstates[A],
    caution = FALSE
  )
})


attractors.df <- rbindlist(
  lapply(seq_along(bin_attractors), function(A) {
    
    x <- bin_attractors[[A]]
    
    # skip NULL elements
    if (is.null(x)) return(NULL)
    
    # convert binary matrix/array to data.table
    x <- as.data.table(x)
    
    # assign component names
    setnames(x, components)
    
    # add attractor and state columns
    x[, attractor := paste0("attractor", A)]
    x[, state := paste0("state", seq_len(.N))]
    
    # put metadata columns first
    setcolorder(x, c("attractor", "state", components))
    
    x
  }),
  fill = TRUE
)

# attractors.df <- as.data.frame(matrix(NA, 1, 1 + length(states$stateInfo$genes), ))
# colnames(attractors.df) <- c("attractors", states$stateInfo$genes)
# 
# 
# n_attractors <- length(states$attractors)
# m <- 1
# for (i in 1:n_attractors) {
#   A1 <- states$attractors[[i]]$involvedStates
#   if (length(A1) == 1) {
#     attractors.df[m, 1] <- paste0("attractor ", i)
#     attractors.df[m, 2:ncol(attractors.df)] <- dec2bin(A1, length(states$stateInfo$genes))
#     
#     m <- m + 1
#   } else {
#     for (j in 1:length(A1)) {
#       A1b <- states$attractors[[i]]$involvedStates[j]
#       attractors.df[m, 1] <- paste0("attractor ", i)
#       attractors.df[m, 2:ncol(attractors.df)] <- dec2bin(A1b, length(states$stateInfo$genes))
#       
#       m <- m + 1
#     }
#   }
#   print(i)
#   flush.console()
# }
# attractors.df

write.csv(attractors.df, file = "arctic/cld3nation_attractors_bool.csv")

nodes<-names(attractors.df)

goals<-nodes[c(5,6,12,27,28,43,44,59,61,63)]
goalsIS<-nodes[c(5,6,12,59)]
goalsFO<-nodes[c(12,43,44,63)]
goalsGL<-nodes[c(12,27,28,61)]

desired_attractors<-attractors.df[rowSums(attractors.df[, goals, drop = FALSE] == 1) == length(goals), ]
desired_attractorsIS<-attractors.df[rowSums(attractors.df[, goalsIS, drop = FALSE] == 1) == length(goalsIS), ]
desired_attractorsFO<-attractors.df[rowSums(attractors.df[, goalsFO, drop = FALSE] == 1) == length(goalsFO), ]
desired_attractorsGL<-attractors.df[rowSums(attractors.df[, goalsGL, drop = FALSE] == 1) == length(goalsGL), ]

