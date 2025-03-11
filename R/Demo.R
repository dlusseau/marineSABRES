########################################################################
####                                                        ############
#### Demo workflow marineSABRES T5.3 network analyses       ############
####                                                        ############
########################################################################
####
####   davlu
####   bmjv
####   October 2024
####
########################################################################


library(igraph)
library(data.table)
library(BoolNet)

source('~/marineSABRES/R/utils.R')

####################################################################################

# Define paths
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/Macaronesia demo/"

# Load data
macaronesia <- fread(paste0(folder,"macaronesia_isa5_May2024.csv")) # One example

###################################################
#### Load and graph

macaronesia.SES <- data.load(df = macaronesia, folder = folder, graph.name = "Macaronesia_cld_graph", graph.title = "Macaronesia SES (CLD)")

#########################################################
##### Qualitative (signed) analysis

# INPUT
# The SES matrix from the data.load function
# The folder where output should be saved
# Filenames for the .csv and graph

# In one function:
macaronesia.qual <- qualitative.analyses(SES.mat = macaronesia.SES, 
                                         folder, 
                                         filename.boolean.csv = "Macaronesia_boolean",
                                         filename.boolean.graph = "Macaronesia_boolean2")

# Stepwise:

laplacian <- SES.laplacian(SES.mat, from = "rows") 

filename.boolean.csv = "Macaronesia_boolean"
boolean.df <- boolean.file.creation(SES.mat, folder, filename.boolean.csv)
macaronesia_boolean <- loadNetwork(paste0(folder, filename.boolean.csv, ".csv"))

macaronesia_bool_ana <- boolean.analyses(macaronesia_boolean, folder, filename.boolean.graph)



attractor1 <- data.frame(name = colnames(macaronesia.SES), boolean = as.numeric(macaronesia.qual$boolean.results[[4]][[1]]))


###############################################################################
#### quantitative projection

# TODO: This can be improved including the targets
# The targets are function of stakeholders: they 

# Use new Hamiltonian Bayesian models - likelihood-free functions?
# We now go from random sampling - random forest. This could be faster if 
# Try an dmak mcmc-abc working

#####################################################################################

# INPUT
# The SES matrix from the data.load function
# The folder where output should be saved
# Filenames for the .csv and graph

greed = 500

iter = 500
filename.simu.figure = "macaronesia_weighted_projection" 
title.simu = "Macaronesia SES (CLD)"
filename.PR = "Macaronesia_participation_ratio"
filename.greedy.res = "macaronesia_greedy_simulation.Rdata"
indicators = c('Food Provision','MPA biodiversity','Residents')


quantitative.analyses <- function(SES.mat, greed = 100, iter = 500, folder, filename.simu.figure = "macaronesia_weighted_projection", 
                                  title.simu = "Macaronesia SES (CLD)",filename.PR = "Macaronesia_participation_ratio",
                                  filename.greedy.res = "macaronesia_greedy_simulation.Rdata", indicators = c('Food Provision','MPA biodiversity','Residents'),
                                  
                                  ...){
  
  macaronesia.sim <- SES.simulate(SES.mat = macaronesia.SES, iter = iter, save.fig = T, folder, fig.filename = filename.simu.figure, 
                                  fig.title = title.simu)
  
  PR <- participation_ratio(SES = macaronesia.SES, folder = folder, filename = filename.PR , title = title.simu) 
  
  macaronesia.state.shift <- state.shift(mat = macaronesia.SES, greed, iter, type = 'uniform', folder, file = 'filename.greedy.res')
  
  
  
  
  
  sim.outcomes.macaronesia <- RF.prep(state.shift.res = macaronesia.state.shift, targets = c('Food Provision','MPA biodiversity','Residents'))

  forest.macaronesia <- random.forest(sim.outcomes = sim.outcomes.macaronesia, 
                                      ntree = 500, folder, 
                                      filename = "macaronesia_forest.Rdata")
  
  gc()

  importance_frame.macaronesia <- random.forest.res(forest.macaronesia, folder, 
                                                    filename1 = "macaronesia_importance.Rdata", 
                                                    filename2 = "macaronesia_randomforest_variable_importance")

}

# .... Introduce measures .....


##### introduce measures = new link

affected <- row.names(macaronesia.SES)[2] # choose by end

indicators <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit")]
indicators <- indicators[indicators %in% row.names(tuscany.SES)]

tuscany.alt <- simulate.measure(macaronesia.SES, "tourism throttling", affected, indicators, lower = -1, upper = 0)


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








##### this is going to be a bit different, we are going to try to be memory efficient and only retain the added row and col
library(reshape2)
greed <- 1000000
iter <- 500
mat <- tuscany.SES
starting.value <- c(starting.value, 0) # we start without measures
measure <- "tourism throttling"


# Duplicated again....
targets <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit" | tuscany_element$Description == "Ecosystem Service")]
## for some reason elements are not the same as the graphs
targets <- targets[targets %in% row.names(tuscany.SES)]
affected <- row.names(tuscany.SES)[2] # choose by end

# Run this to figure out what it is exactly
indicators <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit")]
indicators <- indicators[indicators %in% row.names(tuscany.SES)]


# Can I make a state shift function that can work with either/or situation with measures?


tol <- 0.000001
state.sim.bin <- state.sim
state.sim.bin[abs(state.sim.bin) < tol] <- 0
state.sim.bin <- sign(state.sim.bin)


desirable.outcomes <- which(colSums(state.sim.bin[match(targets, row.names(state.sim.bin)), ]) == length(targets))
mat.works <- mat.sim.df[, c(1, 2, (desirable.outcomes + 2))]

##################################################################

sim.ana.df <- t(mat.sim.df[, 3:(greed + 2)])
colnames(sim.ana.df) <- apply(mat.sim.df[, 1:2], 1, function(x) paste0(levels(focus.mat.df$from)[x[1]], " to ", levels(focus.mat.df$to)[x[2]]))
sim.ana.df <- as.data.frame(sim.ana.df)

sim.ana.df$outcomes <- 0
sim.ana.df$outcomes[desirable.outcomes] <- 1

sim.ana.df$outcomes <- factor(sim.ana.df$outcomes)
select <- colSums(sim.ana.df[, 1:42])
which(select != 0)

sim.ana.df <- sim.ana.df[, c(which(select != 0), 43)]

library(GGally)


ggpairs(subset(sim.ana.df, outcomes == 1),
        columns = 1:5
)


ggpairs(sim.ana.df,
        columns = 1:5,
        aes(
          color = outcomes,
          alpha = 0.5
        )
)

# Random Forest

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename <- "tuscany_forest_measure.Rdata"

forest.tuscany_measure <- random.forest(data = sim.ana.df, ntree = 500, folder, filename)


gc()

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename1 <- "tuscany_measure_importance.Rdata"
filename2 <- "Tuscany_measure_randomforest_variable_importance"

importance_frame.tuscany_measure <- random.forest.res(forest.tuscany_measure, folder, filename1, filename2)







# .... ABC ..... This is still in trial and should not be considered for now

library(EasyABC)

# Make a list of x uniform distributions
length.vec <- prod(dim(mat))
my_prior <- rep(list(c("unif",0,1)),length.vec) # Create N 

indicators <- c("MPA biodiversity", "Food Provision") # It somehow doesn't work with only one target
set.seed(1)

sum_stat_obs <- rep(1, length(indicators))
tol.lev <- 0.000001 #indicating the proportion of simulations retained nearest the targeted summary statistics.

# It works with uniform distribution, there's an issue with the selection process
ABC_rej <- ABC_rejection(model = state.shift.ABC, prior = my_prior, nb_simul = 10000) #At least this is faster than the greedy approach
ABC_rej # Takes 320 sec
filename <- 'ABC_rej'
save(ABC_rej, file = paste0(folder, filename, ".RData"))

# RF can be applied to this

# Run long greedy simulation and check output format + test RF code
# Then try and better identify or change the selected summary statistics: can I play more around with those?


ABC_mxmc <- ABC_mcmc(method = 'Marjoram',model = state.shift.ABC, prior = my_prior, 
                     summary_stat_target = sum_stat_obs) #At least this is faster than the greedy approach
# The proposal jumps outside of the prior distribution too often - consider using the option 'inside_prior=FALSE' or enlarging the prior distribution


# Try with abc package
library(abc)

rej <- abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.1, method="rejection")
# Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.


tolerance <- c(0.01,0.0001)
BC_Beaumont <- ABC_sequential(method = "Beaumont", model = state.shift.ABC,
                              prior = my_prior, nb_simul = 1000, summary_stat_target = sum_stat_obs,
                              tolerance_tab = tolerance, inside_prior = F)
BC_Beaumont # Stats are good, but all weights are NaN

# Explore results
dim(BC_Beaumont$param)
summary(BC_Beaumont$stats)

# How to report these results .....
plot(BC_Beaumont)


ABC_rej <- ABC_rejection(model = state.shift.ABC, prior = my_prior, nb_simul = 1000,
                         summary_stat_target = sum_stat_obs, tol = 0.54, 
                         use_seed = F, progress_bar = F) # Doesn't work with the tolerance being very low or zero. Find a better way
ABC_rej

# How to red the output?
ABC_rej$param #The model parameters used in the model simulations.
ABC_rej$stats # The summary statistics obtained at the end of the model simulations.
ABC_rej$weights #	The weights of the different model simulations. In the standard rejection scheme, all model simulations have the same weights.
ABC_rej$stats_normalization #The standard deviation of the summary statistics across the model simulations.
ABC_rej$nsim # The number of model simulations performed.
ABC_rej$nrec # The number of retained simulations (if targeted summary statistics are provided).

# It doesn't work with tol.lev of 0, then nothing is retained




