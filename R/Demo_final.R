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

source('~/marineSABRES/R/utils_final.R')

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

#####################################################################################

# INPUT
# The SES matrix from the data.load function
# The folder where output should be saved
# Filenames for the .csv and graph

macaronesia.quant <- quantitative.analyses(SES.mat = macaronesia.SES, folder = folder, 
                                           filename.simu.figure = "macaronesia_weighted_projection", 
                                           title.simu = "Macaronesia SES (CLD)",
                                           filename.PR = "Macaronesia_participation_ratio",
                                           filename.greedy.res = "macaronesia_greedy_simulation.Rdata", 
                                           indicators = c('Food Provision','MPA biodiversity','Residents'),
                                           file.RF = "macaronesia_forest.Rdata", 
                                           filename1.RF = "macaronesia_importance.Rdata", 
                                           filename2.RF = "macaronesia_randomforest_variable_importance")



# Stepwise (see also utils_final.R l.529:

macaronesia.sim <- SES.simulate(SES.mat = macaronesia.SES, iter = iter, save.fig = T, folder, fig.filename = filename.simu.figure, 
                                fig.title = title.simu)

PR <- participation_ratio(SES.mat = macaronesia.SES, folder = folder, filename = filename.PR , title = title.simu) 

macaronesia.state.shift <- state.shift(SES.mat = macaronesia.SES, greed, iter, type = 'uniform', folder, file = filename.greedy.res)

sim.outcomes.macaronesia <- RF.prep(state.shift.res = macaronesia.state.shift, targets = c('Food Provision','MPA biodiversity','Residents'))

forest.macaronesia <- random.forest(sim.outcomes = sim.outcomes.macaronesia, 
                                    ntree = 500, folder, 
                                    file.RF = "macaronesia_forest.Rdata")
gc()

importance_frame.macaronesia <- random.forest.res(forest.macaronesia, folder, 
                                                  filename1 = "macaronesia_importance.Rdata", 
                                                  filename2 = "macaronesia_randomforest_variable_importance")



# .... Introduce measures .....


macaronesia.measure <- state.shift.measure(measure = '...', affected = '...', SES.mat = macaronesia.SES, greed, iter, 
                    folder, file.measure = 'greedy.res.measure')





##### introduce measures = new link

affected <- row.names(macaronesia.SES)[2] # choose by end

indicators <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit")]
indicators <- indicators[indicators %in% row.names(tuscany.SES)]

tuscany.alt <- simulate.measure(macaronesia.SES, "tourism throttling", affected, indicators, lower = -1, upper = 0)


##### this is going to be a bit different, we are going to try to be memory efficient and only retain the added row and col
library(reshape2)
greed <- 1000000
iter <- 500
mat <- tuscany.SES
starting.value <- c(starting.value, 0) # we start without measures
measure <- "tourism throttling"















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










