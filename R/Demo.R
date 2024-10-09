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


####################################################################################

# TODO:
# Where did the Macaronesia and Atlantic graphs go?
# Ended with the definition of outcomes and desirable outcomes for Macaronesia at random Forest function. Based on what is this defined?

## Macaronesia


# Output folder
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/Macaronesia demo/"

# Load data
macaronesia <- fread(paste0(folder,"macaronesia_Oct2024.csv"))

macaronesia$weight <- 1
macaronesia$weight[macaronesia$strength == "Medium Positive"] <- .5
macaronesia$weight[macaronesia$strength == "Medium Negative"] <- (-.5)
macaronesia$weight[macaronesia$strength == "Strong Negative"] <- (-1)
macaronesia$sign <- sign(macaronesia$weight)


macaronesia.SES <- make_matrix(from = macaronesia$from, to = macaronesia$to, weight = macaronesia$weight)

macaronesia.laplacian <- SES.laplacian(macaronesia.SES, from = "rows")


###################################################
#### Graph

macaronesia.graph <- plot.SES(macaronesia.SES,
                          folder = folder,
                          filename = "Macaronesia_cld_graph",
                          title = "Macaronesia SES (CLD)",
                          label.cex = 1.1,
                          vsize = 30
)

# To double check

#########################################################
##### Boolean

filename <- "Macaronesia_boolean"

boolean_file_creation(macaronesia.SES, folder, filename)


macaronesia_boolean <- loadNetwork(paste0(folder, filename, ".csv"))

filename <- "Macaronesia_boolean2"
macaronesia_bool_ana <- boolean_analyses(macaronesia_boolean, folder, filename)

attractor1 <- data.frame(name = colnames(macaronesia.SES), boolean = as.numeric(macaronesia_bool_ana[[4]][[1]]))
# All positive

###############################################################################
#### quantitative projection

iter <- 1000

macaronesia.sim <- SES.simulate(macaronesia.SES, 1000, folder, filename = "macaronesia_weighted_projection", 
                                title = "Macaronesia SES (CLD)")

sim.melt <- melt(macaronesia.sim)

fillcol <- glasbey.colors(nrow(macaronesia.sim))
names(fillcol) <- unique(sim.melt$Var1)

ggplot(subset(sim.melt, Var2 > 50 & Var2 < 70), aes(x = (Var2), y = (value), colour = (Var1))) +
  geom_path(size = 2) +
  labs(colour = "Element", y = "Progress", x = "(time)") +
  scale_colour_manual(values = fillcol) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none")

####################################################################################################
##############################
### current SES participation ratio

filename <- "Macaronesia_participation_ratio"
title <- "Macaronesia (CLD)"
tuscany.PR <- participation_ratio(macaronesia.SES, folder = folder, filename = filename, title = title) 


###############################################################################################
### Greedy approach

# TODO: This can be improved including the targets

starting.value <- runif(nrow(mat), 0, 1) # Create random starting values between 0 and 1

greed <- 1000#000 Take less for  now befor double checking everything
iter <- 500
mat <- macaronesia.SES
starting.value <- runif(nrow(mat), 0, 1)

file <- "macaronesia_greedy_simulation.Rdata"

macaronesia.state.shift <- state.shift(greed, iter, mat, starting.value, folder, file)


tol <- 0.000001

# Put this into function as well?
macaronesia.state.sim.bin <- macaronesia.state.shift
macaronesia.state.sim.bin[abs(macaronesia.state.sim.bin) < tol] <- 0
macaronesia.state.sim.bin <- sign(macaronesia.state.sim.bin)

# Why do we calculate these desirable outcomes??
targets <- unique(macaronesia$to) # Double check this with David

desirable.outcomes <- which(colSums(macaronesia.state.sim.bin[match(targets, row.names(macaronesia.state.sim.bin)), ]) == length(targets))
mat.works <- mat.sim.df[, c(1, 2, (desirable.outcomes + 2))]


##################################################################
# Analyse with randomForest

# To test
sim.ana <- macaronesia.state.shift$mat.sim.df[which(macaronesia.state.shift$mat.m$value != 0),]
sim.ana.df <- t(sim.ana[, 3:(greed + 2)])
colnames(sim.ana.df) <- apply(sim.ana[, 1:2], 1, function(x) paste0(row.names(macaronesia.state.shift$mat)[x[1]], " to ", row.names(macaronesia.state.shift$mat)[x[2]]))
sim.ana.df <- as.data.frame(sim.ana.df)

sim.ana.df$outcomes <- 0
sim.ana.df$outcomes[desirable.outcomes] <- 1

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename <- "tuscany_forest.Rdata"

forest.tuscany <- random.forest(data = sim.ana.df, ntree = 500, folder, filename)


gc()

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename1 <- "tuscany_importance.Rdata"
filename2 <- "Tuscany_randomforest_variable_importance"

importance_frame.tuscany <- random.forest.res(forest.tuscany, folder, filename1, filename2)


########################################
########################################
########################################
########################################

###
##### introduce measures

library(igraph)
library(data.table)


tuscany <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD.csv")
tuscany_element <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD_element.csv")
tuscany_element <- tuscany_element[, 1:2]

unique(tuscany$strength)

tuscany$weight <- 1
tuscany$weight[tuscany$strength == "Medium Positive"] <- .5
tuscany$weight[tuscany$strength == "Medium Negative"] <- (-.5)
tuscany$weight[tuscany$strength == "Strong Negative"] <- (-1)
tuscany$sign <- sign(tuscany$weight)

tuscany.SES <- make_matrix(from = tuscany$From, to = tuscany$To, weight = tuscany$weight)


targets <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit" | tuscany_element$Description == "Ecosystem Service")]
## for some reason elements are not the same as the graphs

# This is already duplicated
targets <- targets[targets %in% row.names(tuscany.SES)]

affected <- row.names(tuscany.SES)[2] # choose by end

indicators <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit")]
indicators <- indicators[indicators %in% row.names(tuscany.SES)]

tuscany.alt <- simulate.measure(tuscany.SES, "tourism throttling", affected, indicators, lower = -1, upper = 0)






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






####################################################################################


## Tuscany


# Output folder
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"

# Load data

tuscany <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD.csv")
tuscany_element <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD_element.csv")
tuscany_element <- tuscany_element[, 1:2]

unique(tuscany$strength)

tuscany$weight <- 1
tuscany$weight[tuscany$strength == "Medium Positive"] <- .5
tuscany$weight[tuscany$strength == "Medium Negative"] <- (-.5)
tuscany$weight[tuscany$strength == "Strong Negative"] <- (-1)
tuscany$sign <- sign(tuscany$weight)

tuscany.SES <- make_matrix(from = tuscany$From, to = tuscany$To, weight = tuscany$weight)

tuscany.laplacian <- SES.laplacian(tuscany.SES, from = "rows")

###################################################
#### graph

tuscany.graph <- plot.SES(tuscany.SES,
  folder = folder,
  filename = "tuscany_cld_graph",
  title = "Tuscany SES (CLD)",
  label.cex = 1.1,
  vsize = 30
)

#########################################################
##### boolean

filename <- "Tuscany_boolean"

boolean_file_creation(tuscany.SES, folder, filename)


Tuscany_boolean <- loadNetwork(paste0(folder, filename, ".csv"))


filename <- "Tuscany_boolean2"

tuscany_bool_ana <- boolean_analyses(Tuscany_boolean, folder, filename)

attractor1 <- data.frame(name = colnames(tuscany.SES), boolean = as.numeric(tuscany_bool_ana[[4]][[1]]))


###############################################################################
#### quantitative projection

iter <- 1000

tuscany.sim <- SES.simulate(tuscany.SES, 1000, folder, filename = "tuscany_weighted_projection", title = "Tuscany SES (CLD)")

sim.melt <- melt(tuscany.sim)

fillcol <- glasbey.colors(nrow(tuscany.sim))
names(fillcol) <- unique(sim.melt$Var1)

ggplot(subset(sim.melt, Var2 > 50 & Var2 < 70), aes(x = (Var2), y = (value), colour = (Var1))) +
  geom_path(size = 2) +
  labs(colour = "Element", y = "Progress", x = "(time)") +
  scale_colour_manual(values = fillcol) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none")

####################################################################################################
##############################
### current SES participation ratio

filename <- "Tuscany_participation_ratio"
title <- "Tuscany (CLD)"
tuscany.PR <- participation_ratio(tuscany.SES, folder = folder, filename = filename, title = title) # Politics has the broader reach


###############################################################################################

targets <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit" | tuscany_element$Description == "Ecosystem Service")]
## for some reason elements are not the same as the graphs

targets <- targets[targets %in% row.names(tuscany.SES)]

#### what weight is important?

starting.value <- runif(nrow(mat), 0, 1) # Create random starting values between 0 and 1


greed <- 1000000
iter <- 500
mat <- tuscany.SES
starting.value <- runif(nrow(mat), 0, 1)
folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
file <- "tuscany_greedy_simulation.Rdata"

tuscany.state.shift <- state.shift(greed, iter, mat, starting.value, folder)

# Put this into function as well?
state.sim.bin <- state.sim
state.sim.bin[abs(state.sim.bin) < tol] <- 0
state.sim.bin <- sign(state.sim.bin)

desirable.outcomes <- which(colSums(state.sim.bin[match(targets, row.names(state.sim.bin)), ]) == length(targets))
mat.works <- mat.sim.df[, c(1, 2, (desirable.outcomes + 2))]

save(tuscany.state.shift, file = "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/tuscany_greedy_simulation.Rdata")


##################################################################
# Analyse with randomForest

# To test
sim.ana <- tuscany.state.shift$mat.sim.df[which(tuscany.state.shiftmat.m$value != 0), ]
sim.ana.df <- t(sim.ana[, 3:(greed + 2)])
colnames(sim.ana.df) <- apply(sim.ana[, 1:2], 1, function(x) paste0(row.names(mat)[x[1]], " to ", row.names(mat)[x[2]]))
sim.ana.df <- as.data.frame(sim.ana.df)

sim.ana.df$outcomes <- 0
sim.ana.df$outcomes[desirable.outcomes] <- 1

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename <- "tuscany_forest.Rdata"

forest.tuscany <- random.forest(data = sim.ana.df, ntree = 500, folder, filename)


gc()

folder <- "C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/outputs for Madeira GA Oct 2024/"
filename1 <- "tuscany_importance.Rdata"
filename2 <- "Tuscany_randomforest_variable_importance"

importance_frame.tuscany <- random.forest.res(forest.tuscany, folder, filename1, filename2)


########################################
########################################
########################################
########################################

###
##### introduce measures

library(igraph)
library(data.table)


tuscany <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD.csv")
tuscany_element <- fread("C:/Users/bmjv/OneDrive - Danmarks Tekniske Universitet/SABRES/loopanalyses/the graphs/Tuscany_Jul_2024_CLD_element.csv")
tuscany_element <- tuscany_element[, 1:2]

unique(tuscany$strength)

tuscany$weight <- 1
tuscany$weight[tuscany$strength == "Medium Positive"] <- .5
tuscany$weight[tuscany$strength == "Medium Negative"] <- (-.5)
tuscany$weight[tuscany$strength == "Strong Negative"] <- (-1)
tuscany$sign <- sign(tuscany$weight)

tuscany.SES <- make_matrix(from = tuscany$From, to = tuscany$To, weight = tuscany$weight)


targets <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit" | tuscany_element$Description == "Ecosystem Service")]
## for some reason elements are not the same as the graphs

# This is already duplicated
targets <- targets[targets %in% row.names(tuscany.SES)]

affected <- row.names(tuscany.SES)[2] # choose by end

indicators <- tuscany_element$Label[which(tuscany_element$Description == "Good and Benefit")]
indicators <- indicators[indicators %in% row.names(tuscany.SES)]

tuscany.alt <- simulate.measure(tuscany.SES, "tourism throttling", affected, indicators, lower = -1, upper = 0)






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
