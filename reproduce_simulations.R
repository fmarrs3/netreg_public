# Tyler H. McCormick, Bailey K. Fosdick, and Frank W. Marrs
# 
# Reproduce simulations as in "Standard errors for regression on relational data with exchangeable errors"
#   writes out coefficients of fit, variance estimates, and coverage thereof, along with plots in paper
#


rm(list = ls())          # clear environment


################
###  Inputs  ###
################

wd <- "FILEPATH"     # working directory (location of function file)
write_dir <- "results_sims"    # name/location to write results.
N.test <- c(20,40,80)   # size of networks for which to simulate. WARNING: Networks over n=50 may take awhile to simulate. Paper was N in 10,20,40,80,160,320 
X.range <- 1:10   # random X matrices to test per simulation, paper was 1:100



########################
###  Function calls  ###
########################

setwd(wd)     # set working directory (location of files)
dir.create(write_dir, showWarnings = F)   # create output directory 
source('./function_file.R')   # load all functions


# Function to reproduce simulation results in paper
reproduce_simulation(N.test, X.range, write_dir)

# Function to reproduce simulation plots in paper
plot_coverage(N.test, X.range, write_dir)


