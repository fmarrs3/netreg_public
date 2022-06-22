# Frank W. Marrs, Bailey K. Fosdick, and Tyler H. McCormick
# 
# Reproduce simulations as in "Regression of exchangeable relational arrays"
#   writes out coefficients of fit, variance estimates, and coverage thereof, along with plots in paper
#


rm(list = ls())          # clear environment
gc()


################
###  Inputs  ###
################
wd <- "FILEPATH"     # working directory (location of function file)
N.test <- c(20,40,80)   # size of networks for which to simulate. WARNING: Networks over n=50 may take awhile to simulate. Paper was N in 10,20,40,80,160,320 
X.range <- 1:10   # random X matrices to test per simulation, paper was 1:100



########################
###  Function calls  ###
########################
write_dir <- "results_sims"    # name/location to write results.
setwd(wd)     # set working directory (location of files)
source('./function_file.R')   # load all functions


# Function to reproduce simulation results in paper
reproduce_simulation(N.test, X.range, write_dir)

# Function to reproduce simulation plots in paper
plot_coverage(N.test, X.range, write_dir)


