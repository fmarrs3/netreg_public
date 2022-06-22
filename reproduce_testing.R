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
n <- 50  # number of actors in simulated data
nerr <- 200  # Number of data sets simulated; 1,000 used in paper but takes a while to run


###################
###  Test sims  ###
###################
write_dir <- "./results_testing"   # look here for results
setwd(wd)
source("./function_file.R")

# Run simulations and write out results
reproduce_testing_sim(n, nerr=nerr, write_dir=write_dir)


