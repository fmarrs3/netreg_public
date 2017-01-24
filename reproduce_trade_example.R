# Tyler H. McCormick, Bailey K. Fosdick, and Frank W. Marrs
# 
# Reproduce trade example results as in "Standard errors for regression on relational data with exchangeable errors"
#   writes out coefficients of fit, standard error estimates, and plots thereof from the paper
#


rm(list = ls())          # clear environment


################
###  Inputs  ###
################
wd <- "FILEPATH"     # working directory (location of function file)
tmax <- 4     # last time period to consider; doing all 20 takes some time to run
write_dir <- paste0("./results_trade_t", tmax)    # name/location to write results


########################
###  Function calls  ###
########################
setwd(wd)     # set working directory (location of files)
dir.create(write_dir, showWarnings = F)   # create output directory 
source('./function_file.R')   # load all functions


# Function to reproduce trade results in paper
reproduce_trade(tmax, mattype="EE", write_dir)

# Function to plot trade results as in paper
plot_beta_hoff_comparison(tmax, write_dir)

