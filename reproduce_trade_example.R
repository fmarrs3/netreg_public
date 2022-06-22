# Frank W. Marrs, Bailey K. Fosdick, and Tyler H. McCormick
# 
# Reproduce trade example results as in "Regression of exchangeable relational arrays"
#   writes out and plots prediction performance from the paper
#   coefficients of fit, standard error estimates
#


rm(list = ls())          # clear environment
gc()


################
###  Inputs  ###
################
wd <- "FILEPATH"     # working directory (location of function file)
tmax <- 6     # last time period to consider; doing all 20 takes some time to run



########################
###  Function calls  ###
########################
setwd(wd)     # set working directory (location of files)
source('./function_file.R')   # load all functions
source('./bayes_mcmc_function.R')   # load Westveld and Hoff comparison code

# Run fits for all time periods under consideration, starting with 4
for(t in 4:tmax){
  # Function to reproduce proposed exchangeable fit
  write_dir <- paste0("./results_trade_gee_t", t)    # name/location to write results
  reproduce_trade(t, write_dir)

  # Run mixed effects Bayesian code
  write_dir <- paste0("./results_trade_bayes_t", t)    # name/location to write results
  fit.bayes(t, write_dir)

  cat("##############\n")
  cat("t=", t, "\n")
}

# Reproduce prediction plot
plot_r2_predictions(tmax)
  




