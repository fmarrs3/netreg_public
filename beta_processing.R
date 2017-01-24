# Post-process Westveld and Hoff (2011) MCMC output

setwd('FILEPATH')
# read in the posterior samples
beta <- read.table('Beta.txt')

# compute posterior median and 95% credible interval
median.beta <- apply(beta, 2, median)
LB.beta <- apply(beta,2,quantile,.025)
UB.beta <- apply(beta,2,quantile,.975)

# organize output into matrix (8 betas X 20 time periods)
median.beta <- matrix(median.beta, 8, 20)
UB.beta <- matrix(UB.beta, 8, 20)
LB.beta <- matrix(LB.beta, 8, 20)

# save ouput
save(median.beta, UB.beta, LB.beta, file='./Betas_hoff.RData')



