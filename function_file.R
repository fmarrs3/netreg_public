# Tyler H. McCormick, Bailey K. Fosdick, and Frank W. Marrs
# 
# functions in support of reproduce_simulations.R and reproduce_trade_example.R
#


library(Matrix)
library(MASS)


########################
###  Main Functions  ###
########################

# Reproduce simulations in paper
reproduce_simulation <- function(N.test, X.range=1:100, write_dir="./results")
{
  # Replicate results from McCormick et. al. 
  # N.test is a vector of network sizes
  # X.test is the number of X matrices to test per simulation
  # All outputs are written out to data files
  #
  
  
  Model.test <- c('iid', 'ble', 'rand')          # iid, sr, ble, rand, data generation type
  f.e <- .25    # proportion of total variance of iid/measurement error
  s.e.max <- sqrt(3)                 # KNOB FOR S/N RATIO:  iid standard deviation (maximum, when f.e = 1)
  Rep.test <- 1000            # Number of replications to run
  rho <- .5                   # send/recieve correlation coefficient
  n.z <- 2                   # number of dimensions in simulated latent space

  
  for (n.loop in 1:length(N.test)) {
    n.tot <- N.test[n.loop]     # number of members
    d.tot <- n.tot*(n.tot-1)    # number of dyads in directed setting
    
    cat("\n\n ************** n =", n.tot, "nodes ************** \n")
    
    # Generate dyad indices for shared membership, no overlaps
    n.list <- node.set(n.tot)
    
    for (m.loop in 1:length(Model.test)) {
      model <- Model.test[m.loop] 
      
      cat("\n *************** Model", model, "*************** \n")
      
      for (x.loop in X.range) { 
        
        if(x.loop == 1){
          cat("Working on X:  1  ")
        } else { cat(x.loop, " ")}
        
        
        # source(file.path(getwd(),'160216_run_file_bigN.R'))
        set.seed(x.loop)
        X.1 <- Generate.X.bv(n.tot)
        p.tot <- ncol(X.1)
        beta.1 <- rep(1, p.tot)
        # EY.1 <- X.1 %*% beta.1
        # XtXinv.1 <- chol2inv(chol(crossprod(X.1,)))      # sandwich estimator "bread"
        # h.ii <- rep(0,d.tot)
        # for (i in 1:d.tot) {
        #   h.ii[i] <- X.1[i,] %*% XtXinv.1 %*% (X.1[i,])
        # }
          
        beta <- data.frame(matrix(NA, Rep.test, ncol(X.1)))
        names(beta) <- paste0("b", 1:ncol(X.1))
        se <- data.frame(matrix(NA, Rep.test, ncol(X.1)*3))
        names(se) <- c(sapply(c("se_HC", "se_DC", "se_E"), paste0, 1:4))
        # COV.pd <- data.frame(matrix(rep(0,4),ncol=2,nrow=2))
        # names(COV.pd) <- c("singular","negdef")
        # rownames(COV.pd) <- c("E","DC")
          for (r.loop in 1:Rep.test) {
            
            # Generate frame of simulated data
            data.sim <- Generate.data(X.1, beta.1, model=model, f.in=f.e, rho=rho, s.e.max=s.e.max, n.z=n.z)
            
            # Fit models and Save out data
            result  <- fit.data(data.sim$Y, X.1, n.list)
            beta[r.loop, ] <- result$beta
            se[r.loop, ] <- sqrt(c(diag(result$var_HC), diag(result$var_DC), diag(result$var_E)))
            # if(result$E.sing==TRUE){COV.pd[1,1] <- COV.pd[1,1] + 1}
            # if(result$E.negdef==TRUE){COV.pd[1,2] <- COV.pd[1,2] + 1}
            # if(result$DC.sing==TRUE){COV.pd[2,1] <- COV.pd[2,1] + 1}
            # if(result$DC.negdef==TRUE){COV.pd[2,2] <- COV.pd[2,2] + 1}
          }  # end of rep loop
        
          # filename <- file.path(write_dir, paste0(model, "_n", n.tot, "_x", x.loop, "_COV_pd.txt"))
          # write.table(COV.pd/Rep.test,filename,sep='\t')
          # 
          filename <- file.path(write_dir, paste0(model, "_n", n.tot, "_x", x.loop, "_beta.txt"))
          write.table(beta, filename, row.names=F, sep='\t')

          filename <- file.path(write_dir, paste0(model, "_n", n.tot, "_x", x.loop, "_se.txt"))
          write.table(se, filename, row.names=F, sep='\t')
        
      }}}
  
  
}  

# Reproduce simulation coverage plots in paper
plot_coverage <- function(N.test, X.range, write_dir="./results", make_legend=FALSE)
{
  blw=1.25
  cex.pch=3/4
  
  p.tot <- 4
  coverage.prob <- array(0, c(length(N.test), length(X.range), 3*p.tot))
  Model.test <- c('iid', 'ble', 'rand')          # iid, sr, ble, rand, data generation type

  for (m.loop in 1:length(Model.test)) {
    model <- Model.test[m.loop]  
    
    for (n.loop in 1:length(N.test)) {
      n.tot <- N.test[n.loop]     # number of members
      d.tot <- n.tot*(n.tot-1)    # number of dyads in directed setting
      
      for (x.loop in X.range) {
        filename <- file.path(write_dir, paste0(model, "_n", n.tot, "_x", x.loop, "_beta.txt"))
        beta.hat <- read.table(file=filename, header=T)
        if(ncol(beta.hat) != p.tot){ "Error in number of columns of beta estimates"}
        Rep.test <- nrow(beta.hat)
        beta.true <- outer(rep(1,Rep.test), rep(1,p.tot))
        
        filename <- file.path(write_dir, paste0(model, "_n", n.tot, "_x", x.loop, "_se.txt"))
        se.hat <- read.table(file=filename, header=T)

        
        # Post process empirical variances and MSEs
        for (i in 1:p.tot) {
          i.names <- paste0(c("se_HC", "se_DC", "se_E"), i)
          i.indices <- sapply(i.names, function(x) which(names(se.hat) == x))
          UB.beta <- beta.hat[,i] + 1.96*se.hat[, i.names]
          LB.beta <- beta.hat[,i] - 1.96*se.hat[, i.names]
          coverage.ind <- (LB.beta <= beta.true[,i]) * (UB.beta >= beta.true[,i])  # = 1 when both true, i.e. beta is in CI
          coverage.prob[n.loop, x.loop, i.indices] <- apply(coverage.ind, 2, mean, na.rm=T)   #*100
          
        }}  # end of p, X
      
      # Write out coverage probability
      table.out <- data.frame(coverage.prob[n.loop,,]) 
      names(table.out) <- c(sapply(c("cp_HC", "cp_DC", "cp_E"), paste0, 1:4))
      write.table(coverage.prob, file.path(write_dir, paste0("coverage_", model, "_n", n.tot, ".txt")), row.names=F, sep='\t')
  
    }  # end of N
    
    
    legend.loc <- rep("topright",6);  #legend.loc[c(1:3,6)]  <- 'bottomright'
    col.vec <- c('#d7191c', '#8da0cb', '#fdae61')
    nn <- 5
  
    # for (i in 1:p.tot){
    #   y.name <- 'Coverage'
    #   filename <- paste0("CP", i, "_", model,".pdf")
    #   pdf( file=file.path(write_dir, filename), width = 4, height = 4.4)
    #   par(mar = c(3, 3, .2, .2), mgp=c(1.8,.5,0), cex.lab=1.3, cex.axis=1.2) 
    #   
    #   blw <- 2   # box line width  
    #   
    #   i.names <- paste0(c("se_HC", "se_DC", "se_E"), i)
    #   seq.plots <- sapply(i.names, function(x) which(names(se.hat) == x))
    #   y.range <- c(.2, 1.00)
    #   
    #   boxplot(t(coverage.prob[,,seq.plots[3]]), ylab=y.name, xlab="Number of nodes", names=N.test, border=col.vec[1], boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,  
    #           ylim=y.range, at=(nn*(1:length(N.test))-(nn - 1)), xlim=c(0,(nn*length(N.test))-1), pch=1, xaxt="n" )
    #   boxplot(t(coverage.prob[,,seq.plots[2]]), names=N.test, add=T, boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,
    #           border=col.vec[2], at=(nn*(1:length(N.test))-(nn - 2)),  pch=2)
    #   boxplot(t(coverage.prob[,,seq.plots[1]]), add=T, border=col.vec[3], boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,
    #           at=(nn*(1:length(N.test))- (nn - 3) ),xaxt="n", pch=3)
    #   abline(h=.95,col='black',lty=2, lwd=1.5)
    #   
    #   if (i == 4) { 
    #     legend('bottomleft', c("Exchangeable","Dyadic Clustering", "Heteroskedasticity-Consistent", 'True 95%'), lty=c(1,1,1,2), lwd=3*c(1,1,1,1),
    #            col = c(col.vec, 'black'), 
    #            bty='n', horiz=F, cex = 1.15)
    #   }
    #   
    #   dev.off()
    # }  # end of p
    
    
    # box_quant <- apply(coverage.prob, c(1,3), quantile, probs=c(.025,.1,.5,.9,.975))   # array of quantiles x N x error type
    box_quant = coverage.prob
    
  ## Quantile plots
    for (i in 2:p.tot){
      if(i == 2){
        y.name <- 'Coverage'
        draw.y.axis <- NULL
      } else {
        y.name <- ''
        draw.y.axis <- "n"
      }
      # y.name <- 'Coverage'
      filename <- paste0("CP", i, "_", model,"_quantiles.pdf")
      pdf( file=file.path(write_dir, filename), width = 5, height = 5.5)

      par(mar = c(3, 3, .2, .2), mgp=c(1.8,.5,0), cex.lab=1.5, cex.axis=1.3) 
      
      
      
      i.names <- paste0(c("se_HC", "se_DC", "se_E"), i)
      seq.plots <- sapply(i.names, function(x) which(names(se.hat) == x))
      y.range <- c(.2, 1.00)
      
      plot(NA, NA, xlim=c(0,(nn*length(N.test))-1), ylim=y.range, 
           yaxt=draw.y.axis,
           xlab="Number of nodes", xaxt="n", ylab=y.name)
      sapply(c(seq(.2,1,by=.2), .9), function(x) abline(h=x, col="gray80", lwd=.8))
      
      boxplot(t(box_quant[,,seq.plots[3]]), names=N.test, border=col.vec[1],
              yaxt=draw.y.axis, cex=cex.pch,
              boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,  
              ylim=y.range, at=(nn*(1:length(N.test))-(nn - 1)), xlim=c(0,(nn*length(N.test))-1), pch=20, xaxt="n", add=T)
      boxplot(t(box_quant[,,seq.plots[2]]), names=N.test, add=T, boxlwd=blw, 
              yaxt=draw.y.axis, cex=cex.pch,
              whisklty=2, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,
              border=col.vec[2], at=(nn*(1:length(N.test))-(nn - 2)),  pch=17)
      boxplot(t(box_quant[,,seq.plots[1]]), add=T, border=col.vec[3], 
              yaxt=draw.y.axis, cex=cex.pch,
              boxlwd=blw, whisklty=3, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw,
              at=(nn*(1:length(N.test))- (nn - 3) ),xaxt="n", pch=15)
      abline(h=.95,col='black',lty=2, lwd=.75)
      
      if (i == 4 & make_legend) { 
        legend('bottomleft', c("Exchangeable","Dyadic Clustering", "Heteroskedasticity-Consistent", 'True 95%'), lty=c(1,1,1,2), lwd=3*c(1,1,1,1),
               col = c(col.vec, 'black'), 
               bty='n', horiz=F, cex = 1.15)
      }
      
      
      
      dev.off()
    }  # end of p
    
    
    } # end of model
  }

# Reproduce trade example in paper
reproduce_trade <- function(tmax, mattype, write_dir="./results", beta_start="")
{
  
  data <- read_trade_data(tmax)
  
  output <- GEE.est(data$Y, data$X, n.tot=max(data$node.mat), t.max=max(data$t.vec), type=mattype, write_dir=write_dir, beta_start=beta_start, node.mat=data$node.mat)
  
  # GEE.est <- function(Y.in, X.in, n.tot, t.max, tol.in=1e-6, directed=T, type="EE", write_dir=getwd(), beta_start="", node.mat=NA) 
    
}

# # Reproduce plots as in paper
# plot_beta_hoff_comparison <- function(tmax, write_dir='./results', loop="max")
# {
#   # wd <- getwd()
#   # 
#   # setwd(wd)
#   # write_dir <- "./results"
#   
#   tmp.hoff <- as.matrix(read.table(file.path(getwd(), 'beta_hoff.txt'), header=F))
#   beta.hoff <- t(matrix(tmp.hoff, 20, 8))
#   
#   tmp.ols <- as.matrix(read.table(file.path(write_dir, 'beta_ols.txt'), header=F))
#   beta.ols <- t(matrix(tmp.ols, tmax, 8))
#   
#   files <- list.files(write_dir, pattern= paste0("(beta_)+.*(", "t", tmax, ".txt)+"))  # find last simulation file
#   if(loop == "max"){
#     tmp.EE <- as.matrix(read.table(file.path(write_dir, files[length(files)]), header=F))
#   } else {
#     tmp.EE <- as.matrix(read.table(file.path(write_dir, files[grep(paste0("loop", loop, "_"), files)]), header=F))
#   }
#   beta.EE <- t(matrix(tmp.EE, tmax, 8))
#   
#   
#   
#   # Plots
#   #with lines
#   names <- c('intercept', 'log gdp of exporter', 'log gdp of importer', 'log distance',
#              'policy of exporter', 'policy of importer', 
#              'cooperation in conflict', 'policy interaction')
#   for(i in 1:8){
#     pdf(file=file.path(write_dir, paste0("beta_comparison_t", tmax, "_", i,'.pdf')), width=4, height=4)
#     par(mar = c(3, 3, 3, .7), mgp=c(2,.5,0), cex.lab=1.3, cex.axis=1.2)  
#     y.range <- range( c( beta.EE[i,], beta.ols[i,], beta.hoff[i,]) )
#     plot(1:tmax, beta.hoff[i, 1:tmax], xlab='', ylab='Estimated coefficient', main=names[i], pch=19, ylim=y.range,type='l',xlim=c(.7, tmax + .2),xaxt="n")
#     if(tmax == 20){
#       axis(side=1,at=c(0,5,10,15,20),labels=c("1980","1985","1990","1995","2000"),las=2)
#     } else {
#       axis(side=1,at=seq(0, tmax, by=5),labels=c("1980","1985","1990","1995","2000")[1:(1+floor(tmax/5))]  ,las=2)
#     }
#     for(ttt in 1:tmax){abline(v=ttt,lty=3,lwd=.8,col='gray')}
#     points(1:tmax, beta.ols[i, 1:tmax], pch=17, col="green", lwd=1.5,type='l',lty=1)
#     points(1:tmax, beta.ols[i, 1:tmax], pch=17, col="green", lwd=2,type='p',cex=1.3)
#     #points(1:20, beta.gee.homo[i,], pch=2, col="red", lwd=2,type='l')
#     #points(1:20, beta.gee.hetero[i,], pch=4, col="purple", lwd=2,type='l')
#     points(1:tmax, beta.EE[i, 1:tmax], pch=18, col="orange", lwd=1.5,type='l',lty=1)
#     points(1:tmax, beta.EE[i, 1:tmax], pch=18, col="orange", lwd=2,type='p',cex=1.3)
#     points(1:tmax, beta.hoff[i, 1:tmax], pch=19,lwd=1.5,type='l',lty=1)
#     points(1:tmax, beta.hoff[i, 1:tmax], pch=19,lwd=2,type='p')
#     #if(i==8){legend('topleft', legend=c('Hoff', 'OLS', 'GEE Homo E/I', 'GEE Hetero E/I', 'GEE Homo E/E'), cex=.9,
#     #pch=c(19, 0, 2, 4, 6), col=c("black", "green", "red", "purple", "orange"), bty='n')}
#     if(i==8){legend('topleft', legend=c('GEE Homoskedastic E/E', 'Westveld & Hoff 2011','OLS'), cex=1.2,
#                     pch=c(18, 19, 17), col=c("orange", "black", "green"), lty=c(1,1,1),bty='n')}
#     dev.off()
#   }
# }

# Function to reproduce coefficient plots from trade example in paper
plot_beta_hoff_comparison <- function(tmax, write_dir='./results_t20_test_160904_EE')
{
  dir.create(write_dir, showWarnings=F)
  
  whfile <- 'Betas_hoff.RData'
  cflag <- whfile %in% list.files()
  if(cflag){
    load(file=whfile)
    beta.hoff <- median.beta
  } 
  tmp.ols <- as.matrix(read.table(file.path(write_dir, 'beta_ols.txt'), header=F))
  beta.ols <- t(matrix(tmp.ols, tmax, 8))
  
  # tmp.ols.se <- as.matrix(read.table(file.path(write_dir, 'beta_olsSE.txt'), header=F))
  tmp.ols.se <- as.matrix(read.table(file.path(write_dir, paste0('se_dc_t', tmax, '.txt')), header=F))
  beta.ols.se <- t(matrix(tmp.ols.se, tmax, 8))
  
  files <- list.files(write_dir, pattern= paste0("(beta_)+.*(", "t", tmax, ".txt)+"))  # find last simulation file
  tmp.EE <- as.matrix(read.table(file.path(write_dir, files[length(files)]), header=F))
  beta.EE <- t(matrix(tmp.EE, tmax, 8))
  se.EE.tmp<-as.matrix(read.table(file.path(write_dir, paste0("se_t", tmax,".txt")), header=F))
  se.EE<-t(matrix(se.EE.tmp,tmax,8))
  
  
  # Plots *with lines*
  names <- c('intercept', 'log gdp of exporter', 'log gdp of importer', 'log distance',
             'policy of exporter', 'policy of importer', 
             'cooperation in conflict', 'policy interaction')
  for(i in 1:8){
    pdf(file=file.path(write_dir, paste0("beta_comparison_t", tmax, "_", i,'.pdf')), width=4.5, height=3) #scratch2 width=8, height=6
    par(mar = c(2.6, 2.1, 1.2, .7), mgp=c(1.3,.5,0), cex.lab=1.3, cex.axis=1.2)
    rangevec <- c( beta.EE[i,]-1.96*se.EE[i,], beta.ols[i,]-1.96*beta.ols.se[i,], beta.EE[i,]+1.96*se.EE[i,], beta.ols[i,]+1.96*beta.ols.se[i,])
    if(cflag){ rangevec <- c(rangevec, LB.beta[i,], UB.beta[i,])}
    y.range <- range( rangevec )
    if(i%%2==1){plot(1:tmax, rep(NA,tmax), xlab='', ylab='Estimated coefficient', main=names[i], pch=19, ylim=y.range,type='l',xlim=c(.7, tmax + .2),xaxt="n",cex.axis=.95,cex.lab=.95,font.main=1)
    }else{plot(1:tmax, rep(NA,tmax), xlab='', ylab='', main=names[i], pch=19, ylim=y.range,type='l',xlim=c(.7, tmax + .2),xaxt="n",cex.axis=.95,cex.lab=.95,font.main=1)}
    if(tmax == 20){ #beta.hoff[i, 1:tmax]
      axis(side=1,at=c(0,5,10,15,20),labels=c("1980","1985","1990","1995","2000"),las=2,cex.axis=.95)
    } else {
      axis(side=1,at=seq(0, tmax, by=5),labels=c("1980","1985","1990","1995","2000")[1:(1+floor(tmax/5))],las=2,cex=.95)
    }
    eps <- .1
    ##'#d7191c', '#8da0cb', '#fdae61'
    
    # for(ttt in 1:tmax){abline(v=ttt,lty=3,lwd=.2,col='gray')}
    
    ## OLS
    abline(h=0,lty=2,lwd=.7)
    points((1:tmax)-eps, beta.ols[i, 1:tmax], pch=17, col='#8da0cb', lwd=1.5,type='l',lty=1)
    points((1:tmax)-eps, beta.ols[i, 1:tmax], pch=17, col='#8da0cb', lwd=2,type='p',cex=1.3)
    segments(x0=(1:tmax)-eps, y0=beta.ols[i, 1:tmax]-(1.96*beta.ols.se[i,1:tmax]),x1=(1:tmax)-eps,y1=beta.ols[i, 1:tmax]+(1.96*beta.ols.se[i,1:tmax]), pch=17, col='#8da0cb', lwd=1.5,lty=1)
    
    #points(1:20, beta.gee.homo[i,], pch=2, col="red", lwd=2,type='l')
    #points(1:20, beta.gee.hetero[i,], pch=4, col="purple", lwd=2,type='l')
    
    # Westveld and Hoff
    if(cflag){
      points((1:tmax)+eps, beta.hoff[i, 1:tmax], pch=19,lwd=1.5,type='l',lty=1,col='#fdae61')
      points((1:tmax)+eps, beta.hoff[i, 1:tmax], pch=19,lwd=2,type='p',col='#fdae61')
      segments(x0=(1:tmax)+eps, y0=LB.beta[i, 1:tmax],x1=(1:tmax)+eps,y1=UB.beta[i,1:tmax], pch=19, col='#fdae61', lwd=1.5,lty=1)
      #if(i==8){legend('topleft', legend=c('Hoff', 'OLS', 'GEE Homo E/I', 'GEE Hetero E/I', 'GEE Homo E/E'), cex=.9,
      #pch=c(19, 0, 2, 4, 6), col=c("black", "green", "red", "purple", "orange"), bty='n')}
    }
    
    # EE
    points(1:tmax, beta.EE[i, 1:tmax], pch=18, col='#d7191c', lwd=1.5,type='l',lty=1)
    segments(x0=1:tmax, y0=beta.EE[i, 1:tmax]-(1.96*se.EE[i,1:tmax]),x1=1:tmax,y1=beta.EE[i, 1:tmax]+(1.96*se.EE[i,1:tmax]), pch=18, col='#d7191c', lwd=1.5,lty=1)
    points(1:tmax, beta.EE[i, 1:tmax], pch=18, col='#d7191c', lwd=2,type='p',cex=1.3)
    if(i==7){legend('topright', legend=c(expression(paste("GEE, Figure 6(a) ",Omega)), 'Westveld & Hoff 2011','OLS, w/ DC stand. err.'), cex=1,
                    pch=c(18, 19, 17), col=c('#d7191c', '#fdae61', '#8da0cb'), lty=c(1,1,1),bty='n')}
    
    dev.off()
  }
}


#######################
###  Key Functions  ###
#######################

# Perform fit on data and calculate variances for simulations in paper
fit.data <- function(Y.1, X.1, n.list, type='ols', data='continuous') 
{
  # Fit data and calculate variances (4 different ways)
  # 
  
  d.tot <- length(Y.1)
  n.tot <- floor((1+sqrt(1+4*d.tot))/2)
  p.tot <- ncol(X.1)
  
  out.length <- p.tot + 4*p.tot*(p.tot+1)/2 + 2
  data.vec.out <- rep(0,out.length)    # beta and four vectorized 
  
  XtXinv.1 <- chol2inv(chol(crossprod(X.1,)))      # sandwich estimator "bread"
  h.ii <- rep(0,d.tot)
  for (i in 1:d.tot) {
    h.ii[i] <- X.1[i,] %*% XtXinv.1 %*% (X.1[i,])
  }
  
  
  if (data == 'continuous') {
    # continuous
    if (type == 'ols'){
      # ols
      beta.hat <- XtXinv.1 %*% t(X.1) %*% Y.1
      e.1 <- Y.1 - X.1 %*% beta.hat
      e.BM <- e.1/sqrt(1-h.ii)
      e.HC <- sqrt(as.vector(e.1^2/(1-h.ii)))
      X.e.HC <- X.1*tcrossprod(e.HC,rep(1,p.tot))
      
      V.OLS.tmp <- XtXinv.1 * sum(e.1^2)/(d.tot - p.tot)
      meat.E <- meat.S(n.list, X.1, e.1) 
      # E.sing <- meat.E$sing_flag
      # E.negdef <- meat.E$negdef_flag
      
      V.E.tmp <- sandwich.var(XtXinv.1, meat.E$M)
      meat.DC <- meat.U(n.list, X.1, e.1)
      
      # DC.sing <- meat.DC$sing_flag
      # DC.negdef <- meat.DC$negdef_flag
      V.DC.tmp <- sandwich.var(XtXinv.1, meat.DC$M)
      V.HC.tmp <- sandwich.var(XtXinv.1, crossprod(X.e.HC))
      
      # make positive variances
      V1 <- make.positive.var(V.E.tmp)
      V2 <- make.positive.var(V.DC.tmp)
      
      V.E.tmp <- V1$V
      V.DC.tmp <- V2$V
      
      
    } else {
      
      stop("Only OLS implemented")
      
    }
    
    
    
  }else if (data == 'binary') {
    
    stop("Binary not yet implemented")
    
    
  } else {
    stop('Invalid data type. Choose continuous or binary')
  }
  
  
  return(list(beta = beta.hat, var_HC = V.HC.tmp, var_DC = V.DC.tmp, var_E = V.E.tmp, var_OLS = V.OLS.tmp, flag_E = V1$flag, flag_DC = V2$flag))
              #, E.sing = E.sing, E.negdef=E.negdef,DC.sing=DC.sing,DC.negdef=DC.negdef))
}

# Perform GEE estimate of coefficients as in trade example
GEE.est <- function(Y.in, X.in, n.tot, t.max, tol.in=1e-6, directed=T, type="EE", write_dir=getwd(), beta_start="", node.mat=1) 
{
  # Calculate GEE estimate of regressors beta for continuous data
  # Return beta, residuals, weight matrix, # iterations, objective function Q, and actual tolerance
  # 
  
  
  t1 <- proc.time()[3]
  
  
  if(nchar(beta_start)>0){
    cat("Starting coefficients lcoated here: ", beta_start, "\n")
  }
  
  
  count.in <- 0
  d.tot <- n.tot*(n.tot-1)/2*(1 + directed)
  
  XX <- chol2inv(chol(crossprod(X.in,)))
  
  # Overlapping pair indicators
  S.list <- Sigma.ind(n.tot, directed, sta_flag=F)
  W.ind <- Sigma.ind(n.tot, directed, sta_flag=T)   # 7 parameter version
  ilist <- lapply(S.list, find_ones)
  d.tot <- nrow(S.list[[1]])
  
  # Initial fit
  beta.0 <- XX %*% t(X.in) %*% Y.in
  se.ols <- sqrt(diag(XX * sum((Y.in - X.in %*% beta.0)^2)/(nrow(X.in) - ncol(X.in))))
  write.table(beta.0, file.path(write_dir, "beta_ols.txt"), row.names=F, col.names=F)
  write.table(se.ols, file.path(write_dir, "beta_olsSE.txt"), row.names=F, col.names=F)
  cat("OLS fit complete \n")
  
  eols <- Y.in - X.in %*% beta.0
  meat_DC <- meatDC2(node.mat, X.in, eols)
  Vdc <- XX %*% meat_DC %*% XX
  Vdc <- make.positive.var(Vdc)
  se_dc <- sqrt(diag(Vdc$V))
  write.table(se_dc, file.path(write_dir, paste0("se_dc_t", tmax,".txt")), row.names=F, col.names=F)
  cat("DC SE complete \n")
  
  # Start looping
  if(nchar(beta_start) == 0){
    # if(type != "STA") {
    beta.gls <- beta.0
  } else {
    
    beta.gls <- as.vector(read.table(beta_start, header=F)[,1])
    # # For stationary fit, first do EE fit!!
    # cat('************************************ \n')
    # cat("Fitting GEE with EE/6a before doing STA/6b \n")
    # temp <- GEE.est(Y.in, X.in, n.tot, t.max, tol.in, directed, type="EE", write_dir)
    # beta.gls <- temp$beta
    # rm(temp)
    # gc()
    # cat("Finished preliminary EE/6a fit \n")
    # cat('************************************ \n\n')
  }
  
  e.test <- Y.in - X.in %*% beta.gls
  Q.old <- Q.new <- sum(e.test^2)
  delta.loop <- tol.in*100
  
  out <- data.frame(matrix(c(count.in, Q.new), ncol=2))
  names(out) <-  c("iteration", "criterion")
  write.table(out, file.path(write_dir, "GEE_criterion.txt"), row.names = F, sep="\t", quote=F)
  
  
  
  while (abs(delta.loop) > tol.in) {
    count.in <- count.in + 1
    
    params <- calculate_matrix_params(ilist, e.test, t.max, type)  # unique parameters in Var(Y)
    invParams <- calculate_parameter_inverse(n.tot, t.max, params, type)  # calculate inverse parameters of Var(Y)^{-1}
    
    # Following if-statement calculates GEE criterion, beta, etc. based on type of matrix
    # Avoid building large inverse matrix, W = Var(Y)^{-1} 
    # find bread = (X^T W X)^{-1} and X^T W Y WITHOUT building W
    #
    if(type == "EI"){
      Wmini <- Reduce("+", lapply(1:6, function(x) invParams$p[x]*S.list[[x]]))
      XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,],
                                                              Wmini %*% X.in[1:d.tot + (x-1)*d.tot,]) ))
      bread <- solve(XWX)
      
      XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,],
                                                              Wmini %*% Y.in[1:d.tot + (x-1)*d.tot]) ))
      
      beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      eWe <- lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot],
                                                  Wmini %*% e.test[1:d.tot + (x-1)*d.tot]) )
      Q.new <- Reduce("+", eWe)
      
      
    } else if (type == "EE"){
      A <- Reduce("+", lapply(1:6, function(x) invParams$p[x]*S.list[[x]]))
      B <- Reduce("+", lapply(1:6, function(x) invParams$p[6 + x]*S.list[[x]]))
      
      ABX <- lapply(1:tmax, function(x) A %*% X.in[1:d.tot + (x-1)*d.tot,]
                    + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) X.in[1:d.tot + (t-1)*d.tot,])) )
      XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], ABX[[x]])))
      bread <- solve(XWX)
      rm(ABX)
      
      ABY <- lapply(1:tmax, function(x) A %*% Y.in[1:d.tot + (x-1)*d.tot]
                    + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) Y.in[1:d.tot + (t-1)*d.tot])) )
      XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], ABY[[x]])))
      rm(ABY)
      
      beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      
      ABe <- lapply(1:tmax, function(x) A %*% e.test[1:d.tot + (x-1)*d.tot]
                    + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) e.test[1:d.tot + (t-1)*d.tot])) )
      Q.new <- Reduce("+", lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot], ABe[[x]])))
      rm(ABe)
      gc()
      
    } else if (type == "STA"){
      
      blocklist <- build_blocklist(tmax, type=invParams$type)
      Psi_blocks <- matrix(NA, tmax, tmax)
      for(i in 1:length(blocklist)){
        Psi_blocks[matrix(blocklist[[i]], ncol=2)] <- i
      }
      
      Alist <- lapply(1:length(blocklist), function(y) Reduce("+", lapply(1:7, function(x) invParams$p[x + 7*(y-1)]*W.ind[[x]])))   # list of all blocks
      
      WX <- lapply(1:tmax, function(x) Reduce("+", 
                                              lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% X.in[1:d.tot + (t-1)*d.tot, ])) )   # list entries are block columns in WX
      XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], WX[[x]]))  )
      bread <- solve(XWX)
      rm(WX)
      
      
      WY <- lapply(1:tmax, function(x) Reduce("+", 
                                              lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% Y.in[1:d.tot + (t-1)*d.tot])) )   # list entries are block columns in WX
      XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], WY[[x]])) )
      rm(WY)
      
      beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      
      We <- lapply(1:tmax, function(x) Reduce("+", 
                                              lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% e.test[1:d.tot + (t-1)*d.tot])) )
      Q.new <- Reduce("+", lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot], We[[x]])))
      rm(We)
      gc()
      
    } else {
      stop("Still need to implement smart inversion for HET")
    }
    
    
    
    write.table(cbind(count.in, as.numeric(Q.new)), file.path(write_dir, "GEE_criterion.txt"), row.names = F, sep="\t", append=T, col.names=F)  # write out criterion
    
    cat('GEE iteration: \t\t', count.in, '\n')
    cat('Change in criterion: \t', abs(delta.loop), '\n')
    cat('Elapsed time: \t\t', proc.time()[3]- t1, 'sec \n')
    write.table(beta.gls, file=file.path(write_dir, paste('beta_loop',count.in,'_t', tmax, '.txt',sep='')), row.names=F, col.names=F, sep='\t')
    warnings()
    cat('************************************ \n')
    
    delta.loop <- as.numeric(Q.old - Q.new)
    
    # Update for next loop
    Q.old <- Q.new # new baseline
    
  }
  
  write.table(sqrt(diag(bread)), file=file.path(write_dir, paste('se_t', tmax, '.txt',sep='')), row.names=F, col.names=F, sep='\t')
  
  cat("\n GEE estimation complete \n")
  output <- list(beta.gls, e.test, bread, invParams=invParams$p, count.in, delta.loop)
  names(output) <- c('beta','resid','bread','W', 'n', 'tol')
  return(output)
}




##########################
###  Helper Functions  ###
##########################

# Generate node sets of various overlapping dyad pairs
node.set <- function(n.tot)
{  
  # Generate node sets of various overlapping dyad pairs
  # Return list of each set of nodes, with null for first set (diagonal)
  
  nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                   rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
  
  nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <-  c()
  
  for (i in 1:n.tot){ 
    
    # ij,ji
    if (i<n.tot){
      c1 <- rep(i,(n.tot-i))
      #       c2 <- (1:n.tot)[-i]
      c2 <- ((i+1):n.tot)
      
      nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
    }
    
    # ij,il  ;  ij,kj
    c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
    c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
    c3 <- c1
    c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
    c4 <- c4.mat[lower.tri(diag(n.tot-1))]
    
    nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
    nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
    
    # ij,jl and ij,ki
    nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3), cbind(c2,c1,c3,c4))   
  }
  return(list(n1=nodes.1,n2=nodes.2,n3=nodes.3,n4=nodes.4,n5=nodes.5))
}

# Generate covariate matrix used in simulation
Generate.X.bv <- function(n.in) 
{
  # Generate X matrix of a given size for bias-variance paper
  # 1st column 1's
  # 2nd column random 0,1 with prob = 1/2
  # 3rd column N(0,1)
  d.in <- n.in*(n.in-1)     # number of dyads
  X.out <- matrix(1,d.in,4)
  
  # Generate correlated categorical variable
  sex.i <- sample(c(0,1), size=n.in, replace=T)
  sum.sex <- sum(sex.i)                         
  i.rand <- sample(1:n.in, size=1)
  sex.i[i.rand] <- sex.i[i.rand] - 1*(sum.sex==n.in) + 1*(sum.sex==0)  # control for all 1's or all 0's
  sex.mat <- outer(sex.i,sex.i,"==")*1
  X.out[,2] <- as.vector(sex.mat)[-seq(1,n.in^2,n.in+1)]   # unfold without diagonal entries
  
  # Generate correlated continuous variable
  X.rand <- rnorm(n.in,0,1)
  X.mat <- outer(X.rand,rep(1,n.in)) - outer(rep(1,n.in),X.rand)
  X.out[,3] <- as.vector(X.mat)[-seq(1,n.in^2,n.in+1)]   # unfold without diagonal entries
  
  # Random continuous covariate
  X.out[,4]  <- rnorm(d.in,0,1)
  return(X.out)
}

# Generate random errors and Y for simulations in paper
Generate.data <- function(X, beta, model='iid', data='continuous', f.in=.25, rho=.5, s.e.max=sqrt(3), n.z=2)
{
  # Generate data from E[Y], beta, return Y and error terms
  # Model must be 'iid', 'sr', 'ble', or 'rand'
  # normalize to approximate same total variance for each model of given size
  #
  
  EY.1 <- X %*% beta
  d.tot <- length(EY.1)
  n.tot <- (1+sqrt(1+4*d.tot))/2
  
  # Generate errors 
  if (model == 'iid'){  
    # generate iid errors
    eps.1 <- rnorm(d.tot,0,s.e.max)
    error.1 <- 0
    
  } else if (model == 'sr') {   
    # generate send/recieve errors, a la Aronow
    
    # Calculate desired variances based on SNR
    rho.ab <- rho
    se.2 <- f.in*s.e.max^2
    s.e <- sqrt(se.2)
    psi.ab <- sqrt(2)
    sb.2 <- (s.e.max^2 - se.2)/3
    s.b <- sqrt(sb.2)
    s.a <- psi.ab*s.b
    
    # Complete error term
    eps.1 <- rnorm(d.tot,0,s.e)
    
    sigma.ab <- matrix(c(s.a^2, rho.ab*s.a*s.b, rho.ab*s.a*s.b, s.b^2 ),2,2)
    ab.mat <- mvrnorm(n.tot,c(0,0), sigma.ab)        # generate ab.i
    ab.ij <- tcrossprod(ab.mat[,1],rep(1,n.tot)) + tcrossprod(rep(1,n.tot),ab.mat[,2])
    ab.1 <- as.vector(ab.ij)[-seq(1,n.tot^2,n.tot+1)]   # unfold without diagonal
    
    error.1  <- ab.1
    
    
  } else if (model == 'hsr') {   
    # generate send/recieve errors, a la Aronow
    
    stop("hsr deprecated")
    
    
   
    
  } else if (model == 'ble') { 
    # bilinear effects a la Hoff
    
    # Calculate desired variances based on SNR
    rho.ab <- rho
    se.2 <- f.in*s.e.max^2
    s.e <- sqrt(se.2)
    psi.ab <- sqrt(2)
    psi.z <- 1
    psi.g <- 1
#     sb.2 <- (s.e.max^2 - se.2)/(1 + psi.ab^2 + psi.z^4*n.z + psi.g^2)
#     s.b <- sqrt(sb.2)
    b.e <- (1 + psi.ab^2 + psi.g^2)
    a.e <- n.z*psi.z^4
    c.e <- (s.e^2 - s.e.max^2)
    
    s.b.2 <- (-b.e +sqrt(b.e^2 - 4*a.e*c.e))/2/a.e
    s.b <- sqrt(s.b.2)

    s.a <- psi.ab*s.b
    s.z <- psi.z*s.b
    s.g <- psi.g*s.b
    
    # Generate complete error term
    eps.1 <- rnorm(d.tot,0,s.e)
    
    sigma.ab <- matrix(c(s.a^2, rho.ab*s.a*s.b, rho.ab*s.a*s.b, s.b^2 ),2,2)
    ab.mat <- mvrnorm(n.tot,c(0,0), sigma.ab)        # generate ab.i
    ab.ij <- tcrossprod(ab.mat[,1],rep(1,n.tot)) + tcrossprod(rep(1,n.tot),ab.mat[,2])
    ab.1 <- as.vector(ab.ij)[-seq(1,n.tot^2,n.tot+1)]   # unfold without diagonal
    
    g.ij = matrix(rnorm(n.tot^2,0,s.g),n.tot,n.tot)   # random symmetric effects
    g.ij = forceSymmetric(g.ij)   # force symmetry
    g.1 <- as.vector(g.ij)[-seq(1,n.tot^2,n.tot+1)]   # unfold without diagonal
    
    z.i = mvrnorm(n.tot, rep(0,n.z), s.z^2*diag(n.z))
    z.mat <- tcrossprod(z.i,z.i)
    z.1 <- as.vector(z.mat)[-seq(1,n.tot^2,n.tot+1)]
    
    error.1  <- ab.1 + z.1 + g.1 
    
    
  } else if (model == 'hble') { 
    
    stop("hble deprecated")

    
  } else if (model == 'rand') {  
    # random locations intended to break unstructured estimate
    n.mid <- floor(n.tot/2)      # middle # of dyads   
    rho.ab <- rho
    se.2 <- f.in*s.e.max^2
    s.e <- sqrt(se.2)
    sc.2 <- (s.e.max^2 - se.2)*n.tot/( n.mid  )
    s.c <- sqrt(sc.2)
    e.c.mat <- matrix(0,n.tot,n.tot)
    e.c.mat[1:n.mid,1:n.mid] <- rnorm(1,0,s.c)
    
    e.c.1 <- as.vector(e.c.mat)[-seq(1,n.tot^2,n.tot+1)]   # unfold without diagonal
    
    # Create error terms  
    eps.1 <- rnorm(d.tot,0,s.e)
    error.1 <- e.c.1
    
  } else if (model == 'rand2') { 
    
    stop("rand2 deprecated")
  
    
  } else {
    stop('Invalid model type')
  }
  
  
  if (data == 'continuous'){  
    Y.1 <- EY.1 + error.1 + eps.1  # Y data
    data.out <- data.frame(cbind(Y.1, X, error.1))
    
  } else if (data == 'binary'){
    # p.gen <- 1/(1+exp(-EY.1 - error.1));        # generating probability with errors  
    # Y.1 <- rbinom(d.tot, size=1, prob=p.gen) 
    # data.out <- data.frame(cbind(Y.1,X.1,error.1))
    stop("Binary not yet implemented")
    
  } else {
    stop('Invalid data type. Choose continuous or binary.')
  }
  
  names(data.out) <- c('Y','X1','X2','error')
  
  return(data.out)
}

# Calculate E meat
meat.S <- function(node.list,X.1,e.1)
{
  # Build meat for structured/exchangeable sandwich standard error estimator
  # Takes in design matrix X.1, residuals e.1, and list of overlapping nodes 
  
  # sizes
  d <- nrow(X.1)
  n <- floor((1+sqrt(1+4*d))/2)
  p <- ncol(X.1)
  
  # Estimates
  phi <- rep(0,5)
  phi[1] <- sum(e.1^2)/d
  meat.1 <- phi[1]*crossprod(X.1)
  
  e.mat <- mat.net(e.1)
  meat.2 <- matrix(0,p,p)
  for (i in 2:length(node.list)){
    nodes.tmp <- node.list[[i]]
    phi[i] <- mean(e.mat[ nodes.tmp[,1:2]]*e.mat[nodes.tmp[,3:4]])
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n)
    meat.2 <- meat.2 + crossprod(X.1[d1,], X.1[d2,])*phi[i]  + crossprod(X.1[d2,], X.1[d1,])*phi[i]
  }
  
  # params <- c(phi,0)
  # S.list <- Sigma.ind(n, directed=TRUE)    # list of indicator matrices, 6 parameter  
  # C = Reduce("+", lapply(1:length(S.list), function(x) params[x]*S.list[[x]]))
  # if(rankMatrix(C,method = "qr")[1]<nrow(C)){sing_flag<-TRUE}else{sing_flag<-FALSE}
  # if(det(C)<0 && sing_flag==FALSE){negdef_flag<-TRUE}else{negdef_flag<-FALSE}
  
  meat.out <- meat.2 + meat.1
  # return(list(M=meat.out,sing_flag = sing_flag, negdef_flag = negdef_flag))
  return(list(M=meat.out,sing_flag = NA, negdef_flag = NA))
  
}

# Calculate DC meat
meat.U <- function(node.list,X.1,e.1)
{
  # Build meat for unstructured/DC sandwich standard error estimator
  # Takes in design matrix X.1, residuals e.1, and list of overlapping nodes 
  
  # sizes
  d <- nrow(X.1)
  n <- floor((1+sqrt(1+4*d))/2)
  p <- ncol(X.1)
  
  # Estimates
  X.e <- X.1*tcrossprod(e.1, rep(1,p))
  meat.1 <- crossprod(X.e)
  
  meat.2 <- matrix(0,p,p)
  for (i in 2:length(node.list)){
    nodes.tmp <- node.list[[i]]
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n)
    meat.2 <- meat.2 + crossprod(X.e[d1,], X.e[d2,]) + crossprod(X.e[d2,], X.e[d1,])
  }
  
  # S.list <- Sigma.ind(n, directed=TRUE)    # list of indicator matrices, 6 parameter  
  # C_DC = outer(c(e.1),c(e.1),'*')
  # C_DC <- C_DC*(1-S.list[[6]])
  # if(rankMatrix(C_DC,method = "qr")[1]<nrow(C_DC)){sing_flag<-TRUE}else{sing_flag<-FALSE}
  # if(det(C_DC)<0 && sing_flag==FALSE){negdef_flag<-TRUE}else{negdef_flag<-FALSE}
  
  meat.out <- meat.2 + meat.1
  # return(list(M=meat.out,sing_flag=sing_flag,negdef_flag=negdef_flag))
  return(list(M=meat.out,sing_flag=NA,negdef_flag=NA))
}

# Calculate DC meat in a cleaner manner with less memory overhead
#  use rowrange to parallelize if desired
meatDC2 <- function(nodes,X,e,rowrange=1:nrow(X))
{
  if(nrow(nodes) != nrow(X)){stop("Columns in nodes must match columns in x")}
  if(ncol(nodes) < 2){stop("Nodes must be # dyads x at least 2 cols")}
  nodes <- as.matrix(nodes)
  
  p <- ncol(X)
  Xe <- sweep(X, 1, e, "*")
  meat <- matrix(0,p,p)  # initialize
  
  cat("%-age done with DC meat: \n")
  # may want dimensionality flag here, not sure, i.e. loop rather than lapply
  for(d in rowrange){
    ij <- nodes[d,]  # dyads
    rows <- which(nodes[,1] %in% ij | nodes[,2] %in% ij)    # overlapping dyads with ij
    temp <- lapply(rows, function(x) tcrossprod(Xe[d,], Xe[x,]))
    meat <- meat + Reduce("+", temp)   # increment meat by ij^th "cluster"
    if(d%%1000 == 0){
      cat(100*d/nrow(X), " ")
    }
  }
  
  return(meat)
}

# Replace negative eigenvalues with zeros in variance matrix
make.positive.var <- function(V.test)
{
  # Make variance matrix positive definite by zeroing out negative eigenvalues
  eig.V <- eigen(as.matrix(V.test), symmetric=T)
  eig.flag <- sum(eig.V$values < 0)                # creat flag if correction made
  
  if (eig.flag>0) { 
    eig.V$values[which(eig.V$values<0)] <- 0         # zero out negative eigenvalues 
    V.corrected <- eig.V$vectors %*% diag(eig.V$values) %*% t(eig.V$vectors)
  } else { 
    V.corrected <- V.test
  }
  
  return(list(V=V.corrected, flag=1*(eig.flag>0) ))
}

# Read in trade data for trade example in paper
read_trade_data <- function(tmax, read_dir=getwd())
{
  
  data <- read.csv(file.path(read_dir, 'Trade.csv'), sep=' ')
  data <- data[,-1]   # remove dummy column of row names
  
  
  # Now perform on range of time periods
  data <- data[order(data$t, data$i),]
  i.time <- data$t<=tmax
  
  
  # Set parameters, limiting to max.time
  Y <- data$ltrade[i.time]
  X <- cbind(1, data$lgdp.exp, data$lgdp.imp, data$ldist, data$pty.exp, data$pty.imp, 
             data$cc, data$pty.exp*data$pty.imp)[i.time,]
  node.mat <- cbind(data$i, data$j)[i.time,]
  t.vec <- data$t[i.time]
  
  n.t <- length(unique(t.vec))
  I.time <- t(apply(as.matrix(t.vec), 1, rep, times=n.t))
  I.time <- sweep(I.time,2,1:n.t,'==')       # indicator of time period 1:20
  X.time <- NULL
  for (i in 1:ncol(X)){
    X.tmp <- sweep(I.time,1,X[,i],'*')
    X.time <- cbind(X.time, X.tmp)
  }
  
  rm(X.tmp, I.time)
  
  return(list(X=X.time, Y=Y, node.mat=node.mat, t.vec=t.vec))
}

# Vectorize a network matrix (without diagonal)
vec.net <- function(A, directed=T)
{
  # convert matrix A to vector
  # unfold along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  n <- nrow(A)
  if (directed==T){ 
    vec.out <- as.vector(A)[-seq(1,n^2,n+1)] 
  } else if (directed==F) {
    vec.out <- A[lower.tri(A)]
  } else {stop('Logical input only for directed input')}
  return(vec.out)
}

# Matricize a network vector (without diagonal)
mat.net <- function(V, directed=T)
{
  # convert vector V to a matrix
  # build along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  if (directed==T){   
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+4*d))/2)
    
    A <- matrix(1:n^2,n,n)
    ind <- vec.net(A, directed)            # indices in matrix
    Mat.out <- matrix(0,n,n)
    Mat.out[ind] <- V
  } else if (directed==F){
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+8*d))/2)
    
    Mat.out <- matrix(0,n,n)
    Mat.out[lower.tri(Mat.out)] <- V
    Mat.out <- Mat.out + t(Mat.out)
  } else {stop('Logical only for directed input')}
  return(Mat.out)
}

# Generate node pairs for complete network
node.gen <- function(n, directed=T)
{
  # Generate pairs of nodes for n.code given network size and directed vs. undirected
  # for columnwise unfolding of matrix
  
  if (directed==T){
    mat.out <- c()
    for (i in 1:n){
      mat.out <- rbind(mat.out, cbind((1:n)[-i], rep(i,n-1)))
    }
  } else if (directed==F){
    mat.out <- c()
    for (i in 1:(n-1)){
      mat.out <- rbind(mat.out, cbind( (i+1):n, rep(i,n-i)))
    }
  } else (stop('Choose T or F logical for directed'))
  
  return(mat.out)
}

# Find all possible combinations of elements in two vectors, or all combinations of all elements in one without repeats
combine <- function(v1, v2=NA)
{
  # Find all possible combinations of elements in two vectors
  # No pairs (a,a)
  # return Nx2 matrix of all pairs
  # If v2 = NULL, do all possible combinations of elements in v1
  #
  # Make sure don't have repeats in v1 or v2, not coded to remove repeats
  #
  v1 <- sort(v1)
  l1 <- length(v1)
  
  if(is.na(v2[1])){  # if null v2 or not input
    c1 <- rep(v1, times=(l1-1):0)    # first column
    v.mat <- outer(v1,rep(1,l1))
    c2 <- v.mat[lower.tri(v.mat)]    # second column
    m.out <- cbind(c1,c2) 
  } else {    # all combinations of v1 and v2
    l2 <- length(v2)
    v2 <- sort(v2)
    c1 <- rep(v1, each=l2)    # first column
    c2 <- rep(v2, times=l1)   # second column
    m.out <- cbind(c1,c2)
    remove <- which(m.out[,1]==m.out[,2], arr.ind=T)
    if(length(remove) > 0){
      m.out <- m.out[-remove,]   # remove any duplicate rows
    }
  }
  
  return(m.out)
}

## Dyad map from nodes i,j --> dyad d
dyad <- function(i.in, j.in, n.tot, directed=T) 
{
  # transform n indices to dyad indices
  if (directed==T){ 
    dyad.result <- ((i.in-1) + (j.in-1)*(n.tot-1) + (j.in > i.in)) * (i.in != j.in)
  } else if (directed==F) {
    #     if (j.in>i.in){stop('Only lower triangular matrix for undirected network')}
    ij <- cbind(i.in, j.in)
    ij <- t(apply(ij, 1, sort, decreasing=T))
    j.1 <- ij[,2]    # want invariance to ordering
    i.1 <- ij[,1]    # want invariance to ordering
    dyad.result <- ( i.1 - 1 -.5*j.1^2 + (n.tot - 1/2)*j.1 - n.tot + 1  ) * (i.1 != j.1)
  } else {stop('Logical T/F only for directed input')}
  return(dyad.result)
}

# Multiply bread and meat for sandwich variance
sandwich.var <- function(bread, meat) 
{
  # Calculate estimate of variance using input bread and meat matrices
  V.out <- bread %*% meat %*% bread
  
  return(V.out)
}

# Generate list indicator matrix of overlapping dyads
Sigma.ind <- function(n.tot, directed, sta_flag=F)
{  
  # Generate dyad sets of various overlapping dyad pairs
  # Return list of each set of dyads
  # Map dyads from matricization numbering to actual rows through d.map
  #    may include zeros
  
  
  if (directed==T){  
    nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                     rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
    
    nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <- nodes.6 <- c()
    
    for (i in 1:n.tot){ 
      
      # ij,ji
      if (i<n.tot){
        c1 <- rep(i,(n.tot-i))
        #       c2 <- (1:n.tot)[-i]
        c2 <- ((i+1):n.tot)
        
        nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
      }
      
      # ij,il  ;  ij,kj
      c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
      c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
      c3 <- c1
      c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
      c4 <- c4.mat[lower.tri(diag(n.tot-1))]
      
      nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
      nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
      
      # ij,jl and ij,ki
      nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3), cbind(c2,c1,c3,c4))   
      n.list.tmp <- list(n1=nodes.1, n2=nodes.2, n3=nodes.3, n4=nodes.4, n5=nodes.5)
    }
    
  } else if (directed==F) {
    nodes.1 <- nodes.3 <- nodes.4 <- nodes.5 <- c()
    
    rows <- outer(1:n.tot, rep(1,n.tot))
    cols <- outer(rep(1,n.tot), 1:n.tot)
    r.1 <- rows[lower.tri(rows)]
    c.1 <- cols[lower.tri(cols)]
    
    nodes.1 <- cbind(r.1,c.1,r.1,c.1)    # diagonal/variance
    for (i in 1:n.tot){
      if(i>2){ 
        c.1 <- rep(1:(i-2), (i-2):1)
        c.mat <- outer(1:(i-1), rep(1,(i-1)))
        c.2 <- c.mat[lower.tri(c.mat)]
        nodes.3 <- rbind(nodes.3, cbind(rep(i, length(c.1)),c.1,rep(i, length(c.1)),c.2))
        #         nodes.5 <- rbind(nodes.5, cbind(rep(i, length(c.1)),c.1,rep(i, length(c.1)),c.2))
      }
      if (i < (n.tot - 1)){
        r.1 <- rep( (i + 1):(n.tot-1), (n.tot-i-1):1 ) 
        r.mat <- outer((1+i):n.tot, rep(1,(n.tot - i)))
        r.2 <- r.mat[lower.tri(r.mat)]
        nodes.4 <- rbind(nodes.4, cbind(r.1, rep(i, length(r.1)), r.2, rep(i, length(r.1))))
        nodes.5 <- rbind(nodes.5, cbind(r.1, rep(i, length(r.1)), r.2, r.1))
      }
    }
    n.list.tmp <- list(n1=nodes.1, n2=rbind(nodes.3,nodes.4,nodes.5))
  } else {stop('Logical only for directed input')}
  
  
  # Form dyad list 
  d.tot <- n.tot*(n.tot-1)/2*(1 + directed)
  S.list <- vector('list', length(n.list.tmp))
  d.list <- vector('list', length(n.list.tmp))
  names(S.list) <- paste('S',1:length(n.list.tmp),sep='')
  for (i in 1:length(n.list.tmp)){
    nodes.tmp <- n.list.tmp[[i]]
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n.tot, directed)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n.tot, directed)
    #     mat <- cbind(d1, d2)
    #     mat <- apply(mat,1,sort)
    L1 <- (d1<=d2)
    L2 <- d2<d1
    s1 <- d1*L1 + d2*L2
    s2 <- d1*L2 + d2*L1
    #     S.list[[i]] <- cbind(d1,d2)
    #     S.list[[i]] <- sparseMatrix(i=mat[,1], j=mat[,2], x=1, dims=c(d.tot, d.tot))  #, symmetric=T)
    S.list[[i]] <- sparseMatrix(i=s1, j=s2, x=1, dims=c(d.tot, d.tot), symmetric=T)
    d.list[[i]] <- cbind(s1,s2)
  }
  
  if(sta_flag == T & directed==T){   # change listing of 5th overlap 
    tmp <- matrix(1, n.tot, n.tot)
    diag(tmp) <- 0
    nodemat <- which(tmp == 1, arr.ind=T)   # columnwise unfolding dyad list
    i1 <- which(S.list[[5]]==1, arr.ind=T)
    nodes.5 <- cbind(nodemat[i1[,1], ], nodemat[i1[,2], ])
    i5 <- which(nodes.5[,1] == nodes.5[,4])
    S.list[[5]] <- S.list[[6]] <- sparseMatrix(i=1, j=1, x=0, dims=c(d.tot, d.tot), symmetric=T)
    S.list[[5]][i1[i5, ]] <- 1
    S.list[[6]][i1[-i5, ]] <- 1
    
  } else if (sta_flag == T & directed==F){
    stop("stationary listing not impememented for undirected")
  }
  
  S.list[[length(S.list) + 1]] <- 1*(Reduce("+", S.list) == 0)   # all the rest (zeros)
  
  return(S.list)
}

# Build intermediate C(phi,n) matrix in inversion of Exchangeable variance matrix
build_phi_matrix<- function(n, phi, sta_flag=F)
{
  # Build intermediate phi/p matrix in inversion of E/E matrix
  # A(n, phi) %*% p = c(1,rep(0, 5))
  #

  if(sta_flag == F){
    A <- matrix(NA, 6, 6)
    A[1,] <- c(phi[1], phi[2], (n-2)*phi[3], (n-2)*phi[4], 2*(n-2)*phi[5], (n-3)*(n-2)*phi[6])
    A[2,] <- c(phi[2], phi[1], (n-2)*phi[5], (n-2)*phi[5], (n-2)*phi[3] + (n-2)*phi[4], (n-3)*(n-2)*phi[6])
    A[3,] <- c(phi[3], phi[5], phi[1] + (n-3)*phi[3], phi[5] + (n-3)*phi[6], phi[2] + phi[4] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[4,] <- c(phi[4], phi[5], phi[5] + (n-3)*phi[6], phi[1] + (n-3)*phi[4], phi[2] + phi[3] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
    A[5,] <- c(phi[5], phi[4], phi[2] + (n-3)*phi[5], phi[3] + (n-3)*phi[6], phi[1] + phi[5] + (n-3)*phi[4] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
      # A[5,] <- c(phi[5], phi[3], phi[4] + (n-3)*phi[6], phi[2] + (n-3)*phi[5], phi[1] + phi[5] + (n-3)*phi[3] + (n-3)*phi[6],
                 # (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[6,] <- c(phi[6], phi[6], phi[4] + phi[5] + (n-4)*phi[6], phi[3] + phi[5] + (n-4)*phi[6], phi[3] + phi[4] + 2*phi[5] + 2*(n-4)*phi[6],
               phi[1] + phi[2] + (n-4)*(phi[3] + phi[4] + 2*phi[5] + (n-5)*phi[6]))

  } else {
    A <- matrix(NA, 7, 7)
    A[1,] <- c(phi[1], phi[2], (n-2)*phi[3], (n-2)*phi[4], (n-2)*phi[6], (n-2)*phi[5], (n-3)*(n-2)*phi[7])
    
    A[2,] <- c(phi[2], phi[1], (n-2)*phi[6], (n-2)*phi[5], (n-2)*phi[3], (n-2)*phi[4], (n-3)*(n-2)*phi[7])
    
    A[3,] <- c(phi[3],   phi[5],   phi[1] + (n-3)*phi[3], 
               phi[6] + (n-3)*phi[7],   phi[4] + (n-3)*phi[7],   phi[2] + (n-3)*phi[5],
               (n-3)*(phi[4] + phi[6] + (n-4)*phi[7]))
    
    A[4,] <- c(phi[4],   phi[6],   phi[5] + (n-3)*phi[7], 
               phi[1] + (n-3)*phi[4],   phi[2] + (n-3)*phi[6],   phi[3] + (n-3)*phi[7],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[7]))
    
    A[5,] <- c(phi[5],   phi[3],  phi[4] + (n-3)*phi[7], 
               phi[2] + (n-3)*phi[5],   phi[1] + (n-3)*phi[3],   phi[6] + (n-3)*phi[7],
               (n-3)*(phi[4] + phi[6] + (n-4)*phi[7]))
    
    A[6,] <- c(phi[6],   phi[4],   phi[2] + (n-3)*phi[6], 
               phi[3] + (n-3)*phi[7],   phi[5] + (n-3)*phi[7],   phi[1] +  (n-3)*phi[4],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[7]))
    
    A[7,] <- c(phi[7],   phi[7],   phi[4] + phi[5] + (n-4)*phi[7], 
               phi[3] + phi[6] + (n-4)*phi[7],    phi[4] + phi[5] + (n-4)*phi[7],    phi[3] + phi[6] + (n-4)*phi[7],
               phi[1] + phi[2] + (n-4)*(phi[3] + phi[4] + phi[5] + phi[6] + (n-5)*phi[7]))
  }
  return(A)
}

# Build intermediate C(phi,n) matrix in inversion of Exchangeable variance matrix, all 7 equations for six-parameter phi
build_phi_matrix76 <- function(n, phi)
{
  # Build intermediate phi/p matrix in inversion of E/E matrix
  # A(n, phi) %*% p = c(1,rep(0, 5))
  #
  
    A <- matrix(NA, 7, 6)
    A[1,] <- c(phi[1], phi[2], (n-2)*phi[3], (n-2)*phi[4], 2*(n-2)*phi[5], (n-3)*(n-2)*phi[6])
    A[2,] <- c(phi[2], phi[1], (n-2)*phi[5], (n-2)*phi[5], (n-2)*phi[3] + (n-2)*phi[4], (n-3)*(n-2)*phi[6])
    A[3,] <- c(phi[3], phi[5], phi[1] + (n-3)*phi[3], phi[5] + (n-3)*phi[6], phi[2] + phi[4] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[4,] <- c(phi[4], phi[5], phi[5] + (n-3)*phi[6], phi[1] + (n-3)*phi[4], phi[2] + phi[3] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
    A[6,] <- c(phi[5], phi[4], phi[2] + (n-3)*phi[5], phi[3] + (n-3)*phi[6], phi[1] + phi[5] + (n-3)*phi[4] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
    A[5,] <- c(phi[5], phi[3], phi[4] + (n-3)*phi[6], phi[2] + (n-3)*phi[5], phi[1] + phi[5] + (n-3)*phi[3] + (n-3)*phi[6],
                (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[7,] <- c(phi[6], phi[6], phi[4] + phi[5] + (n-4)*phi[6], phi[3] + phi[5] + (n-4)*phi[6], phi[3] + phi[4] + 2*phi[5] + 2*(n-4)*phi[6],
               phi[1] + phi[2] + (n-4)*(phi[3] + phi[4] + 2*phi[5] + (n-5)*phi[6]))

  return(A)
}

# find entries of ones in a symmetric matrix, including diagonal
find_ones <- function(x, k=1)
{   # function to find entries of ones in a symmetric matrix, including diagonal
  # x <- matrix(x, nrow=nrow(x))
  # x[upper.tri(x, diag=F)] <- 0
  return(which(x == k, arr.ind=T))
}

param_est_single_ilist <- function(ilist, blockmat, e, d.tot)
{
  # Given matrix of time blocks and a particular exchangeable parameter set (within each block),
  # calculate a single parameter/phi. ASSUMES NO MISSING DATA
  #
  # e, vector double, residuals
  # ilist, matrix of 2 columns, list of i^th parameter entries (within block)
  # blockmat, matrix of 2 columns, time blocks to consider
  # 
  
  partial_sum <- sapply(1:nrow(blockmat), function(x) 
    sum(e[(blockmat[x,1]-1)*d.tot + ilist[,1]]*e[(blockmat[x,2]-1)*d.tot + ilist[,2]]) )
  
  return(sum(partial_sum)/(nrow(ilist)*nrow(blockmat)))
}

# calculate parameter estimates for different types of matrices, i.e. 6a, 6b, etc.
calculate_matrix_params <- function(ilist, e, tmax, type="EI")
{ # calculate parameter estimates for different types of matrices
  # EI is 6a with \Omega_2 = 0
  # EE is 6a
  # STA is 6b
  # HET is 6c
  # 
  d.tot <- length(e)/tmax
  
  ilist <- ilist[1:(length(ilist) - 1)]   # remove zero-entry list
  blocklist <- build_blocklist(tmax, type)
  
  # Calculate parameters for first 5 entries or first 2 (directed vs. undirected)
  params <- c(sapply(blocklist, function(x) sapply(ilist, param_est_single_ilist, x, e, d.tot)))
  
  # Add in zeros for 6th parameters
  param_mat <- matrix(params, nrow=length(ilist))
  param_mat <- rbind(param_mat, 0)
  params <- c(param_mat)
  
  return(params)
}

# build an exchangeable matrix given the parameter vector
build_matrix_from_params_kronecker <- function(S.list, params, n, tmax, type="EI")
{  # build a matrix given the parameter vector
  # EI is 6a with \Omega_2 = 0
  # EE is 6a
  # STA is 6b
  # HET is 6c
  # 
  
  # S.list <- S.list[1:(length(ilist) - 1)]   # remove zero-entry list
  
  d <- n*(n-1)
  dt <- tmax*d
  if(params[6] == 0  | type=="EI"){
    V <- sparseMatrix(i=1,j=1,x=0, dims=c(dt, dt))    # initialize final matrix
    inmat <- sparseMatrix(i=1,j=1,x=0, dims=c(d, d))
  } else {
    V <- matrix(0,dt,dt)    # initialize final matrix
    inmat <- matrix(0,d,d)
  }

  if(type == "STAINV" | type == "HET" | type == "HET2"){
    pmax <- 7
    S.list <- Sigma.ind(n, directed = T, sta_flag = T)
  } else {pmax <- 6}
  
  blocklist <- build_blocklist(tmax, type)  # find time blocks
  if(length(params) != pmax*length(blocklist)){  stop("Build_matrix: Wrong length of parameter vector for matrix type")}
  
  for(b in 1:length(blocklist)){
    pin <- params[1:pmax + pmax*(b - 1)]
    inmat <- Reduce("+", lapply(1:length(S.list), function(x) pin[x]*S.list[[x]]))
    
    # Kronecker things into place
    Kmat <- matrix(0, tmax, tmax)
    Kmat[matrix(blocklist[[b]], ncol=2)] <- 1
    if(type != "HET2"){
      Kmat[matrix(blocklist[[b]][,c(2,1)], ncol=2)] <- 1   # symmetric
    }
    # print(Kmat)
    V <- V + kronecker(Kmat, inmat)
  
  # sapply(1:nrow(y), function(row) V[(y[row, 1]-1)*d + 1:d, (y[row, 2]-1)*d + 1:d] <- inmat)
  }
  
  return(V)
}

# build list of time blocks that are correlated based on the maximum time and type of temporal model
build_blocklist <- function(tmax, type="EI")
{  # build list of time blocks that are correlated based on the maximum time and type of temporal model
  # EI is 6a with \Omega_2 = 0
  # EE is 6a
  # STA is 6b
  # HET is 6c
  # 
  if(type == "EI"){
    blocklist <- vector("list", 1)    # EE matrix
    blocklist[[1]] <- cbind(1:tmax, 1:tmax)   # on-diagonal
  } else if (type == "EE"){
    blocklist <- vector("list", 2)    # EE matrix
    blocklist[[1]] <- cbind(1:tmax, 1:tmax)   # on-diagonal
    blocklist[[2]] <- combine(1:tmax)  # off-diagonal
  } else if (type == "STA"){
    blocklist <- lapply(1:tmax, function(x) cbind(1:(tmax - x + 1), x:tmax))   # stationary block mat
  } else if (type == "STAINV") {
    Psi_blocks <- diag(c(1:(tmax/2),ceiling(tmax/2):1))  # block matrix of Omega^{-1} for STA
    Psi_blocks[lower.tri(Psi_blocks)] <- max(Psi_blocks) +1:(tmax*(tmax-1)/2)
    Psi_blocks <- Psi_blocks + t(Psi_blocks*lower.tri(Psi_blocks, diag=F))
    blocklist <- lapply(1:max(Psi_blocks), function(x) which(Psi_blocks == x, arr.ind=T))
  } else if (type == "HET") {
    blockmat <- rbind(cbind(1:tmax, 1:tmax), combine(1:tmax))
    blocklist <- lapply(1:nrow(blockmat), function(x) matrix(blockmat[x, ], ncol=2))   # same as stationary
  } else if (type == "HET2") {
    blockmat <- rbind(cbind(1:tmax, 1:tmax), combine(1:tmax), combine(1:tmax)[,c(2,1)])
    blocklist <- lapply(1:nrow(blockmat), function(x) matrix(blockmat[x, ], ncol=2))   # same as stationary
  } else {
    stop("Invalid type")
  }
  
  return(blocklist)
}

# Invert matrix parameters based on inputs.
calculate_parameter_inverse <- function(n, tmax, params, type="EI")
{
  # Invert matrix parameters based on inputs. 
  #
  # Returns: 
  #   p, vector double, is parameters of inverse matrix, length 6*# unique blocks
  #   type, the type of matrix that should be built using build_matrix function      
  #
  # Args:
  #   n, integer, is network size (nodes)
  #   t, integer, is number of ``time'' periods
  #   phi, vector double, is parameters STA matrix to invert.
  
  if(tmax==1 | type=="EI"){ 
    A <- build_phi_matrix(n, params[1:6])
    return(list(p = solve(A, c(1, rep(0,5))), type="EI"))
  }    # single time period solution
  
  else if (type == "EE" & tmax!=1) {
    A <- build_phi_matrix(n, params[1:6])
    B <- build_phi_matrix(n, params[7:12])
    C <- matrix(NA, 12, 12)
    C[1:6, 1:6] <- A
    C[1:6, 7:12] <- (tmax-1)*B
    C[7:12, 1:6] <- B
    C[7:12, 7:12] <- A + (tmax-2)*B
    
    return(list(p=solve(C, c(1, rep(0,11))), type="EE"))
    
  } else if (type == "STA" & tmax!=1) {  
    # Stationary, needs 7-parameter solution
    nt <- tmax*(tmax - 1)/2 + ceiling(tmax/2)   # number of unique blocks in inverse
    
    numcol <- nt*7 - tmax   # number of unique parameters in inverse
    
    # Stationary, needs 7-parameter solution
    pmat <- matrix(params, 6, tmax)   # blow up parameters to 7 of them
    pmat <- rbind(pmat, pmat[6,])
    pmat[6,] <- pmat[5,]
    
    blocklist <- build_blocklist(tmax, "HET2")
    nblocks <- length(blocklist)
    blockmat <- matrix(unlist(blocklist), ncol=2, byrow=2)
    
    Cmat <- matrix(0, 7*nblocks, 7*nblocks)
    Dvec <- rep(0, 7*nblocks)
    
    Alist7 <- lapply(1:tmax, function(x) build_phi_matrix(n, as.vector(pmat[,x]), sta_flag=T))   # find all phi matrices
    
    Oblocks <- diag(tmax)
    Oblocks <- abs(row(Oblocks) - col(Oblocks)) + 1  # block matrix of Omega's for STA
    
    Psi_blocks <- diag(1:tmax)  # block matrix of Omega^{-1} for STA
    Psi_blocks[blockmat] <- 1:nblocks
    
    rc <- blockmat  # each row is row in O and column in Psi to multiply and get an equation
    
    # pair_mat <- NULL
    for(i in 1:nrow(rc)){
      x <- rc[i,1]
      y <- rc[i,2]
      Dvec[1 + (i-1)*7] <- x == y   # if along diagonal, should be Identity, o.w. zero
      
      for(j in 1:tmax){
        
        Cmat[(i - 1)*7 + 1:7, (Psi_blocks[j, y]-1)*7 + 1:7] <- Alist7[[ Oblocks[x, j] ]]
        
      }
    }
    Dvec <- Dvec*1
    
    
    return(list(p=solve(Cmat, Dvec), type="HET2"))

  } else if (type == "HET" & tmax!=1) {
      stop("HET inverse not yet implemented")
    
  }
  
    
}


