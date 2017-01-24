# netreg_public
This directory contains code to reproduce the simulation and trade data figures from ``Standard errors for regression on relational data with exchangeable errors’’ by Tyler H. McCormick, Bailey K. Fosdick, and Frank W. Marrs. 

## Original Code
`reproduce_simulations.R`  -  Run this to reproduce plots and data from simulation study. 

`reproduce_trade_example.R`  -  Run this to reproduce plots and data from trade example. 

`function_file.R`  -  The above scripts call this file with supporting functions.

`beta_processing.R`  -  Script used to process the results of Westveld and Hoff’s code.


## Supporting data from Westveld and Hoff

Before fitting the trade data or running Westveld and Hoff’s files, it was necessary to make slight modifications to the country and column labels in the data file. Downloading and modifying the source data is necessary to reproduce the trade example plots. To add Westveld and Hoff's estimates to the trade example plots, it is necessary to run their code and save the results. 

#### Preparing the data
1) Download the trade data and R scripts here:  https://projecteuclid.org/euclid.aoas/1310562208#supplemental 

2) Open `Trade.csv` using Microsoft Excel or another file editor.  Although this is
labeled as a comma separated file (csv), it is space separated.

3) Spaces in the country names causes problems when loading the file in R so 
perform the following find and replace operations.
 - Find: “Costa Rica”; Replace: “CostaRica”
 - Find: “El Salvador”; Replace: “ElSalvador”
 - Find: “New Zealand”; Replace: “NewZealand”
 - Find: “United Kingdom”; Replace: “UnitedKingdom”
 - Find: “Trinidad and Tobago”; Replace: “TrinidadandTobago”
 - Find: “Egypt, Arab Rep.”; Replace: “ArabRepEgypt”
 - Find: “Korea, Rep.”; Replace: “RepKorea”
 - Find: “United States”; Replace: “UnitedStates”
 
In addition, the column names must be separated by spaces.  This can be fixed by performing
the following find and replace operation.
 - Find: “i.j.t.exp.imp.ltrade.lgdp.exp.lgdp.imp.ldist.pty.exp.pty.imp.cc”;
    Replace: “a i j t exp imp ltrade lgdp.exp lgdp.imp ldist pty.exp pty.imp cc”

#### Running Westveld and Hoff's code
1) Run script `LSRTradeGRW.R`. We found line 25 that reads

`data <- read.table("Trade.txt”)`

needed to be replaced by

`data <- read.table(“Trade.csv”,header=TRUE,row.names=1)`

to ensure the script ran appropriately.  

2) The posterior samples generated were post-processed using the script `beta_processing.R`,
which creates an Rdata file containing the posterior mean of the regression coefficients and lower and upper bounds for a 95% credible interval.
