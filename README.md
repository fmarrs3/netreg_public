# netreg_public
This directory contains code to reproduce the simulation and trade data figures from ``Regression of exchangeable relational arrays’’ by Frank W. Marrs, Bailey K. Fosdick, and Tyler H. McCormick.

## Original Code
`reproduce_simulations.R`  -  Reproduce plots and data from simulation study. 

`reproduce_trade_example.R`  -  Reproduce plots and data from trade example. 

`reproduce_testing.R`  -  Reproduce testing simulations from supplementary material. 

`function_file.R`  -  The above scripts call this file with supporting functions.

`bayes_mcmc_function.R`  - Function wrapped around Bayesian MCMC code from Westveld and Hoff (2011), see link below.


## Supporting data from Westveld and Hoff

Before fitting the trade data or running Westveld and Hoff’s files, it was necessary to make slight modifications to the country and column labels in the data file. Downloading and modifying the source data is necessary to reproduce the trade example plots. 

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
 
4) The column names should be separated by spaces.  This can be fixed by performing the following find and replace operation.
 - Find: “i.j.t.exp.imp.ltrade.lgdp.exp.lgdp.imp.ldist.pty.exp.pty.imp.cc”;
    Replace: “a i j t exp imp ltrade lgdp.exp lgdp.imp ldist pty.exp pty.imp cc”

5) Place the modified Trade.csv file into the directory with the code

