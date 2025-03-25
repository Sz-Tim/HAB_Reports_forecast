# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Run operational scripts



# setup -------------------------------------------------------------------
library(tidyverse)
logfile <- paste0("./out/logs/daily/refresh_", format(today(), "%F"), ".log")
cat("Starting refresh at", format(now(), "%F %T"), "\n", file=logfile)


# update datasets ---------------------------------------------------------
cat("  Running ./code/4a_opDataUpdates.R", "\n", file=logfile, append=T)
source("./code/4a_opDataUpdates.R")


# generate candidate predictions ------------------------------------------
cat("  Running ./code/4b_opPredCandidates.R", "\n", file=logfile, append=T)
source("./code/4b_opPredCandidates.R")


# generate ensemble predictions -------------------------------------------
cat("  Running ./code/4c_opPredEnsemble.R", "\n", file=logfile, append=T)
source("./code/4c_opPredEnsemble.R")


# calculate updated validation metrics ------------------------------------
cat("  Running ./code/4d_opValidation.R", "\n", file=logfile, append=T)
source("./code/4d_opValidation.R")


# create new plots and tables ---------------------------------------------
cat("  Running ./code/4e_opViz.R", "\n", file=logfile, append=T)
source("./code/4e_opViz.R")


# publish to HABreports ---------------------------------------------------
cat("  Running ./code/4f_opPublish.R", "\n", file=logfile, append=T)
source("./code/4f_opPublish.R")


# success -----------------------------------------------------------------
cat("Finished refresh at", format(now(), "%F %T"), "\n", file=logfile, append=T)
