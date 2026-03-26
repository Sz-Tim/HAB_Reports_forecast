# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Run operational scripts



# setup -------------------------------------------------------------------
library(tidyverse)
logfile <- paste0("./out/logs/daily/refresh_", format(today(), "%F"), ".log")

tryCatch(
  {
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
    
    # prepare objects for shiny app -------------------------------------------
    cat("  Running ./code/4d_opShinyPrep.R", "\n", file=logfile, append=T)
    source("./code/4d_opShinyPrep.R")
    
    # success -----------------------------------------------------------------
    cat("Finished refresh at", format(now(), "%F %T"), "\n", file=logfile, append=T)
    
  }, error=function(err.msg) {
    write(toString(err.msg), logfile, append=T)
  }
)

