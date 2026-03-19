# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Train candidate models


# Sys.sleep(3600*0.2)
# setup -------------------------------------------------------------------
library(tidyverse)
library(glue)
library(tidymodels)
library(nnet)
library(randomForest)
library(glmnet)
library(baguette)
library(earth)
library(bonsai)
library(lightgbm)
library(brms)
library(bayesian)
library(future)
library(butcher)
library(habforecastr)

options(future.globals.maxSize=5000*1024^2)

write_covSet_seeds <- F
n_covSets <- 15 # number of covariate subsets
prop_covs <- 0.25 # proportion of covariates included per subset
train_prop <- 0.75 # proportion of years, withholding rest to train ensemble

# read target info
target_sets <- c("hab", "tox", "habfish")[-3]
targ_exclude <- c("AZP", "YTX")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

# generate covariate set info
if(write_covSet_seeds) {
  set.seed(111)
  covSet.df <- expand_grid(y=targ_i$abbr,
                           id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
    mutate(seed=sample(1000, n()),
           prop_covs=prop_covs) |>
    arrange(y, id)
  write_csv(covSet.df, "data/covSet_df.csv")  
} else {
  covSet.df <- read_csv("data/covSet_df.csv")
}

candidates <- c("Ridge", 
                "MARS", 
                "NN",
                "RF",
                "lgbm",
                "HB")

for(m in (1:length(candidates))) {
  for(i in (1:nrow(covSet.df))) {
    try({
      fit_covSet(y_i=targ_i,
                 run_type="0_init",
                 covSet=covSet.df[i,],
                 mod=candidates[m],
                 train_prop=train_prop, 
                 nTuneVal=ifelse(candidates[m] %in% c("ENet", "Ridge", "MARS"), 1e3, 1e2),
                 ncores=30,
                 responses=c(alert="alert"))
    })
  }
}

