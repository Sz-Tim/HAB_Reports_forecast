# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Train candidate models


# Sys.sleep(3600*24)
# setup -------------------------------------------------------------------
library(tidyverse)
library(glue)
library(tidymodels)
library(nnet)
library(randomForest)
library(glmnet)
library(xgboost)
library(earth)
library(bonsai)
library(lightgbm)
library(brms)
library(bayesian)
library(future)
library(butcher)
library(habforecastr)

options(future.globals.maxSize=5000*1024^2)


n_covSets <- 15 # number of covariate subsets
prop_covs <- 0.25 # proportion of covariates included per subset
train_prop <- 0.75 # proportion of years, withholding rest to train ensemble


# Hold over: to be replaced in future
use_og_seeds <- TRUE
if(use_og_seeds) {
  y_i_HABTOX <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="hab"),
                   read_csv("data/i_tox.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="tox")) |>
    filter(! abbr %in% c("AZP", "YTX", "Prli"))
  y_i_HABFISH <- read_csv("data/i_fish.csv", show_col_types=F) |> 
    arrange(abbr) |> mutate(type="hab-fish")
  
  
  set.seed(1)
  covSet.df_HABTOX <- expand_grid(y=y_i_HABTOX$abbr,
                           id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
    mutate(seed=sample(1000, n()),
           prop_covs=prop_covs) |>
    arrange(y, id)
  
  set.seed(1)
  covSet.df_HABFISH <- expand_grid(y=y_i_HABFISH$abbr,
                           id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
    mutate(seed=sample(1000, n()),
           prop_covs=prop_covs) |>
    arrange(y, id)
  
  targ_i <- bind_rows(y_i_HABTOX, y_i_HABFISH) |>
    arrange(type, abbr)
  covSet.df <- bind_rows(covSet.df_HABTOX, covSet.df_HABFISH) |>
    arrange(y, id)
  write_csv(covSet.df, "data/covSet_df.csv")  
} else {
  # read target info
  target_sets <- c("hab", "tox", "habfish")
  targ_exclude <- c("AZP", "YTX", "Prli")
  targ_i <- map_dfr(target_sets, 
                    ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                      mutate(type=.x)) |>
    filter(! abbr %in% targ_exclude) |>
    arrange(type, abbr)
  
  # generate covariate set info
  set.seed(1)
  covSet.df <- expand_grid(y=targ_i$abbr,
                           id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
    mutate(seed=sample(1000, n()),
           prop_covs=prop_covs) |>
    arrange(y, id)
  write_csv(covSet.df, "data/covSet_df.csv")  
}



candidates <- c("Ridge", 
                "MARS", 
                "NN",
                "RF",
                "Boost",
                "lgbm",
                "HB")

for(m in 1:length(candidates)) {
  for(i in 1:nrow(covSet.df)) {
    fit_covSet(y_i=targ_i,
               run_type="0_init",
               covSet=covSet.df[i,],
               mod=candidates[m],
               train_prop=train_prop, 
               nTuneVal=ifelse(candidates[m] %in% c("Ridge", "MARS"), 1e3, 1e2),
               ncores=20,
               responses=c(alert="alert"))
  }
}

