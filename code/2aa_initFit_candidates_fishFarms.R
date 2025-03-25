# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Train candidate models: Fish farms



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

y_i <- read_csv("data/i_fish.csv", show_col_types=F) |> 
  arrange(abbr) |> mutate(type="hab-fish") #|>
  # filter(abbr %in% c("Alex", 
  #                    # "Assp", # n(A1) == 1
  #                    "Cesp", 
  #                    "Chco", 
  #                    # "Chhsp", 
  #                    "Gysp", 
  #                    "Kasp")) 
  #                    # "Ppsp" 
  #                    # "Rhsp"))

n_covSets <- 15 # number of covariate subsets
prop_covs <- 0.3 # proportion of covariates included per subset
set.seed(1)
covSet.df <- expand_grid(y=y_i$abbr,
                         id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
  mutate(seed=sample(1000, n()),
         prop_covs=prop_covs) |>
  arrange(y, id)
write_csv(covSet.df, "data/covSet_fish.csv")

candidates <- c("Ridge", 
                "MARS", 
                "NN",
                "RF",
                "Boost",
                "lgbm",
                "HB")[c(4)]



# fit models --------------------------------------------------------------

covSet.df <- covSet.df |> arrange(id, y)
# covSet.df <- covSet.df |> filter(id=="d01")
for(m in 1:length(candidates)) {
  for(i in 1:nrow(covSet.df)) {
    fit_covSet(y_i=y_i,
               run_type="0_init_smote_10",
               covSet=covSet.df[i,],
               mod=candidates[m],
               train_prop=0.75,
               nTuneVal=ifelse(candidates[m] %in% c("Ridge", "MARS"), 1e3, 1e2),
               ncores=20,
               responses=c(alert="alert"),
               rebalance_thresh=0.1)
  }
}

