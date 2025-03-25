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


y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="tox")) |>
  filter(! abbr %in% c("AZP", "YTX", "Prli"))


set.seed(1)
covSet.df <- expand_grid(y=y_i$abbr,
                         id=paste0("d", str_pad(1:n_covSets, 2, "left", "0"))) |>
  mutate(seed=sample(1000, n()),
         prop_covs=prop_covs) |>
  arrange(y, id)
write_csv(covSet.df, "data/covSet_hab_tox.csv")

candidates <- c("Ridge", 
                "MARS", 
                "NN",
                "RF",
                "Boost",
                "lgbm",
                "HB")[7]

todo_df <- list("Psse"=c(11, 12, 13), 
                "PSP"=c(5, 7, 8, 9, 10, 11, 12, 13), 
                "DSP"=c(4, 5, 7, 8, 9, 10, 11, 12, 13, 14), 
                "Pssp"=c(3, 10, 12)) |>
  imap_dfr(~paste0("d", str_pad(.x, 2, 'left', '0')) |>
             as_tibble() |>
             mutate(y=.y) |> 
             rename(id=value))
# covSet.df <- filter(covSet.df, y %in% y_i$abbr[c(1, 3, 6, 7)])
# covSet.df <- filter(covSet.df, y %in% y_i$abbr[6])
covSet.df <- inner_join(covSet.df, todo_df[24,], by=c("y", "id"))
for(m in 1:length(candidates)) {
  for(i in rev(1:nrow(covSet.df))) {
    # if(m==1 & i < 68) next
    fit_covSet(y_i=y_i,
               run_type="0_init",
               covSet=covSet.df[i,],
               mod=candidates[m],
               train_prop=train_prop, 
               nTuneVal=ifelse(candidates[m] %in% c("Ridge", "MARS"), 1e3, 1e2),
               ncores=20,
               responses=c(alert="alert"))
  }
}

