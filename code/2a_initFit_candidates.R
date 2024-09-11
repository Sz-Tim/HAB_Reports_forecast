# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Model periodic re-fits


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
walk(dir("code/fn", ".R", full.names=T), source)

y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="tox")) |>
  filter(! abbr %in% c("AZP", "YTX", "Prli"))

covSet.df <- expand_grid(y=y_i$abbr,
                         Avg=c(0,1), 
                         Xf=c(0,1),
                         XN=c(0,1),
                         Del=c(0,1)) |>
  group_by(y) |>
  mutate(id=paste0("d", str_pad(row_number(), 2, "left", "0")),
         f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}")) |>
  ungroup() |>
  arrange(y, id) 

candidates <- c("Ridge", 
                "MARS", 
                "NN",
                "RF",
                "Boost",
                "lgbm",
                "HBL")
m <- 1 # model type to run


for(i in 1:nrow(covSet.df)) {
  fit_covSet(y_i=y_i, 
             run_type="0_init",
             covSet=covSet.df[i,], 
             mod=candidates[m], 
             test_prop=0.75, 
             nTuneVal=100, 
             ncores=40, 
             responses=c(alert="alert"))
}

