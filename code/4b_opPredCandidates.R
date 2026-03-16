# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update candidate forecasts




# setup -------------------------------------------------------------------
library(tidyverse)
library(glue)
library(tidymodels)
library(nnet)
library(randomForest)
library(glmnet)
library(earth)
library(bonsai)
library(lightgbm)
library(brms)
library(bayesian)
library(future)
library(doFuture)
library(habforecastr)
library(butcher)


set.seed(1)
ncores <- 40
responses <- c(alert="alert")

target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

covSet.df <- read_csv("data/covSet_df.csv")



# Model predictions -------------------------------------------------------

if(.Platform$OS.type=="unix") {
  plan(multicore, workers=ncores)
} else {
  plan(multisession, workers=ncores)
}

foreach(i=1:nrow(covSet.df), .options.future=list(seed=TRUE), .errorhandling="remove") %dofuture% {
  # covariate set / response info
  id <- covSet.df$id[i]
  y.i <- covSet.df$y[i]
  y_i.i <- targ_i |> filter(abbr==y.i)
  
  cat(y.i, id, file=glue("out/logs/fcst/{y.i}_{id}.log"))
  
  # directories
  data.dir <- glue("data/2_new/compiled/")
  base.dir <- glue("out/0_init/")
  fit.dir <- glue("{base.dir}/model_fits/{id}/")
  out.dir <- glue("out/1_forecast/compiled/")
  dir.create(glue("{out.dir}/{id}/"), showWarnings=F, recursive=T)
    
  # load datasets
  d.y <- dirf(data.dir, glue("{y.i}_{id}_dy")) |> last() |> readRDS()
  dPCA.y <- dirf(data.dir, glue("{y.i}_{id}_dPCAy")) |> last() |> readRDS()
    
  # generate all forecast predictions
  fcst.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
  saveRDS(fcst.ls, glue("{out.dir}/{id}/{y.i}_fcst_{format(max(fcst.ls$alert$date), '%F')}.rds"))
  gc()
}

plan(sequential)
