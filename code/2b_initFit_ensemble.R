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
library(doFuture)
source("code/00_fn.R")

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
  mutate(id=row_number(),
         f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}")) |>
  ungroup() |>
  arrange(y, id) 

ncores <- 3
base.dir <- "out/0_init" 
responses <- c(alert="alert")




# Model predictions -------------------------------------------------------

if(.Platform$OS.type=="unix") {
  plan(multisession, workers=ncores)
} else {
  plan(multicore, workers=ncores)
}

foreach(i=1:nrow(covSet.df)) %dofuture% {
  
  # covariate set / response info
  f <- covSet.df$f[i]
  d <- covSet.df$id[i]
  y <- covSet.df$y[i]
  
  # directories
  fit.dir <- glue("{base.dir}/model_fits/{f}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{f}/")
  
  # load datasets
  d.y <- readRDS(glue("{base.dir}/{y}_{covSet}_dy_testPct-{test_prop}.rds")) 
  dPCA.y <- readRDS(glue("{base.dir}/{y}_{covSet}_dPCAy_testPct-{test_prop}.rds"))
  
  # generate all fitted values
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, dPCA.y$train, .x, fit.dir, y_i.i))
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))

  # generate all out-of-sample predictions
  oos.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  gc()
}


plan(sequential); gc()



# Compile -----------------------------------------------------------------

ens.dir <- glue("{base.dir}/ensembles/")

for(i in 1:nrow(y_i)) {
  
  y_i.i <- y_i[i]
  set.seed(1003)
  
  # . ensemble --------------------------------------------------------------
  
  fit.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y}_fit_ls.rds"), recursive=T))
  oos.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y}_oos_ls.rds"), recursive=T))
  cv.ls <- list(alert=full_join(
    merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y}_.*_HBL_CV"), recursive=T), CV="HB"),
    merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y}_.*_CV.rds"), recursive=T), CV="ML"),
    by=c("y", "obsid"))
  )
  saveRDS(cv.ls, glue("{base.dir}/compiled/{y}_cv.rds"))
  
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  saveRDS(wt.ls, glue("{base.dir}/compiled/{y}_wt.rds"))
  
  fit.ls <- map(responses, ~fit_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
  fit.ls <- map(responses, ~fit_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e4))
  saveRDS(fit.ls, glue("{base.dir}/compiled/{y}_fit.rds"))
  
  oos.ls <- map(responses, ~fit_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
  oos.ls <- map(responses, ~fit_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  saveRDS(oos.ls, glue("{base.dir}/compiled/{y}_oos.rds"))
  
  
  # . null ------------------------------------------------------------------
  
  fit.ls <- readRDS(glue("{base.dir}/compiled/{y}_fit.rds"))
  oos.ls <- readRDS(glue("{base.dir}/compiled/{y}_oos.rds"))

  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls,
                 ~left_join(.x |> mutate(yday=yday(date)), .y$yday.df) |> select(-yday)) |>
    map2(.x=_, fit.ls, ~bind_cols(.x, .y |> select(contains("nullGrand")) |> slice_head(n=1)))

  saveRDS(fit.ls, glue("{base.dir}/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("{base.dir}/compiled/{y}_oos.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("{base.dir}/compiled/{y}_null.rds"))
  
}

plan(sequential); gc()










