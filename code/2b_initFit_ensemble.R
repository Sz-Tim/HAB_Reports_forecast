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
dir("code/fn", ".R", full.names=T) |> walk(source)

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

ncores <- 3
run_type <- "0_init" 
test_prop <- 0.75
responses <- c(alert="alert")




# Model predictions -------------------------------------------------------

if(.Platform$OS.type=="unix") {
  plan(multisession, workers=ncores)
} else {
  plan(multicore, workers=ncores)
}

foreach(i=1:nrow(covSet.df)) %dofuture% {
  
  # covariate set / response info
  run_type <- "0_init" 
  test_prop <- 0.75
  f <- covSet.df$f[i]
  id <- covSet.df$id[i]
  y.i <- covSet.df$y[i]
  y_i.i <- y_i |> filter(abbr==y.i)
  
  # directories
  data.dir <- glue("data/{run_type}/")
  base.dir <- glue("out/{run_type}/")
  fit.dir <- glue("{base.dir}/model_fits/{f}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{f}/")
  
  # load datasets
  d.y <- readRDS(glue("{data.dir}/compiled/{y.i}_{id}_dy_testPct-{test_prop}.rds"))
  dPCA.y <- readRDS(glue("{data.dir}/compiled/{y.i}_{id}_dPCAy_testPct-{test_prop}.rds"))
  
  # generate all fitted values
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, dPCA.y$train, .x, fit.dir, y_i.i))
  saveRDS(fit.ls, glue("{out.dir}/{y.i}_fit_ls.rds"))

  # generate all out-of-sample predictions
  oos.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
  saveRDS(oos.ls, glue("{out.dir}/{y.i}_oos_ls.rds"))
  gc()
}


plan(sequential); gc()



# Compile -----------------------------------------------------------------

ens.dir <- glue("{base.dir}/ensembles/")

for(i in 1:nrow(y_i)) {
  
  y_i.i <- y_i[i,]
  set.seed(1003)
  
  # . ensemble --------------------------------------------------------------
  
  fit.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y.i}_fit_ls.rds"), recursive=T))
  oos.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y.i}_oos_ls.rds"), recursive=T))
  cv.ls <- list(alert=full_join(
    merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y.i}_.*_HB_CV"), recursive=T), CV="HB"),
    merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y.i}_.*_CV.rds"), recursive=T), CV="ML"),
    by=c("y", "obsid"))
  )
  saveRDS(cv.ls, glue("{base.dir}/compiled/{y.i}_cv.rds"))
  
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  saveRDS(wt.ls, glue("{base.dir}/compiled/{y.i}_wt.rds"))
  
  fit.ls$alert$year <- year(fit.ls$alert$date)
  oos.ls$alert$year <- year(oos.ls$alert$date)
  cv.ls$alert$year <- year(cv.ls$alert$date)
  fit.ls <- map(responses, ~fit_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
  fit.ls <- map(responses, ~fit_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e2))
  fit.ls <- map(responses, ~fit_ensemble(fit.ls, cv.ls, .x, y_i.i, "RF_fit", ens.dir, 10))
  saveRDS(fit.ls, glue("{base.dir}/compiled/{y.i}_fit.rds"))
  
  oos.ls <- map(responses, ~fit_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
  oos.ls <- map(responses, ~fit_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  oos.ls <- map(responses, ~fit_ensemble(oos.ls, cv.ls, .x, y_i.i, "RF_oos", ens.dir))
  saveRDS(oos.ls, glue("{base.dir}/compiled/{y.i}_oos.rds"))
  
  
  # . null ------------------------------------------------------------------
  
  fit.ls <- readRDS(glue("{base.dir}/compiled/{y.i}_fit.rds"))
  oos.ls <- readRDS(glue("{base.dir}/compiled/{y.i}_oos.rds"))

  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls,
                 ~left_join(.x |> mutate(yday=yday(date)), .y$yday.df) |> select(-yday)) |>
    map2(.x=_, fit.ls, ~bind_cols(.x, .y |> select(contains("nullGrand")) |> slice_head(n=1)))

  saveRDS(fit.ls, glue("{base.dir}/compiled/{y.i}_fit.rds"))
  saveRDS(oos.ls, glue("{base.dir}/compiled/{y.i}_oos.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("{base.dir}/compiled/{y.i}_null.rds"))
  
}

plan(sequential); gc()










