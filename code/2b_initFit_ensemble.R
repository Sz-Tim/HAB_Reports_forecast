# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Train ensemble models



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
library(habforecastr)
library(butcher)


set.seed(1)
fit_ens_with_testing_data <- FALSE 
ncores <- 20
train_prop <- 1
run_type <- "0_init_fish" 
responses <- c(alert="alert")

target_sets <- c("hab", "tox", "fish")[3]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

covSet.df <- read_csv("data/covSet_hab_tox.csv")
covSet.df <- read_csv("data/covSet_fish.csv")



# Model predictions -------------------------------------------------------

if(.Platform$OS.type=="unix") {
  plan(multicore, workers=ncores)
} else {
  plan(multisession, workers=ncores)
}

# foreach(i=1:nrow(covSet.df), .options.future=list(seed=TRUE)) %dofuture% {
for(i in 1:nrow(covSet.df)) {
  try({
  # covariate set / response info
  id <- covSet.df$id[i]
  y.i <- covSet.df$y[i]
  y_i.i <- targ_i |> filter(abbr==y.i)
  
  # directories
  dummy <- paste(run_type, train_prop)
  dir.create(paste0("out/", run_type, "/compiled/logs/"), showWarnings=F, recursive=T)
  cat(y.i, id, file=paste0("out/", run_type, "/compiled/logs/", y.i, "_", id, ".log"))
  data.dir <- glue("data/{run_type}/")
  base.dir <- glue("out/{run_type}/")
  fit.dir <- glue("{base.dir}/model_fits/{id}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{id}/")
  
  # load datasets
  if(!file.exists(glue("{data.dir}/compiled/{y.i}_{id}_dy_testPct-{train_prop}.rds"))) {
    next
  }
  d.y <- readRDS(glue("{data.dir}/compiled/{y.i}_{id}_dy_testPct-{train_prop}.rds"))
  dPCA.y <- readRDS(glue("{data.dir}/compiled/{y.i}_{id}_dPCAy_testPct-{train_prop}.rds"))
  
  # generate all fitted values
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, dPCA.y$train, .x, fit.dir, y_i.i))
  saveRDS(fit.ls, glue("{out.dir}/{y.i}_fit_ls.rds"))

  if(train_prop < 1) {
    # generate all out-of-sample predictions
    oos.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
    saveRDS(oos.ls, glue("{out.dir}/{y.i}_oos_ls.rds"))
  }
  gc()
  })
}

plan(sequential); gc()



# Compile -----------------------------------------------------------------

base.dir <- glue("out/{run_type}")
ens.dir <- glue("{base.dir}/ensembles/")

for(i in 1:nrow(targ_i)) {
  try({
    y_i.i <- targ_i[i,]
    y.i <- targ_i$abbr[i]
    set.seed(1003)
    
    # . ensemble --------------------------------------------------------------
    if(length(dirf(glue("{base.dir}/compiled"), glue("{y.i}_fit_ls.rds"), recursive=T))==0) {
      next
    }
    fit.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y.i}_fit_ls.rds"), recursive=T))
    fit.ls$alert <- fit.ls$alert |> select(-ends_with(".x"), -ends_with(".y"))
    
    if(length(dirf(glue("{base.dir}/model_fits"), glue("{y.i}_.*_HB_CV"), recursive=T)) == 0) {
      HB_CV <- fit.ls$alert[, c("y", "obsid", "siteid", "date", "alert")]
    } else {
      HB_CV <- merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y.i}_.*_HB_CV"), recursive=T), CV="HB")
    }
    cv.ls <- list(alert=full_join(
      HB_CV,
      merge_pred_dfs(dirf(glue("{base.dir}/model_fits"), glue("{y.i}_.*_CV.rds"), recursive=T), CV="ML"),
      by=c("y", "obsid"))
    )
    wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
    
    fit.ls$alert$year <- year(fit.ls$alert$date)
    cv.ls$alert$year <- year(cv.ls$alert$date)
    
    na_cols <- names(which(colSums(is.na(cv.ls$alert)) > 0))
    if(length(na_cols) > 0) {
      fit.ls$alert <- fit.ls$alert |> select(-any_of(na_cols))
      cv.ls$alert <- cv.ls$alert |> select(-any_of(na_cols))
    }
    if(train_prop < 1) {
      oos.ls <- merge_pred_dfs(dirf(glue("{base.dir}/compiled"), glue("{y.i}_oos_ls.rds"), recursive=T))
      oos.ls$alert <- oos.ls$alert |> select(-ends_with(".x"), -ends_with(".y"))
      oos.ls$alert$year <- year(oos.ls$alert$date)
      if(length(na_cols) > 0) {
        oos.ls$alert <- oos.ls$alert |> select(-any_of(na_cols))
      }
    }
    if(fit_ens_with_testing_data) {
      # fit ensemble using testing data; validate with future data
      oos.cv <- list(alert=oos.ls$alert |> select(any_of(names(cv.ls$alert))))
      if(sum(oos.cv$alert$alert=="A1") < 10) {
        next
      }
      oos.ls <- map(responses, ~fit_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
      oos.ls <- map(responses, ~fit_ensemble(oos.ls, oos.cv, .x, y_i.i, "GLM_fit", ens.dir, 1e3))
    } else {
      # fit ensembles using training data; validate with testing data
      fit.ls <- map(responses, ~fit_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
      fit.ls <- map(responses, ~fit_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e3))
      # fit.ls <- map(responses, ~fit_ensemble(fit.ls, cv.ls, .x, y_i.i, "RF_fit", ens.dir, 1e3))
      if(train_prop < 1) {
        oos.ls <- map(responses, ~fit_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
        oos.ls <- map(responses, ~fit_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
        # oos.ls <- map(responses, ~fit_ensemble(oos.ls, cv.ls, .x, y_i.i, "RF_oos", ens.dir))
      }
    }
    
    
    # . null ------------------------------------------------------------------
    null.ls <- map(responses, ~calc_null(fit.ls, .x))
    fit.ls <- map(null.ls, ~.x$obs.df)
    
    saveRDS(wt.ls, glue("{base.dir}/compiled/{y.i}_wt.rds"))
    saveRDS(cv.ls, glue("{base.dir}/compiled/{y.i}_cv.rds"))
    saveRDS(fit.ls, glue("{base.dir}/compiled/{y.i}_fit.rds"))
    saveRDS(map(null.ls, ~.x$yday.df), glue("{base.dir}/compiled/{y.i}_null.rds"))
    if(train_prop < 1) {
      oos.ls <- map2(oos.ls, null.ls,
                     ~left_join(.x |> mutate(yday=yday(date)), .y$yday.df) |> select(-yday)) |>
        map2(.x=_, fit.ls, ~bind_cols(.x, .y |> select(contains("nullGrand")) |> slice_head(n=1)))
      saveRDS(oos.ls, glue("{base.dir}/compiled/{y.i}_oos.rds"))
    }
    
    
  })
}

plan(sequential); gc()










