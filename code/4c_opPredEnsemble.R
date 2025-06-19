# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update ensemble forecasts




# setup -------------------------------------------------------------------
library(tidyverse)
library(glue)
library(tidymodels)
library(glmnet)
library(future)
library(doFuture)
library(habforecastr)
library(butcher)


set.seed(1)
ncores <- 9
responses <- c(alert="alert")

target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

covSet.df <- read_csv("data/covSet_hab_tox.csv")



base.dir <- glue("out/0_init")
ens.dir <- glue("{base.dir}/ensembles/")

if(.Platform$OS.type=="unix") {
  plan(multicore, workers=ncores)
} else {
  plan(multisession, workers=ncores)
}

foreach(i=1:nrow(targ_i), .options.future=list(seed=TRUE), .errorhandling="remove") %dofuture% {
  try({
    y_i.i <- targ_i[i,]
    y.i <- targ_i$abbr[i]
    set.seed(1003)
    
    # load latest dataset for each yi_dXX
    f_fcst <- dirf("out/1_forecast/compiled", y.i, recursive=T)
    fcst_maxDate <- str_sub(f_fcst, -14, -5) |> ymd() |> max()
    f_fcst <- grep(fcst_maxDate, f_fcst, value=T)
    fcst.ls <- merge_pred_dfs(f_fcst)
    fcst.ls$alert <- fcst.ls$alert |> select(-ends_with(".x"), -ends_with(".y"))
    fcst.ls$alert$year <- year(fcst.ls$alert$date)
    cv.ls <- readRDS(glue("{base.dir}/compiled/{y.i}_cv.rds"))
    
    fcst.ls <- map(responses, ~fit_ensemble(fcst.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
    saveRDS(fcst.ls, glue("out/1_forecast/compiled/{y.i}_fcst_{fcst_maxDate}.rds"))
    
  })
}

plan(sequential); gc()
