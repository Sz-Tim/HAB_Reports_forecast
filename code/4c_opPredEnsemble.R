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
ncores <- 10
responses <- c(alert="alert")

target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

covSet.df <- read_csv("data/covSet_df.csv")
history_f <- "out/1_forecast/compiled/fcst_history_df.rds"
if(file.exists(history_f)) {
  fcst_history_df <- readRDS(history_f)
} else {
  fcst_history_df <- tibble(y=character(), 
                            obsid=numeric(),
                            siteid=numeric(),
                            date_generated=as_date(character()),
                            date_forecast=as_date(character()),
                            
                            ensGLM2_alert_A1=numeric())
}


base.dir <- glue("out/0_init")
ens.dir <- glue("{base.dir}/ensembles/")

if(.Platform$OS.type=="unix") {
  plan(multicore, workers=ncores)
} else {
  plan(multisession, workers=ncores)
}

for(i in 1:nrow(targ_i)) {
  y_i.i <- targ_i[i,]
  y.i <- targ_i$abbr[i]
  set.seed(1003)
  
  # load latest dataset for each yi_dXX
  f_fcst <- dirf("out/1_forecast/compiled", y.i, recursive=T) |>
    grep("d[0-9]", x=_, value=T)
  fcst_maxDate <- str_sub(f_fcst, -14, -5) |> ymd() |> max()
  f_fcst <- grep(fcst_maxDate, f_fcst, value=T)
  fcst.ls <- merge_pred_dfs(f_fcst)
  fcst.ls$alert <- fcst.ls$alert |> select(-ends_with(".x"), -ends_with(".y"))
  fcst.ls$alert$year <- year(fcst.ls$alert$date)
  cv.ls <- readRDS(glue("{base.dir}/compiled/{y.i}_cv.rds"))
  
  fcst.ls <- map(responses, ~fit_ensemble(fcst.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  saveRDS(fcst.ls, glue("out/1_forecast/compiled/{y.i}_fcst_{fcst_maxDate}.rds"))
  fcst_history_df <- bind_rows(
    fcst_history_df,
    fcst.ls$alert |> 
      rename(date_forecast=date) |>
      mutate(date_generated=today()) |>
      select(y, obsid, siteid, date_generated, date_forecast, ensGLM2_alert_A1)
  )
}

plan(sequential); gc()
saveRDS(fcst_history_df, history_f)
