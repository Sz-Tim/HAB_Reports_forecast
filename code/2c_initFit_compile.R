# HAB Forecasting in Scotland
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
library(habforecastr)
library(butcher)


base.dir <- "out/0_init"
target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(glue("data/i_{.x}.csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM", "ensGLM2",
                         "ensHB", "ensRF", "ensRF2", 
                         "HB1", "Ridge", "MARS", "NN", 
                         "RF", "Boost", "lgbm"),
                labels=c("Null[0]", "Null[Date]", "Null[auto]", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ensemble", "Ensemble2", 
                         "Ens-HB", "Ens-RF", "Ens-RF2", 
                         "HB", "Ridge", "MARS", "NN",
                         "RF", "XGB", "lgbm"))
mod_cols <- c(rep("grey", 3), "grey30",
              rep("grey40", 7),
              "#1f78b4", "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a", "#6a3d9a") |>
  setNames(mod_i$labels)
d_i <- tibble(f=dir(glue("{base.dir}/model_fits")))

fit.ls <- dirf(glue("{base.dir}/compiled"), "_fit.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
oos.ls <- dirf(glue("{base.dir}/compiled"), "_oos.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
# spatTime.ls <- dirf(glue("{base.dir}/compiled"), "_spatTime.rds") |>
#   map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)

fit.ls$alert_L <- fit.ls$alert |>
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |>
  mutate(model=str_split_fixed(run, "_", 3)[,1],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model),
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |>
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(d_i$f, "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensRF", "ensRF2",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  mutate(covSet=factor(covSet, levels=c(d_i$f,
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensRF", "ensRF2",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
         # y=factor(y, levels=targ_i$abbr)) |>
  arrange(y, run, siteid, date) |>
  group_by(y, run, siteid) |>
  mutate(prevAlert=lag(alert)) |>
  ungroup()

oos.ls$alert_L <- oos.ls$alert |>
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |>
  mutate(model=str_split_fixed(run, "_", 3)[,1],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model),
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |>
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(d_i$f, "ens", "ensLogitMn", "ensGLM", "ensGLM2",  "ensRF", "ensRF2",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  mutate(covSet=factor(covSet, levels=c(d_i$f,
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2",  "ensRF", "ensRF2",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
         # y=factor(y, levels=targ_i$abbr)) |>
  arrange(y, run, siteid, date) |>
  group_by(y, run, siteid) |>
  mutate(prevAlert=lag(alert)) |>
  ungroup()



# Threshold analysis ------------------------------------------------------

gc()
constituents <- fit.ls$alert_L |> count(y, run)
opt.F1 <- opt.mcc <- vector("list", nrow(constituents))
for(i in seq_along(opt.mcc)) {
  m_i <- constituents$run[i]
  fit_i <- fit.ls$alert_L |> 
    filter(run==constituents$run[i] & y==constituents$y[i]) |> 
    filter(!is.na(prevAlert)) |>
    select(y, model, PCA, covSet, obsid, siteid, date, alert, prevAlert, prA1)
  thresh.fit <- compute_thresholds2(fit_i, 
                                   0, 1, 0.01,
                                   # 0.001, 0.9, 0.0025, 
                                   byPrevAlert=!grepl("Null", m_i), 
                                   cores=50)
  if(grepl("null|ens", m_i)) {
    saveRDS(thresh.fit, 
            glue("out/0_init/metric_thresh/", 
                 "{constituents$y[i]}_", 
                 "{str_replace(str_remove(m_i, '_alert_A1'), '\\\\.', '_')}",
                 "_thresh.rds"))
  }
  if(grepl("null", m_i)) {
    opt.F1[[i]] <- thresh.fit |> filter(!is.na(F1)) |>
      group_by(y, model, PCA, covSet) |>
      arrange(desc(F1)) |> slice_head(n=1) |> ungroup() |>
      select(y, model, PCA, covSet, thresh, F1, precision, recall) |>
      rename(optF1=thresh, F1_fit=F1, F1_precision=precision, F1_recall=recall)
    opt.mcc[[i]] <- thresh.fit |> filter(!is.na(mcc)) |>
      group_by(y, model, PCA, covSet) |>
      arrange(desc(mcc)) |> slice_head(n=1) |> ungroup() |>
      select(y, model, PCA, covSet, thresh, mcc) |>
      rename(optMCC=thresh)
  } else {
    opt.F1[[i]] <- thresh.fit |> filter(!is.na(F1)) |>
      group_by(y, model, PCA, covSet, prevAlert) |>
      arrange(desc(F1)) |> slice_head(n=1) |> ungroup() |>
      select(y, model, PCA, covSet, thresh, F1, precision, recall, prevAlert) |>
      rename(optF1=thresh, F1_fit=F1, F1_precision=precision, F1_recall=recall)
    opt.mcc[[i]] <- thresh.fit |> filter(!is.na(mcc)) |>
      group_by(y, model, PCA, covSet, prevAlert) |>
      arrange(desc(mcc)) |> slice_head(n=1) |> ungroup() |>
      select(y, model, PCA, covSet, thresh, mcc, prevAlert) |>
      rename(optMCC=thresh)
  }
  gc()
  cat("Finished", i, "of", length(opt.F1), "\n")
}
m_null <- grep("null", constituents$run)
m_mods <- grep("null", constituents$run, invert=T)
opt.F1 <- list(
  do.call('rbind', opt.F1[m_mods]),
  do.call('rbind', opt.F1[m_null]) |>
    mutate(optF1=if_else(model=="Null[0]", 0.99, optF1))
)
opt.mcc <- list(
  do.call('rbind', opt.mcc[m_mods]),
  do.call('rbind', opt.mcc[m_null]) |>
    mutate(optMCC=if_else(model=="Null[0]", 0.99, optMCC))
)

fit.ls$alert_L <- bind_rows(
  fit.ls$alert_L |>
    filter(!grepl("Null", model)) |>
    left_join(opt.F1[[1]] |> select(-starts_with("F1"))) |>
    left_join(opt.mcc[[1]]),
  fit.ls$alert_L |>
    filter(grepl("Null", model)) |>
    left_join(opt.F1[[2]] |> select(-starts_with("F1"))) |>
    left_join(opt.mcc[[2]])
) |>
  mutate(predF1=factor(if_else(prA1 > optF1, "A1", "A0"), levels=c("A0", "A1")),
    predMCC=factor(if_else(prA1 > optMCC, "A1", "A0"), levels=c("A0", "A1")))

oos.ls$alert_L <- bind_rows(
  oos.ls$alert_L |>
    filter(!grepl("Null", model)) |>
    left_join(opt.F1[[1]] |> select(-F1_fit)) |>
    left_join(opt.mcc[[1]]),
  oos.ls$alert_L |>
    filter(grepl("Null", model)) |>
    left_join(opt.F1[[2]] |> select(-F1_fit)) |>
    left_join(opt.mcc[[2]])
) |>
  mutate(predF1=factor(if_else(prA1 > optF1, "A1", "A0"), levels=c("A0", "A1")),
         predMCC=factor(if_else(prA1 > optMCC, "A1", "A0"), levels=c("A0", "A1")))


saveRDS(fit.ls$alert_L, "out/clean/out_fit.rds")
saveRDS(oos.ls$alert_L, "out/clean/out_oos.rds")





# performance metrics -----------------------------------------------------


library(kerneval)
schoenr_alt <- function(d1, d2, a = NULL, b = NULL)
{
  
  if(sum(!is.na(d1)) < 2 | sum(!is.na(d2)) < 2) {
    return(NA_real_)
  }
  if (min(d1$x) > max(d2$x) | max(d1$x) < min(d2$x)) {
    return(0)
  }
  kerneval::schoenr(d1, d2, a, b)
}
density_alt <- function(x) {
  x_ <- x[!is.na(x)]
  if(length(x_) < 2) {
    NA_real_
  } else {
    density(x_)
  }
}



# ALL performance ---------------------------------------------------------

# all_df <- bind_rows(
#   fit.ls$alert_L |>
#   filter(model=="Ensemble") |>
#   mutate(week=floor_date(date, "week")),
# oos.ls$alert_L |>
#   filter(model=="Ensemble") |>
#   mutate(week=floor_date(date, "week"))
# )

all_ens <- #bind_rows(
  fit.ls$alert_L |>
    filter(grepl("Null|Ensemble", model)) |>
    mutate(week=floor_date(date, "week"))#,
  # oos.ls$alert_L |>
  #   filter(grepl("Null|Ensemble", model)) |>
  #   mutate(week=floor_date(date, "week"))
  #)
saveRDS(all_ens, "out/clean/out_all_ens.rds")

performance_all_ens <- all_ens |>
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(all_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_all_ens, "out/clean/performance_all_ens.rds")


performance_site_ens <- all_ens |>
  select(y, siteid, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, siteid) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, siteid) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(all_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, siteid, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, siteid, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, siteid) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_site_ens, "out/clean/performance_site_ens.rds")


performance_week_ens <- all_ens |>
  select(y, week, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, week) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, week) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(all_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, week, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, week) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(all_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, week, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, week) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_week_ens, "out/clean/performance_week_ens.rds")




# OOS performance ---------------------------------------------------------

oos_ens <- oos.ls$alert_L |>
  filter(grepl("Null|Ensemble", model)) |>
  mutate(week=floor_date(date, "week"))
saveRDS(oos_ens, "out/clean/out_oos_ens.rds")

performance_all_ens_oos <- oos_ens |>
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_all_ens_oos, "out/clean/performance_all_ens_oos.rds")


performance_site_ens_oos <- oos_ens |>
  select(y, siteid, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, siteid) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, siteid) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, siteid, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, siteid, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, siteid) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_site_ens_oos, "out/clean/performance_site_ens_oos.rds")


performance_week_ens_oos <- oos_ens |>
  select(y, week, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, week) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, week) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, week, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, week) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, week, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, week) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_week_ens_oos, "out/clean/performance_week_ens_oos.rds")

oos_ens <- oos_ens |>
  mutate(month=month(date))
performance_month_ens_oos <- oos_ens |>
  select(y, month, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, month) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, month) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_ens |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, month, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, month) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_ens |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, month, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, month) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_month_ens_oos, "out/clean/performance_month_ens_oos.rds")





# OOS all models ----------------------------------------------------------

oos_allMod <- oos.ls$alert_L |>
  mutate(week=floor_date(date, "week"))
saveRDS(oos_allMod, "out/clean/out_oos_allMod.rds")

performance_all_allMod_oos <- oos_allMod |>
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_allMod |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_all_allMod_oos, "out/clean/performance_all_allMod_oos.rds")


performance_site_allMod_oos <- oos_allMod |>
  select(y, siteid, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, siteid) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, siteid) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_allMod |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, siteid, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, siteid, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, siteid) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_site_allMod_oos, "out/clean/performance_site_allMod_oos.rds")


performance_week_allMod_oos <- oos_allMod |>
  select(y, week, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, week) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, week) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_allMod |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, week, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, week) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, week, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, week) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_week_allMod_oos, "out/clean/performance_week_allMod_oos.rds")

oos_allMod <- oos_allMod |>
  mutate(month=month(date))
performance_month_allMod_oos <- oos_allMod |>
  select(y, month, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, month) |>
  nest(dat=c(prA1, alert)) |>
  mutate(N=map_int(dat, nrow),
         nAlert=map_int(dat, ~sum(.x$alert=="A1")),
         AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, month) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  ungroup() |>
  bind_rows(oos_allMod |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate),
                     .estimate=if_else(covSet=="nullGrand", 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, month, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
                        N=n(),
                        nAlert=sum(!is.na(A1))) |>
              group_by(y, month) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D") |>
              ungroup()) |>
  bind_rows(oos_allMod |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model) |>
              mutate(N=n(),
                     nAlert=sum(alert=="A1")) |>
              ungroup() |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, month, N, nAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ_trunc") |>
              ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, month) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_month_allMod_oos, "out/clean/performance_month_allMod_oos.rds")

