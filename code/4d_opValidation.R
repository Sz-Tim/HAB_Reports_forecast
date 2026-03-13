# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update validation




# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(habforecastr)
library(yardstick)


# monitoring targets
target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(paste0("data/i_", .x, ".csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr) |>
  mutate(plotGroups=c(1, 2, 4, 5, 3, 3, 3, 3, 2, 1)) |>
  mutate(targ_ordered=factor(abbr,
                             levels=c("Alsp", "PSP", "Disp", "DSP",
                                      "Pssp", "Psde", "Psse", "ASP",  
                                      "Prli", "Kami"),
                             labels=c("Alexandrium", "PSTs",
                                      "Dinophysis", "DSTs (OA/DTXs/PTXs)",
                                      paste("Pseudo-nitzschia", c("spp.", "del.", "ser.")), "DA", 
                                      "Prorcentrum lima", "Karenia mikimotoi"))) |>
  arrange(targ_ordered)

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
d_i <- tibble(f=dir("out/0_init/model_fits"))


# monitoring locations
site_i <- map_dfr(target_sets,
                  ~readRDS(paste0("data/site_", .x, "_df.rds")) |>
                    mutate(type=.x))

fcst.ls <- dirf(glue("out/1_forecast/compiled"), "_fcst_") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)

fcst.ls$alert_L <- fcst.ls$alert |>
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
saveRDS(fcst.ls$alert_L, 
        glue("out/1_forecast/fcst_all_{format(max(fcst.ls$alert$date), '%F')}.rds"))
saveRDS(fcst.ls$alert_L |> filter(model=="Ensemble"),
        glue("out/1_forecast/fcst_ens_{format(max(fcst.ls$alert$date), '%F')}.rds"))



# calculate validation ----------------------------------------------------

MCC_df <- readRDS("out/clean/out_oos_ens.rds") |>
  filter(!is.na(optMCC)) |>
  group_by(model, y, prevAlert) |>
  slice_head(n=1) |>
  ungroup() |>
  select(model, y, prevAlert, optMCC, prA1) 
F1_df <- readRDS("out/clean/out_oos_ens.rds") |>
  filter(!is.na(optF1)) |>
  group_by(model, y, prevAlert) |>
  slice_head(n=1) |>
  ungroup() |>
  select(model, y, prevAlert, optF1, prA1)
ens_oos <- readRDS("out/clean/out_oos_ens.rds")
ens_fcst <- dirf("out/1_forecast", "fcst_all_.*rds") |>
  sort() |> last() |> readRDS() |>
  filter(grepl("Ensemble", model)) |>
  bind_rows(readRDS("out/clean/out_oos_ens.rds") |> 
              filter(grepl("Null", model)) |> 
              select(-matches("MCC|mcc|F1"))) |>
  inner_join(MCC_df |> select(-prA1), by=join_by(model, y, prevAlert)) |>
  inner_join(F1_df |> select(-prA1), by=join_by(model, y, prevAlert)) |>
  mutate(predMCC=factor(if_else(prA1 > optMCC, "A1", "A0"), levels=c("A0", "A1")),
         week=floor_date(date, "week"))
saveRDS(ens_fcst, "out/clean/out_ens_fcst.rds")

oos_combined <- bind_rows(ens_oos, ens_fcst) |>
  mutate(month=month(date))

performance_all_ens_fcst <- oos_combined |>
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
  bind_rows(oos_combined |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             select(y, covSet, PCA, model, obsid, alert, prA1) %>%
  #             filter(!is.na(prA1)) |>
  #             pivot_wider(names_from="alert", values_from="prA1") |>
  #             group_by(y, model, PCA, covSet) |>
  #             summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
  #                       N=n(),
  #                       nAlert=sum(!is.na(A1))) |>
  #             group_by(y) |>
  #             mutate(rank=min_rank(.estimate),
  #                    .metric="Schoener's D") |>
  #             ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, model) |>
  #             mutate(N=n(),
  #                    nAlert=sum(alert=="A1")) |>
  #             ungroup() |>
  #             mutate(prA1=if_else(prA1==0, 1e-5, prA1),
  #                    prA1=if_else(prA1==1, 1-1e-5, prA1),
  #                    alert=as.numeric(alert=="A1")) |>
  #             calc_R2(type="vz", y, N, nAlert) |>
  #             rename(.estimate=R2) |>
  #             mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
  #             na.omit() |>
  #             group_by(y) |>
  #             mutate(rank=min_rank(desc(.estimate)),
  #                      .metric="R2-VZ_trunc") |>
  #             ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_all_ens_fcst, "out/clean/performance_all_ens_fcst.rds")
saveRDS(performance_all_ens_fcst, "out/2_shiny/validation_df.rds")


performance_site_ens_fcst <- oos_combined |>
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
  bind_rows(oos_combined |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, siteid, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, siteid) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             select(y, siteid, covSet, PCA, model, obsid, date, alert, prA1) |>
  #             filter(!is.na(prA1)) |>
  #             pivot_wider(names_from="alert", values_from="prA1") |>
  #             group_by(y, siteid, model, PCA, covSet) |>
  #             summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
  #                       N=n(),
  #                       nAlert=sum(!is.na(A1))) |>
  #             group_by(y, siteid) |>
  #             mutate(rank=min_rank(.estimate),
  #                    .metric="Schoener's D") |>
  #             ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, siteid, model) |>
  #             mutate(N=n(),
  #                    nAlert=sum(alert=="A1")) |>
  #             ungroup() |>
  #             mutate(prA1=if_else(prA1==0, 1e-5, prA1),
  #                    prA1=if_else(prA1==1, 1-1e-5, prA1),
  #                    alert=as.numeric(alert=="A1")) |>
  #             calc_R2(type="vz", y, siteid, N, nAlert) |>
  #             rename(.estimate=R2) |>
  #             mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
  #             na.omit() |>
  #             group_by(y, siteid) |>
  #             mutate(rank=min_rank(desc(.estimate)),
  #                    .metric="R2-VZ_trunc") |>
  #             ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, siteid) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) |>
  left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
  left_join(site_i |> select(siteid, sin, lon, lat, type))
saveRDS(performance_site_ens_fcst, "out/clean/performance_site_ens_fcst.rds")
saveRDS(performance_site_ens_fcst, "out/2_shiny/validation_sin_df.rds")


performance_week_ens_fcst <- oos_combined |>
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
  bind_rows(oos_combined |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, week, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, week) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             select(y, week, covSet, PCA, model, obsid, date, alert, prA1) %>%
  #             filter(!is.na(prA1)) |>
  #             pivot_wider(names_from="alert", values_from="prA1") |>
  #             group_by(y, week, model, PCA, covSet) |>
  #             summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
  #                       N=n(),
  #                       nAlert=sum(!is.na(A1))) |>
  #             group_by(y, week) |>
  #             mutate(rank=min_rank(.estimate),
  #                    .metric="Schoener's D") |>
  #             ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, week, model) |>
  #             mutate(N=n(),
  #                    nAlert=sum(alert=="A1")) |>
  #             ungroup() |>
  #             mutate(prA1=if_else(prA1==0, 1e-5, prA1),
  #                    prA1=if_else(prA1==1, 1-1e-5, prA1),
  #                    alert=as.numeric(alert=="A1")) |>
  #             calc_R2(type="vz", y, week, N, nAlert) |>
  #             rename(.estimate=R2) |>
  #             mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
  #             na.omit() |>
  #             group_by(y, week) |>
  #             mutate(rank=min_rank(desc(.estimate)),
  #                    .metric="R2-VZ_trunc") |>
  #             ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, week) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate)))) 
saveRDS(performance_week_ens_fcst, "out/clean/performance_week_ens_fcst.rds")
saveRDS(performance_week_ens_fcst, "out/2_shiny/validation_week_df.rds")

performance_month_ens_fcst <- oos_combined |>
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
  bind_rows(oos_combined |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=roc_auc_vec(prA1, truth=alert, event_level="second"),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=mcc_vec(truth=alert, estimate=predMCC),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              ungroup()) |>
  bind_rows(oos_combined |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, month, model, PCA, covSet) |>
              summarise(.estimate=f_meas_vec(truth=alert, estimate=predF1),
                        N=n(),
                        nAlert=sum(alert=="A1")) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, month) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             select(y, month, covSet, PCA, model, obsid, date, alert, prA1) %>%
  #             filter(!is.na(prA1)) |>
  #             pivot_wider(names_from="alert", values_from="prA1") |>
  #             group_by(y, month, model, PCA, covSet) |>
  #             summarise(.estimate=schoenr_alt(density_alt(A0), density_alt(A1)),
  #                       N=n(),
  #                       nAlert=sum(!is.na(A1))) |>
  #             group_by(y, month) |>
  #             mutate(rank=min_rank(.estimate),
  #                    .metric="Schoener's D") |>
  #             ungroup()) |>
  # bind_rows(oos_combined |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, month, model) |>
  #             mutate(N=n(),
  #                    nAlert=sum(alert=="A1")) |>
  #             ungroup() |>
  #             mutate(prA1=if_else(prA1==0, 1e-5, prA1),
  #                    prA1=if_else(prA1==1, 1-1e-5, prA1),
  #                    alert=as.numeric(alert=="A1")) |>
  #             calc_R2(type="vz", y, month, N, nAlert) |>
  #             rename(.estimate=R2) |>
  #             mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
  #             na.omit() |>
  #             group_by(y, month) |>
  #             mutate(rank=min_rank(desc(.estimate)),
  #                    .metric="R2-VZ_trunc") |>
  #             ungroup()) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric, month) |>
  mutate(score_v_null=.estimate - first(.estimate),
         skill=if_else(.metric=="Schoener's D",
                       score_v_null/(0 - first(.estimate)),
                       score_v_null/(1 - first(.estimate))),
         date_std=ymd(paste0("2024-", month, "-01")))
saveRDS(performance_month_ens_fcst, "out/clean/performance_month_ens_fcst.rds")
saveRDS(performance_month_ens_fcst, "out/2_shiny/validation_month_df.rds")






# merge for shiny ---------------------------------------------------------

# update observation dataframe
obs_old <- bind_rows(
  readRDS("data/1_current/hab_obs.rds"),
  readRDS("data/1_current/tox_obs.rds")
) |>
  select(sin, date, siteid, y, N, lnN, tl, alert) |>
  mutate(y_site_date=paste(y, siteid, date, sep="_"),
         src="old")
obs_new <- bind_rows(
  readRDS("data/2_new/hab_obs.rds"),
  readRDS("data/2_new/tox_obs.rds")
) |>
  select(sin, date, siteid, y, N, lnN, tl, alert) |>
  mutate(y_site_date=paste(y, siteid, date, sep="_"),
         src="new")
obs_ysd <- unique(obs_old$y_site_date)
obs_combined <- obs_old |>
  bind_rows(obs_new |> filter(! y_site_date %in% obs_ysd)) |>
  left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
  left_join(site_i |> select(siteid, sin, lon, lat, type)) |>
  mutate(week=floor_date(date, "week"),
         tl=factor(tl, levels=c("TL0", "TL1", "TL2", "TL3"), ordered=T))
saveRDS(obs_combined, "out/2_shiny/obs_df.rds")


# update forecast dataframe
ens_all <- readRDS("out/clean/out_oos_ens.rds") |>
  mutate(y_site_date=paste(y, siteid, date, sep="_"),
         src="old")
ens_fcst <- readRDS("out/clean/out_ens_fcst.rds") |>
  mutate(y_site_date=paste(y, siteid, date, sep="_"),
         src="new")

fcst_ysd <- unique(filter(ens_fcst, model=="Ensemble")$y_site_date)
all_ysd <- unique(filter(ens_all, model=="Ensemble")$y_site_date)

ens_combined <- ens_all |>
  bind_rows(ens_fcst |> 
              filter(y_site_date %in% fcst_ysd) |>
              filter(! y_site_date %in% all_ysd)) |>
  group_by(y, siteid) |>
  mutate(N=n(), nYrs=n_distinct(year(date))) |>
  ungroup() |> filter(N >= 30 & nYrs > 1) |>
  arrange(y, siteid, date) |>
  select(y, siteid, date, week, obsid, model, prevAlert, alert, prA1) |>
  left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
  left_join(site_i |> select(siteid, sin, lon, lat, type)) |>
  mutate(week=floor_date(date, "week"))
saveRDS(ens_combined, "out/2_shiny/fcst_df.rds")


# update observation dataframe



# output and validation ---------------------------------------------------

# # point-wise predictions
# all_df <- readRDS("out/clean/temp/out_all_ens.rds") |>
#   filter(model=="Ensemble")
# # oos_df <- readRDS(paste0(dirs$proj, "out/clean/temp/out_ens_fcst.rds"))
# # all_df <- bind_rows(fit_df, oos_df) |>
# all_df <- all_df |>
#   group_by(y, siteid) |>
#   mutate(N=n(), nYrs=n_distinct(year(date))) |>
#   ungroup() |> filter(N >= 30 & nYrs > 1) |>
#   arrange(y, siteid, date) |>
#   select(y, siteid, date, week, obsid, prevAlert, alert, prA1) |>
#   left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
#   left_join(site_i |> select(siteid, sin, lon, lat, type), by=join_by(type, siteid))
# saveRDS(all_df, "out/clean/all_df.rds")
# 
# obs_df <- bind_rows(
#   readRDS("data/0_init/data_hab_all.rds") |>
#     mutate(type="hab"),
#   readRDS("data/0_init/data_tox_all.rds") |>
#     mutate(type="tox")
# ) |>
#   select(type, y, obsid, sin, date, N, lnN, tl, alert) |>
#   filter(sin %in% unique(all_df$sin)) |>
#   # inner_join(all_df |> select(y, obsid), by=join_by(y, obsid)) |>
#   mutate(week=floor_date(date, unit="week"),
#          tl=factor(tl, levels=paste0("TL", 0:3), ordered=T)) |>
#   # complete(type, y, sin, week, fill=list(N=0, lnN=0, tl="TL0", alert="A0")) |>
#   group_by(type, y, sin, week) |>
#   summarise(lnN=mean(lnN), tl=max(tl)) |>
#   group_by(y) |>
#   mutate(lnN_rel=lnN/max(lnN)) |>
#   ungroup()
# saveRDS(obs_df, "out/clean/obs_df.rds")
# 
# # skill scores: Overall
# validation_df <- readRDS("out/clean/temp/performance_all_ens_fcst.rds") |>
#   filter(model %in% c("Null[0]", "Null[Date]", "Ensemble")) |> 
#   filter(.metric %in% c("MCC", "PR-AUC", "R2-VZ_trunc", "ROC-AUC", "Schoener's D")) |>
#   mutate(.metric=if_else(.metric=="R2-VZ_trunc", "R2-VZ", .metric)) |>
#   arrange(y, model, covSet, PCA) |>
#   group_by(y, .metric) |>
#   mutate(null_to_perfect=case_when(.metric=="Schoener's D" ~ first(.estimate),
#                                    .default=1-first(.estimate)),
#          m_to_null=case_when(.metric=="Schoener's D" ~ first(.estimate) - .estimate,
#                              .default=.estimate - first(.estimate)),
#          skill=m_to_null/null_to_perfect) |>
#   ungroup() |>
#   mutate(modType=case_when(model=="Null[0]" ~ model,
#                            model=="Null[Date]" ~ model,
#                            model=="Ensemble" ~ "Ensemble",
#                            .default="Constituent"),
#          modType=factor(modType, levels=c("Null[0]", "Null[Date]", "Constituent", "Ensemble")),
#          .metric=factor(.metric, 
#                         levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
#                         labels=c("AUC['ROC']", "AUC['PR']", "MCC", "R['VZ']^2", "D['Overlap']"))) |>
#   group_by(y, .metric, modType) |>
#   mutate(modNum=row_number()) |>
#   ungroup() |>
#   select(y, modType, modNum, .metric, .estimate, skill)
# saveRDS(validation_df, "out/clean/validation_df.rds")
# 
# # skill scores: By site
# validation_sin_df <- readRDS("out/clean/temp/performance_site_ens_fcst.rds") |>
#   filter(model %in% c("Null[0]", "Null[Date]", "Ensemble")) |> 
#   filter(.metric %in% c("MCC", "PR-AUC", "R2-VZ_trunc", "ROC-AUC", "Schoener's D")) |>
#   filter(nAlert > 0) |>
#   mutate(.metric=if_else(.metric=="R2-VZ_trunc", "R2-VZ", .metric)) |>
#   arrange(y, model, covSet, PCA) |>
#   group_by(y, .metric, siteid) |>
#   mutate(null_to_perfect=case_when(.metric=="Schoener's D" ~ first(.estimate),
#                                    .default=1-first(.estimate)),
#          m_to_null=case_when(.metric=="Schoener's D" ~ first(.estimate) - .estimate,
#                              .default=.estimate - first(.estimate)),
#          skill=m_to_null/null_to_perfect) |>
#   ungroup() |>
#   mutate(modType=case_when(model=="Null[0]" ~ model,
#                            model=="Null[Date]" ~ model,
#                            model=="Ensemble" ~ "Ensemble",
#                            .default="Constituent"),
#          modType=factor(modType, levels=c("Null[0]", "Null[Date]", "Constituent", "Ensemble")),
#          .metric=factor(.metric, 
#                         levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
#                         labels=c("AUC['ROC']", "AUC['PR']", "MCC", "R['VZ']^2", "D['Overlap']"))) |>
#   group_by(y, .metric, modType, siteid) |>
#   mutate(modNum=row_number()) |>
#   ungroup() |>
#   select(siteid, y, modType, modNum, .metric, .estimate, skill) |>
#   left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
#   left_join(site_i |> select(siteid, sin, lon, lat, type), by=join_by(type, siteid)) |>
#   filter(sin %in% unique(all_df$sin))
# saveRDS(validation_sin_df, "out/clean/validation_sin_df.rds")
# 
# # skill scores: By month
# validation_month_df <- readRDS("out/clean/temp/performance_month_ens_fcst.rds") |>
#   filter(model %in% c("Null[0]", "Null[Date]", "Ensemble")) |> 
#   filter(.metric %in% c("MCC", "PR-AUC", "R2-VZ_trunc", "ROC-AUC", "Schoener's D")) |>
#   filter(nAlert > 0) |>
#   mutate(.metric=if_else(.metric=="R2-VZ_trunc", "R2-VZ", .metric)) |>
#   arrange(y, model, covSet, PCA) |>
#   group_by(y, .metric, month) |>
#   mutate(null_to_perfect=case_when(.metric=="Schoener's D" ~ first(.estimate),
#                                    .default=1-first(.estimate)),
#          m_to_null=case_when(.metric=="Schoener's D" ~ first(.estimate) - .estimate,
#                              .default=.estimate - first(.estimate)),
#          skill=m_to_null/null_to_perfect) |>
#   ungroup() |>
#   mutate(modType=case_when(model=="Null[0]" ~ model,
#                            model=="Null[Date]" ~ model,
#                            model=="Ensemble" ~ "Ensemble",
#                            .default="Constituent"),
#          modType=factor(modType, levels=c("Null[0]", "Null[Date]", "Constituent", "Ensemble")),
#          .metric=factor(.metric, 
#                         levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
#                         labels=c("AUC['ROC']", "AUC['PR']", "MCC", "R['VZ']^2", "D['Overlap']"))) |>
#   group_by(y, .metric, modType, month) |>
#   mutate(modNum=row_number()) |>
#   ungroup() |>
#   select(month, y, modType, modNum, .metric, .estimate, skill) |>
#   mutate(date_std=ymd(paste0("2024-", month, "-01"))) |>
#   full_join(targ_i |> select(abbr, plotGroups), by=join_by(y==abbr))
# saveRDS(validation_month_df, "out/clean/validation_month_df.rds")

