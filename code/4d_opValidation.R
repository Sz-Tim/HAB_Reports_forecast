# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update validation




# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(habforecastr)


# monitoring targets
target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(paste0("data/i_", .x, ".csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr) |>
  mutate(plotGroups=c(1, 2, 4, 3, 3, 3, 3, 2, 1)) |>
  mutate(targ_ordered=factor(abbr,
                             levels=c("Alsp", "PSP", "Disp", "DSP",
                                      "Pssp", "Psde", "Psse", "ASP",  
                                      "Kami"),
                             labels=c("Alexandrium", "PSTs",
                                      "Dinophysis", "DSTs (OA/DTXs/PTXs)",
                                      paste("Pseudo-nitzschia", c("spp.", "del.", "ser.")), "DA", 
                                      "Karenia mikimotoi"))) |>
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

# output and validation ---------------------------------------------------

# # point-wise predictions
# all_df <- readRDS("out/clean/temp/out_all_ens.rds") |>
#   filter(model=="Ensemble")
# # oos_df <- readRDS(paste0(dirs$proj, "out/clean/temp/out_oos_ens.rds"))
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
# validation_df <- readRDS("out/clean/temp/performance_all_ens_oos.rds") |>
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
# validation_sin_df <- readRDS("out/clean/temp/performance_site_ens_oos.rds") |>
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
# validation_month_df <- readRDS("out/clean/temp/performance_month_ens_oos.rds") |>
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

