# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model analysis




# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "brms", "furrr")
lapply(pkgs, library, character.only=T)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv") |> arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv") |> arrange(abbr) |> mutate(type="tox"))
y_resp <- c("Alsp", "PSP", "Disp", "DSP", "Pssp", "Psde", "Psse", "ASP", "Kami")

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM", "ensGLM2",
                         "ensHB", "ensRF", "ensRF2", 
                         "HBL1", "Ridge", "MARS", "NN", 
                         "RF", "Boost"),
                labels=c("Null[0]", "Null[Date]", "Null[auto]", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ensemble", "Ensemble2", 
                         "Ens-HB", "Ens-RF", "Ens-RF2", 
                         "HB", "Ridge", "MARS", "NN",
                         "RF", "XGB"))
mod_cols <- c(rep("grey", 3), "grey30",
              rep("grey40", 7),
              "#1f78b4", "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a") |>
  setNames(mod_i$labels)
d_i <- tibble(f=dir(glue("{base.dir}/model_fits"))) |>
  mutate(Regional=grepl("Avg1", f),
         EnvInt=grepl("Xf1", f),
         AutoInt=grepl("XN1", f),
         EnvDelta=grepl("Del1", f),
         covSet=paste0("d", str_split_fixed(f, "-", 2)[,1])) |>
  arrange(Regional, EnvDelta, AutoInt, EnvInt) |>
  mutate(covSet_reorder=factor(covSet, levels=unique(covSet)))

fit.ls <- dirf(glue("{base.dir}/compiled"), "_fit.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
oos.ls <- dirf(glue("{base.dir}/compiled"), "_oos.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
# spatTime.ls <- dirf(glue("{base.dir}/compiled"), "_spatTime.rds") |>
#   map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)

fit.ls$alert_L <- fit.ls$alert |>
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |>
  mutate(model=str_split_fixed(run, "_", 3)[,1],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model),
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |>
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  left_join(d_i |> select(-f)) |>
  mutate(covSet=factor(covSet, levels=c(levels(d_i$covSet_reorder),
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect")),
         y=factor(y, levels=y_resp))

oos.ls$alert_L <- oos.ls$alert |>
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |>
  mutate(model=str_split_fixed(run, "_", 3)[,1],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model),
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |>
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  left_join(d_i |> select(-f)) |>
  mutate(covSet=factor(covSet, levels=c(levels(d_i$covSet_reorder),
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect")),
         y=factor(y, levels=y_resp))

# spatTime.ls$alert_L <- spatTime.ls$alert |>
#   select(-starts_with("d")) |>
#   pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |>
#   mutate(model=str_split_fixed(run, "_", 3)[,1],
#          PCA=grepl("PCA", model),
#          covSet=str_split_fixed(model, "\\.", 2)[,1],
#          model=if_else(grepl("^d", model),
#                        str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
#                        str_remove(model, "PCA."))) |>
#   mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
#          covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
#                                         "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
#   mutate(covSet=factor(covSet, levels=c(levels(d_i$covSet_reorder),
#                                         "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
#                                         "null4wk", "nullAuto", "nullGrand", "perfect")),
#          y=factor(y, levels=y_resp))
# 
# saveRDS(spatTime.ls$alert_L, "out/clean/out_spatTime.rds")




# Threshold analysis ------------------------------------------------------

gc()
opt.F1 <- opt.mcc <- vector("list", n_distinct(fit.ls$alert_L$model))
for(i in seq_along(opt.F1)) {
  m_i <- unique(fit.ls$alert_L$model)[i]
  fit_i <- fit.ls$alert_L |> filter(model==m_i)
  thresh.fit <- compute_thresholds(fit_i, 
                                   0.001, 0.9, 0.0025, 
                                   byPrevAlert=!grepl("Null", m_i), 
                                   cores=10)
  if(grepl("Null", m_i)) {
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
m_null <- grep("Null", unique(fit.ls$alert_L$model))
m_mods <- grep("Null", unique(fit.ls$alert_L$model), invert=T)
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
rank.df <- oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  bind_rows(oos.ls$alert_L |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              roc_auc(prA1, truth=alert, event_level="second") |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              mcc(truth=alert, estimate=predMCC) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=schoenr(density(A0, na.rm=T), density(A1, na.rm=T))) |>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D")) |>
  # bind_rows(oos.ls$alert_L |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, model, PCA, covSet) |>
  #             f_meas(predF1, truth=alert, beta=1, event_level="second") |>
  #             group_by(y) |>
  #             mutate(rank=min_rank(desc(.estimate)),
  #                    .metric="F1") |>
  #             select(-.estimator)) |>
  # bind_rows(oos.ls$alert_L |>
  #             filter(!grepl("perfect|auto", model)) |>
  #             group_by(y, model, PCA, covSet) |>
#             kap(predMCC, truth=alert, event_level="second") |>
#             group_by(y) |>
#             mutate(rank=min_rank(desc(.estimate)),
#                    .metric="Kappa (MCC opt)") |>
#             select(-.estimator)) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(predF1=="A1"))|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="Precision: TP/(TP+FP) (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(alert=="A1"))|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="Recall: TP/(TP+FN) (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A1" & alert=="A0")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="FPR (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A1" & alert=="A1")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(desc(.estimate)),
#                    .metric="TPR (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A0" & alert=="A1")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="FNR (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predF1=="A0" & alert=="A0")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(desc(.estimate)),
#                    .metric="TNR (F1)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/sum(predMCC=="A1"))|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="Precision: TP/(TP+FP) (MCC)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/sum(alert=="A1"))|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="Recall: TP/(TP+FN) (MCC)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A1" & alert=="A0")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="FPR (MCC)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(desc(.estimate)),
#                    .metric="TPR (MCC)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A0" & alert=="A1")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="FNR (MCC)")) |>
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=sum(predMCC=="A0" & alert=="A0")/n())|>
#             group_by(y) |>
#             mutate(rank=min_rank(desc(.estimate)),
#                    .metric="TNR (MCC)")) |>
bind_rows(oos.ls$alert_L |>
            filter(!grepl("perfect|auto", model)) |>
            mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                   prA1=if_else(prA1==1, 1-1e-5, prA1),
                   alert=as.numeric(alert=="A1")) |>
            calc_R2(type="vz", y) |>
            rename(.estimate=R2) |>
            mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
            na.omit() |>
            group_by(y) |>
            mutate(rank=min_rank(desc(.estimate)),
                   .metric="R2-VZ_trunc")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y) |>
              rename(.estimate=R2) |>
              na.omit() |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="mf", y) |>
              rename(.estimate=R2) |>
              na.omit() |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-MF"))

saveRDS(rank.df, "out/clean/rank_oos.rds")



rankPrev.df <-  oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prevAlert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y, prevAlert) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, prevAlert) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  bind_rows(oos.ls$alert_L |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              roc_auc(prA1, truth=alert, event_level="second") |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prevAlert, prA1) %>%
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=schoenr(density(A0, na.rm=T), density(A1, na.rm=T))) |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              f_meas(predF1, truth=alert, beta=1, event_level="second") |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              mcc(truth=alert, estimate=predMCC) |>
              mutate(.estimate=if_else(is.na(.estimate), 0, .estimate)) |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="MCC") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(predF1=="A1"))|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Precision: TP/(TP+FP) (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(alert=="A1"))|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Recall: TP/(TP+FN) (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A0")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FPR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TPR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A0" & alert=="A1")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FNR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predF1=="A0" & alert=="A0")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TNR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/sum(predMCC=="A1"))|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Precision: TP/(TP+FP) (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/sum(alert=="A1"))|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Recall: TP/(TP+FN) (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A1" & alert=="A0")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FPR (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A1" & alert=="A1")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TPR (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A0" & alert=="A1")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FNR (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet, prevAlert) |>
              summarise(.estimate=sum(predMCC=="A0" & alert=="A0")/n())|>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TNR (MCC)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y, prevAlert) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="mf", y, prevAlert) |>
              rename(.estimate=R2) |>
              na.omit() |>
              group_by(y, prevAlert) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-MF"))

saveRDS(rankPrev.df, "out/clean/rankPrevA_oos.rds")








# variable importance -----------------------------------------------------


col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
col_resp <- c("lnN", "tl", "alert")
col_cmems <- readRDS("data/cmems_vars.rds")
col_wrf <- readRDS("data/wrf_vars.rds")
y_levels <- c("Alsp", "PSP", "Disp", "DSP", "Pssp", "Psde", "Psse", "ASP", "Kami")

varTypes <- list(
  spacetime=c("ydayCos", "ydaySin", "ydaySinXydayCos",
              "latz", "lonz", "lonzXlatz"),
  autoreg=c("lnNWt1", "lnNAvg1", "prAlertAvg1", "alert1A1",
            "lnNWt2", "lnNAvg2", "prAlertAvg2", "alert2A1"),
  prevYr=c("lnNPrevYr", "lnNAvgPrevYr", "prAlertPrevYr", "prAlertAvgPrevYr"),
  fetch="fetch",
  CMEMS=col_cmems, 
  WRF=col_wrf,
  UV_interact=c(
    paste("UWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
    paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
    paste("UWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
    paste("VWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
    paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
    paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
    paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
    paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
  ),
  hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
)
varTypes$auto_interact <- paste("lnNWt1", 
                                c(varTypes$autoreg[-1], varTypes$prevYr, 
                                  varTypes$fetch, varTypes$CMEMS, varTypes$WRF), 
                                sep="X")

varType_df <- imap_dfr(varTypes, ~tibble(VariableType=.y, Variable=.x)) |>
  mutate(varTypeClean=factor(VariableType, 
                             levels=names(varTypes),
                             labels=c("Space time", "Autoregression",
                                      "Previous year", "Fetch", "Biogeochemistry",
                                      "Weather", "Wind interaction",
                                      "HAB", "Autoregression interaction"))) |>
  mutate(varTypeClean=lvls_reorder(varTypeClean, c(1,2,9,3:8)))

ML.f <- dirf(glue("{base.dir}/model_fits"), "(Ridge|MARS|NN|RF|Boost)_vi.rds", recursive=T)
HBL.f <- dirf(glue("{base.dir}/model_fits"), "HBL1_vi.rds", recursive=T)

vi.df1 <- bind_rows(
  ML.f |> map_dfr(~readRDS(.x) |> mutate(f=.x) |>
                    select(Variable, Importance, f)),
  HBL.f |>
    map_dfr(~readRDS(.x) |> mutate(f=.x)) |>
    filter(Variable != "Intercept") |>
    rename(Importance=Estimate) |>
    select(Variable, Importance, f)
) |>
  mutate(covSet=paste0("d", str_split_fixed(str_split_fixed(f, "fits/", 2)[,2], "-", 2)[,1]),
         y=str_split_fixed(str_split_fixed(f, "/vi/", 2)[,2], "_alert", 2)[,1],
         model=str_split_fixed(str_split_fixed(f, "alert_", 2)[,2], "_vi", 2)[,1]) |>
  left_join(varType_df)



og.dir <- base.dir
enet.f <- dirf(glue("{og.dir}/model_fits"), "_Ridge.rds", recursive=T) |>
  grep("1[3-6]-", x=_, invert=T, value=T)
rf.f <- dirf(glue("{og.dir}/model_fits"), "_RF.rds", recursive=T) |>
  grep("1[3-6]-", x=_, invert=T, value=T) 
xgb.f <- dirf(glue("{og.dir}/model_fits"), "_Boost.rds", recursive=T) |>
  grep("1[3-6]-", x=_, invert=T, value=T)
HB.f <- dirf(glue("{og.dir}/model_fits"), "_HBL1.rds", recursive=T) |>
  grep("1[3-6]-", x=_, invert=T, value=T)
HB.f_s <- dirf(glue("{og.dir}/model_fits"), "_HBL1.rds", recursive=T) 

vi.df2 <- bind_rows(
  enet.f |> 
    map_dfr(~readRDS(.x) |> tidy() |> mutate(f=.x)) |> 
    filter(term!="(Intercept)") |>
    rename(Importance=estimate, Variable=term) |> select(-penalty), 
  rf.f |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> vip::vi(scale=T) |> mutate(f=.x)),
  xgb.f |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> vip::vi(scale=T) |> mutate(f=.x)),
  HB.f |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> fixef(probs=0.5) |> 
              as_tibble(rownames="Variable") |> mutate(f=.x) |>
              filter(Variable != "Intercept") |>
              rename(Importance=Estimate) |> select(-Est.Error, -Q50)),
  HB.f_s |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> summary() %>%
              magrittr::extract2("splines") |>
              as_tibble(rownames="Variable") |>
              mutate(Variable=c("ydaySinXydayCos", "lonzXlatz"),
                     f=.x) |>
              rename(Importance=Estimate) |>
              select(Variable, f, Importance))
) |>
  mutate(covSet=paste0("d", str_split_fixed(str_split_fixed(f, "fits/", 2)[,2], "-", 2)[,1]),
         y=str_split_fixed(str_split_fixed(f, "Del[0-1]/", 2)[,2], "_alert", 2)[,1],
         model=str_remove(str_split_fixed(f, "alert_", 2)[,2], ".rds")) |>
  left_join(varType_df)

wt.ls <- dirf(glue("{base.dir}/compiled"), "_wt.rds") |> 
  map_dfr(~readRDS(.x)$alert |> mutate(f=.x)) |>
  rename(run=model) |>
  filter(!grepl("PCA", run)) |>
  mutate(covSet=str_split_fixed(run, "\\.", 2)[,1],
         model=str_remove(str_split_fixed(run, "\\.", 2)[,2], "_alert_A1"),
         y=str_split_fixed(str_split_fixed(f, "compiled/", 2)[,2], "_wt", 2)[,1]) |>
  select(y, covSet, model, wt)

vi.df <- bind_rows(vi.df1, vi.df2) |>
  group_by(f) |>
  mutate(Importance=abs(Importance)/max(abs(Importance))*100) |>
  ungroup() |>
  left_join(wt.ls) |>
  mutate(ImpWt=Importance*wt) |>
  mutate(y=factor(y, levels=y_levels))
saveRDS(vi.df, "out/clean/vi_df.rds")

vi_pt_df <- vi.df |>
  group_by(y, Variable, VariableType, varTypeClean) |>
  summarise(mnImp=mean(ImpWt, na.rm=T)) |>
  group_by(y) |>
  mutate(mnImp=mnImp/max(mnImp)*100,
         rank=min_rank(desc(mnImp)),
         top=(rank)/max(rank) <= 0.1) |>
  ungroup()

vi_rng_df <- vi_pt_df |> group_by(y, VariableType, varTypeClean) |>
  reframe(mnImp=c(0, max(mnImp)))

vi_pt_df |>
  ggplot(aes(mnImp, y, colour=y)) + 
  geom_line(data=vi_rng_df, linewidth=0.25) +
  geom_point(aes(shape=top, size=top, alpha=top)) + 
  scale_colour_brewer(type="qual", palette="Paired") + 
  scale_shape_manual(values=c(1, 19), guide="none") +
  scale_size_manual(values=c(0.75, 1.75), guide="none") +
  scale_alpha_manual(values=c(0.3, 1), guide="none") +
  scale_x_continuous("Ensemble Relative Importance (%)", breaks=c(0,50,100)) +
  facet_grid(varTypeClean~., scales="free_y", space="free_y", switch="y",
             labeller=label_wrap_gen()) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position=c(0.9,0.2),
        legend.title=element_blank())
