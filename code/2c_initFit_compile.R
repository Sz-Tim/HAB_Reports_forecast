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
         covSet=factor(covSet, levels=c(d_i$f, "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  mutate(covSet=factor(covSet, levels=c(d_i$f,
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
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
         covSet=factor(covSet, levels=c(d_i$f, "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  mutate(covSet=factor(covSet, levels=c(d_i$f,
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", "ensHB",
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
         # y=factor(y, levels=targ_i$abbr)) |>
  arrange(y, run, siteid, date) |>
  group_by(y, run, siteid) |>
  mutate(prevAlert=lag(alert)) |>
  ungroup()



# Threshold analysis ------------------------------------------------------

gc()
opt.F1 <- opt.mcc <- vector("list", n_distinct(fit.ls$alert_L$model))
for(i in seq_along(opt.F1)) {
  m_i <- unique(fit.ls$alert_L$model)[i]
  fit_i <- fit.ls$alert_L |> filter(model==m_i)
  thresh.fit <- compute_thresholds(fit_i, 
                                   0.001, 0.9, 0.025,
                                   # 0.001, 0.9, 0.0025, 
                                   byPrevAlert=!grepl("Null", m_i), 
                                   cores=12)
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
# bind_rows(oos.ls$alert_L |>
#             filter(!grepl("perfect|auto", model)) |>
#             select(y, covSet, PCA, model, obsid, alert, prA1) %>%
#             filter(!is.na(prA1)) |>
#             pivot_wider(names_from="alert", values_from="prA1") |>
#             group_by(y, model, PCA, covSet) |>
#             summarise(.estimate=schoenr(density(A0, na.rm=T), density(A1, na.rm=T))) |>
#             group_by(y) |>
#             mutate(rank=min_rank(.estimate),
#                    .metric="Schoener's D")) |>
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
                   .metric="R2-VZ_trunc")) 

saveRDS(rank.df, "out/clean/rank_oos.rds")


















ggplot(rank.df, aes(rank, model)) + 
  geom_boxplot() + 
  facet_grid(y~.metric)

ggplot(rank.df, aes(.estimate, model)) + 
  geom_boxplot() + 
  facet_grid(y~.metric, scales="free_x")





















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
library(butcher)
library(habforecastr)


train_prop <- 0.75
y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="tox")) |>
  filter(! abbr %in% c("AZP", "YTX", "Prli"))
d_ids <- paste0("d", str_pad(1:15, 2, "left", "0"))

cv_out <- vector("list", nrow(y_i))

for(i in 1:nrow(y_i)) {
  out_ls <- vector("list", 15)
  for(d in seq_along(d_ids)) {
    obs_df <- bind_rows(
      readRDS(glue("data/0_init/compiled/{y_i$abbr[i]}_{d_ids[d]}_dy_testPct-{train_prop}.rds"))$train$alert |>
        select(obsid, y, date, siteid, alert),
      readRDS(glue("data/0_init/compiled/{y_i$abbr[i]}_{d_ids[d]}_dy_testPct-{train_prop}.rds"))$test$alert |>
        select(obsid, y, date, siteid, alert))
    
    cv_df <- dirf(glue("out/0_init/model_fits/{d_ids[d]}/cv/"), y_i$abbr[i]) |>
      map(readRDS) |> 
      reduce(full_join, by=join_by(y, obsid)) |>
      rename_with(.cols=ends_with("_A1"), .fn=~paste0(d_ids[d], ".", .x))
    
    out_ls[[d]] <- inner_join(obs_df, cv_df, by=join_by(y, obsid))
  }
  
  cv_out[[i]] <- reduce(out_ls, full_join, by=join_by(y, obsid, date, siteid, alert))
}

out_df <- reduce(cv_out, bind_rows)

out_df |> select(ends_with("_A1")) |> as.matrix() |> cor(use="pairwise") |> image(col=viridis::turbo(10))

hist(c(out_df |> select(ends_with("_A1")) |> as.matrix() |> cor(use="pairwise")))

out_long <- out_df |> 
  pivot_longer(ends_with("_A1"))

out_long |> 
  ggplot(aes(value, fill=alert)) + 
  geom_density(alpha=0.5) +
  facet_wrap(~name, scales="free_y")

out_long |>
  group_by(date, siteid) |>
  summarise(mn=mean(value),
            lo=quantile(value, probs=0.05),
            hi=quantile(value, probs=0.95)) |>
  ggplot(aes(date, mn)) + 
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.25, colour=NA) +
  geom_line() + 
  facet_wrap(~siteid)

out_long |>
  group_by(y, date, siteid) |>
  summarise(value=mean(value, na.rm=T)) |>
  group_by(y, date) |>
  summarise(mn=mean(value),
            lo=quantile(value, probs=0.05),
            hi=quantile(value, probs=0.95)) |>
  ungroup() |>
  mutate(year=year(date)) |> 
  ggplot(aes(date, mn, ymin=lo, ymax=hi)) + 
  # geom_linerange(alpha=0.5) +
  geom_point() +
  facet_wrap(~y)


out_long |>
  filter(grepl("d01", name)) |>
  ggplot(aes(date, value, group=date)) + 
  geom_boxplot() + 
  facet_wrap(~siteid)




oos.ls$alert |> 
  summarise(across(ends_with("_A1"), ~average_precision_vec(alert, .x, event_level="second"))) |> 
  pivot_longer(ends_with("_A1")) |> 
  arrange(desc(value)) |> 
  mutate(name=factor(name, levels=unique(name))) |> 
  ggplot(aes(name, value)) + 
  geom_point() + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))


oos.ls$alert |> 
  summarise(across(ends_with("_A1"), 
                   ~average_precision_vec(alert, .x, event_level="second"))) |> 
  pivot_longer(ends_with("_A1")) |> 
  arrange(desc(value)) |> 
  mutate(d=str_split_fixed(str_split_fixed(name, "\\.", 2)[,1], "_", 2)[,1],
         mod=str_remove(name, "_alert_A1") |>
           str_remove("d[0-9][0-9]\\.") |>
           str_remove("PCA\\.")) |> 
  ggplot(aes(mod, value)) + 
  geom_point(shape=1, alpha=0.5) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert |> 
  summarise(across(ends_with("_A1"), 
                   ~average_precision_vec(alert, .x, event_level="second"))) |> 
  pivot_longer(ends_with("_A1")) |> 
  arrange(desc(value)) |> 
  mutate(d=str_split_fixed(str_split_fixed(name, "\\.", 2)[,1], "_", 2)[,1],
         mod=str_remove(name, "_alert_A1") |>
           str_remove("d[0-9][0-9]\\.") |>
           str_remove("PCA\\."),
         PCA=grepl("PCA", name)) |> 
  ggplot(aes(d, value, group=paste(mod, PCA), colour=mod, shape=PCA)) + 
  geom_point() +
  geom_line() +
  scale_colour_viridis_d(option="turbo", begin=0.1) + 
  scale_shape_manual(values=c(1, 19)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

