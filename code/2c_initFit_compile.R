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
