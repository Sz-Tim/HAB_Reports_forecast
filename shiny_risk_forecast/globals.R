# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Shiny App: global libraries and objects



# setup -------------------------------------------------------------------

library(tidyverse)
library(shiny)
library(rprojroot)
library(sf)
library(cowplot)
library(ggdist)
library(scales)
theme_set(theme_classic())

dirs <- list(proj=find_rstudio_root_file(),
             shiny=whereami::thisfile() |> dirname())

scotland_sf <- st_read(paste0(dirs$proj, "/data/northAtlantic_footprint.gpkg"), 
                       quiet=TRUE) |>
  st_crop(xmin=45000, xmax=490000, 
          ymin=500000, ymax=1230000)

# monitoring targets
target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(paste0(dirs$proj, "/data/i_", .x, ".csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr)

# monitoring locations
site_i <- map_dfr(target_sets,
                  ~readRDS(paste0(dirs$proj, "/data/site_", .x, "_df.rds")) |>
                    mutate(type=.x))
site_sf <- site_i |>
  st_as_sf(coords=c("lon", "lat"), crs=27700, remove=FALSE)


# output and validation ---------------------------------------------------

# point-wise predictions
fit_df <- readRDS(paste0(dirs$proj, "/out/clean/out_fit_ens.rds"))
oos_df <- readRDS(paste0(dirs$proj, "/out/clean/out_oos_ens.rds"))
all_df <- bind_rows(fit_df, oos_df) |>
  group_by(y, siteid) |>
  mutate(N=n(), nYrs=n_distinct(year(date))) |>
  ungroup() |> filter(N >= 30 & nYrs > 1) |>
  arrange(y, siteid, date) |>
  select(y, siteid, date, obsid, prevAlert, alert, prA1) |>
  mutate(week=floor_date(date, unit="week")) |>
  left_join(targ_i |> select(abbr, type), by=join_by(y==abbr)) |>
  left_join(site_i |> select(siteid, sin, lon, lat, type), by=join_by(type, siteid))

obs_df <- bind_rows(
  readRDS(paste0(dirs$proj, "/data/0_init/old/data_hab_all.rds")) |>
    mutate(type="hab"),
  readRDS(paste0(dirs$proj, "/data/0_init/old/data_tox_all.rds")) |>
    mutate(type="tox")
) |>
  select(type, y, obsid, sin, date, N, lnN, tl, alert) |>
  filter(sin %in% unique(all_df$sin)) |>
  # inner_join(all_df |> select(y, obsid), by=join_by(y, obsid)) |>
  complete(type, y, sin, date, fill=list(N=0, lnN=0, tl="TL0", alert="A0")) |>
  mutate(week=floor_date(date, unit="week")) |>
  group_by(type, y, sin, week) |>
  summarise(lnN=mean(lnN)) |>
  group_by(y) |>
  mutate(lnN_rel=lnN/max(lnN)) |>
  ungroup()

# skill scores: Overall
validation_df <- readRDS(paste0(dirs$proj, "/out/clean/rank_oos.rds")) |>
  filter(model %in% c("Null[0]", "Null[Date]",
                      "HB", "Ridge", "MARS", "NN", "RF", "XGB", 
                      "Ensemble2")) |> 
  filter(.metric %in% c("MCC", "PR-AUC", "R2-VZ_trunc", "ROC-AUC", "Schoener's D")) |>
  mutate(.metric=if_else(.metric=="R2-VZ_trunc", "R2-VZ", .metric)) |>
  arrange(y, model, covSet, PCA) |>
  group_by(y, .metric) |>
  mutate(null_to_perfect=case_when(.metric=="Schoener's D" ~ first(.estimate),
                                   .default=1-first(.estimate)),
         m_to_null=case_when(.metric=="Schoener's D" ~ first(.estimate) - .estimate,
                             .default=.estimate - first(.estimate)),
         skill=m_to_null/null_to_perfect) |>
  ungroup() |>
  mutate(modType=case_when(model=="Null[0]" ~ model,
                           model=="Null[Date]" ~ model,
                           model=="Ensemble2" ~ "Ensemble",
                           .default="Constituent"),
         modType=factor(modType, levels=c("Null[0]", "Null[Date]", "Constituent", "Ensemble")),
         .metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("AUC['ROC']", "AUC['PR']", "MCC", "R['VZ']^2", "D['Overlap']"))) |>
  group_by(y, .metric, modType) |>
  mutate(modNum=row_number()) |>
  ungroup() |>
  select(y, modType, modNum, .metric, .estimate, skill)

# skill scores: By site


# skill scores: By date

