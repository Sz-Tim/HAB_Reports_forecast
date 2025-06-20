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
library(beeswarm)
library(plotly)
library(leaflet)
library(scales)
theme_set(theme_classic())

dirs <- list(proj=find_rstudio_root_file(),
             shiny=whereami::thisfile() |> dirname())
source(paste0(dirs$shiny, "/fn.R"))

scotland_sf <- st_read(paste0(dirs$proj, "/data/northAtlantic_footprint.gpkg"), 
                       quiet=TRUE) |>
  st_crop(xmin=45000, xmax=490000, 
          ymin=500000, ymax=1230000)
scot_bbox <- st_bbox(scotland_sf)

# monitoring targets
target_sets <- c("hab", "tox", "fish")[1:2]
targ_exclude <- c("AZP", "YTX", "Prli")
targ_i <- map_dfr(target_sets, 
                  ~read_csv(paste0(dirs$proj, "/data/i_", .x, ".csv"), show_col_types=F) |>
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
  arrange(targ_ordered) |>
  mutate(col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
               "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
               "#cab2d6"),
         fig_short=factor(fig_short, levels=unique(fig_short)))

tl_i <- map_dfr(target_sets, 
                  ~read_csv(paste0(dirs$proj, "/data/tl_thresholds_", .x, ".csv"), show_col_types=F) |>
                    mutate(type=.x)) |>
  filter(! abbr %in% targ_exclude) |>
  arrange(type, abbr) |>
  mutate(targ_ordered=factor(abbr,
                             levels=c("Alsp", "PSP", "Disp", "DSP",
                                      "Pssp", "Psde", "Psse", "ASP",  
                                      "Kami"),
                             labels=c("Alexandrium", "PSTs",
                                      "Dinophysis", "DSTs (OA/DTXs/PTXs)",
                                      paste("Pseudo-nitzschia", c("spp.", "del.", "ser.")), "DA", 
                                      "Karenia mikimotoi"))) |>
  arrange(targ_ordered) |>
  select(abbr, targ_ordered, tl, units, alert, min_ge) |>
  filter(!is.na(tl) & tl != "TL0") |>
  mutate(min_lnN=log1p(min_ge)) |>
  group_by(abbr, tl) |>
  slice_head(n=1)

# monitoring locations
site_i <- map_dfr(target_sets,
                  ~readRDS(paste0(dirs$proj, "/data/site_", .x, "_df.rds")) |>
                    mutate(type=.x))
site_sf <- site_i |>
  st_as_sf(coords=c("lon", "lat"), crs=27700, remove=FALSE)
site_wgs <- site_sf |>
  st_transform(4326) |>
  sevcheck::add_lonlat(drop_geom=T) |>
  group_by(sin) |>
  arrange(type) |>
  slice_head(n=1) |>
  ungroup()


# output and validation ---------------------------------------------------

# point-wise predictions
all_df <- readRDS(paste0(dirs$proj, "/out/clean/all_df.rds"))
# observations
obs_df <- readRDS(paste0(dirs$proj, "/out/clean/obs_df.rds"))
# skill scores: Overall
validation_df <- readRDS(paste0(dirs$proj, "/out/clean/validation_df.rds"))
# skill scores: By site
validation_sin_df <- readRDS(paste0(dirs$proj, "/out/clean/validation_sin_df.rds"))
# skill scores: By month
validation_month_df <- readRDS(paste0(dirs$proj, "/out/clean/validation_month_df.rds"))
# list of SINs by year for plotting
sin_yr <- obs_df |>
  mutate(year=year(week)) |>
  select(type, year, sin) |>
  group_by(year, type, sin) |>
  slice_head(n=1) |>
  group_by(year) |>
  group_split() 
sin_yr <- imap(unique(year(obs_df$week)),
               ~unique(filter(obs_df, year(week)==.x)$sin))


