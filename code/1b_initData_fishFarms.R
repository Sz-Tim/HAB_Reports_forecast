# HAB Forecast
# www.habreports.org
# Tim Szewczyk
# Initial dataset compilation


# setup -------------------------------------------------------------------

library(terra)
library(gdistance)
library(tidyverse)
library(glue)
library(lubridate)
library(ncdf4)
library(sf)
library(jsonlite)
library(WeStCOMS)
library(sevcheck)
library(habforecastr)

nDays_avg <- 14
dateStart <- "2015-01-01"
UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)
urls <- c(mowi="mowi_counts",
          mowi_sites="mowi_sites",
          ssf="ssf_counts",
          ssf_sites="ssf_sites") |>
  map(~glue("http://www.habreports.org/dbdatastuff/{.x}"))

# hab_i <- read_csv("data/i_hab.csv")
# tox_i <- read_csv("data/i_tox.csv")
fish_i <- read_csv("data/i_fish.csv")
tl_fish <- get_tl_info(fish_i) |>
  select(abbr, min_ge, A, alert, tl)





# sampling locations and dates --------------------------------------------

fish_sites <- bind_rows(
  read_and_clean_sites(urls$mowi_sites, dateStart),
  read_and_clean_sites(urls$ssf_sites, dateStart))
fish.df <- read_and_clean_fish(urls$mowi, urls$ssf, fish_i, fish_sites, dateStart) |>
  mutate(siteid=as.numeric(factor(sin)))
fish.df |> select(-lon, -lat) |>
  saveRDS("data/0_init/fish_df.rds")
site_fish.df <- fish.df |>
  select(siteid, sin, lon, lat) |>
  group_by(siteid) |> slice_head(n=1) |> ungroup()
saveRDS(site_fish.df, "data/site_fish_df.rds")



# CMEMS -------------------------------------------------------------------

# Atlantic-European North West Shelf Biogeochemistry Analysis and Forecast
# Reanalysis: 1993-present, but updated every 6 months
#  - https://doi.org/10.48670/moi-00058
# Analysis & Forecast: 2019-present, updated daily with 6-day forecast
#  - https://doi.org/10.48670/moi-00056

cmems_i <- expand_grid(
  var=c("chl", # Mass concentration of chlorophyll a
        "no3", # Mole concentration of nitrate
        "o2", # Mole concentration of dissolved molecular oxygen
        "ph", # Sea water ph reported on total scale
        "phyc", # Mole concentration of phytoplankton expressed as carbon
        "po4" # Mole concentration of phosphate
  ),
  source=c("Reanalysis", "AnalysisForecast")) |>
  mutate(server=if_else(source=="Reanalysis", "my.cmems-du.eu", "nrt.cmems-du.eu"), 
         doi=glue("https://doi.org/10.48670/moi-0005{if_else(source=='Reanalysis', 8, 6)}"),
         ID=glue("cmems_mod_nws_bgc-{var}_", 
                 "{if_else(source=='Reanalysis', 'my', 'anfc')}_7km-3D_P1D-m"),
         ID_toolbox=glue("cmems_mod_nws_bgc", 
                         "{if_else(source=='Reanalysis', paste0('-', var, '_my_7km'), '_anfc_0.027deg')}-3D_P1D-m"))
write_csv(cmems_i, "data/cmems_i.csv")

fish.df <- readRDS("data/0_init/fish_df.rds")
get_CMEMS(userid=NULL, pw=NULL, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=nDays_avg, dateRng=range(fish.df$date), 
          out.dir="data/00_env/cmems/",
          toolbox=TRUE)

cmems_LU <- readRDS(dir("data/00_env/cmems/", "coords.*rds", full.names=T)[1]) 
cmems.f <- dir("data/00_env/cmems", "cmems.*.rds", full.names=T)
cmems.ls <- map(cmems.f, ~readRDS(.x)) 
cmems.df <- cmems.ls[[1]] |> mutate(chl=log1p(chl))
cmems.df$no3 <- log1p(cmems.ls[[2]]$no3)
cmems.df$o2 <- log1p(cmems.ls[[3]]$o2)
cmems.df$ph <- cmems.ls[[4]]$ph
cmems.df$phyc <- log1p(cmems.ls[[5]]$phyc)
cmems.df$po4 <- log1p(cmems.ls[[6]]$po4)
rm(cmems.ls)
saveRDS(cmems.df, glue("data/0_init/cmems_end_{max(cmems.df$date)}.rds"))



# WeStCOMS-WRF ------------------------------------------------------------

# There are three different resolutions / domains, from coarse but expansive
# (d01) to finer but restricted (d03). The getWRF() function nests all three, 
# selecting the highest resolution. Setting wrf.dir="https" will download the
# files from the public SAMS THREDDS server.
#  - Wind speed
#  - Wind direction
#  - Shortwave radiation
#  - Precipitation
#  - Sea surface temperature (not used due to extensive NAs)

fish.df <- readRDS("data/0_init/fish_df.rds")
wrf.dir <- ifelse(.Platform$OS.type=="unix",
                  "https",#"/media/archiver/common/sa01da-work/WRF/Archive/",
                  "E:/hydroOut/WRF/Archive/")
wrf.out <- "data/00_env/wrf/"
get_WRF(wrf.dir=wrf.dir, nDays_buffer=nDays_avg, 
        dateRng=range(fish.df$date), 
        out.dir=wrf.out)
wrf.df <- aggregate_WRF(wrf.out)
saveRDS(wrf.df, glue("data/0_init/wrf_end_{max(wrf.df$date)}.rds"))



# pairwise distances ------------------------------------------------------

site_fish.df <- readRDS("data/site_fish_df.rds")
path.ls <- get_shortestPaths(ocean.path="data/ScotlandOcean_footprint.tif", 
                             site.df=site_fish.df, 
                             site_savePath="data/site_fish_df.rds")
write_csv(path.ls$dist.df, "data/site_fish_pairwise_distances.csv")
path.ls <- list(dist.df=read_csv("data/site_fish_pairwise_distances.csv"))
path.ls$dist.df |>
  bind_rows(path.ls$dist.df |> rename(destinations=origins, origins=destinations)) |>
  complete(origins, destinations, fill=list(distances=0)) |>
  filter(distances < 100e3) |> 
  dplyr::select(-distances) |> 
  group_by(origins) |> 
  nest(data=destinations) |>
  mutate(dest_c=c(data[[1]])) |> 
  dplyr::select(-data) |>
  ungroup() |>
  saveRDS("data/site_fish_neighbors_100km.rds")



# fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1

site_fish.df <- readRDS("data/site_fish_df.rds")
site_fish.df <- site_fish.df |>
  get_fetch("data/log10_eu200m1a.tif") |>
  get_openBearing("data/northAtlantic_footprint.gpkg", buffer=200e3)
saveRDS(site_fish.df, "data/site_fish_df.rds")



# autoregressive terms ----------------------------------------------------

# Observed densities and bloom states for focal HAB taxa and toxins
# * = [t]
# *1 = [t-1]
# *2 = [t-2]
# N: reported density or concentration
# lnN: ln(N + 1)
# tl: HABReports traffic light for density or concentration N
# alert: HABReports action (0_none, 1_warn, 2_alert)
# lnNAvg: regional average of lnN within previous week
# prAlertAvg: regional average of alerting sites within previous week

fish.df <- calc_y_features(readRDS("data/0_init/fish_df.rds"),
                           fish_i, tl_fish,
                           readRDS("data/site_fish_neighbors_100km.rds"))
saveRDS(fish.df, "data/0_init/fish_obs.rds")



# site buffers ------------------------------------------------------------

# Buffers for averaging environmental conditions
# CMEMS is coarse and does not resolve lochs, so some sites are not covered
# The buffer size was determined by trial and error to reduce averaging while
# ensuring all sites are represented.

site_fish.df <- readRDS("data/site_fish_df.rds")
site_fish.sf <- site_fish.df |> 
  st_as_sf(coords=c("lon", "lat"), crs=27700) |>
  select(-sin) |>
  st_buffer(dist=100e3) |>
  split_to_NSEW()
st_write(site_fish.sf, "data/site_fish_sf.gpkg", append=F)



# extract sites -----------------------------------------------------------

# . CMEMS  site:date ------------------------------------------------------
cmems_i <- list(all=c("chl", "no3", "o2", "ph", "phyc", "po4"))
cmems.df <- readRDS(glue("data/0_init/cmems_end_{max(fish.df$date)+nDays_avg}.rds"))
cmems.sf <- readRDS("data/00_env/cmems/coords_chl.rds") |>
  mutate(date=first(cmems.df$date)) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)

site_fish.df <- readRDS("data/site_fish_df.rds") |> select(-any_of("cmems_id"))
site_fish.df <- site_fish.df |> find_nearest_feature_id(cmems.sf, "cmems_id")
saveRDS(site_fish.df, "data/site_fish_df.rds")
cmems.site_fish <- extract_env_pts(site_fish.df, cmems_i$all, cmems.df |> mutate(version=1), cmems_id, "cmems_id")
saveRDS(cmems.site_fish, "data/0_init/cmems_sitePt_fish.rds")



# . CMEMS  buffer:date ----------------------------------------------------

site.buffer_fish <- st_read("data/site_fish_sf.gpkg") |>
  find_buffer_intersect_ids(cmems.sf, "cmems_id")
cmems.buffer_fish <- extract_env_buffers(site.buffer_fish, cmems_i, cmems.df, "cmems_id")
saveRDS(cmems.buffer_fish, "data/0_init/cmems_siteBufferNSEW_fish.rds")



# . WRF  site:date --------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/00_env/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/00_env/wrf", glue("domain_d0._{.x}")), readRDS) |>
                      arrange(res, i) |>
                      mutate(wrf_id=row_number()) |>
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/0_init/", "wrf_end_.*rds"))) 

# HABs
site_fish.df <- readRDS("data/site_fish_df.rds") |> select(-starts_with("wrf_id"))
site_fish.df <- map(wrf_versions, ~site_fish.df |> find_nearest_feature_id(.x, "wrf_id")) |>
  reduce(full_join, by=names(site_fish.df), suffix=paste0(".", seq_along(wrf_versions)))
saveRDS(site_fish.df, "data/site_fish_df.rds")
site_fish.versions <- grep("wrf_id", names(site_fish.df), value=T)
wrf.site_fish <- extract_env_pts(site_fish.df, wrf_i$all, wrf.df, wrf_id, site_fish.versions)
saveRDS(wrf.site_fish, "data/0_init/wrf_sitePt_fish.rds")



# . WRF  buffer:date ------------------------------------------------------
# HABs
site.buffer_fish <- map(wrf_versions, 
                       ~st_read("data/site_fish_sf.gpkg") |> 
                         find_buffer_intersect_ids(.x, "wrf_id")) |>
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_fish <- extract_env_buffers(site.buffer_fish, wrf_i, wrf.df, paste0("wrf_id.", 1:2))
saveRDS(wrf.buffer_fish, "data/0_init/wrf_siteBufferNSEW_fish.rds")



# yday averages -----------------------------------------------------------

for(i in c("fish")) {
  # CMEMS points
  calc_ydayAvg(readRDS(glue("data/0_init/cmems_sitePt_{i}.rds")),
               glue("data/1_current/ydayAvg_cmems_sitePt_{i}.rds"),
               cmems_id, version, yday)
  # CMEMS buffers
  calc_ydayAvg(readRDS(glue("data/0_init/cmems_siteBufferNSEW_{i}.rds")),
               glue("data/1_current/ydayAvg_cmems_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
  # WRF points
  calc_ydayAvg(readRDS(glue("data/0_init/wrf_sitePt_{i}.rds")),
               glue("data/1_current/ydayAvg_wrf_sitePt_{i}.rds"),
               wrf_id, version, yday)
  # WRF buffers
  calc_ydayAvg(readRDS(glue("data/0_init/wrf_siteBufferNSEW_{i}.rds")),
               glue("data/1_current/ydayAvg_wrf_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
}



# compile -----------------------------------------------------------------

fish.ls <- load_datasets("0_init", "fish", "1_current")
saveRDS(fish.ls$compiled, "data/0_init/data_fish_all.rds")

# variable names
grep("cmems_id|date|siteid|version", 
     unique(c(names(fish.ls$cmems.pt), 
              names(fish.ls$cmems.buf |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")))),
     value=T, invert=T) |> 
  sort() |>
  saveRDS("data/cmems_vars.rds")
grep("wrf_id|date|siteid|version", 
     unique(c(names(fish.ls$wrf.pt), 
              names(fish.ls$wrf.buf |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")))),
     value=T, invert=T) |> 
  sort() |>
  saveRDS("data/wrf_vars.rds")

obs_end <- list(
  fish=max(fish.ls$obs$date),
  cmems=max(ymd(str_sub(dir("data/0_init", "cmems_end"), 11, 20))),
  wrf=max(ymd(str_sub(dir("data/0_init", "wrf_end"), 9, 18)))
)
saveRDS(obs_end, "data/0_init/obs_end.rds")

write_to_current <- T
if(write_to_current) {
  file.copy(dirf("data/0_init/", "fish_df"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_obs.rds"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_sitePt_"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_siteBufferNSEW_"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "data_.*_all.rds"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "obs_end"), "data/1_current/", overwrite=T)
}

