# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Initial dataset compilation


# setup -------------------------------------------------------------------
library(terra)
library(raster)
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
urls <- c(fsa="fsa_counts",
          fsa_sites="fsa_sites", 
          cefas="cefas_counts",
          cefas_sites="cefas_sites",
          mowi="mowi_counts",
          mowi_sites="mowi_sites",
          ssf="ssf_counts",
          ssf_sites="ssf_sites") |>
  map(~glue("http://www.habreports.org/dbdatastuff/{.x}"))
saveRDS(urls, "data/habreports_urls.rds")

target_sets <- c("hab", "tox", "fish")[1:2]
targ_i <- map(target_sets, ~read_csv(glue("data/i_{.x}.csv"))) |>
  set_names(target_sets)
targ_tl <- list(
  "hab"=read_csv("data/tl_thresholds_hab.csv") |>
    filter(min_ge != -99) |>
    group_by(abbr) |>
    mutate(alert=case_when(is.na(alert)~0,
                           alert=="warn"~1,
                           alert=="alert"~2),
           A=paste0("A", as.numeric(alert>0)),
           tl=factor(tl)) |>
    ungroup() |>
    select(abbr, min_ge, A, alert, tl),
  "tox"=read_csv("data/tl_thresholds_tox.csv") |>
    filter(min_ge != -99) |>
    group_by(abbr) |>
    mutate(alert=case_when(is.na(alert)~0,
                           alert=="warn"~1,
                           alert=="alert"~2),
           A=paste0("A", as.numeric(alert>0)),
           tl=factor(tl)) |>
    ungroup() |>
    select(abbr, min_ge, A, alert, tl) |>
    mutate(A=if_else(abbr %in% c("ASP", "AZAs", "YTXs") & as.numeric(tl)>1, "A1", A),
           alert=if_else(abbr %in% c("ASP", "AZAs", "YTXs") & as.numeric(tl)>1, 2, alert))
)
if("fish" %in% target_sets) {
  targ_tl$fish <- get_tl_info(targ_i$fish) |>
    select(abbr, min_ge, A, alert, tl)
}
if(!is.null(targ_tl$fish)) {
  kasp_rows <- which(targ_tl$fish$abbr=="Kasp" & targ_tl$fish$min_ge < 500000)
  targ_tl$fish$A[kasp_rows] <- "A0"
  targ_tl$fish$alert[kasp_rows] <- c(0, 1)
  targ_tl$fish$tl[kasp_rows] <- c("TL0", "TL2")
}



# sampling locations and dates --------------------------------------------

for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  # read and clean monitoring sites
  if(i == "fish") {
    sites <- read_and_clean_sites(urls$mowi_sites, dateStart) |>
      bind_rows(read_and_clean_sites(urls$ssf_sites, dateStart))
  } else {
    sites <- read_and_clean_sites(urls[[glue("{iSrc}_sites")]], dateStart)
  }
  # read and clean monitoring data
  dat.df <- read_and_clean_monitoring_data(i, urls, targ_i, sites, dateStart) |>
    mutate(siteid=as.numeric(factor(sin)))
  dat.df |> 
    select(-lon, -lat) |> 
    saveRDS(glue("data/0_init/{iSrc}_df.rds"))
  # save sites that align with actual data
  if(i == "fish") {
    site_df <- dat.df |>
      select(siteid, sin, lon, lat) |>
      group_by(siteid) |> slice_head(n=1) |> ungroup()
  } else {
    site_df <- dat.df |>
      select(siteid, sin, site, area, farm_species, lon, lat) |>
      group_by(siteid) |> slice_head(n=1) |> ungroup()
  }
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
}



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
         ID_toolbox=ID) # Access keeps changing........
         # ID_toolbox=glue("cmems_mod_nws_bgc", 
                         # "{if_else(source=='Reanalysis', paste0('-', var, '_my_7km'), '_anfc_0.027deg')}-3D_P1D-m"))
write_csv(cmems_i, "data/cmems_i.csv")

fsa.df <- readRDS("data/0_init/fsa_df.rds")
cefas.df <- readRDS("data/0_init/cefas_df.rds")
# fish.df <- readRDS("data/0_init/fish_df.rds")
get_CMEMS(userid=NULL, pw=NULL, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=nDays_avg, 
          dateRng=range(c(fsa.df$date, 
                          cefas.df$date)), 
                          # fish.df$date)), 
          out.dir="data/00_env/cmems/",
          toolbox=TRUE)

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

fsa.df <- readRDS("data/0_init/fsa_df.rds")
cefas.df <- readRDS("data/0_init/cefas_df.rds")
# fish.df <- readRDS("data/0_init/fish_df.rds")
wrf.dir <- ifelse(.Platform$OS.type=="unix",
                  "https",#"/media/archiver/common/sa01da-work/WRF/Archive/",
                  "E:/hydroOut/WRF/Archive/")
wrf.out <- "data/00_env/wrf/"
get_WRF(wrf.dir=wrf.dir, nDays_buffer=nDays_avg, 
        dateRng=c(ymd("2016-01-01"), ymd("2025-12-31")),
        # dateRng=range(c(fsa.df$date,
        #                 cefas.df$date)),
                        # fish.df$date)), 
        out.dir=wrf.out)
wrf.df <- aggregate_WRF(wrf.out, refreshStart="2015-05-13")
saveRDS(wrf.df, glue("data/0_init/wrf_end_{max(wrf.df$date)}.rds"))



# pairwise distances ------------------------------------------------------

for(i in target_sets) {
  site_df <- readRDS(glue("data/site_{i}_df.rds"))
  path.ls <- get_shortestPaths(ocean.path="data/ScotlandOcean_footprint.tif", 
                               site.df=site_df, 
                               site_savePath=glue("data/site_{i}_df.rds"))
  write_csv(path.ls$dist.df, glue("data/site_{i}_pairwise_distances.csv"))
  path.ls <- list(dist.df=read_csv(glue("data/site_{i}_pairwise_distances.csv")))
  path.ls$dist.df |>
    bind_rows(path.ls$dist.df |> 
                rename(destinations=origins, origins=destinations)) |>
    bind_rows(tibble(origins=1:nrow(site_df), 
                     destinations=1:nrow(site_df), 
                     distances=0)) |>
    filter(distances < 100e3) |> 
    dplyr::select(-distances) |> 
    group_by(origins) |> 
    nest(data=destinations) |>
    mutate(dest_c=c(data[[1]])) |> 
    dplyr::select(-data) |>
    ungroup() |>
    saveRDS(glue("data/site_{i}_neighbors_100km.rds"))
}



# fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1

for(i in target_sets) {
  site_df <- readRDS(glue("data/site_{i}_df.rds"))
  site_df <- site_df |>
    get_fetch("data/log10_eu200m1a.tif") |>
    get_openBearing("data/northAtlantic_footprint.gpkg", buffer=200e3)
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
}



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

for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  y.df <- calc_y_features(
    readRDS(glue("data/0_init/{iSrc}_df.rds")), 
    targ_i[[i]], targ_tl[[i]],
    readRDS(glue("data/site_{i}_neighbors_100km.rds"))
  )
  saveRDS(y.df, glue("data/0_init/{i}_obs.rds"))
}



# site buffers ------------------------------------------------------------

# Buffers for averaging environmental conditions
# CMEMS is coarse and does not resolve lochs, so some sites are not covered
# The buffer size was determined by trial and error to reduce averaging while
# ensuring all sites are represented.

for(i in target_sets) {
  site_df <- readRDS(glue("data/site_{i}_df.rds"))
  site_sf <- site_df |> 
    st_as_sf(coords=c("lon", "lat"), crs=27700) |>
    select(-sin) |>
    st_buffer(dist=100e3) |>
    split_to_NSEW()
  st_write(site_sf, glue("data/site_{i}_sf.gpkg"), append=F)
}



# HAB status for toxins ---------------------------------------------------

# Calculate average HAB densities surrounding each cefas site

habAvg_tox.df <- summarise_hab_states(
  site_tox.sf=st_read("data/site_tox_sf.gpkg") |>
    group_by(siteid) |> summarise(), 
  site_hab.sf=readRDS("data/site_hab_df.rds") |> 
    select(siteid, lon, lat) |> st_as_sf(coords=c("lon", "lat"), crs=27700), 
  tox.obs=readRDS("data/0_init/cefas_df.rds") |> select(obsid, siteid, date), 
  hab.df=readRDS("data/0_init/hab_obs.rds")
)
saveRDS(habAvg_tox.df, "data/0_init/tox_habAvg.rds")



# extract sites -----------------------------------------------------------

# . CMEMS -----------------------------------------------------------------
cmems_i <- list(all=c("chl", "no3", "o2", "ph", "phyc", "po4"))
cmems.df <- readRDS(last(dirf("data/0_init", "cmems_end.*rds")))
cmems.sf <- readRDS("data/00_env/cmems/coords_chl.rds") |>
  mutate(date=first(cmems.df$date)) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)

for(i in target_sets) {
  # filter cmems.df to dates needed
  dateMin_i <- min(readRDS(glue("data/0_init/{i}_obs.rds"))$date)
  cmems.df_i <- cmems.df |> filter(date >= (dateMin_i - 365))
  # find site point locations
  site_df <- readRDS(glue("data/site_{i}_df.rds")) |> select(-any_of("cmems_id"))
  site_df <- site_df |> find_nearest_feature_id(cmems.sf, "cmems_id")
  saveRDS(site_df, glue("data/site_{i}_df_CMEMSUPDATE.rds"))
  # extract point environment
  cmems.site <- extract_env_pts(site_df, cmems_i$all, 
                                cmems.df_i |> mutate(version=1), 
                                cmems_id, "cmems_id")
  saveRDS(cmems.site, glue("data/0_init/cmems_sitePt_{i}.rds"))
  rm(cmems.site); gc()
  # find site buffer locations
  site.buffer <- st_read(glue("data/site_{i}_sf.gpkg")) |>
    find_buffer_intersect_ids(cmems.sf, "cmems_id")
  # extract buffer environment
  cmems.buffer <- extract_env_buffers(site.buffer, cmems_i, 
                                      cmems.df_i, "cmems_id")
  saveRDS(cmems.buffer, glue("data/0_init/cmems_siteBufferNSEW_{i}.rds"))
  rm(cmems.buffer); gc()
}



# . WRF -------------------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/00_env/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/00_env/wrf", glue("domain_d0._{.x}")), readRDS) |>
                      arrange(res, i) |>
                      mutate(wrf_id=row_number()) |>
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/0_init/", "wrf_end_.*rds"))) 

for(i in target_sets) {
  # filter wrf.df to dates needed
  dateMin_i <- min(readRDS(glue("data/0_init/{i}_obs.rds"))$date)
  wrf.df_i <- wrf.df |> filter(date >= (dateMin_i - 365))
  # find site point locations
  site_df <- readRDS(glue("data/site_{i}_df.rds")) |> select(-starts_with("wrf_id"))
  site_df <- map(wrf_versions, ~site_df |> find_nearest_feature_id(.x, "wrf_id")) |>
    reduce(full_join, by=names(site_df), suffix=paste0(".", seq_along(wrf_versions)))
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
  # extract point environment
  site.versions <- grep("wrf_id", names(site_df), value=T)
  wrf.site <- extract_env_pts(site_df, wrf_i$all, wrf.df_i, wrf_id, site.versions)
  saveRDS(wrf.site, glue("data/0_init/wrf_sitePt_{i}.rds"))
  rm(wrf.site); rm(site_df); gc()
  # find site buffer locations
  site.buffer <- map(wrf_versions, 
                     ~st_read(glue("data/site_{i}_sf.gpkg")) |> 
                       find_buffer_intersect_ids(.x, "wrf_id")) |>
    reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
  # extract buffer environment
  wrf.buffer <- extract_env_buffers(site.buffer, wrf_i, 
                                    wrf.df_i, paste0("wrf_id.", 1:2))
  saveRDS(wrf.buffer, glue("data/0_init/wrf_siteBufferNSEW_{i}.rds"))
  rm(wrf.df_i); rm(wrf.buffer); rm(site.buffer); gc()
}



# yday averages -----------------------------------------------------------

for(i in target_sets) {
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

# combine datasets
dat.ls <- map(target_sets, ~load_datasets("0_init", .x, "1_current")) |> 
  set_names(target_sets)
iwalk(dat.ls, ~saveRDS(.x$compiled, glue("data/0_init/data_{.y}_all.rds")))

# identify all cmems/wrf variable names
grep("cmems_id|date|siteid|version",
     c(unlist(map(dat.ls, ~.x$cmems.pt |> names())),
       unlist(map(dat.ls, ~.x$cmems.buf |>
                  pivot_wider(names_from="quadrant",
                              values_from=-(1:3),
                              names_sep="Dir") |>
                  names()))) |>
       unique(),
     value=T, invert=T) |> 
  sort() |>
  saveRDS("data/cmems_vars.rds")
grep("wrf_id|date|siteid|version",
     c(unlist(map(dat.ls, ~.x$wrf.pt |> names())),
       unlist(map(dat.ls, ~.x$wrf.buf |>
                  pivot_wider(names_from="quadrant",
                              values_from=-(1:3),
                              names_sep="Dir") |>
                  names()))) |>
       unique(),
     value=T, invert=T) |> 
  sort() |>
  saveRDS("data/wrf_vars.rds")

# record max date for each dataset
obs_end <- map(dat.ls, ~max(.x$obs$date)) |>
  c(list(cmems=max(ymd(str_sub(dir("data/0_init", "cmems_end"), 11, 20))),
         wrf=max(ymd(str_sub(dir("data/0_init", "wrf_end"), 9, 18)))))
saveRDS(obs_end, "data/0_init/obs_end.rds")

# copy files from 0_init/ to 1_current/
write_to_current <- T
if(write_to_current) {
  file.copy(dirf("data/0_init/", "fsa_df"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "cefas_df"), "data/1_current/", overwrite=T)
  # file.copy(dirf("data/0_init/", "fish_df"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_obs.rds"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_habAvg.rds"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_sitePt_"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "_siteBufferNSEW_"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "data_.*_all.rds"), "data/1_current/", overwrite=T)
  file.copy(dirf("data/0_init/", "obs_end"), "data/1_current/", overwrite=T)
}

