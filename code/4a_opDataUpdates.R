# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update data




# setup -------------------------------------------------------------------
library(tidyverse)
library(habforecastr)
library(glue)
library(future)
library(sf)
library(jsonlite)

UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)

nDays_replace <- 14 # number of days to replace from previous dataset
urls <- readRDS("data/habreports_urls.rds")
old_end <- readRDS("data/1_current/obs_end.rds") |>
  map(~ymd(.x)-nDays_replace)

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



# monitoring data ---------------------------------------------------------

# new_sins <- map(target_sets, ~FALSE) |> set_names(target_sets)
for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  # read and clean monitoring sites
  sites_old <- readRDS(glue("data/site_{i}_df.rds"))
  if(i == "fish") {
    sites <- read_and_clean_sites(urls$mowi_sites, "2015-01-01") |>
      bind_rows(read_and_clean_sites(urls$ssf_sites, "2015-01-01")) |>
      filter(date > old_end[[i]])
  } else {
    sites <- read_and_clean_sites(urls[[glue("{iSrc}_sites")]], "2015-01-01") |>
      filter(date > old_end[[i]])
  }
  if(any(! unique(sites$sin) %in% sites_old$sin)) {
    # new_sins[[i]] <- TRUE
    cat("Warning: New SIN found for", iSrc, "! Models must be re-trained to include!\n")
  }
  # read and clean monitoring data
  dat.df <- read_and_clean_monitoring_data(i, urls, targ_i, sites, old_end[[i]]) |>
    inner_join(sites_old |> select(sin, siteid)) |>
    arrange(siteid, sin, date)
  dat.df |> 
    select(-lon, -lat) |> 
    saveRDS(glue("data/2_new/{iSrc}_df.rds"))
}




# CMEMS -------------------------------------------------------------------

cmems_i <- read_csv("data/cmems_i.csv") |> 
  filter(source=="AnalysisForecast") |>
  mutate(ID_toolbox=ID) # Access methods keep changing!
get_CMEMS(userid=NULL, pw=NULL, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=0, 
          dateRng=c(old_end$cmems-nDays_replace, today()+nDays_replace), 
          out.dir="data/00_env/test/",
          toolbox=TRUE,
          init_LU=readRDS(dir("data/00_env/cmems/", "coords.*rds", full.names=T)[1]))

# cmems_LU <- readRDS(dir("data/00_env/cmems/", "coords.*rds", full.names=T)[1]) 
cmems.f <- dir("data/00_env/test", "cmems.*.rds", full.names=T)
cmems.ls <- map(cmems.f, ~readRDS(.x)) 
cmems.df <- cmems.ls[[1]] |> mutate(chl=log1p(chl))
cmems.df$no3 <- log1p(cmems.ls[[2]]$no3)
cmems.df$o2 <- log1p(cmems.ls[[3]]$o2)
cmems.df$ph <- cmems.ls[[4]]$ph
cmems.df$phyc <- log1p(cmems.ls[[5]]$phyc)
cmems.df$po4 <- log1p(cmems.ls[[6]]$po4)
rm(cmems.ls)
saveRDS(cmems.df, glue("data/2_new/cmems_end_{max(cmems.df$date)}.rds"))
rm(cmems.df); gc()


# WRF ---------------------------------------------------------------------

wrf.dir <- ifelse(.Platform$OS.type=="unix",
                  "https",#"/media/archiver/common/sa01da-work/WRF/Archive/",
                  "E:/hydroOut/WRF/Archive/")
wrf.out <- "data/00_env/wrf/"

get_WRF(wrf.dir=wrf.dir, nDays_buffer=0, 
        dateRng=c(old_end$wrf, today()), 
        out.dir=wrf.out)

# use forecast when hindcast is unavailable
latest_wrf <- dir(wrf.out, "wrf_.*_d01.rds") |> sort() |> last() |> 
  str_sub(5, 14) |> ymd()
get_WRF(wrf.dir=wrf.dir, nDays_buffer=0, 
        dateRng=c(latest_wrf, today()+nDays_replace), 
        out.dir=wrf.out, forecast=T)

# TODO: Still designed for initial preparation
# Need to interpolate within each day instead - currently fills in using
# climatology by cell
wrf.df <- aggregate_WRF(wrf.out, refreshStart=old_end$wrf-nDays_replace, ncores=20)
saveRDS(wrf.df, glue("data/2_new/wrf_end_{max(wrf.df$date)}.rds"))



# autoregressive terms ----------------------------------------------------

max_env_date <- min(
  ymd(str_sub(last(dirf("data/2_new/", "cmems_end_.*rds")), -14, -5)),
  ymd(str_sub(last(dirf("data/2_new/", "wrf_end_.*rds")), -14, -5))
) 
for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  new_obs_df <- readRDS(glue("data/2_new/{iSrc}_df.rds"))
  min_new_date <- if_else(nrow(new_obs_df)>0, min(new_obs_df$date), today())
  if(file.exists("out/1_forecast/compiled/fcst_history_df.rds")) {
    last_obsid <- (readRDS("out/1_forecast/compiled/fcst_history_df.rds") |>
                     filter(y %in% targ_i[[i]]$abbr) |>
                     slice_max(obsid))$obsid[1]
  } else {
    last_obsid <- 1e6
  }
  combined_obs_df <- bind_rows(
    readRDS(glue("data/1_current/{iSrc}_df.rds")) |>
      # include previous year for prevYr calculations -> easy place to make more efficient...
      filter(between(date, ymd(paste(year(min_new_date)-1, "-01-01")), min_new_date)),
    new_obs_df)
  # This is hacky and repetitive but ok for now...
  days_to_forecast <- seq(today(), max_env_date, by=1)
  y_ls <- vector("list", length(days_to_forecast))
  for(j in seq_along(days_to_forecast)) {
    forecastDays_df <- combined_obs_df |>
      group_by(siteid) |>
      slice_head(n=1) |>
      ungroup() |>
      select(sin, site, area, farm_species, siteid) |>
      mutate(date=days_to_forecast[j]) |> 
      mutate(obsid=last_obsid + row_number())
    y_features_df <- calc_y_features(
      bind_rows(combined_obs_df, forecastDays_df), 
      targ_i[[i]], targ_tl[[i]],
      readRDS(glue("data/site_{i}_neighbors_100km.rds")),
      forecastStart=days_to_forecast[j]
    )
    if(j == 1) {
      y_ls[[j]] <- y_features_df |>  filter(date >= min_new_date) 
    } else {
      y_ls[[j]] <- y_features_df |>  filter(date >= first(days_to_forecast))
    }
    last_obsid <- max(forecastDays_df$obsid)
  }
  
  saveRDS(bind_rows(y_ls), glue("data/2_new/{i}_obs.rds"))
}



# HAB status for toxins ---------------------------------------------------

# TODO: This will fail if there is no new toxin data -- add ifelse like above
# TODO: Returns NaN if no hab data within the time frame... replace with 0? Or 
# just exclude? As-is, this means forecasts are not possible.
# Calculate average HAB densities surrounding each cefas site
if(nrow(readRDS("data/2_new/cefas_df.rds"))==0) {
  min_new_date_tox <- old_end$tox + nDays_replace
} else {
  min_new_date_tox <- min(readRDS("data/2_new/cefas_df.rds")$date)
}
hab.df <- bind_rows(
  readRDS("data/1_current/hab_obs.rds") |>
    filter(date > min_new_date_tox - 7*12),
  readRDS("data/2_new/hab_obs.rds")
) |>
  group_by(date, siteid, y) |>
  summarise(across(where(is.numeric), mean, na.rm=T), 
            across(where(is.factor) | where(is.character), first)) |>
  ungroup() 
tox.df <- bind_rows(
  readRDS(glue("data/1_current/cefas_df.rds")) |>
    filter(between(date, min_new_date_tox - 7*8, min_new_date_tox)) |>
    group_by(siteid) |>
    slice_max(date, n=3) |>
    ungroup() |>
    select(obsid, siteid, date),
  readRDS("data/2_new/tox_obs.rds") |> 
    filter(date >= min_new_date_tox) |> 
    group_by(obsid) |>
    slice_head(n=1) |>
    ungroup() |>
    select(obsid, siteid, date)
  )
habAvg_tox.df <- summarise_hab_states(
  site_tox.sf=st_read("data/site_tox_sf.gpkg") |>
    group_by(siteid) |> summarise(), 
  site_hab.sf=readRDS("data/site_hab_df.rds") |> 
    select(siteid, lon, lat) |> st_as_sf(coords=c("lon", "lat"), crs=27700), 
  tox.obs=tox.df, 
  hab.df=hab.df
) |>
  mutate(across(matches("lnNAvg|prA"), ~if_else(is.na(.x), 0, .x)))
saveRDS(habAvg_tox.df, "data/2_new/tox_habAvg.rds")



# extract sites -----------------------------------------------------------

# . CMEMS -----------------------------------------------------------------
cmems_i <- list(all=c("chl", "no3", "o2", "ph", "phyc", "po4"))
cmems.df <- readRDS(last(dirf("data/2_new", "cmems_end.*rds")))
cmems.sf <- readRDS("data/00_env/cmems/coords_chl.rds") |>
  mutate(date=first(cmems.df$date)) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)

for(i in target_sets) {
  # filter cmems.df to dates needed
  dateMin_i <- min(readRDS(glue("data/2_new/{i}_obs.rds"))$date)
  cmems.df_i <- cmems.df |> filter(date >= (dateMin_i - 365))
  site_df <- readRDS(glue("data/site_{i}_df.rds"))
  # extract point environment
  cmems.site <- extract_env_pts(site_df, cmems_i$all, 
                                cmems.df_i |> mutate(version=1), 
                                cmems_id, "cmems_id") |>
    drop_na() # removes first dates with no Wk, delta 
  saveRDS(cmems.site, glue("data/2_new/cmems_sitePt_{i}.rds"))
  # find site buffer locations
  site.buffer <- st_read(glue("data/site_{i}_sf.gpkg")) |>
    find_buffer_intersect_ids(cmems.sf, "cmems_id")
  # extract buffer environment
  cmems.buffer <- extract_env_buffers(site.buffer, cmems_i, 
                                      cmems.df_i, "cmems_id") |>
    drop_na() # removes first dates with no Wk, delta 
  saveRDS(cmems.buffer, glue("data/2_new/cmems_siteBufferNSEW_{i}.rds"))
}



# . WRF -------------------------------------------------------------------
# TODO: Is this currently accurate? Better to:
# 1) convert to interpolated rasters
# 2) apply a land mask
# 3) calculate zonal means
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/00_env/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/00_env/wrf", glue("domain_d0._{.x}")), readRDS) |>
                      arrange(res, i) |>
                      mutate(wrf_id=row_number()) |>
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/2_new/", "wrf_end_.*rds"))) 

for(i in target_sets) {
  # filter wrf.df to dates needed
  dateMin_i <- min(readRDS(glue("data/2_new/{i}_obs.rds"))$date)
  wrf.df_i <- wrf.df |> filter(date >= (dateMin_i - 365))
  # extract point environment
  site_df <- readRDS(glue("data/site_{i}_df.rds"))
  site.versions <- grep("wrf_id", names(site_df), value=T)
  wrf.site <- extract_env_pts(site_df, wrf_i$all, wrf.df_i, wrf_id, site.versions)
  saveRDS(wrf.site, glue("data/2_new/wrf_sitePt_{i}.rds"))
  # find site buffer locations
  site.buffer <- map(wrf_versions, 
                     ~st_read(glue("data/site_{i}_sf.gpkg")) |> 
                       find_buffer_intersect_ids(.x, "wrf_id")) |>
    reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
  # extract buffer environment
  wrf.buffer <- extract_env_buffers(site.buffer, wrf_i, 
                                    wrf.df_i, paste0("wrf_id.", 1:2))
  saveRDS(wrf.buffer, glue("data/2_new/wrf_siteBufferNSEW_{i}.rds"))
}



# yday averages -----------------------------------------------------------

for(i in target_sets) {
  # CMEMS points
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current/cmems_sitePt_{i}.rds")) |>
                           drop_na(),
                         readRDS(glue("data/2_new/cmems_sitePt_{i}.rds")) |>
                           drop_na()),
               glue("data/2_new/ydayAvg_cmems_sitePt_{i}.rds"),
               cmems_id, version, yday)
  # CMEMS buffers
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current/cmems_siteBufferNSEW_{i}.rds")) |>
                           drop_na(),
                         readRDS(glue("data/2_new/cmems_siteBufferNSEW_{i}.rds")) |>
                           drop_na()),
               glue("data/2_new/ydayAvg_cmems_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
  # WRF points
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current/wrf_sitePt_{i}.rds")) |>
                           drop_na(),
                         readRDS(glue("data/2_new/wrf_sitePt_{i}.rds")) |>
                           drop_na()),
               glue("data/2_new/ydayAvg_wrf_sitePt_{i}.rds"),
               wrf_id, version, yday)
  # WRF buffers
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current/wrf_siteBufferNSEW_{i}.rds")) |>
                           drop_na(),
                         readRDS(glue("data/2_new/wrf_siteBufferNSEW_{i}.rds")) |>
                           drop_na()),
               glue("data/2_new/ydayAvg_wrf_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
}


# compile -----------------------------------------------------------------

# combine datasets
dat.ls <- map(target_sets, ~load_datasets("2_new", .x, "2_new")) |> 
  set_names(target_sets)
iwalk(dat.ls, ~saveRDS(.x$compiled, glue("data/2_new/data_{.y}_all.rds")))

# record max date for each dataset
obs_end <- list(
  hab=max(dat.ls$hab$fsa$date),
  tox=max(dat.ls$tox$cefas$date),
  habF=max(dat.ls$hab$obs$date),
  toxF=max(dat.ls$tox$obs$date),
  cmems=max(ymd(str_sub(dir("data/2_new", "cmems_end"), -14, -5))),
  wrf=max(ymd(str_sub(dir("data/2_new", "wrf_end"), -14, -5)))
)
# Correct if any had no new observations
date_check <- map_lgl(obs_end, ~is.infinite(.x) | is.na(.x) | is.nan(.x))
if(any(date_check)) {
  obs_end_old <- readRDS("data/1_current/obs_end.rds")
  obs_end[date_check] <- obs_end_old[date_check]
}
saveRDS(obs_end, "data/2_new/obs_end.rds")



# apply recipe ------------------------------------------------------------

covSet.df <- read_csv("data/covSet_df.csv")
for(i in 1:nrow(covSet.df)) {
  make_forecast_data(imap_dfr(targ_i, ~.x |> mutate(type=.y)), 
                     covSet.df[i,], "data/2_new/")
}



# replace previous data ---------------------------------------------------

# copy files from 2_new/ to 1_current/
write_to_current <- T
if(write_to_current) {
  fnames <- dir("data/2_new") |> grep("_end|compiled|yday", x=_, invert=T, value=T)
  for(f in fnames) {
    current_df <- readRDS(glue("data/1_current/{f}"))
    new_df <- readRDS(glue("data/2_new/{f}"))
    if(any(grepl(".x$", c(names(current_df), names(new_df))))) {
      cat(f, " has issues :( \n")
    }
    bind_rows(current_df |> filter(date < min(new_df$date)),
              new_df) |>
      saveRDS(glue("data/1_current/{f}"))
  }
  fnames <- dir("data/2_new") |> grep("obs_end|yday", x=_, value=T)
  for(f in fnames) {
    file.copy(glue("data/2_new/{f}"), glue("data/1_current/{f}"), overwrite=T)
  }
}
