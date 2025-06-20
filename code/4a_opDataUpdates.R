# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update data




# setup -------------------------------------------------------------------
library(tidyverse)
library(habforecastr)
library(glue)
library(sf)
library(jsonlite)

UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)

nDays_replace <- 14 # number of days to replace from previous dataset
urls <- readRDS("data/habreports_urls.rds")
old_end <- readRDS("data/1_current_new/obs_end.rds") |>
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

cmems_i <- read_csv("data/cmems_i.csv")
get_CMEMS(userid=NULL, pw=NULL, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=0, 
          dateRng=c(old_end$cmems, today()+nDays_replace), 
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
saveRDS(cmems.df, glue("data/2_new/cmems_end_{max(cmems.df$date)}.rds"))



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

wrf.df <- aggregate_WRF(wrf.out, refreshStart=old_end$wrf)
saveRDS(wrf.df, glue("data/2_new/wrf_end_{max(wrf.df$date)}.rds"))



# autoregressive terms ----------------------------------------------------

for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  min_new_date <- min(readRDS(glue("data/2_new/{iSrc}_df.rds"))$date)
  y.df <- calc_y_features(
    bind_rows(readRDS(glue("data/1_current_new/{iSrc}_df.rds")) |>
                filter(date < min_new_date) |>
                group_by(siteid) |>
                slice_max(date, n=3),
              readRDS(glue("data/2_new/{iSrc}_df.rds"))), 
    targ_i[[i]], targ_tl[[i]],
    readRDS(glue("data/site_{i}_neighbors_100km.rds"))
  )
  saveRDS(y.df, glue("data/2_new/{i}_obs.rds"))
}



# HAB status for toxins ---------------------------------------------------

# Calculate average HAB densities surrounding each cefas site
habAvg_tox.df <- summarise_hab_states(
  site_tox.sf=st_read("data/site_tox_sf.gpkg") |>
    group_by(siteid) |> summarise(), 
  site_hab.sf=readRDS("data/site_hab_df.rds") |> 
    select(siteid, lon, lat) |> st_as_sf(coords=c("lon", "lat"), crs=27700), 
  tox.obs=readRDS("data/2_new/cefas_df.rds") |> select(obsid, siteid, date), 
  hab.df=readRDS("data/2_new/hab_obs.rds")
)
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
  # find site point locations
  site_df <- readRDS(glue("data/site_{i}_df.rds")) |> select(-any_of("cmems_id"))
  site_df <- site_df |> find_nearest_feature_id(cmems.sf, "cmems_id")
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
  # extract point environment
  cmems.site <- extract_env_pts(site_df, cmems_i$all, 
                                cmems.df_i |> mutate(version=1), 
                                cmems_id, "cmems_id")
  saveRDS(cmems.site, glue("data/2_new/cmems_sitePt_{i}.rds"))
  # find site buffer locations
  site.buffer <- st_read(glue("data/site_{i}_sf.gpkg")) |>
    find_buffer_intersect_ids(cmems.sf, "cmems_id")
  # extract buffer environment
  cmems.buffer <- extract_env_buffers(site.buffer, cmems_i, 
                                      cmems.df_i, "cmems_id")
  saveRDS(cmems.buffer, glue("data/2_new/cmems_siteBufferNSEW_{i}.rds"))
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
wrf.df <- readRDS(last(dirf("data/2_new/", "wrf_end_.*rds"))) 

for(i in target_sets) {
  # filter wrf.df to dates needed
  dateMin_i <- min(readRDS(glue("data/2_new/{i}_obs.rds"))$date)
  wrf.df_i <- wrf.df |> filter(date >= (dateMin_i - 365))
  # find site point locations
  site_df <- readRDS(glue("data/site_{i}_df.rds")) |> select(-starts_with("wrf_id"))
  site_df <- map(wrf_versions, ~site_df |> find_nearest_feature_id(.x, "wrf_id")) |>
    reduce(full_join, by=names(site_df), suffix=paste0(".", seq_along(wrf_versions)))
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
  # extract point environment
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
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current_new/cmems_sitePt_{i}.rds")),
                         readRDS(glue("data/2_new/cmems_sitePt_{i}.rds"))),
               glue("data/2_new/ydayAvg_cmems_sitePt_{i}.rds"),
               cmems_id, version, yday)
  # CMEMS buffers
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current_new/cmems_siteBufferNSEW_{i}.rds")),
                         readRDS(glue("data/2_new/cmems_siteBufferNSEW_{i}.rds"))),
               glue("data/2_new/ydayAvg_cmems_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
  # WRF points
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current_new/wrf_sitePt_{i}.rds")),
                         readRDS(glue("data/2_new/wrf_sitePt_{i}.rds")),),
               glue("data/2_new/ydayAvg_wrf_sitePt_{i}.rds"),
               wrf_id, version, yday)
  # WRF buffers
  calc_ydayAvg(bind_rows(readRDS(glue("data/1_current_new/wrf_siteBufferNSEW_{i}.rds")),
                         readRDS(glue("data/2_new/wrf_siteBufferNSEW_{i}.rds"))),
               glue("data/2_new/ydayAvg_wrf_siteBufferNSEW_{i}.rds"),
               siteid, quadrant, yday)
}


# compile -----------------------------------------------------------------

# combine datasets
dat.ls <- map(target_sets, ~load_datasets("2_new", .x, "2_new")) |> 
  set_names(target_sets)
iwalk(dat.ls, ~saveRDS(.x$compiled, glue("data/2_new/data_{.y}_all.rds")))

# record max date for each dataset
obs_end <- map(dat.ls, ~max(.x$obs$date)) |>
  c(list(cmems=max(ymd(str_sub(dir("data/2_new", "cmems_end"), 11, 20))),
         wrf=max(ymd(str_sub(dir("data/2_new", "wrf_end"), 9, 18)))))
saveRDS(obs_end, "data/2_new/obs_end.rds")



# apply recipe ------------------------------------------------------------

covSet.df <- read_csv("data/covSet_hab_tox.csv")
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
    current_df <- readRDS(glue("data/1_current_new/{f}"))
    new_df <- readRDS(glue("data/2_new/{f}"))
    bind_rows(current_df |> filter(date < min(new_df$date)),
              new_df) |>
      saveRDS(glue("data/1_current_new/{f}"))
  }
  fnames <- dir("data/2_new") |> grep("obs_end|yday", x=_, value=T)
  for(f in fnames) {
    file.copy(glue("data/2_new/{f}"), glue("data/1_current_new/{f}"), overwrite=T) 
  }
}
