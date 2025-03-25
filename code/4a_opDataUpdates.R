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

nDays_replace <- 10 # number of days to replace from previous dataset
nDays_fcst <- 7 # number of days to forecast ahead from today
urls <- readRDS("data/habreports_urls.rds")
old_end <- readRDS("data/1_current/obs_end.rds") |>
  map(~ymd(.x)-nDays_replace)

target_sets <- c("hab", "tox", "habfish")
targ_i <- map(target_sets, ~read_csv(glue("data/i_{.x}.csv"))) |>
  set_names(target_sets)

dat_old_dir <- "data/1_current"
dat_old <- list(site_hab="data/site_hab_df.rds", 
                obs_hab=glue("{dat_old_dir}/hab_obs.rds"), 
                cmems.pt_hab=glue("{dat_old_dir}/cmems_sitePt_hab.rds"), 
                cmems.buf_hab=glue("{dat_old_dir}/cmems_siteBufferNSEW_hab.rds"), 
                wrf.pt_hab=glue("{dat_old_dir}/wrf_sitePt_hab.rds"), 
                wrf.buf_hab=glue("{dat_old_dir}/wrf_siteBufferNSEW_hab.rds"), 
                site_tox="data/site_tox_df.rds", 
                obs_tox=glue("{dat_old_dir}/tox_obs.rds"), 
                cmems.pt_tox=glue("{dat_old_dir}/cmems_sitePt_tox.rds"), 
                cmems.buf_tox=glue("{dat_old_dir}/cmems_siteBufferNSEW_tox.rds"), 
                wrf.pt_tox=glue("{dat_old_dir}/wrf_sitePt_tox.rds"), 
                wrf.buf_tox=glue("{dat_old_dir}/wrf_siteBufferNSEW_tox.rds"), 
                fsa=glue("{dat_old_dir}/fsa_df.rds"), 
                cefas=glue("{dat_old_dir}/cefas_df.rds")) |>
  map(readRDS)

dat_new <- imap(dat_old, ~NULL)

# monitoring data ---------------------------------------------------------

for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "habfish"="fish")
  # read and clean monitoring sites
  if(i == "habfish") {
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
  site_df <- dat.df |>
    select(siteid, sin, lon, lat) |>
    group_by(siteid) |> slice_head(n=1) |> ungroup()
  saveRDS(site_df, glue("data/site_{i}_df.rds"))
}




# CMEMS -------------------------------------------------------------------





# WRF ---------------------------------------------------------------------





# autoregressive terms ----------------------------------------------------






# HAB status for toxins ---------------------------------------------------





# extract sites -----------------------------------------------------------


# CMEMS  site:date --------------------------------------------------------


# CMEMS  buffer:date ------------------------------------------------------


# WRF  site:date ----------------------------------------------------------


# WRF  buffer:date --------------------------------------------------------




# yday averages -----------------------------------------------------------




# compile -----------------------------------------------------------------


