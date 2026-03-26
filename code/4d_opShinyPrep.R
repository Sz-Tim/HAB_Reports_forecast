# Project: HAB Reports Forecast
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Operational forecast: Update dataset for habreports




# setup -------------------------------------------------------------------
library(tidyverse)
library(jsonlite)
library(glue)
library(habforecastr)

urls <- readRDS("data/habreports_urls.rds")

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





# random functions....? ---------------------------------------------------

read_and_clean_sites2 <- function(url_sites, dateStart="2000-01-01") {
  library(tidyverse)
  paste0(url_sites, "?fromdate=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(east < 7e5 &
             north < 125e4 &
             !(east==0 & north==0) &
             sin != "-99" &
             !is.na(fromdate) & !is.na(todate) &
             fromdate != todate) |>
    mutate(fromdate=lubridate::date(fromdate), todate=lubridate::date(todate)) |>
    rowwise() |>
    mutate(date=list(seq(fromdate, todate, by=1))) |>
    ungroup() |>
    arrange(sin, fromdate) |>
    select(sin, site, area, farm_species, east, north, date) |>
    unnest(date) |>
    arrange(sin, date) |>
    group_by(sin, date) |>
    slice_head(n=1) |>
    ungroup()
}

read_and_clean_fsa2 <- function(url_fsa, hab_i, sites, dateStart="2000-01-01") {
  library(tidyverse)
  paste0(url_fsa, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(hab_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, site, area, farm_species, date, easting, northing, all_of(hab_i$full)) |>
    left_join(sites |> select(sin, east, north, date), by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    rename(all_of(setNames(hab_i$full, hab_i$abbr))) |>
    select(obsid, lon, lat, sin, site, area, farm_species, date, all_of(hab_i$abbr)) |>
    arrange(sin, date)
}


read_and_clean_cefas2 <- function(url_cefas, tox_i, sites, dateStart="2000-01-01") {
  library(tidyverse)
  paste0(url_cefas, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected) & sin != "-99") |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(tox_i$full), ~if_else(.x == -99, NA_real_, .x)),
           across(any_of(tox_i$full), ~if_else(.x < 0, 0, .x))) |>
    group_by(sin, date) |> slice_head(n=1) |> ungroup() |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, site, area, farm_species, date, easting, northing, all_of(tox_i$full)) |>
    left_join(sites |> select(sin, east, north, date), by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    filter(lat > 500000) |>
    rename(all_of(setNames(tox_i$full, tox_i$abbr))) |>
    select(obsid, lon, lat, sin, site, area, farm_species, date, all_of(tox_i$abbr)) |>
    arrange(sin, date)
}


read_and_clean_fish2 <- function(url_mowi, url_ssf, fish_i, sites, dateStart="2000-01-01") {
  library(tidyverse)
  bind_rows(url(paste0(url_mowi, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble(),
            url(paste0(url_ssf, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble()) |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(fish_i$full), ~if_else(is.na(.x) | .x < 0, 0, .x))) |>
    rowwise() |>
    mutate(TOTAL=sum(c_across(any_of(fish_i$full)))) |>
    ungroup() |>
    filter(TOTAL > 0) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, any_of(fish_i$full)) |>
    mutate(easting=if_else(easting==0 & northing==0, NA_real_, easting),
           northing=if_else(easting==0 & northing==0, NA_real_, northing)) |>
    left_join(sites |> select(sin, east, north, date), by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east, na.rm=T), lat=median(north, na.rm=T)) |> ungroup() |>
    filter(!is.na(lon)) |>
    rename(any_of(setNames(fish_i$full, fish_i$abbr))) |>
    select(obsid, lon, lat, sin, date, any_of(fish_i$abbr)) |>
    arrange(sin, date) |>
    filter(sin!=0) |>
    group_by(date, sin) |>
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
    ungroup()
}

read_and_clean_monitoring_data2 <- function(i, urls, targ_i, sites, dateStart="2000-01-01") {
  if(i == "hab") {
    dat.df <- read_and_clean_fsa2(urls$fsa, targ_i$hab, sites, dateStart)
  } else if(i == "tox") {
    dat.df <- read_and_clean_cefas2(urls$cefas, targ_i$tox, sites, dateStart)
  } else if(i == "habfish") {
    dat.df <- read_and_clean_fish2(urls$mowi, urls$ssf, targ_i$fish, sites, dateStart)
  }
  return(dat.df)
}





# compile updated data ----------------------------------------------------

for(i in target_sets) {
  iSrc <- switch(i, 
                 "hab"="fsa",
                 "tox"="cefas",
                 "fish"="fish")
  # read and clean monitoring sites
  if(i == "fish") {
    sites <- read_and_clean_sites2(urls$mowi_sites) |>
      bind_rows(read_and_clean_sites2(urls$ssf_sites))
  } else {
    sites <- read_and_clean_sites2(urls[[glue("{iSrc}_sites")]])
  }
  # read and clean monitoring data
  dat.df <- read_and_clean_monitoring_data2(i, urls, targ_i, sites) |>
    group_by(sin, date) |>
    slice_head(n=1) |>
    ungroup() |>
    pivot_longer(any_of(targ_i[[i]]$abbr), names_to="y", values_to="N") |>
    filter((!is.na(N))) |>
    mutate(lnN=log1p(N)) |>
    get_trafficLights(N, targ_tl[[i]]) |>
    arrange(y, sin, date) |>
    mutate(year=year(date),
           week=floor_date(date, "week"),
           date_std=ymd("2000-01-01") + yday(date)-1,
           type=i)
  dat.df |> 
    select(type, sin, date, year, week, date_std, y, N, lnN, tl)  |> 
    saveRDS(glue("shiny_risk_forecast/{iSrc}_df.rds"))
  # save sites that align with actual data
  site_df <- dat.df |>
    select(sin, site, area, farm_species, lon, lat) |>
    group_by(sin) |> slice_head(n=1) |> ungroup()
  saveRDS(site_df, glue("shiny_risk_forecast/site_{i}_df.rds"))
}


obs_df <- bind_rows(
  readRDS("shiny_risk_forecast/fsa_df.rds"),
  readRDS("shiny_risk_forecast/cefas_df.rds")
)
saveRDS(obs_df, "shiny_risk_forecast/obs_df.rds")
site_df <- bind_rows(
  readRDS("shiny_risk_forecast/site_hab_df.rds") |> mutate(type="hab"),
  readRDS("shiny_risk_forecast/site_tox_df.rds") |> mutate(type="tox")
)
saveRDS(site_df, "shiny_risk_forecast/site_allObs_df.rds")

site_i <- bind_rows(
  readRDS("data/site_hab_df.rds") |> mutate(type="hab"),
  readRDS("data/site_tox_df.rds") |> mutate(type="tox")
) |>
  select(type, siteid, sin, site, area, farm_species, lon, lat)

fcst_df <- readRDS("out/1_forecast/compiled/fcst_history_df.rds") |>
  slice_max(date_generated, by=c(y, siteid, date_forecast)) |>
  mutate(type=if_else(y %in% targ_i$hab$abbr, "hab", "tox")) |>
  rename(prA1=ensGLM2_alert_A1) |>
  left_join(site_i, by=join_by(type, siteid)) |>
  mutate(week=floor_date(date_forecast, "week"))
saveRDS(fcst_df, "shiny_risk_forecast/fcst_df.rds")
