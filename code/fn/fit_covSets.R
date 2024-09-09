


fit_covSet <- function(y_i, base.dir, covSet, mod, split_test=TRUE,
                       nTuneVal=2, prior_strength=1, ncores=4, 
                       responses=c(alert="alert")) {
 
  # covariate set / response info
  f <- covSet$f
  d <- covSet$d
  y <- covSet$y
  
  # directories
  fit.dir <- glue("{base.dir}/model_fits/{f}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{f}/")
  dir.create(ens.dir, recursive=T, showWarnings=F)
  dir.create(cv.dir, recursive=T, showWarnings=F)
  dir.create(out.dir, recursive=T, showWarnings=F)
  dir.create(glue("{fit.dir}/vi/"), recursive=T, showWarnings=F)
  
  # covariate type columns
  col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
  col_resp <- c("lnN", "tl", "alert")
  col_cmems <- readRDS("data/cmems_vars.rds")
  col_wrf <- readRDS("data/wrf_vars.rds")
  
  # All possible covariates
  all_covs <- list(
    spacetime=c("yday", "lat", "lon"),
    main=c(
      "fetch",
      "lnNWt1", "lnNAvg1", "prAlertAvg1", "alert1A1",
      "lnNWt2", "lnNAvg2", "prAlertAvg2", "alert2A1",
      "lnNPrevYr", "lnNAvgPrevYr", "prAlertPrevYr", "prAlertAvgPrevYr",
      col_cmems, col_wrf
    ),
    interact=c(
      paste("UWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
      paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
      paste("UWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
      paste("VWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
      paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
      paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
      paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
      paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
    ),
    hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
  )
  all_covs$interact <- c(all_covs$interact,
                         paste("lnNWt1", c(all_covs$main[-2]), sep="X"))
  # Identify excluded covariates
  covs_exclude <- get_excluded_cov_regex(covSet)
  
  # testing/training splits
  obs.ls <- map_dfr(dirf(base.dir, "data_.*_all.rds"), readRDS) |>
    filter(y==y) |>
    select(all_of(col_metadata), all_of(col_resp),
           "alert1", "alert2", any_of(unname(unlist(all_covs)))) |>
    mutate(across(starts_with("alert"), ~factor(.x)),
           across(starts_with("tl"), ~factor(.x, ordered=T))) |>
    group_by(obsid) |>
    slice_head(n=1) |>
    ungroup() |>
    select(where(~any(!is.na(.x)))) |>
    drop_na()
  
  set.seed(1003)
  if(test_prop > 0) {
    obs.split <- group_initial_split(obs.ls, group=year, prop=test_prop)
    saveRDS(obs.split, glue("{base.dir}/{y}_{covSet}_dataSplit.rds"))
    obs.train <- training(obs.split)
    obs.test <- testing(obs.split)
  } else {
    obs.train <- obs.ls
  }
  
  # prepare dataset
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  prepPCA.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude, TRUE))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)))
  dPCA.y <- list(train=map(prepPCA.ls, ~bake(.x, obs.train)))
  if(test_prop > 0) {
    d.y$test <- map(prep.ls, ~bake(.x, obs.test))
    dPCA.y$test <- map(prepPCA.ls, ~bake(.x, obs.test))
  }
  saveRDS(d.y, glue("{base.dir}/{y}_{covSet}_dy_testPct-{test_prop}.rds"))
  saveRDS(dPCA.y, glue("{base.dir}/{y}_{covSet}_dPCAy_testPct-{test_prop}.rds"))
  covs <- filter_corr_covs(all_covs, d.y, covs_exclude)
  covsPCA <- names(dPCA.y$train[[1]] |> select(starts_with("PC")))
  
  # formulas
  form.ls <- map(
    responses,
    ~list(ML=formula(glue("{.x} ~ .")),
          ML_PCA=formula(glue("{.x} ~ .")),
          HBL=make_HB_formula(.x, c(covs$main, covs$interact)),
          HBL_PCA=make_HB_formula(.x, covsPCA),
          HB_vars=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                               "{paste(unlist(covs), collapse='+')}")),
          HB_vars_PCA=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                                   "{paste(covsPCA, collapse='+')}"))
    )
  )
  
  cat("Starting", covSet, "for", y, ":", as.character(Sys.time()), "\n",
      file=glue("out/logs/d{d}_{y}_{mod}.log"))
  
  for(r in responses) {
    set.seed(1003)
    
    if(mod!="HB") { 
      # ML models
      if(.Platform$OS.type=="unix") {
        plan(multisession, workers=ncores)
      } else {
        plan(multicore, workers=ncores)
      }
      opts <- vfold_cv(d.y$train[[r]], strata=r)
      tunes <- list(nTuneVal) |> set_names(mod) 
    } else {
      # HB models
      priStr <- switch(
        prior_strength,
        "1"=list(r1=0.2, r2=2, hs1=0.5, hs2=0.6, b=0.75, de=0.3, i=1),
        "2"=list(r1=0.1, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i=2),
        "3"=list(r1=0.1, r2=1, hs1=5, hs2=0.5, b=0.5, de=0.05, i=3)
      )
      opts <- list(
        iter=1500,
        warmup=1000,
        chains=3,
        cores=ncores,
        ctrl=list(adapt_delta=0.8, max_treedepth=10),
        prior_i=priStr$i
      )
      tunes <- map(
        responses, 
        ~list(HBL=make_HB_priors(priStr, "HBL", .x, covs),
              HBL_PCA=make_HB_priors(priStr, "HBL", .x, covsPCA, PCA=T))
      )
    }
    
    # fit non-PCA covariates
    fit_candidate(mod, r, form.ls, d.y$train, opts, tunes, fit.dir, y)
    
    # fit PCA covariates
    set.seed(1003)
    if(mod != "HB") {
      opts <- vfold_cv(dPCA.y$train[[r]], strata=r)
    }
    fit_candidate(mod, r, form.ls, dPCA.y$train, opts, tunes, fit.dir, y, "_PCA")
  }
  plan(sequential)
  
  cat("Finished", covSet, "for", y, ":", as.character(Sys.time()), "\n")
}
