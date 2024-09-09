
#' Create recipe and prepare using training data
#'
#' @param train.df 
#' @param response 
#' @param dimReduce 
#'
#' @return
#' @export
#'
#' @examples
prep_recipe <- function(train.df, response, covsExclude="NA", dimReduce=FALSE) {
  respExclude <- grep(response, c("lnN", "tl", "alert"), value=T, invert=T)
  pred_vars <- names(train.df)
  include_UVX <- !grepl("Xfetch", covsExclude)
  include_lnNX <- !grepl("lnNWt1X", covsExclude)
  rec <- recipe(train.df) |>
    step_mutate(prevAlert=alert1, role="ID") |>
    update_role(all_of(pred_vars), new_role="predictor") |>
    update_role(all_of(response), new_role="outcome") |>
    update_role(obsid, y, date, siteid, year, new_role="ID") |>
    update_role(lon, lat, yday, new_role="RE")  |>
    step_select(-any_of(respExclude)) |> 
    step_dummy(all_factor_predictors()) |>
    step_logit(starts_with("prAlert"), offset=0.01) |>
    step_logit(ends_with("A1"), offset=0.01) |>
    step_mutate_at(lon, lat, fn=list(z=~.)) |>
    step_interact(terms=~lonz:latz, sep="X")
  if(include_UVX) {
    rec <- rec |>
      step_interact(terms=~UWk:fetch:matches("Dir[EW]"), sep="X") |>
      step_interact(terms=~VWk:fetch:matches("Dir[NS]"), sep="X") |>
      step_interact(terms=~UWk:fetch:matches("Dir[NS]"), sep="X") |>
      step_interact(terms=~VWk:fetch:matches("Dir[EW]"), sep="X")
  }
  if(include_lnNX) {
    rec <- rec |>
      step_interact(terms=~lnNWt1:all_predictors(), sep="X") 
  }
  rec <- rec |>
    step_YeoJohnson(all_predictors()) |>
    step_normalize(all_predictors()) |>
    step_harmonic(yday, frequency=1, cycle_size=365, keep_original_cols=T) |>
    step_rename(ydayCos=yday_cos_1, ydaySin=yday_sin_1) |>
    step_interact(terms=~ydaySin:ydayCos, sep="X") |>
    step_bs(lonz, deg_free=tune()) |>
    step_bs(latz, deg_free=tune()) |>
    step_bs(lonzXlatz, deg_free=tune()) |>
    step_rename_at(contains("_"), fn=~str_remove_all(.x, "_")) |>
    step_select(-matches(covsExclude))
  if(dimReduce) {
    rec <- rec |>
      step_pca(all_predictors(), threshold=tune())
  }
  rec |>
    prep(training=train.df)
}








#' Filter list of covariates following recipe thinning
#'
#' @param all_covs 
#' @param data.y 
#' @param covsExclude 
#'
#' @return
#' @export
#'
#' @examples
filter_corr_covs <- function(all_covs, data.y, covsExclude="NA") {
  uncorr_covs <- unique(unlist(map(data.y, ~unlist(map(.x, names)))))
  if(any(grepl("^PC", uncorr_covs))) {
    covs_incl <- list(main=grep("^PC", uncorr_covs, value=T),
                      interact=NULL,
                      nonHB=NULL)
  } else {
    uncorr_covs <- unique(c(uncorr_covs, "yday", "ydayCos", "ydaySin", "lon", "lat"))
    covs_incl <- map(all_covs, ~.x[.x %in% uncorr_covs])
  }
  
  covs_incl <- map(covs_incl, ~grep(covsExclude, .x, value=T, invert=T))
  return(covs_incl)
}







#' Get regex for excluded covariates given a covariate set
#'
#' @param covSet 
#'
#' @return
#' @export
#'
#' @examples
get_excluded_cov_regex <- function(covSet) {
  covSet.df <- expand_grid(Avg=c(0,1), 
                           Xf=c(0,1),
                           XN=c(0,1),
                           Del=c(0,1)) |>
    mutate(id=row_number(),
           f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}"),
           exclude=glue("Dt|", 
                        "{ifelse(Avg==0, 'Avg|', 'NA|')}",
                        "{ifelse(Xf==0, 'Xfetch|', 'NA|')}",
                        "{ifelse(XN==0, 'lnNWt1X|', 'NA|')}",
                        "{ifelse(Del==0, 'Delta', 'NA')}"))
  filter(covSet.df, f==covSet)$exclude
}
