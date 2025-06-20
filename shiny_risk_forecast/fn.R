

update_map_range <- function(usr_in, map_ranges, type="brush", bbox=NULL) {
  if(!is.null(usr_in)) {
    if(type=="brush") {
      map_ranges$x <- c(round(usr_in$xmin), round(usr_in$xmax))
      map_ranges$y <- c(round(usr_in$ymin), round(usr_in$ymax))
    }
    if(type=="dblClick") {
      map_ranges$x <- bbox[c("xmin", "xmax")]
      map_ranges$y <- bbox[c("ymin", "ymax")]
    }
  }
}



update_sin_click <- function(usr_in, sin_reactive, site_sf, sin_avail) {
  if(!is.null(usr_in)) {
    res <- nearPoints(site_sf |>
                        st_drop_geometry() |> 
                        filter(sin %in% sin_avail),
                      usr_in, "lon", "lat", maxpoints=1)
    if(nrow(res) > 0) {
      res$sin[1]
    }
  }
}