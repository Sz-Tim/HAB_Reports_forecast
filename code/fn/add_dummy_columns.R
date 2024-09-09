#' Add Dt dummy variables
#' 
#' Should only be necessary for initial model fitting since there were some
#' where Dt variables were included but then excluded within the recipe. The
#' bake() function still expects them initially even though they are excluded.
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
add_dummy_columns <- function(df) {
  dummy_cols <- c("chlAvgDtDirE", "chlAvgDtDirN", "chlAvgDtDirS", "chlAvgDtDirW", "chlDt", 
                  "kdAvgDtDirE", "kdAvgDtDirN", "kdAvgDtDirS", "kdAvgDtDirW", "kdDt", 
                  "no3AvgDtDirE", "no3AvgDtDirN", "no3AvgDtDirS", "no3AvgDtDirW", "no3Dt", 
                  "o2AvgDtDirE", "o2AvgDtDirN", "o2AvgDtDirS", "o2AvgDtDirW", "o2Dt",
                  "phAvgDtDirE", "phAvgDtDirN", "phAvgDtDirS", "phAvgDtDirW", "phDt", 
                  "phycAvgDtDirE", "phycAvgDtDirN", "phycAvgDtDirS", "phycAvgDtDirW", "phycDt", 
                  "po4AvgDtDirE", "po4AvgDtDirN", "po4AvgDtDirS", "po4AvgDtDirW", "po4Dt", 
                  "ppAvgDtDirE", "ppAvgDtDirN", "ppAvgDtDirS", "ppAvgDtDirW", "ppDt", 
                  "PrecipAvgDtDirE", "PrecipAvgDtDirN", "PrecipAvgDtDirS", "PrecipAvgDtDirW", "PrecipDt",
                  "ShortwaveAvgDtDirE", "ShortwaveAvgDtDirN", "ShortwaveAvgDtDirS", "ShortwaveAvgDtDirW", "ShortwaveDt", 
                  "sstAvgDtDirE", "sstAvgDtDirN", "sstAvgDtDirS", "sstAvgDtDirW", "sstDt",
                  "UAvgDtDirE", "UAvgDtDirN", "UAvgDtDirS", "UAvgDtDirW", "UDt", 
                  "UVAvgDtDirE", "UVAvgDtDirN", "UVAvgDtDirS", "UVAvgDtDirW", "UVDt", 
                  "VAvgDtDirE", "VAvgDtDirN", "VAvgDtDirS", "VAvgDtDirW", "VDt")
  df |>
    bind_cols(map_dfr(dummy_cols, ~tibble(x=.x, y=0)) |> pivot_wider(names_from=x, values_from=y))
}
