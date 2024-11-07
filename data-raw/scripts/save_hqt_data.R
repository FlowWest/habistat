hqt <- readRDS("data-raw/source/hqt/hqt_valley_lowland.Rds") |>
  sf::st_transform("+proj=longlat +datum=WGS84 +no_defs") |>
  st_zm() |>
  st_make_valid() |>
  st_union() |>
  st_cast("POLYGON") |>  #
  st_boundary()  #

usethis::use_data(hqt, overwrite = TRUE)
