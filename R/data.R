#' @name flowline_geom
#' @title NHD Flowline (ComID) Geometries
#' @examples
#' flowline_geom |> pillar::glimpse()
"flowline_geom"

#' @name flowline_geom_proj
#' @title NHD Flowline (ComID) Geometries - Projected
#' @examples
#' flowline_geom_proj |> pillar::glimpse()
"flowline_geom_proj"

#' @name flowline_attr
#' @title NHD Flowline (ComID) Attributes
#' @examples
#' flowline_attr |> pillar::glimpse()
"flowline_attr"

#' @name flowline_spawning_attr
#' @title NHD Flowline (ComID) Spawning Suitability Flags
#' @examples
#' flowline_spawning_attr |> pillar::glimpse()
"flowline_spawning_attr"

#' @name wua_hydraulic
#' @title WUA/LF by ComID by Flow (training dataset)
#' @examples
#' wua_hydraulic |> pillar::glimpse()
"wua_hydraulic"

#' @name wua_hydraulic_interp
#' @title WUA/LF by ComID by Flow (training dataset, interpolated)
#' @examples
#' wua_hydraulic_interp |> pillar::glimpse()
"wua_hydraulic"

#' @name wua_predicted
#' @title WUA/LF by ComID by Flow (predicted by statistical model)
#' @examples
#' wua_predicted |> pillar::glimpse()
"wua_predicted"

#' @name cv_mainstems
#' @title CVPIA Mainstems (geometry)
#' @examples
#' cv_mainstems |> ggplot2::ggplot() + ggplot2::geom_sf()
"cv_mainstems"

#' @name cv_watersheds
#' @title CVPIA Watersheds (geometry)
#' @examples
#' cv_watersheds |> ggplot2::ggplot() + ggplot2::geom_sf()
"cv_watersheds"

#' @name wua_predicted_cv_mainstems
#' @title WUA/LF by CV Mainstem by Flow (predicted by statistical model)
#' @examples
#' wua_predicted_cv_mainstems |> pillar::glimpse()
"wua_predicted_cv_mainstems"

#' @name wua_predicted_cv_watersheds
#' @title WUA/LF by CV Watershed by Flow (predicted by statistical model)
#' @examples
#' wua_predicted_cv_watersheds |> pillar::glimpse()
"wua_predicted_cv_watersheds"

#' @name cv_mainstems_flow_xw
#' @title Flow Area Crosswalk for CV Mainstems
"cv_mainstems_flow_xw"

#' @name cv_watersheds_flow_xw
#' @title Flow Area Crosswalk for CV Watersheds
"cv_watersheds_flow_xw"

#' @name streamgage_attr
#' @title CDEC Mainstem Streamgage Attributes
"streamgage_attr"

#' @name streamgage_geom
#' @title CDEC Mainstem Streamgage Point Geometries
"streamgage_geom"

#' @name streamgage_duration_rating_curves
#' @title CDEC Mainstem Streamgage Duration Rating Curves
"streamgage_duration_rating_curves"

#' @name streamgage_da_attr
#' @title CDEC Mainstem Streamgage Drainage Area Attributes
"streamgage_da_attr"

#' @name suitability_hsi_vel_spawning
#' @title Velocity suitability criteria for spawning by run
"suitability_hsi_vel_spawning"

#' @name hqt
#' @title HQT Valley Lowland Polygon
"hqt"
