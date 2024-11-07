#' Import tabular SRH-2D results
#'
#' @param csv_tbl A data frame or tibble containing the file path to each CSV containing SRH-2D model results in the standard SRH-2D export format (columns for vertex id [optional], X coordinate, Y coordinate, Z coordinate, WSE, depth, velocity X component, velocity Y component, velocity magnitude, froude number, and shear stress). There should be a numeric column called `flow_cfs` and a character column called `filename` containing valid file paths to the CSV corresponding to each flow.
#' @param units A string indicating units of `ft` or `m`
#' @param units_xy A string indicating horizontal coordinate units of `ft` or `m`. Defaults to `units` if not separately specified.
#' @param units_z A string indicating vertical coordinate units of `ft` or `m`. Defaults to `units` if not separately specified.
#' @param has_vid_col Boolean indicating whether or not there is already a vertex id (`vid`) column in the CSV. If `FALSE`, one will be generated based on the row number.
#'
#' @return A `tbl_df` data frame containing all SRH-2D model results at all flows with standard column names, one row per vertex per flow.
#' @md
#' @examples
#'
#' @export
vector_import_srh2d <- function(csv_tbl,
                                units = "ft",
                                units_xy = units,
                                units_z = units,
                                has_vid_col = FALSE) {

  conv_fact_xy <- case_when(units_xy == "ft" ~ 1,
                            units_xy == "m" ~ 3.28084)

  conv_fact_z <- case_when(units_z == "ft" ~ 1,
                           units_z == "m" ~ 3.28084)


  cols <- c("x"="n", "y"="n", "z"="n",
            "wse"="n", "wdepth"="n",
            "vel_x"="n", "vel_y"="n", "vel_mag"="n",
            "froude"="n", "stress"="n")

  if (has_vid_col) {
    cols <- c("vid"="n", cols)
  }

  csv_tbl |>
    mutate(tbl = map(filename, function(x) {
      message("import SRH-2D results from ", x)
      csv_in <- read_csv(x, col_names=names(cols), col_types=cols, skip=1)
      if (has_vid_col) {
        return(csv_in)
      } else {
        return(csv_in |> mutate(vid = row_number()))
      }
    })) |>
    unnest(tbl) |>
    mutate(across(everything(), function(n) if_else(n==-999, NA, n))) |>
    mutate(across(c("x", "y", "vel_x", "vel_y", "vel_mag"), function(n) n * conv_fact_xy)) |>
    mutate(across(c("z", "wse", "wdepth"), function(n) n * conv_fact_z))
}

#' Vector Habitat Suitability Index (HSI) function: HQT
#'
#' @param d a vector of depths in feet
#' @param v a vector of velocities in feet per second
#'
#' @return A vector of habitat suitability index values between 0 and 1.
#' @md
#' @examples
#'
#' @export
vector_dvhsi_hqt <- function(d, v) {
  if_else(d>1.0 & d<=3.28 & v>0 & v<=1.5, 1, 0)
}

#' Vector Habitat Suitability Index (HSI) function: LYR
#'
#' @param d a vector of depths in feet
#' @param v a vector of velocities in feet per second
#'
#' @return A vector of habitat suitability index values between 0 and 1.
#' @md
#' @examples
#'
#' @export
vector_dvhsi_lyr <- function(d, v) {
  if_else(d>0.5 & d<=5.2 & v>0 & v<=4.0, 1, 0)
}

#' Vector Habitat Suitability Index (HSI) function: Fall-Run Chinook Spawning
#'
#' @param d a vector of depths in feet
#' @param v a vector of velocities in feet per second
#' @param run identifier for the salmon run of interest (one of `fall`, `late_fall`, `winter`, `spring`) or `max` to use the outer envelope of all runs. Defaults to `max`.
#' @param strict specify whether both criteria must strictly be met (minimum operation) or whether the factors are combined using a geometric mean
#' @return A vector of habitat suitability index values between 0 and 1.
#' @md
#' @examples
#'
#' @export
vector_dvhsi_spawning <- function(d, v, run="max", strict=T) {
  # if_else(d>0 & d<=6.56 & v>0 & v<=5.9, 1, 0)

  fct_dep <-
    case_when(run == "winter" ~ if_else(d>1 & d<=6.56, 1, 0),
              TRUE ~            if_else(d>0 & d<=6.56, 1, 0))

  fct_vel <-
    approx(x = habistat::suitability_hsi_vel_spawning[[run]]$vel_ft_s,
           y = habistat::suitability_hsi_vel_spawning[[run]]$hsi,
           xout = v, yleft = 0, yright = 0)$y
  if (strict) {
    # must strictly meet both criteria
    return(pmin(fct_dep, fct_vel))
  } else {
    return(sqrt(fct_dep * fct_vel)) # geometric mean
  }

}

#' Calculate HSI for each row of SRH-2D results
#'
#' @param data A `tbl_df` of the format returned by `vector_import_srh2d()`
#' @param hsi_func A function that takes a depth (ft) vector and velocity (ft/s) vector and returns a vector of habitat suitability index values between 0 and 1. Function arguments should be `d` and `v` in US units. Defaults to the `vector_dvhsi_hqt` function in this package.
#'
#' @return A `tbl_df` data frame containing the calculated habitat suitability indices by flow (`flow_cfs`) and vertex (`vid`) along with the associated depths and velocities.
#' @md
#' @examples
#'
#' @export
vector_calculate_hsi <- function(data, hsi_func = vector_dvhsi_hqt) {

  data |>
    select(flow_cfs, vid, wdepth, vel_mag) |>
    mutate(inundated = wdepth > 0,
           hsi = hsi_func(wdepth, vel_mag))

}

#' Prep SRH-2D mesh, assigning groups and calculating areas
#'
#' @param mesh An `sf` polygon feature collection with one polygon for each vertex (computation point) identified by a `vid` attribute. These may be Voronoi (aka Thiessen) polygons generated based on vertex points in a tool such as sf (`st_voronoi`), QGIS (Voronoi Polygons), or ArcGIS (Create Thiessen Polygons).
#' @param group_polygons An `sf` polygon feature collection with one polygon for each analysis reach. These may have been delineated previously by splitting the hydraulic model domain into subsections.
#' @param mask_polygons Optional: An `sf` polygon feature collection defining areas to be masked out of the HSI calculation. These areas will still be included in the total inundated area.
#' @param .group_var The unquoted name of the attribute used to identify groups in the `group_polygons` layer. Defaults to `comid`
#' @param .id_var The unquoted name of the attribute used to identify vertices in the `mesh` layer. Defaults to `vid`
#'
#' @return A `tbl_df` with one row per vertex (computation point), including the vertex id (e.g. `vid`), the associated group id (e.g. `comid`), and the area of the polygon in square feet.
#' @md
#' @examples
#'
#' @export
vector_prep_mesh <- function(mesh,
                             group_polygons = NULL,
                             mask_polygons = NULL,
                             .group_var = comid,
                             .id_var = vid) {

  groups <- group_polygons |>
    group_by({{.group_var}}) |>
    summarize() |>
    st_union(by_feature = TRUE)

  result <- mesh |>
    mutate(area_ft2 = st_area(geometry) |> units::set_units("ft2") |> units::drop_units()) |>
    st_centroid() |>
    st_intersection(groups) |>
    transmute({{.id_var}},
              {{.group_var}},
              area_ft2)

  if(!is.null(mask_polygons)) {

    mask_polygons <-
      mask_polygons |>
      st_sf() |>
      transmute(mask_val = 0)

    result <- result |>
      st_join(mask_polygons, largest = T, left = T) |>
      mutate(mask_val = coalesce(mask_val, 1)) |>
      st_drop_geometry()

  } else {

    result <- result |>
      mutate(mask_val = 1) |>
      st_drop_geometry()

  }

  return(result)

}

#' Summarize habitat suitability for SRH-2D results
#'
#' @param mesh_tbl A `tbl_df` with one row per vertex, containing group IDs and mesh cell areas, as created with `vector_prep_mesh()`
#' @param hsi_tbl A `tbl_df` with one row per vertex per flow, containing depths, velocities, and habitat suitability indices, as created with `vector_calculate_hsi()`
#' @param .group_var The unquoted name of the attribute used to identify groups in the `group_polygons` layer. Defaults to `comid`
#' @param .id_var The unquoted name of the attribute used to identify vertices in the `mesh` layer. Defaults to `vid`
#'
#' @return A `tbl_df` data frame with one row per group (`comid`) per flow (`flow_cfs`), containing the total inundated area (`area_tot`) in square feet, the suitable habitat area (`area_wua`) in square feet, and the ratio of suitable to total area (`area_pct`).
#' @md
#' @examples
#'
#' @export
vector_summarize_hsi <- function(mesh_tbl, hsi_tbl, .group_var = comid, .id_var=vid) {
  mesh_tbl |>
    inner_join(hsi_tbl, by=join_by({{.id_var}}), relationship="one-to-many") |>
    group_by({{.group_var}}, flow_cfs) |>
    summarize(area_tot = sum(if_else(wdepth > 0, area_ft2, 0)),
              area_wua = sum(area_ft2 * hsi * mask_val)) |>
    mutate(area_pct = area_wua / area_tot)
}
