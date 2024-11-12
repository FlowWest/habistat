#' Create Flow-To-Suitable-Area-Curve for Individual Reach (unscaled)
#'
#' @param reach Integer `comid` identifier for reach
#' @param habitat_type Either "rearing" or "spawning"
#' @param units Desired units for output, either "ft" for square ft per linear ft, or "ac" for total acres
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and suitable habitat (`habitat`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_fsa_reach(reach = 7978069)
#'
habitat_fsa_reach <- function(reach, habitat_type = "rearing", units = "ft") {

  mode_unit <-
    if (units %in% c("ft", "feet", "ft2/lf")) "ft" else
      if(units %in% c("ac", "acres")) "ac"

  fsa <-
    switch(mode_unit,
           "ft" =
             habistat::wua_predicted |>
             filter(habitat == habitat_type) |>
             filter(comid == reach) |>
             select(flow_cfs,
                    habitat = wua_per_lf_pred),
           "ac" =
             habistat::wua_predicted |>
             filter(habitat == habitat_type) |>
             filter(comid == reach) |>
             transmute(flow_cfs,
                       habitat = wua_per_lf_pred * reach_length_ft / 43560))

  return(drop_na(fsa))

}

#' Create Flow-To-Suitable-Area Curve for Individual Reach (scaled by multiplier)
#'
#' @param reach Integer `comid` identifier for reach
#' @param multiplier Numeric multiplier to scale flow
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param units Desired units for output, either "ft" for square ft per linear ft, or "ac" for total acres
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and suitable habitat (`habitat`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_fsa_reach_scaled(reach = 7978069, multiplier = 0.8)
#'
habitat_fsa_reach_scaled <- function(reach, multiplier, habitat_type = "rearing", units = "ft") {

  fsa <-
    habistat::habitat_fsa_reach(reach,
                                habitat_type = habitat_type,
                                units = units)

  if (nrow(fsa) > 0) {

    return(tibble(flow_cfs = fsa$flow_cfs,
                  habitat = approx(x = fsa$flow_cfs * multiplier,
                                   y = fsa$habitat,
                                   xout = fsa$flow_cfs,
                                   rule = 2:2,
                                   method = "linear",
                                   na.rm = F)$y))

  } else {

    return(tibble(flow_cfs = numeric(0),
                  habitat = numeric(0)))

  }

}

#' Create Flow-To-Suitable-Area Curve for Reach, Mainstem, or Watershed
#'
#' @param reach Integer `comid` identifier for reach (must provide either reach, mainstem, OR watershed)
#' @param mainstem String name for mainstem (refer to `cv_mainstems` dataset for options)
#' @param watershed String name for watershed (refer to `cv_watersheds` dataset for options)
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param units Either "ft" for square ft per linear ft, or "ac" for total acres
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and suitable habitat (`habitat`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_fsa(reach = 7978069, units = "ft")
#'
#' habitat_fsa(mainstem = "Feather River", units = "ac")
#'
habitat_fsa <- function(reach, mainstem, watershed,
                        habitat_type = "rearing", units = "ft") {

  mode_geom <-
    if (!missing(reach)) "reach" else
      if (!missing(mainstem)) "mainstem" else
        if (!missing(watershed)) "watershed"

  mode_unit <-
    if (units %in% c("ft", "feet", "ft2/lf")) "ft" else
      if(units %in% c("ac", "acres")) "ac"

  if (mode_geom == "reach") {

    habitat_fsa_reach(reach = reach,
                      habitat_type = habitat_type,
                      units = units)

  } else if (mode_geom %in% c("mainstem", "watershed")) {

    group_var <-
      switch(mode_geom,
             "mainstem" = "river_cvpia",
             "watershed" = "watershed_level_3")

    flow_xw <-
      switch(mode_geom,
             "mainstem" = habistat::cv_mainstems_flow_xw,
             "watershed" = habistat::cv_watersheds_flow_xw)

    comids <-
      switch(mode_geom,
             "mainstem" =
               habistat::cv_mainstems |>
               filter(river_cvpia == mainstem) |>
               pull(comid),
             "watershed" =
               habistat::cv_watersheds |>
               filter(river_cvpia == mainstem) |>
               pull(comid))

    scaled_predictions <-
      flow_xw |>
      filter(comid %in% comids) |>
      mutate(fsa = map2(comid, multiplier, function(x, y) {
        habistat::habitat_fsa_reach_scaled(reach = x,
                                           multiplier = y,
                                           habitat_type = habitat_type,
                                           units = "ft")})) |> # output in units ft2/ft
      inner_join(habistat::flowline_attr |> select(comid, reach_length_ft),
                 by = join_by(comid)) |>
      unnest(fsa) |>
      mutate(habitat_ft2 = habitat * reach_length_ft)

    switch(mode_unit,
           "ft" =
             scaled_predictions |>
             group_by(flow_cfs) |>
             summarize(habitat = sum(habitat_ft2, na.rm = T) / sum(reach_length_ft),
                       .groups = "drop"),
           "ac" =
             scaled_predictions |>
             group_by(flow_cfs) |>
             summarize(habitat = sum(habitat_ft2, na.rm = T) / 43560,
                       .groups = "drop"))

  }

}

#' Create Duration Suitability Rating Curve for a streamgage
#'
#' @param streamgage Three-letter CDEC identifier for the streamgage
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param run One of "fall" (default), "late fall", "winter", "spring", or "steelhead"
#' @param wy_group Either "Dry" or "Wet" for water year type scenario
#' @param gradient For rearing, either "Valley Foothill" (default) or "Valley Lowland"
#' @param .n_days If TRUE, adds a third column to the output with the number of days inundated (`n_days`) value that was used to generate the duration suitability score. Defaults to FALSE.
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and duration habitat suitability index (`durhsi`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_drc("FSB", habitat_type="rearing", gradient = "Valley Foothill")
#'
habitat_drc <- function(streamgage,
                        habitat_type = "rearing",
                        run = "fall",
                        wy_group = "Dry",
                        gradient = "Valley Foothill",
                        .n_days = FALSE) {

  drc_filtered <-
    habistat::streamgage_duration_rating_curves |>
    filter((station_id == str_to_lower(streamgage)) &
             (run == str_to_lower(run)) &
             (habitat_type == str_to_lower(habitat_type)) &
             (wy_group == str_to_title(wy_group)))

  varname <- case_when(
    habitat_type == "spawning" ~ "durhsi_spawning",
    str_to_title(gradient) == "Valley Lowland" ~ "durhsi_rearing_vl",
    TRUE ~ "durhsi_rearing_vf")

  if (.n_days) {
    drc_selected <-
      drc_filtered$data[[1]] |>
      select(flow_cfs = model_q,
             durhsi = !!sym(varname),
             n_days = avg_max_days_inundated)
  } else {
    drc_selected <-
      drc_filtered$data[[1]] |>
      select(flow_cfs = model_q,
             durhsi = !!sym(varname))
  }

  return(drc_selected)

}

#' Create Duration Suitability Rating Curve for a streamgage (scaled to prediction point by drainage area and precipitation)
#'
#' @param streamgage Three-letter CDEC identifier for the streamgage
#' @param reach Integer `comid` identifier for reach identifying the prediction point (used for scaling)
#' @param scale Either "durhsi" to keep original flow scale and adjust duration HSI values (default), or "flow" to keep original duration HSI scale and adjust flow values.
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param run One of "fall" (default), "late fall", "winter", "spring", or "steelhead"
#' @param wy_group Either "Dry" or "Wet" for water year type scenario
#' @param gradient For rearing, either "Valley Foothill" (default) or "Valley Lowland"
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and duration habitat suitability index (`durhsi`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_drc_weighted("FSB", habitat_type="rearing", gradient = "Valley Foothill", reach = 7978069)
#'
habitat_drc_weighted <- function(streamgage,
                                 reach,
                                 scale = "durhsi",
                                 habitat_type = "rearing",
                                 run = "fall",
                                 wy_group = "Dry",
                                 gradient = "Valley Foothill") {

  gage_attr <-
    habistat::streamgage_da_attr |>
    filter(station_id == str_to_lower(streamgage)) |>
    as.list()

  comid_attr <-
    habistat::flowline_attr |>
    filter(comid == reach) |>
    select(comid, da_area_sq_km, da_ppt_mean_mm) |>
    as.list()

  multiplier <-
    (comid_attr$da_area_sq_km / gage_attr$da_gage) *
    (comid_attr$da_ppt_mean_mm / gage_attr$pc_gage)

  drc_unweighted <-
    habistat::habitat_drc(streamgage = streamgage,
                habitat_type = habitat_type,
                run = run,
                wy_group = wy_group,
                gradient = gradient)

  result <-
    switch(scale,
           "durhsi" =
             tibble(flow_cfs = drc_unweighted$flow_cfs,
                    durhsi = approx(x = drc_unweighted$flow_cfs * multiplier,
                                    xout = drc_unweighted$flow_cfs,
                                    y = drc_unweighted$durhsi,
                                    method = "linear",
                                    rule = 2:2,
                                    na.rm = FALSE)$y),
           "flow" =
             drc_unweighted |>
             mutate(flow_cfs = flow_cfs * multiplier))
  return(result)

}


#' Create Flow-To-Suitable-Area Curve (with Duration HSI applied) for Reach
#'
#' @param reach Integer `comid` identifier for reach
#' @param streamgage Three-letter CDEC identifier for the streamgage
#' @param habitat_type Either "rearing" or "spawning"
#' @param units Desired units for output, either "ft" for square ft per linear ft, or "ac" for total acres
#' @param run One of "fall" (default), "late fall", "winter", "spring", or "steelhead"
#' @param wy_group Either "Dry" or "Wet" for water year type scenario
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and suitable habitat (`habitat`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_fsa_duration_reach(reach = 7978069, streamgage = "FSB")
#'
habitat_fsa_duration_reach <- function(reach,
                                       streamgage,
                                       habitat_type = "rearing",
                                       units = "ft",
                                       run = "fall",
                                       wy_group = "Dry") {

  fsa <-
    habistat::habitat_fsa(reach = reach,
                          habitat_type = habitat_type,
                          units = units)

  grad <- if (habistat::flowline_attr$hqt_gradient_class
              [[which(habistat::flowline_attr$comid == reach)]]
              == "Valley Lowland") "Valley Lowland" else "Valley Foothill"

  drc <-
    habistat::habitat_drc_weighted(streamgage = streamgage,
                                   reach = reach,
                                   habitat_type = habitat_type,
                                   run = run,
                                   wy_group = wy_group,
                                   gradient = grad)

  fsa_weighted <-
    habistat::duration_apply_dhsi_to_fsa_curve(fsa = fsa,
                                               drc = drc,
                                               # var names
                                               fsa_q = flow_cfs,
                                               fsa_wua = habitat,
                                               drc_q = flow_cfs,
                                               drc_dhsi = durhsi)

  return(tibble(flow_cfs = fsa$flow_cfs,
                habitat = approx(x = fsa_weighted$q,
                                 xout = fsa$flow_cfs,
                                 y = fsa_weighted$durwua,
                                 rule = 2:2,
                                 na.rm = F)$y))

}


#' Create Flow-To-Suitable-Area Curve (with Duration HSI applied) for Reach, Mainstem, or Watershed
#'
#' @param reach Integer `comid` identifier for reach (must provide either reach, mainstem, OR watershed)
#' @param mainstem String name for mainstem (refer to `cv_mainstems` dataset for options)
#' @param watershed String name for watershed (refer to `cv_watersheds` dataset for options)
#' @param streamgage Three-letter CDEC identifier for the streamgage
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param units Desired units for output, either "ft" for square ft per linear ft, or "ac" for total acres
#' @param run One of "fall" (default), "late fall", "winter", "spring", or "steelhead"
#' @param wy_group Either "Dry" or "Wet" for water year type scenario
#'
#' @returns A `tibble` with columns for flow (`flow_cfs`) and suitable habitat (`habitat`)
#' @md
#' @export
#'
#' @examples
#'
#' habitat_fsa_duration(reach = 7978069, streamgage = "FSB")
#'
#' habitat_fsa_duration(mainstem = "Feather River", streamgage = "FSB")
#'
habitat_fsa_duration <- function(reach,
                                 mainstem,
                                 watershed,
                                 streamgage,
                                 habitat_type = "rearing",
                                 units = "ft",
                                 run = "fall",
                                 wy_group = "Dry") {

  mode_geom <-
    if (!missing(reach)) "reach" else
      if (!missing(mainstem)) "mainstem" else
        if (!missing(watershed)) "watershed"

  mode_unit <-
    if (units %in% c("ft", "feet", "ft2/lf")) "ft" else
      if(units %in% c("ac", "acres")) "ac"

  if (mode_geom == "reach") {

    habitat_fsa_duration_reach(reach = reach,
                               streamgage = streamgage,
                               habitat_type = habitat_type,
                               units = units,
                               run = run,
                               wy_group = wy_group)

  } else if (mode_geom %in% c("mainstem", "watershed")) {

    group_var <-
      switch(mode_geom,
             "mainstem" = "river_cvpia",
             "watershed" = "watershed_level_3")

    flow_xw <-
      switch(mode_geom,
             "mainstem" = habistat::cv_mainstems_flow_xw,
             "watershed" = habistat::cv_watersheds_flow_xw)

    comids <-
      switch(mode_geom,
             "mainstem" =
               habistat::cv_mainstems |>
               filter(river_cvpia == mainstem) |>
               pull(comid),
             "watershed" =
               habistat::cv_watersheds |>
               filter(river_cvpia == mainstem) |>
               pull(comid))

    scaled_predictions <-
      flow_xw |>
      filter(comid %in% comids) |>
      inner_join(habistat::flowline_attr |>
                   select(comid, hqt_gradient_class),
                 by = join_by(comid)) |>
      mutate(grad = if_else(hqt_gradient_class == "Valley Lowland",
                            "Valley Lowland",
                            "Valley Foothill"),
             .keep = "unused") |>
      mutate(fsa = map2(comid, multiplier, function(x, y) {
        habistat::habitat_fsa_reach_scaled(reach = x,
                                           multiplier = y,
                                           habitat_type = habitat_type,
                                           units = "ft")})) |> # output in units ft2/ft
      mutate(drc = pmap(list(comid, streamgage, grad), function(x, y, z) {
        habistat::habitat_drc_weighted(reach = x,
                                       streamgage = y,
                                       run = run,
                                       wy_group = wy_group,
                                       gradient = z,
                                       scale = "flow")})) |>
      mutate(fsa_x_drc = map2(fsa, drc, function(f, d) {
        fsa_weighted <-
          habistat::duration_apply_dhsi_to_fsa_curve(fsa = f,
                                                     drc = d,
                                                     # var names
                                                     fsa_q = flow_cfs,
                                                     fsa_wua = habitat,
                                                     drc_q = flow_cfs,
                                                     drc_dhsi = durhsi)
        if(nrow(fsa_weighted) > 0) {
          return(tibble(flow_cfs = f$flow_cfs,
                        habitat = approx(x = fsa_weighted$q,
                                         xout = f$flow_cfs,
                                         y = fsa_weighted$durwua,
                                         rule = 2:2,
                                         na.rm = F)$y))
        } else {
          return(tibble(flow_cfs = numeric(0),
                        habitat = numeric(0)))
        }
      })) |>
      inner_join(habistat::flowline_attr |>
                   select(comid, reach_length_ft),
                 by = join_by(comid)) |>
      select(-fsa, -drc) |>
      unnest(fsa_x_drc) |>
      mutate(habitat_ft2 = habitat * reach_length_ft)

    switch(mode_unit,
           "ft" =
             scaled_predictions |>
             group_by(flow_cfs) |>
             summarize(habitat = sum(habitat_ft2, na.rm = T) / sum(reach_length_ft),
                       .groups = "drop"),
           "ac" =
             scaled_predictions |>
             group_by(flow_cfs) |>
             summarize(habitat = sum(habitat_ft2, na.rm = T) / 43560,
                       .groups = "drop"))

  }

}

#' Predict Habitat Area using a Flow-to-Suitable-Area Curve
#'
#' @param x A vector of flows at which to predict suitable habitat values
#' @param reach Integer `comid` identifier for reach (must provide either reach, mainstem, OR watershed)
#' @param mainstem String name for mainstem (refer to `cv_mainstems` dataset for options)
#' @param watershed String name for watershed (refer to `cv_watersheds` dataset for options)
#' @param streamgage (optional) Three-letter CDEC identifier for a streamgage
#' @param habitat_type Either "rearing" (default) or "spawning"
#' @param units Desired units for output, either "ft" for square ft per linear ft, or "ac" for total acres
#' @param run (required if streamgage is provided) One of "fall" (default), "late fall", "winter", "spring", or "steelhead"
#' @param wy_group (required if streamgage is provided) Either "Dry" or "Wet" for water year type scenario
#' @param fsa_x (optional) Vector of flows for manually defined flow-to-suitable-area curve. If provided, overrides reach, mainstem, or watershed lookup.
#' @param fsa_y (optional) Vector of habitat values for manually defined flow-to-suitable-area curve, of same length as `fsa_x`
#'
#' @returns A vector of the same length as `x` containing predicted suitable habitat values
#' @md
#' @export
#'
#' @examples
#'
#' # example non-temporal prediction given comid
#' habitat_predict(c(0, 10, 150, 300, 9999), reach = 7978069)
#'
#' # example temporal prediction given comid and streamgage
#' habitat_predict(c(0, 10, 150, 300, 9999), reach = 7978069, streamgage = "FSB")
#'
#' # example non-temporal prediction given mainstem
#' habitat_predict(c(0, 10, 150, 300, 9999), mainstem = "Feather River")
#'
#' # example temporal prediction given comid and streamgage
#' habitat_predict(c(0, 10, 150, 300, 9999), mainstem = "Feather River", streamgage = "FSB")
#'
#' # example manual prediction given custom provided flow-to-suitable-area curve
#' habitat_predict(c(0, 10, 150, 300, 9999), fsa_x = c(0, 100, 1000), fsa_y = c(47, 470, 4700))
#'
habitat_predict <- function(x,
                            reach, mainstem, watershed,
                            streamgage,
                            fsa_x,
                            fsa_y,
                            habitat_type = "rearing",
                            units = "ft",
                            run = "fall",
                            wy_group = "Dry") {

  mode_geom <-
    if (!missing(reach)) "reach" else
      if (!missing(mainstem)) "mainstem" else
        if (!missing(watershed)) "watershed"

  mode_unit <-
    if (units %in% c("ft", "feet", "ft2/lf")) "ft" else
      if(units %in% c("ac", "acres")) "ac"

  if (isTRUE((!missing(fsa_x) && (!missing(fsa_y))))) {

    fsa <- tibble(flow_cfs = fsa_x,
                  habitat = fsa_y)

  } else if (!missing(streamgage)) {

    fsa <- habistat::habitat_fsa_duration(reach = reach,
                                          mainstem = mainstem,
                                          watershed = watershed,
                                          streamgage = streamgage,
                                          habitat_type = habitat_type,
                                          units = units,
                                          run = run,
                                          wy_group = wy_group)

  } else {

    fsa <- habistat::habitat_fsa(reach = reach,
                                 mainstem = mainstem,
                                 watershed = watershed,
                                 habitat_type = habitat_type,
                                 units = units)

  }

  result <- approx(x = if_else(fsa$flow_cfs > 0, log10(fsa$flow_cfs), 0),
                   y = fsa$habitat,
                   xout = if_else(x > 0, log10(x), 0),
                   method = "linear",
                   rule = 2:2)

  return(result$y)

}
