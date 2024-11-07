library(tidyverse)
library(sf)
library(shiny)
library(shinyjs)
library(leaflet)
library(leafgl)
library(habistat)
library(patchwork)

# reactlog::reactlog_enable()

get_data <- function(...) {
  # use this for loading data from R package such that it can be easily garbage-collected
  # get_data(wua_predicted, package = "habistat")
  # gc()
  e <- new.env()
  name <- data(..., envir = e)[1]
  e[[name]]
}

# DATA IMPORT ------------------------------------------------------------------

predictions <- get_data(wua_predicted, package = "habistat") |>
  #filter(habitat == "rearing") |>
  #mutate(
                           #if_else((habitat=="rearing" & model_bfc),
                           #paste0(model_name, "_ph_bfc_rm"), # post-hoc baseflow channel removal
                           #model_name)) |>
  select(comid, flow_idx, flow_cfs, habitat, model_id = model_name, wua_per_lf_pred, river_cvpia, watershed_level_3, reach_length_ft) |>
  pivot_wider(names_from = model_id,
              values_from = wua_per_lf_pred,
              names_glue = c("wua_per_lf_pred_{model_id}")) |>
  left_join(get_data(wua_hydraulic_interp, package = "habistat") |>
              # filter((habitat=="rearing" & !bfc) | (habitat=="spawning")) |> # only showing the actuals with prior bfc removed for now
              select(habitat, comid, flow_cfs, wua_per_lf_actual = wua_per_lf),
            by = join_by(habitat, comid, flow_cfs)) |>
  left_join(get_data(flowline_attr, package = "habistat") |>
              transmute(comid,
                        chan_width_ft = chan_width_m/0.3048,
                        baseflow_cfs = nf_bfl_dry_cfs),
            by = join_by(comid)) |>
  mutate(wua_acres_pred_SD = wua_per_lf_pred_SD * reach_length_ft / 43560,
         wua_acres_pred_SN = wua_per_lf_pred_SN * reach_length_ft / 43560,
         wua_acres_actual = wua_per_lf_actual * reach_length_ft / 43560)

gc() # garbage collect after loading from habistat package

all_flows <- unique(predictions$flow_cfs)
all_flows_idx <- unique(predictions$flow_idx)

var_names <- c(
          comid = "NHDPlusV2 ComID",
          gnis_name = "GNIS Name",
          river_cvpia = "Habitat Mainstem Name",
          watershed_level_3 = "Watershed Name",
          hqt_gradient_class = "HQT Gradient Class",
          hyd_cls = "UCD Hydrologic Class",
          chan_width_ft = "Channel Width (ft)",
          reach_length_ft = "Reach Length (linear ft)",
          slope = "Reach Slope (ft/ft)",
          sinuosity = "Reach Sinuosity",
          divergence_ratio = "Reach Divergence Fraction",
          loc_bfi = "Reach Baseflow Index",
          loc_pct_clay = "Reach Soil Pct. Clay",
          loc_pct_sand = "Reach Soil Pct. Sand",
          loc_permeability = "Reach Soil Permeability",
          loc_bedrock_depth = "Reach Depth to Bedrock",
          loc_ppt_mean_mm = "Reach Mean Annual Precip. (mm)",
          mtpi30_min = "Reach Topographic Position Index (mTPI)",
          frac_leveed_longitudinal = "Reach Levee Confinement Fraction",
          da_area_sq_km = "DA Area (sq km)",
          da_avg_slope = "DA Avg. Slope (ft/ft)",
          da_k_erodibility = "DA Erodibility",
          mean_ndvi = "DA Mean NDVI",
          da_elev_mean = "DA Mean Elevation",
          da_ppt_mean_mm = "DA Mean Annual Precip. (mm)",
          bf_depth_m = "Bankfull Depth (m)",
          bf_w_d_ratio = "Bankfull Width:Depth Ratio",
          vb_width_transect = "Valley Bottom Width",
          vb_bf_w_ratio = "Valley Bottom:Channel Width Ratio",
          flow_idx = "Flow Identifier",
          flow_cfs = "Flow (cfs)",
          wua_per_lf_pred_SD = "Predicted Suitable Habitat (ft2/ft) (SD Model)",
          wua_per_lf_pred_SN = "Predicted Suitable Habitat (ft2/ft) (SN Model)",
          wua_per_lf_actual = "Actual Suitable Habitat (ft2/ft) ",
          wua_acres_pred_SD = "Predicted Suitable Habitat (acres) (SD Model)",
          wua_acres_pred_SN = "Predicted Suitable Habitat (acres) (SN Model)",
          wua_acres_actual = "Actual Suitable Habitat (acres)",
          baseflow_cfs = "Estimated Baseflow (cfs)",
          spawning_geographic_context = "Spawning: Within Geographic Range",
          spawning_gravel_either_method = "Spawning: Likely Geomorphic Suitability")

attr <- get_data(flowline_attr, package = "habistat") |>
  left_join(get_data(flowline_spawning_attr, package = "habistat"), by=join_by(comid)) |>
  filter(comid %in% predictions$comid) |>
  select(any_of(names(var_names)))
  # stream size filter -- make sure this matches model_cleaned.Rmd
  # filter(((stream_order >= 4) & (da_area_sq_km > 1)) | ((stream_order >= 3) & (da_area_sq_km >= 50))) |>
  # filter(ftype == "StreamRiver")

gc() # garbage collect after loading from habistat package

geom <- get_data(flowline_geom, package = "habistat") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  select(comid, geometry) |>
  inner_join(attr |> transmute(comid, gnis_name), by = join_by(comid)) |>
  mutate(object_id = paste0("comid_", comid))

gc() # garbage collect after loading from habistat package

watersheds <- get_data(cv_watersheds, package = "habistat") |>
  group_by(watershed_level_1, watershed_level_2, watershed_level_3) |> #, range_pisces) |>
  summarize(.groups = "drop") |>
  st_union(by_feature=T) |>
  st_transform("+proj=longlat +datum=NAD83") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  mutate(watershed_id = paste0("watershed_", row_number()),
         watershed_label = paste0(watershed_level_1, " Basin",
                                  if_else(watershed_level_2 != watershed_level_1,
                                          paste0("<br />", watershed_level_2, " Watershed"),
                                          " (Direct Drainage)"),
                                  if_else(watershed_level_3 != watershed_level_2,
                                          paste0("<br />", watershed_level_3, " Subwatershed"),
                                          "")))

cv_watersheds_flow_xw <- get_data(cv_watersheds_flow_xw, package = "habistat")

mainstems <- get_data(cv_mainstems, package = "habistat") |>
  group_by(river_group, river_cvpia) |> #, range_pisces) |>
  summarize(.groups = "drop") |>
  st_union(by_feature=T) |>
  st_transform("+proj=longlat +datum=NAD83") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  mutate(mainstem_id = paste0("mainstem_", row_number()),
         mainstem_label = paste0(river_cvpia))

cv_mainstems_flow_xw <- get_data(cv_mainstems_flow_xw, package = "habistat")

gc()

predictions_watershed <- get_data(wua_predicted_cv_watersheds, package = "habistat") |>
  ungroup() |>
  select(watershed_level_3, flow_idx, flow_cfs, habitat, model_id = model_name, wua_per_lf_pred, wua_acres_pred) |>
  pivot_wider(names_from = model_id,
              values_from = c(wua_per_lf_pred, wua_acres_pred),
              names_glue = c("{.value}_{model_id}"))

predictions_mainstem <- get_data(wua_predicted_cv_mainstems, package = "habistat") |>
  ungroup() |>
#  select(watershed_level_3, flow_idx, flow_cfs, habitat, model_id = model_name, wua_per_lf_pred, wua_acres_pred) |>
  select(river_cvpia, flow_idx, flow_cfs, habitat, model_id = model_name, wua_per_lf_pred, wua_acres_pred) |>
  pivot_wider(names_from = model_id,
              values_from = c(wua_per_lf_pred, wua_acres_pred),
              names_glue = c("{.value}_{model_id}"))

gc() # garbage collect after loading from habistat package

# STREAMGAGES ------------------------------------------------------------------

streamgage_attr <- get_data(streamgage_attr, package = "habistat") |>
  inner_join(get_data(streamgage_da_attr, package = "habistat"), by = join_by(station_id)) |>
  transmute(river_group, river_cvpia, station_id,
            gage_comid, da_gage, pc_gage,
            station_label = glue::glue("{str_to_upper(station_id)}: {str_to_upper(name)} ({min_wy}-{max_wy})"))

streamgage_pts <- get_data(streamgage_geom, package = "habistat") |>
  st_transform("+proj=longlat +datum=NAD83") |>
  st_set_crs("+proj=longlat +datum=WGS84") |>
  inner_join(streamgage_attr, by=join_by(station_id))

streamgages <- streamgage_attr |>
  select(river_cvpia, station_id, station_label) |>
  nest(data = c(station_id, station_label)) |>
  mutate(data = map(data, function(x) deframe(x) |> as.list())) |>
  deframe()

streamgages_by_watershed <- streamgage_attr |>
  inner_join(attr |> select(comid, watershed_level_3), by=join_by(gage_comid == comid)) |>
  select(watershed_level_3, station_id, station_label) |>
  nest(data = c(station_id, station_label)) |>
  mutate(data = map(data, function(x) deframe(x) |> as.list())) |>
  deframe()

hqt_boundary <- get_data(hqt, package = "habistat")

# PLOT STYLES AND PALETTES -----------------------------------------------------

# pal <- function(x, n = 10,
#                 palette_function = viridisLite::cividis,
#                 palette_args = list(direction = -1)) {
#   #palette <- palette_function(n)
#   palette <- do.call(palette_function, c(list(n=n), palette_args))
#   palette[cut(x, n)]
# }

theme_set(theme_minimal())

ihs <- scales::trans_new("ihs",
                         transform = function(x) asinh(x),
                         inverse = function(y) sinh(y),
                         breaks = function(i) scales::breaks_log(n=5, base=10)(pmax(i,0.01)),
                         #minor_breaks = scales::minor_breaks_n(n = 0),
                         domain=c(0, Inf),
                         format = scales::label_comma())

# flow_scale_colors <- c("darkblue", "turquoise", "gold", "darkorange", "darkred", "violetred4", "mediumvioletred"),
flow_scale_colors <- list(rearing = c("#00008B", "#40E0D0", "#FFD700", "#FF8C00", "#8B0000", "#8B2252", "#C71585"),
                          spawning = c("#00008B", "#40E0D0", "#FFD700", "#FF8C00", "#8B0000", "#8B2252", "#C71585"))
flow_scale_breaks <- list(rearing = c(0, 1, 3, 10, 30, 100, 300),
                          spawning = c(0, 20, 40, 60, 80, 100, 120))

pal <- function(x, type = "rearing") {
  if(type == "rearing") {
    breaks_scaled <- scales::rescale(habistat::semiIHS(flow_scale_breaks$rearing))
    values_scaled <- scales::rescale(habistat::semiIHS(x), from = range(habistat::semiIHS(flow_scale_breaks$rearing)))
  } else if(type == "spawning") {
    breaks_scaled <- scales::rescale(flow_scale_breaks$spawning)
    values_scaled <- scales::rescale(x, from = range(flow_scale_breaks$spawning))
  }
  cut(x = values_scaled,
      breaks = c(-Inf, breaks_scaled[2:length(breaks_scaled)], Inf),
      labels = flow_scale_colors[[type]]) |> as.character()
}


add_color_scale <- function(g, type="rearing", ...) {
  g + scale_color_gradientn(name = "WUA per LF",
                            limits = c(flow_scale_breaks[[type]][[1]],
                                       flow_scale_breaks[[type]][[length(flow_scale_breaks[[type]])]]),
                            breaks = flow_scale_breaks[[type]],
                            trans = if (type=="rearing") habistat::trans_semiIHS else NULL,
                            values = if (type=="rearing") scales::rescale(habistat::semiIHS(flow_scale_breaks[[type]])) else NULL,
                            colors = flow_scale_colors[[type]],
                            na.value = "#ffffff00",
                            ...)
}


# before deploying shiny app, install package from github:
# options(timeout=600)
# devtools::install_github("FlowWest/swc-habitat-suitability@26-shiny-app-interface-fixes")

scale_fsa <- function(data, multiplier, .wua_var = selected_wua, .flow_var = flow_cfs) {
  data |>
    expand_grid(scaled = c(FALSE, TRUE)) |>
    mutate({{.flow_var}} := if_else(scaled, {{.flow_var}} * multiplier, {{.flow_var}})) |>
    arrange({{.flow_var}}) |>
    mutate({{.wua_var}} := if_else(scaled, NA, {{.wua_var}})) |>
    mutate({{.wua_var}} := zoo::na.approx({{.wua_var}}, x = {{.flow_var}}, na.rm=F, rule=2)) |>
    filter(scaled) |>
    select(-scaled) |>
    ungroup()
}
