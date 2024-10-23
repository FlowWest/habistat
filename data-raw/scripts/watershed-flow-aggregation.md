Predictor Data Preparation and Consolidation
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-10-23

- [Basic Flow Crosswalk via Drainage Area and
  Precipitation](#basic-flow-crosswalk-via-drainage-area-and-precipitation)
  - [Applying the crosswalk to aggregate model
    results](#applying-the-crosswalk-to-aggregate-model-results)
- [DSMHabitat Comparison](#dsmhabitat-comparison)
- [Applying this method to the duration
  analyis](#applying-this-method-to-the-duration-analyis)
  - [Streamgage drainage area
    attributes](#streamgage-drainage-area-attributes)
  - [Example application to COMID-level flow-to-suitable-area
    curve](#example-application-to-comid-level-flow-to-suitable-area-curve)
  - [Example application to watershed-level flow-to-suitable-area
    curve](#example-application-to-watershed-level-flow-to-suitable-area-curve)

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(sf)
```

    ## Linking to GEOS 3.11.2, GDAL 3.7.2, PROJ 9.3.0; sf_use_s2() is TRUE

``` r
library(habistat)
library(patchwork)
theme_set(theme_minimal())
```

## Basic Flow Crosswalk via Drainage Area and Precipitation

See watershed_delineation.R for dependencies

``` r
make_flow_xw <- function(group_var) {
  habistat::flowline_attr |>
  select({{group_var}}, 
         comid = comid,
         da_reach = da_area_sq_km_tot,
         pc_reach = da_ppt_mean_mm) |>
  drop_na({{group_var}}) |>
  group_by({{group_var}}) |>
  mutate(outlet_comid = comid[which.max(da_reach)]) |>
  mutate(da_outlet = da_reach[which(comid == outlet_comid)],
         pc_outlet = pc_reach[which(comid == outlet_comid)]) |>
  mutate(multiplier = (da_reach / da_outlet) * (pc_reach / pc_outlet)) |>
  arrange({{group_var}}, -multiplier)
}

cv_watersheds_flow_xw <- make_flow_xw(watershed_level_3)
cv_watersheds_flow_xw |> usethis::use_data(overwrite=T)
```

    ## ✔ Setting active project to
    ##   "C:/Users/Skyler/Documents/GitHub/swc-habitat-suitability".

    ## ✔ Saving "cv_watersheds_flow_xw" to "data/cv_watersheds_flow_xw.rda".

    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

``` r
cv_mainstems_flow_xw <- make_flow_xw(river_cvpia)
cv_mainstems_flow_xw |> usethis::use_data(overwrite=T)
```

    ## ✔ Saving "cv_mainstems_flow_xw" to "data/cv_mainstems_flow_xw.rda".
    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

``` r
habistat::flowline_geom_proj |>
  inner_join(cv_watersheds_flow_xw, by=join_by(comid)) |>
  filter(comid %in% habistat::wua_predicted$comid) |>
  ggplot() + 
  geom_sf(aes(color = multiplier)) + 
  geom_sf(data = habistat::flowline_geom_proj |> 
            filter(comid %in% cv_watersheds_flow_xw$outlet_comid) |> 
            st_line_sample(sample=1), 
          aes(color = 1)) + 
  scale_color_viridis_c(direction = -1) +
  ggtitle("Flow Multipliers by Watershed")
```

![](watershed-flow-aggregation_files/figure-gfm/multipliers-1.png)<!-- -->

``` r
habistat::flowline_geom_proj |>
  inner_join(cv_mainstems_flow_xw, by=join_by(comid)) |>
  filter(comid %in% habistat::wua_predicted$comid) |>
  ggplot() + 
  geom_sf(aes(color = multiplier)) + 
  geom_sf(data = habistat::flowline_geom_proj |> 
            filter(comid %in% cv_mainstems_flow_xw$outlet_comid) |> 
            st_line_sample(sample=1), 
          aes(color = 1)) + 
  scale_color_viridis_c(direction = -1) + 
  ggtitle("Flow Multipliers by Mainstem")
```

![](watershed-flow-aggregation_files/figure-gfm/multipliers-2.png)<!-- -->

### Applying the crosswalk to aggregate model results

``` r
wua_predicted_cv_watersheds <- 
  habistat::wua_predicted |>
  inner_join(cv_watersheds_flow_xw, 
             by=join_by(watershed_level_3, comid)) |>
  expand_grid(scaled = c(FALSE, TRUE)) |>
  mutate(flow_cfs = if_else(scaled, flow_cfs * multiplier, flow_cfs)) |>
  arrange(habitat, model_name, comid, flow_cfs) |>
  group_by(habitat, model_name, comid) |>
  mutate(wua_per_lf_pred = if_else(scaled, NA, wua_per_lf_pred)) |>
  mutate(wua_per_lf_pred = zoo::na.approx(wua_per_lf_pred, x = flow_cfs, na.rm=F, rule=2)) |>
  filter(scaled) |>
  select(-scaled) |>
  ungroup() |>
  group_by(habitat, model_name, watershed_level_3, flow_idx) |>
  summarize(wua_per_lf_pred = sum(wua_per_lf_pred * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred = sum(wua_per_lf_pred * reach_length_ft) / 43560, .groups="drop") |>
  inner_join(habistat::wua_predicted |> 
               group_by(flow_idx) |> 
               summarize(flow_cfs = first(flow_cfs), .groups="drop"),
             by = join_by(flow_idx))

wua_predicted_cv_watersheds |> usethis::use_data(overwrite=T)
```

    ## ✔ Saving "wua_predicted_cv_watersheds" to
    ##   "data/wua_predicted_cv_watersheds.rda".

    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

``` r
wua_predicted_cv_mainstems <- 
  habistat::wua_predicted |> 
  inner_join(cv_mainstems_flow_xw, 
             by=join_by(river_cvpia, comid)) |>
  expand_grid(scaled = c(FALSE, TRUE)) |>
  mutate(flow_cfs = if_else(scaled, flow_cfs * multiplier, flow_cfs)) |>
  arrange(habitat, model_name, comid, flow_cfs) |>
  group_by(habitat, model_name, comid) |>
  mutate(wua_per_lf_pred = if_else(scaled, NA, wua_per_lf_pred)) |>
  mutate(wua_per_lf_pred = zoo::na.approx(wua_per_lf_pred, x = flow_cfs, na.rm=F, rule=2)) |>
  filter(scaled) |>
  select(-scaled) |>
  ungroup() |>
  group_by(habitat, model_name, river_cvpia, flow_idx) |>
  summarize(wua_per_lf_pred = sum(wua_per_lf_pred * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred = sum(wua_per_lf_pred * reach_length_ft) / 43560, .groups="drop") |>
  inner_join(habistat::wua_predicted |> 
               group_by(flow_idx) |> 
               summarize(flow_cfs = first(flow_cfs), .groups="drop"),
             by = join_by(flow_idx))

wua_predicted_cv_mainstems |> usethis::use_data(overwrite=T)
```

    ## ✔ Saving "wua_predicted_cv_mainstems" to "data/wua_predicted_cv_mainstems.rda".
    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

``` r
# these are the mainstems grouped for CVPIA model comparison
wua_predicted_cv_mainstems_grouped <- 
  habistat::wua_predicted |> 
  inner_join(cv_mainstems_flow_xw, 
             by=join_by(river_cvpia, comid)) |>
  expand_grid(scaled = c(FALSE, TRUE)) |>
  mutate(flow_cfs = if_else(scaled, flow_cfs * multiplier, flow_cfs)) |>
  arrange(habitat, model_name, comid, flow_cfs) |>
  group_by(habitat, model_name, comid) |>
  mutate(wua_per_lf_pred = if_else(scaled, NA, wua_per_lf_pred)) |>
  mutate(wua_per_lf_pred = zoo::na.approx(wua_per_lf_pred, x = flow_cfs, na.rm=F, rule=2)) |>
  filter(scaled) |>
  select(-scaled) |>
  ungroup() |>
  group_by(habitat, model_name, river_group, flow_idx) |>
  summarize(wua_per_lf_pred = sum(wua_per_lf_pred * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred = sum(wua_per_lf_pred * reach_length_ft) / 43560, .groups="drop") |>
  inner_join(habistat::wua_predicted |> 
               group_by(flow_idx) |> 
               summarize(flow_cfs = first(flow_cfs), .groups="drop"),
             by = join_by(flow_idx))

wua_predicted_cv_mainstems_grouped |> usethis::use_data(overwrite=T)
```

    ## ✔ Saving "wua_predicted_cv_mainstems_grouped" to
    ##   "data/wua_predicted_cv_mainstems_grouped.rda".
    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

Plotting the output

``` r
wua_predicted_cv_watersheds |>
  ggplot(aes(x = flow_cfs)) + 
  facet_wrap(~watershed_level_3, scales="free") + 
  geom_line(aes(y = wua_acres_pred, color = paste(habitat, model_name))) +
  scale_x_log10() +
  theme(legend.position = "top", panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Paired") +
  ylab("Predicted Total Habitat (acres)") + xlab("Flow at Outlet (cfs)") + 
  ggtitle("Watersheds (incl. tributaries)")
```

![](watershed-flow-aggregation_files/figure-gfm/pred-watersheds-new-1.png)<!-- -->

``` r
wua_predicted_cv_mainstems |>
  ggplot(aes(x = flow_cfs)) + 
  facet_wrap(~river_cvpia, scales="free") + 
  geom_line(aes(y = wua_acres_pred, color = paste(habitat, model_name))) +
  scale_x_log10() +
  theme(legend.position = "top", panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Paired") +
  ylab("Predicted Total Habitat (acres)") + xlab("Flow at Outlet (cfs)") + 
  ggtitle("Mainstems")
```

![](watershed-flow-aggregation_files/figure-gfm/pred-watersheds-new-2.png)<!-- -->

``` r
wua_predicted_cv_mainstems_grouped |>
  ggplot(aes(x = flow_cfs)) + 
  facet_wrap(~river_group, scales="free") + 
  geom_line(aes(y = wua_acres_pred, color = paste(habitat, model_name))) +
  scale_x_log10() +
  theme(legend.position = "top", panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Paired") +
  ylab("Predicted Total Habitat (acres)") + xlab("Flow at Outlet (cfs)") + 
  ggtitle("Mainstems (grouped by watershed)")
```

![](watershed-flow-aggregation_files/figure-gfm/pred-watersheds-new-3.png)<!-- -->

Versus the “naive” method

``` r
wua_predicted |>
  filter(!is.na(watershed_level_3)) |>
  group_by(habitat, model_name, watershed_level_3, flow_idx, flow_cfs) |>
  summarize(wua_per_lf_pred = sum(wua_per_lf_pred * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred = sum(wua_per_lf_pred * reach_length_ft) / 43560, .groups="drop") |>
  ggplot(aes(x = flow_cfs)) + 
  facet_wrap(~watershed_level_3, scales="free") + 
  geom_line(aes(y = wua_acres_pred, color = paste(habitat, model_name))) +
  scale_x_log10() +
  theme(legend.position = "top", panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Paired") +
  ylab("Predicted Total Habitat (acres)") + xlab("Flow (cfs)")
```

![](watershed-flow-aggregation_files/figure-gfm/pred-watersheds-old-1.png)<!-- -->

## DSMHabitat Comparison

Moved this over from model_cleaned.Rmd

``` r
mainstems_comid <- 
  read_sf(file.path("/vsizip", here::here("data-raw", "source", "rearing_spatial_data", "nhdplusv2_comid_habitat_xw.shp.zip"))) |>
  janitor::clean_names() |>
  st_zm() |>
  st_transform(st_crs(habistat::flowline_geom_proj)) |>
  mutate(length_ft = st_length(geometry) |> units::set_units("ft") |> units::drop_units()) |>
  filter(str_detect(habitat, "rearing")) |>
  left_join(habistat::flowline_attr |> select(comid, hqt_gradient_class), by=join_by(comid)) |>
  filter(!(watershed %in% c("Sacramento River", "San Joaquin River")))

mainstems <-
  mainstems_comid |>
  group_by(watershed, river, hqt_gradient_class) |>
  summarize() 
```

    ## `summarise()` has grouped output by 'watershed', 'river'. You can override
    ## using the `.groups` argument.

``` r
mainstems_comid |> 
  ggplot() + 
  geom_sf(aes(group=river, color=hqt_gradient_class)) + 
  theme(legend.key.height = unit(12, "point"))
```

![](watershed-flow-aggregation_files/figure-gfm/dsmhabitat-prep-1.png)<!-- -->

``` r
#remotes::install_github("CVPIA-OSC/DSMhabitat")
watersheds <- mainstems |> pull(watershed) |> unique()
watershed_name <- tolower(gsub(pattern = "-| ", replacement = "_", x = watersheds))
watershed_rda_name <- paste(watershed_name, "floodplain", sep = "_")

dsm_habitat_floodplain <- map_df(watershed_rda_name, function(watershed) {
  df <- as.data.frame(do.call(`::`, list(pkg = "DSMhabitat", name = watershed)))
}) |> 
  transmute(river_group = watershed,
            flow_cfs,
            FR_floodplain_m2 = FR_floodplain_acres * 4046.86,
            FR_floodplain_m2_suitable = DSMhabitat::apply_suitability(FR_floodplain_m2),
            FR_floodplain_acres_suitable = FR_floodplain_m2_suitable / 4046.86)

dsm_habitat_instream <- map_df(paste(watershed_name, "instream", sep = "_"), 
                               possibly(function(watershed) {
                                 df <- as.data.frame(do.call(`::`, list(pkg = "DSMhabitat", name = watershed)))
                                 }, otherwise = NULL)) |> 
  transmute(river_group = watershed,
            flow_cfs,
            FR_juv_wua) 

dsm_flows <- bind_rows(dsm_habitat_floodplain, dsm_habitat_instream) |>
  group_by(river_group, flow_cfs) |>
  summarize() |>
  ungroup() |>
  arrange(river_group, flow_cfs)
```

    ## `summarise()` has grouped output by 'river_group'. You can override using the
    ## `.groups` argument.

``` r
dsm_flow_ranges <- 
  dsm_flows |> 
  group_by(river_group) |> 
  summarize(min_flow_cfs = min(flow_cfs), max_flow_cfs = max(flow_cfs))

mainstems_comid |> 
  st_zm() |> 
  filter(comid %in% mainstems_comid$comid) |>
  ggplot() + 
  geom_sf(aes(color=watershed)) + 
  theme(legend.key.height = unit(12, "point"))
```

![](watershed-flow-aggregation_files/figure-gfm/dsmhabitat-join-1.png)<!-- -->

``` r
# combining both floodplain and instream rearing acreages/WUAs for comparison
dsm_habitat_combined <- mainstems |> 
  mutate(length_ft = st_length(geometry) |> units::set_units("ft") |> units::drop_units()) |>
  group_by(watershed) |> 
  summarize(length_ft = sum(length_ft)) |>
  st_drop_geometry() |>
  inner_join(full_join(dsm_habitat_instream, dsm_habitat_floodplain, by=join_by(river_group, flow_cfs)), by=join_by(watershed==river_group)) |>
  group_by(watershed) |>
  arrange(watershed, flow_cfs) |>
  mutate(FR_juv_wua = zoo::na.approx(FR_juv_wua, flow_cfs, na.rm=F),
         FR_floodplain_acres_suitable = zoo::na.approx(FR_floodplain_acres_suitable, flow_cfs, na.rm=F)) |>
  transmute(watershed, flow_cfs, 
            instream_wua_per_lf = coalesce(FR_juv_wua/1000,0),
            instream_suitable_ac = coalesce(FR_juv_wua/1000,0)*length_ft/43560,
            floodplain_wua_per_lf = coalesce(FR_floodplain_acres_suitable,0)*43560/length_ft,
            floodplain_suitable_ac = coalesce(FR_floodplain_acres_suitable,0),
            combined_wua_per_lf = instream_wua_per_lf + floodplain_wua_per_lf,
            combined_suitable_ac =  instream_suitable_ac + floodplain_suitable_ac) |>
  ungroup()
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `FR_floodplain_acres_suitable =
    ##   zoo::na.approx(FR_floodplain_acres_suitable, flow_cfs, na.rm = F)`.
    ## ℹ In group 20: `watershed = "Stanislaus River"`.
    ## Caused by warning in `regularize.values()`:
    ## ! collapsing to unique 'x' values

``` r
dsm_habitat_wua_per_lf <- dsm_habitat_combined |>
  select(watershed, flow_cfs, instream_wua_per_lf, floodplain_wua_per_lf) |>
#  mutate(combined_wua_per_lf = pmax(instream_wua_per_lf, floodplain_wua_per_lf)) |>
  pivot_longer(cols=c(instream_wua_per_lf, floodplain_wua_per_lf)) |>
  mutate(name = paste("DSMhabitat", str_replace(name, "_wua_per_lf", "")),
         value = if_else(value>0, value, NA))

dsm_habitat_suitable_ac <- dsm_habitat_combined |>
  select(watershed, flow_cfs, instream_suitable_ac, floodplain_suitable_ac) |>
#  mutate(combined_suitable_ac = pmax(instream_suitable_ac, floodplain_suitable_ac)) |>
  pivot_longer(cols=c(instream_suitable_ac, floodplain_suitable_ac)) |>
  mutate(value = if_else(value>0, value, NA)) |>
  mutate(name = paste("DSMhabitat", str_replace(name, "_suitable_ac", "")),
         value = if_else(value>0, value, NA))
```

``` r
wua_predicted_cv_mainstems_grouped |>
  filter(habitat == "rearing") |>
  filter(river_group %in% dsm_habitat_wua_per_lf$watershed) |>
  mutate(watershed = river_group) |>
  ggplot(aes(x = flow_cfs)) + 
  facet_wrap(~watershed, scales="free") + 
  geom_line(aes(y = wua_acres_pred, color = paste("habistat", model_name))) +
  geom_line(data=dsm_habitat_suitable_ac |> 
              filter(flow_cfs <= 15000), 
            aes(x = flow_cfs, y = value, color = paste("vs", name))) +
  scale_x_log10() +
  theme(legend.position = "top", panel.grid.minor = element_blank()) + 
  scale_color_brewer(name = "", palette = "Paired") +
  ylab("Predicted Total Habitat (acres)") + xlab("Flow at Outlet (cfs)") + 
    guides(color = guide_legend(nrow = 2))
```

    ## Warning in scale_x_log10(): log-10 transformation introduced infinite values.

    ## Warning: Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

![](watershed-flow-aggregation_files/figure-gfm/pred-watersheds-dsmhabitat-comparison-1.png)<!-- -->

## Applying this method to the duration analyis

- Associate each streamgage with a COMID, drainage area, and DA
  precipitation
- When doing COMID scale duration analysis, use these values to scale
  the x axis of the inundation curve directly
- When doing watershed scale duration analysis
  - pick a streamgage
  - for each selected comid,
    - create a flow-to-suitable-area curve along downscaled flows using
      complete-approx method (above)
    - make a flow-scaled version of the gage’s inundation curve
    - apply the flow-scaled inundation curve to the flow-scaled habitat
      curve
  - consolidate and sum the result curves

### Streamgage drainage area attributes

``` r
streamgage_da_attr <-
  habistat::cv_mainstems |>
  nest(stem_ln_geom = c(comid, habitat, geometry), .by = c(river_group, river_cvpia)) |>
  inner_join(habistat::streamgage_attr, by=join_by(river_group, river_cvpia)) |>
  inner_join(x=habistat::streamgage_geom |> nest(gage_pt_geom = c(geometry)), 
             y=_, 
             by=join_by(station_id)) |>
  mutate(gage_comid = map2_int(gage_pt_geom, stem_ln_geom,
                               function(gage_pt, stem_ln) stem_ln$comid[[st_nearest_feature(gage_pt, stem_ln)]])) |>
  select(station_id, gage_comid) |>
  inner_join(habistat::flowline_attr |>
              select(gage_comid = comid,
                     da_gage = da_area_sq_km,
                     pc_gage = da_ppt_mean_mm), 
             by=join_by(gage_comid)) |>
  glimpse()
```

    ## Rows: 163
    ## Columns: 4
    ## Groups: station_id [163]
    ## $ station_id <chr> "afo", "amf", "amk", "amr", "atc", "bas", "bat", "bcd", "bc…
    ## $ gage_comid <dbl> 948021150, 15024885, 14982590, 2829122, 8019804, 2766098, 2…
    ## $ da_gage    <dbl> 4928.9751, 4833.0171, 493.1172, 2838.5181, 321.4746, 798.25…
    ## $ pc_gage    <dbl> 1349.6419, 1363.6610, 1350.1780, 1019.7141, 1096.5603, 1227…

``` r
streamgage_da_attr |> usethis::use_data(overwrite = T)
```

    ## ✔ Saving "streamgage_da_attr" to "data/streamgage_da_attr.rda".

    ## ☐ Document your data (see <https://r-pkgs.org/data.html>).

### Example application to COMID-level flow-to-suitable-area curve

``` r
selected_comid <- 8020536  
selected_gage <- "dcv"
selected_habitat <- "rearing"
selected_model <- "SD"
selected_run <- "fall"
selected_wy_type <- "Dry"

fsa <-
  habistat::wua_predicted |>
  filter((model_name == selected_model) & 
           (habitat == selected_habitat) &
           (comid == selected_comid)) |>
  glimpse()
```

    ## Rows: 50
    ## Columns: 10
    ## $ comid             <dbl> 8020536, 8020536, 8020536, 8020536, 8020536, 8020536…
    ## $ model_name        <chr> "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD", "SD"…
    ## $ flow_idx          <int> 50, 56, 63, 70, 79, 89, 100, 112, 125, 141, 158, 177…
    ## $ flow_cfs          <dbl> 50.11872, 56.23413, 63.09573, 70.79458, 79.43282, 89…
    ## $ wua_per_lf_pred   <dbl> 0.006930601, 0.008865896, 0.009664865, 0.013621307, …
    ## $ river_group       <fct> Deer Creek, Deer Creek, Deer Creek, Deer Creek, Deer…
    ## $ river_cvpia       <fct> Deer Creek, Deer Creek, Deer Creek, Deer Creek, Deer…
    ## $ watershed_level_3 <fct> Deer Creek, Deer Creek, Deer Creek, Deer Creek, Deer…
    ## $ reach_length_ft   <dbl> 7732.94, 7732.94, 7732.94, 7732.94, 7732.94, 7732.94…
    ## $ habitat           <chr> "rearing", "rearing", "rearing", "rearing", "rearing…

``` r
selected_comid_attr <- 
  habistat::flowline_attr |>
  filter(comid == selected_comid) |>
  as.list()

drc <-
  habistat::streamgage_duration_rating_curves |>
  filter((station_id == selected_gage) & 
           (habitat == selected_habitat) & 
           (run == selected_run) & 
           (wy_group == selected_wy_type)) |>
  unnest(data) 

drc_scaled <- drc |>
  inner_join(streamgage_da_attr, by=join_by(station_id)) |>
  mutate(model_q = model_q * (selected_comid_attr$da_area_sq_km_tot / da_gage) * (selected_comid_attr$da_ppt_mean_mm / pc_gage)) |>
  glimpse()
```

    ## Rows: 150
    ## Columns: 12
    ## $ station_id             <chr> "dcv", "dcv", "dcv", "dcv", "dcv", "dcv", "dcv"…
    ## $ run                    <fct> fall, fall, fall, fall, fall, fall, fall, fall,…
    ## $ habitat                <fct> rearing, rearing, rearing, rearing, rearing, re…
    ## $ wy_group               <fct> Dry, Dry, Dry, Dry, Dry, Dry, Dry, Dry, Dry, Dr…
    ## $ model_q                <dbl> 46.81536, 93.63073, 140.44609, 187.26145, 234.0…
    ## $ avg_max_days_inundated <dbl> 142.0909091, 41.0909091, 14.0909091, 7.1818182,…
    ## $ durhsi_rearing_vl      <dbl> 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66,…
    ## $ durhsi_rearing_vf      <dbl> 0.66, 0.66, 0.66, 1.00, 1.00, 1.00, 1.00, 1.00,…
    ## $ durhsi_spawning        <dbl> 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ gage_comid             <dbl> 8020924, 8020924, 8020924, 8020924, 8020924, 80…
    ## $ da_gage                <dbl> 542.574, 542.574, 542.574, 542.574, 542.574, 54…
    ## $ pc_gage                <dbl> 1484.876, 1484.876, 1484.876, 1484.876, 1484.87…

``` r
habistat::duration_apply_dhsi_to_fsa_curve(
        fsa = fsa,
        fsa_q = flow_cfs,
        fsa_wua = wua_per_lf_pred,
        drc = drc_scaled,
        drc_q = model_q,
        drc_dhsi = durhsi_rearing_vf) |>
        glimpse()
```

    ## Rows: 192
    ## Columns: 10
    ## $ q                       <dbl> 50.11872, 56.23413, 63.09573, 70.79458, 79.432…
    ## $ wua                     <dbl> 0.006930601, 0.008865896, 0.009664865, 0.01362…
    ## $ dhsi                    <dbl> 0.6600000, 0.6600000, 0.6600000, 0.6600000, 0.…
    ## $ marginal_wua            <dbl> 0.0069306007, 0.0019352949, 0.0007989699, 0.00…
    ## $ wua_added               <dbl> 0.0069306007, 0.0019352949, 0.0007989699, 0.00…
    ## $ marginal_durwua_added   <dbl> 0.0045741965, 0.0012772946, 0.0005273201, 0.00…
    ## $ wua_removed             <dbl> 0.000000000, 0.000000000, 0.000000000, 0.00000…
    ## $ marginal_durwua_removed <dbl> 0.000000000, 0.000000000, 0.000000000, 0.00000…
    ## $ marginal_durwua         <dbl> 0.0045741965, 0.0012772946, 0.0005273201, 0.00…
    ## $ durwua                  <dbl> 0.004574196, 0.005851491, 0.006378811, 0.00899…

### Example application to watershed-level flow-to-suitable-area curve

``` r
scale_fsa <- function(data, multiplier, .wua_var = wua_per_lf_pred, .flow_var = flow_cfs) {
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
```

``` r
selected_watershed <- "Deer Creek"
selected_gage <- "dcv"
selected_habitat <- "rearing"
selected_model <- "SD"
selected_run <- "fall"
selected_wy_type <- "Dry"

fsas <-
  habistat::wua_predicted |>
  filter((watershed_level_3 == selected_watershed) &
           (model_name == selected_model) & 
           (habitat == selected_habitat)) |>
  nest(fsa_raw = c(flow_idx, flow_cfs, reach_length_ft, wua_per_lf_pred)) |>
  inner_join(cv_watersheds_flow_xw, by=join_by(watershed_level_3, comid)) |>
  mutate(fsa_scaled = pmap(list(fsa_raw, multiplier), scale_fsa))

drcs <-
  habistat::streamgage_duration_rating_curves |>
  filter((station_id == selected_gage) & 
           (habitat == selected_habitat) & 
           (run == selected_run) & 
           (wy_group == selected_wy_type)) |>
    rename(drc_raw = data) |>
    inner_join(streamgage_da_attr, by=join_by(station_id)) |>
    expand_grid(comid = fsas$comid) |>
    inner_join(cv_watersheds_flow_xw |> 
                 select(comid, watershed_level_3, da_reach, pc_reach), 
               by=join_by(comid)) |>
    mutate(drc_scaled = 
             pmap(list(drc_raw, da_gage, pc_gage), 
                  function(drc_raw, da_gage, pc_gage) {
                        drc_raw |>
                          mutate(model_q = model_q * 
                                   (da_reach / da_gage) * 
                                   (pc_reach / pc_gage))
                    }))
```

    ## Warning: There were 57 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `model_q = model_q * (da_reach/da_gage) * (pc_reach/pc_gage)`.
    ## Caused by warning:
    ## ! There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `model_q = model_q * (da_reach/da_gage) * (pc_reach/pc_gage)`.
    ## Caused by warning in `model_q * (da_reach / da_gage)`:
    ## ! longer object length is not a multiple of shorter object length
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 56 remaining warnings.

``` r
joined <- inner_join(fsas, drcs, by=join_by(comid, watershed_level_3, habitat, da_reach, pc_reach)) |>
  mutate(result = map2(fsa_scaled, drc_scaled, function(fsa, drc) {
    duration_applied <- 
      habistat::duration_apply_dhsi_to_fsa_curve(
        fsa = fsa,
        fsa_q = flow_cfs,
        fsa_wua = wua_per_lf_pred,
        drc = drc,
        drc_q = model_q,
        drc_dhsi = durhsi_rearing_vf)
    fsa |>
      left_join(duration_applied |> 
                   select(flow_cfs = q, wua_per_lf_pred_dur = durwua),
                 by=join_by(flow_cfs)) |>
      mutate(wua_per_lf_pred_dur = zoo::na.approx(wua_per_lf_pred_dur, x = flow_cfs, na.rm=F, rule=2))
  })) |>
  select(-fsa_raw, -drc_raw, -fsa_scaled, -drc_scaled) |>
  unnest(result)

joined |>
  ggplot(aes(x = flow_cfs)) +
  facet_wrap(~paste0(comid,"\n",round(da_reach), " km2"), scales = "free_x") +
  geom_line(aes(y = wua_per_lf_pred, color = "original")) +
  geom_line(aes(y = wua_per_lf_pred_dur, color = "duration")) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Flow (cfs)") + scale_x_continuous(breaks = scales::breaks_extended(n=3)) + ylab("Suitable Habitat Area (ft2) per LF")
```

![](watershed-flow-aggregation_files/figure-gfm/example-watershed-1.png)<!-- -->

``` r
summarized <-
  joined |> 
  group_by(model_name, habitat, run, wy_group,
           watershed_level_3, flow_idx) |>
  summarize(wua_per_lf_pred = sum(wua_per_lf_pred * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred = sum(wua_per_lf_pred * reach_length_ft) / 43560, 
            wua_per_lf_pred_dur = sum(wua_per_lf_pred_dur * reach_length_ft) / sum(reach_length_ft),
            wua_acres_pred_dur = sum(wua_per_lf_pred_dur * reach_length_ft) / 43560, 
            .groups="drop")

summarized |>
  ggplot(aes(x = flow_idx)) +
  geom_line(aes(y = wua_acres_pred, color = "original")) +
  geom_line(aes(y = wua_acres_pred_dur, color = "duration")) +
  ggtitle(glue::glue("{selected_watershed} watershed total calculated using gage {selected_gage}")) +
  scale_x_log10() + theme(panel.grid.minor = element_blank()) +
  xlab("Flow (cfs) at Outlet") + ylab("Suitable Habitat Area (acres)")
```

![](watershed-flow-aggregation_files/figure-gfm/example-watershed-summary-1.png)<!-- -->

``` r
knitr::knit_exit()
```
