library(tidyverse)
library(sf)
library(future)
library(future.apply)
library(furrr)
#source(here::here("data-raw", "scripts", "data-functions.R"))
#source(here::here("data-raw", "scripts", "suitability-functions.R"))
library(habistat)

flowlines <- readRDS(here::here("data-raw", "results", "flowline_geometries_proj.Rds"))

googledrive::drive_auth() # run this before proceeding with rest of script

# STANISLAUS RIVER -------------------------------------------------------------

stan_dir <- here::here("data-raw", "temp", "stanislaus")
dir.create(stan_dir, recursive = TRUE)
drive_file_by_id("1J7Iw-PHdGuy5o-VdilwZzkv02Qn7ajks") |>
  zip::unzip(overwrite = T, exdir = stan_dir)

# SRH2D model domain, dissolved from SRH2D mesh faces using QGIS
stan_domain <-
  st_read(file.path("/vsizip", stan_dir, "StanMesh072313_Domain.shp.zip"), as_tibble=T)

# Model domain manually broken up into corresponding COMIDs in GIS
# TODO: Is there a programmatic way we can do this delineation?
stan_comid <-
  st_read(file.path("/vsizip", stan_dir, "StanMesh072313_Domain_COMID.shp.zip"), as_tibble=T) |>
  st_zm() |>
  janitor::clean_names()

# Mesh vertex points identified by sequential VID number
stan_vertices <- st_read(file.path("/vsizip", stan_dir, "StanMesh072313_Vertices.shp.zip"), as_tibble=T) |>
  mutate(vid = row_number()) |>
  select(vid)

# Thiessen (aka Voronoi) polygons for mesh vertices, generated using QGIS
stan_thiessen <- st_read(file.path("/vsizip", stan_dir, "StanMesh072313_Thiessen.shp.zip"), as_tibble=T) |>
  st_join(stan_vertices, join=st_nearest_feature) |> # row order doesn't match the SRH2D outputs so need to spatial join
  select(vid) |>
  arrange(vid)

stan_csvs <- tribble(~flow_cfs, ~filename,
                       500, "500cfs_072313.csv.gz",
                       750, "750cfs_072413.csv.gz",
                      1000, "1000cfs_072413.csv.gz",
                      1250, "1250cfs_090313.csv.gz",
                      1500, "1500cfs_071013.csv.gz",
                      1750, "1750cfs_090313.csv.gz",
                      2250, "2250cfs_121713.csv.gz",
                      3000, "3000cfs_061913.csv.gz",
                      5000, "5000cfs_071113.csv.gz") |>
  mutate(filename = file.path(stan_dir, filename))

# ### ORIGINAL VERSION
#
# stan_result <-
#   vector_summarize_hsi(stan_mesh_prepped, stan_hsi_by_flow) |>
#   suitability_postprocess(stan_comid, .group_var = comid)
#
# outpath <- here::here("data-raw", "results", "fsa_stan.Rds")
#
# stan_result |> saveRDS(outpath)

### VARIATION FOR REARING - BASEFLOW REMOVAL

stan_hsi_by_flow <-
  stan_csvs |>
  vector_import_srh2d(units = "m") |>
  vector_calculate_hsi(hsi_func = vector_dvhsi_hqt)

# define baseflow mask as 500 cfs (COMID 2819818 observed dry season baseflow is 306 cfs)
stan_baseflow_poly <- stan_thiessen |>
  filter(vid %in% filter(stan_hsi_by_flow, inundated & flow_cfs == 500)$vid) |>
  summarize() |>
  st_union()

stan_mesh_prepped_nbfc <-
  stan_thiessen |>
  vector_prep_mesh(group_polygons = stan_comid,
                   mask_polygons = stan_baseflow_poly,
                   .group_var = comid,
                   .id_var = vid)

stan_result_nbfc <-
  vector_summarize_hsi(stan_mesh_prepped_nbfc, stan_hsi_by_flow) |>
  suitability_postprocess(stan_comid, .group_var = comid)

outpath <- here::here("data-raw", "results", "fsa_stan_nbfc.Rds")

stan_result_nbfc |> saveRDS(outpath)

### VARIATION FOR SPAWNING - ALTERNATIVE HSI

stan_mesh_prepped <-
  stan_thiessen |>
  vector_prep_mesh(group_polygons = stan_comid,
                   .group_var = comid,
                   .id_var = vid)

stan_hsi_by_flow_spawning <-
  stan_csvs |>
  vector_import_srh2d(units = "m") |>
  vector_calculate_hsi(hsi_func = vector_dvhsi_spawning)

stan_result_spawning <-
  vector_summarize_hsi(stan_mesh_prepped, stan_hsi_by_flow_spawning) |>
  suitability_postprocess(stan_comid, .group_var = comid)

outpath <- here::here("data-raw", "results", "fsa_stan_spawning.Rds")

stan_result_spawning |> saveRDS(outpath)


# YUBA RIVER -------------------------------------------------------------

plan(multisession, workers = availableCores() - 1) # split processing by reach across sessions
plan(sequential)

yuba_dir <- here::here("data-raw", "temp", "yuba")
dir.create(yuba_dir, recursive = TRUE)
yuba_csv_results <- drive_file_by_id("13Sw8WRJr7lGElgx2MfzXAYYy0OFQr2U5", dir=yuba_dir)
yuba_csv_results |> archive::archive_extract(dir=yuba_dir)

yuba_comid <- drive_file_by_id("1-xSi142jtNZKQS9-VgXZT1yzTfH1KjyR", vsizip=T) |>
  st_read(as_tibble=T) |>
  st_zm() |>
  select(reach = ReachLYR, comid = COMID)

yuba_reaches <- c("EDR", "TBR", "HR", "DGR", "FR")

yuba_flows <- c(300, 350, 400, 450, 530, 600, 622, 700, 800, 880, 930,
                1000, 1300, 1500, 1700, 2000, 2500, 3000, 4000, 5000, 7500,
                10000, 15000, 21100, 30000, 42200, 84400, 110400)

yuba_srh2d_result_by_reach <-
  expand_grid(reach = yuba_reaches, flow_cfs = yuba_flows) |>
  mutate(filename = file.path(yuba_dir, glue::glue("{reach}_{flow_cfs}_SMS.csv"))) |>
  filter(file.exists(filename)) |>
  nest(tbl_filenames = c(flow_cfs, filename)) |>
  mutate(tbl_results = future_map(tbl_filenames,
                           function(x) vector_import_srh2d(x, units="ft", has_vid_col=T)))

yuba_srh2d_result_by_reach |> saveRDS(file.path(yuba_dir, "yuba_srh2d_result_by_reach.Rds"))

yuba_hsi_by_flow_by_reach <-
  yuba_srh2d_result_by_reach |>
  mutate(tbl_hsi = future_map(tbl_results,
                       function(x) vector_calculate_hsi(x, hsi_func = vector_dvhsi_hqt))) |>
  select(reach, tbl_hsi)

yuba_hsi_by_flow_by_reach |> saveRDS(file.path(yuba_dir, "yuba_hsi_by_flow_by_reach.Rds"))

yuba_datasets <-
  tribble(~reach, ~vertices, ~thiessen,
          "EDR", "1t_oT37auM_25SnU3TLrXa_xCyoV9H9qt", "1s-FhTTDlclwcuu-H-9wIwwC1K7lEijXI",
          "TBR", "1tHAvrhIaw_Gz47K5pL0SpeniEMYnczln", "1E7ucN5WOyynzNjeYTiPiS_ixbaJXA6cW",
          "HR" , "1cByc4uVCf46lwVnvfECmKi3IOsR3T2sU", "1ca5tpyAAe2DzeiVStl4Cl8qEvwI7Jaxx",
          "DGR", "1aSsCAHVwRyquDH_-tA89-7qObJ4eqc-v", "1aCvZAd7YqU5kzQfshmxIBNECWQqRJKoC",
          "FR" , "1teWo-kp-ujtKtHfvoLfZZCKcGsh0twiD", "1lkGlWPqBDHHP8yN8U-lkxKBnkRKmfVgS")

yuba_mesh_data_by_reach <-
  yuba_datasets |>
  mutate(filename_vertices = future_map_chr(vertices,
                                     function(x) drive_file_by_id(x, vsizip=T)),
         filename_thiessen = future_map_chr(thiessen,
                                     function(x) drive_file_by_id(x, vsizip=T))) |>
  mutate(geom_vertices = future_map(filename_vertices,
                             function(x) st_read(x, as_tibble=T)),
         geom_thiessen = future_map(filename_thiessen,
                             function(x) st_read(x, as_tibble=T)))

# ### ORIGINAL VERSION
#
# yuba_output_by_reach <-
#   inner_join(yuba_hsi_by_flow_by_reach, yuba_mesh_prepped_by_reach,
#              by=join_by(reach), relationship="one-to-one") |>
#   mutate(hsi_output = pmap(list(mesh_prepped, tbl_hsi),
#                            function(x, y) vector_summarize_hsi(x, y, .group_var=comid, .id_var=vid))) |>
#   mutate(hsi_final = pmap(list(hsi_output, reach),
#                           function(x, y) suitability_postprocess(x, yuba_comid |> filter(reach==y), .group_var=comid))) |>
#   select(reach, hsi_final)
#
# yuba_output_by_reach |> saveRDS(file.path(yuba_dir, "yuba_output_by_reach.Rds"))
#
# yuba_result <-
#   yuba_output_by_reach |>
#   select(reach, hsi_final) |>
#   unnest(hsi_final)
#
# outpath <- here::here("data-raw", "results", "fsa_yuba.Rds")
#
# yuba_result |> saveRDS(outpath)
#
### VARIATION FOR REARING - BASEFLOW REMOVAL

if(!file.exists(file.path(yuba_dir, "yuba_mesh_prepped_by_reach_nbfc.Rds"))) {

  yuba_baseflow_poly <-
    yuba_mesh_data_by_reach |>
      mutate(result = future_pmap(list(geom_vertices, geom_thiessen, reach),
                                  function(vtx, thp, rch) {
                                    thp |>
                                      st_join(vtx |>
                                                transmute(vid=row_number()),
                                              join=st_contains_properly, left = FALSE) |>
                                      arrange(vid) |>
                                      filter(vid %in% filter(unnest(yuba_hsi_by_flow_by_reach, tbl_hsi), inundated & flow_cfs == 700 & reach == rch)$vid) |>
                                      summarize() |>
                                      st_union()
                                  })) |>
    pull(result) |>
    data.table::rbindlist() |>
    st_sfc(crs = st_crs(yuba_comid)) |>
    st_union()

  yuba_mesh_prepped_by_reach_nbfc <-
    yuba_mesh_data_by_reach |>
    mutate(mesh_prepped = future_pmap(list(geom_vertices, geom_thiessen, reach),
                                      function(vtx, thp, rch) {
                                        thp |>
                                          st_join(vtx |>
                                                    transmute(vid=row_number()),
                                                  join=st_contains_properly, left = FALSE) |> #st_nearest_feature) |>
                                          arrange(vid) |>
                                          vector_prep_mesh(group_polygons = yuba_comid |> filter(reach==rch),
                                                           mask_polygons = yuba_baseflow_poly,
                                                           .group_var = comid,
                                                           .id_var = vid)
                                      })) |>
    select(reach, mesh_prepped)

  yuba_mesh_prepped_by_reach_nbfc |> saveRDS(file.path(yuba_dir, "yuba_mesh_prepped_by_reach_nbfc.Rds"))

} else {

  yuba_mesh_prepped_by_reach_nbfc <- readRDS(file.path(yuba_dir, "yuba_mesh_prepped_by_reach_nbfc.Rds"))

}

yuba_output_by_reach_nbfc <-
  inner_join(yuba_hsi_by_flow_by_reach, yuba_mesh_prepped_by_reach_nbfc,
             by=join_by(reach), relationship="one-to-one") |>
  mutate(hsi_output = pmap(list(mesh_prepped, tbl_hsi),
                           function(x, y) vector_summarize_hsi(x, y, .group_var=comid, .id_var=vid))) |>
  mutate(hsi_final = pmap(list(hsi_output, reach),
                          function(x, y) suitability_postprocess(x, yuba_comid |> filter(reach==y), .group_var=comid))) |>
  select(reach, hsi_final)

yuba_output_by_reach_nbfc |> saveRDS(file.path(yuba_dir, "yuba_output_by_reach_nbfc.Rds"))

yuba_result_nbfc <-
  yuba_output_by_reach_nbfc |>
  select(reach, hsi_final) |>
  unnest(hsi_final)

outpath <- here::here("data-raw", "results", "fsa_yuba_nbfc.Rds")

yuba_result_nbfc |> saveRDS(outpath)

### VARIATION FOR SPAWNING - ALTERNATIVE HSI

if(!file.exists(file.path(yuba_dir, "yuba_mesh_prepped_by_reach.Rds"))) {

  yuba_mesh_prepped_by_reach <-
    yuba_mesh_data_by_reach |>
    mutate(mesh_prepped = future_pmap(list(geom_vertices, geom_thiessen, reach),
                                      function(vtx, thp, rch) {
                                        thp |>
                                          st_join(vtx |>
                                                    transmute(vid=row_number()),
                                                  join=st_contains_properly, left = FALSE) |> #st_nearest_feature) |>
                                          arrange(vid) |>
                                          vector_prep_mesh(group_polygons = yuba_comid |> filter(reach==rch),
                                                           .group_var = comid,
                                                           .id_var = vid)
                                      })) |>
    select(reach, mesh_prepped)

  yuba_mesh_prepped_by_reach |> saveRDS(file.path(yuba_dir, "yuba_mesh_prepped_by_reach.Rds"))

} else {

  yuba_mesh_prepped_by_reach <- readRDS(file.path(yuba_dir, "yuba_mesh_prepped_by_reach.Rds"))

}

yuba_hsi_by_flow_by_reach_spawning <-
  yuba_srh2d_result_by_reach |>
  mutate(tbl_hsi = future_map(tbl_results,
                              function(x) vector_calculate_hsi(x, hsi_func = vector_dvhsi_spawning))) |>
  select(reach, tbl_hsi)

yuba_hsi_by_flow_by_reach_spawning |> saveRDS(file.path(yuba_dir, "yuba_hsi_by_flow_by_reach_spawning.Rds"))

yuba_output_by_reach_spawning <-
  inner_join(yuba_hsi_by_flow_by_reach_spawning, yuba_mesh_prepped_by_reach,
             by=join_by(reach), relationship="one-to-one") |>
  mutate(hsi_output = pmap(list(mesh_prepped, tbl_hsi),
                           function(x, y) vector_summarize_hsi(x, y, .group_var=comid, .id_var=vid))) |>
  mutate(hsi_final = pmap(list(hsi_output, reach),
                          function(x, y) suitability_postprocess(x, yuba_comid |> filter(reach==y), .group_var=comid))) |>
  select(reach, hsi_final)

yuba_output_by_reach_spawning |> saveRDS(file.path(yuba_dir, "yuba_output_by_reach_spawning.Rds"))

yuba_result_spawning <-
  yuba_output_by_reach_spawning |>
  select(reach, hsi_final) |>
  unnest(hsi_final)

outpath <- here::here("data-raw", "results", "fsa_yuba_spawning.Rds")

yuba_result_spawning |> saveRDS(outpath)
