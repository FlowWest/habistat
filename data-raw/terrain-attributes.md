Terrain-Derived Predictors
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-01-24

- [Load DEM](#load-dem)
- [Channel Confinement via mTPI](#channel-confinement-via-mtpi)
  - [Valley Bottom Width via Slope Cutoff
    Method](#valley-bottom-width-via-slope-cutoff-method)
  - [Valley Bottom Width via VBET](#valley-bottom-width-via-vbet)

## Load DEM

Data Source: 10m NHDPlusHR DEM

This input has already been clipped to the AOI, converted from cm to m,
and crs redefined

Future improvement: script should pull straight from the source geotiff
and apply these changes programatically

``` r
dem <- read_stars("nhdplushr/dem_nhdplushr_yuba_meters_v2.tif") |>
  rename(elev = dem_nhdplushr_yuba_meters_v2.tif)

ggplot() + geom_stars(data=dem) + coord_fixed()
```

![](terrain-attributes_files/figure-gfm/import-dem-1.png)<!-- -->

Example calculating slope (can use `stars` or `terra`)

``` r
slope <- dem |> starsExtra::slope()

#slope <- dem |> terra::rast() |> terra::slope() |> st_as_stars()

ggplot() + geom_stars(data=slope) + coord_fixed()
```

![](terrain-attributes_files/figure-gfm/calc-slope-1.png)<!-- -->

## Channel Confinement via mTPI

mTPI as a measure of topographic confinement

Souce: Theobald DM, Harrison-Atlas D, Monahan WB, Albano CM. 2015.
Ecologically-relevant maps of landforms and physiographic diversity for
climate adaptation planning. PLOS ONE. DOI: 10.1371/journal.pone.0143619

``` r
dem_rast <- terra::rast(dem) 
cell_size <- 10 # mean(terra::res(dem_rast))
tpi_90 <- dem_rast - terra::focal(dem_rast, w=90/cell_size, fun="mean")
tpi_270 <- dem_rast - terra::focal(dem_rast, w=270/cell_size, fun="mean")
tpi_810 <- dem_rast - terra::focal(dem_rast, w=810/cell_size, fun="mean")
#tpi_2430 <- dem_rast - terra::focal(dem_rast, w=2430/cell_size, fun="mean")
mtpi <- terra::app(c(tpi_90, tpi_270, tpi_810), mean)
mtpi |> st_as_stars() |> plot()
```

    ## downsample set to 1

![](terrain-attributes_files/figure-gfm/calc-mtpi-1.png)<!-- -->

``` r
flow_buffered <- terra::vect(flowlines) |> terra::buffer(width=50)
vec_mtpi_min <- terra::zonal(mtpi, flow_buffered, "min") |> rename(mtpi_min = mean)

mtpi_comid <- flowlines |>
  select(comid) |>
  mutate(vec_mtpi_min) |>
  st_drop_geometry() |> 
  select(comid, mtpi_min) |> 
  filter(!is.nan(mtpi_min))

flowlines |>
  inner_join(mtpi_comid) |>
  st_zm() |>
  ggplot() + 
  geom_sf(aes(color = mtpi_min), linewidth=1) + 
  scale_color_viridis_c(direction=-1)
```

    ## Joining with `by = join_by(comid)`

![](terrain-attributes_files/figure-gfm/intersect-mtpi-1.png)<!-- -->

### Valley Bottom Width via Slope Cutoff Method

Method described in
<https://watermanagement.ucdavis.edu/download_file/view_inline/144>

``` r
valley_bottom <- function(catchment, flowline, dem) {
  dem |> 
    st_crop(catchment) |> 
    starsExtra::slope() |>
    mutate(valley_bottom = if_else(slope<units::set_units(atan(0.25)*180/pi,"degrees"), 1, NA)) |>
    select(valley_bottom) |>
    st_as_sf(merge=TRUE) |>
    st_filter(st_zm(flowline)) #|> st_as_sfc()
}

valley_bottom_by_comid <- function(x) {
  c <- catchments |> filter(comid==x)
  f <- flowlines |> filter(comid==x)
  valley_bottom(c, f, dem)
}

valley_width_naive <- function(geom_valley_bottom) {
   geom_valley_bottom |> 
    terra::vect() |> 
    terra::width()
}

perpendicular_transect <- function(pair=st_multipoint(), length=numeric()) {
  x1 <- pair[1]
  x2 <- pair[2]
  y1 <- pair[3]
  y2 <- pair[4]
  midpoint_x <- (x1+x2)/2
  midpoint_y <- (y1+y2)/2
  orig_bearing <- atan2((y2-y1), (x2-x1))
  perp_bearing <- orig_bearing + pi/2
  radius <- length/2
  perp_x1 <- midpoint_x + radius*cos(perp_bearing)
  perp_y1 <- midpoint_y + radius*sin(perp_bearing)
  perp_x2 <- midpoint_x - radius*cos(perp_bearing)
  perp_y2 <- midpoint_y - radius*sin(perp_bearing)
  return(st_linestring(c(st_point(c(perp_x1,perp_y1)), 
                         st_point(c(perp_x2, perp_y2)))))
}

perpendicular_transects <- function(ls=st_linestring(), length=numeric()) {
  if ("sf" %in% class(ls)){
    ls <- st_as_sfc(ls, crs=st_crs(ls))
  }
  if ("sfc" %in% class(ls)){
    m <- t(st_zm(ls[[1]]))
  } 
  if ("sfg" %in% class(ls)){
    m <- t(ls)
  }
  point_pairs <- lapply(seq(1, length(m)-3, 2), function(i) c(st_point(c(m[i], m[i+1])), st_point(c(m[i+2], m[i+3]))))
  perp_lines <- lapply(point_pairs, function(pair) perpendicular_transect(pair, length))
  return(st_sfc(perp_lines, crs=st_crs(ls)))
}

valley_width_via_transects <- function(geom_valley_bottom, flowline, max_width) {
  transects <- perpendicular_transects(flowline, length=max_width)
  if("sfg" %in% class(geom_valley_bottom)) {
    geom_valley_bottom <- st_sfc(geom_valley_bottom, crs=st_crs(transects))
  }
  clipped <- st_intersection(transects, geom_valley_bottom)
  return(mean(st_length(clipped)))
}
```

``` r
flowline <- flowlines |> filter(comid==8062593) |> st_zm()
catchment <- catchments |> filter(comid==8062593) |> st_transform(st_crs(flowline))
ggplot() + 
  #geom_stars(data=catchment_confinement(catchment, flowline, dem)) + 
  geom_sf(data=valley_bottom(catchment, flowline, dem)) +
  geom_sf(data=st_zm(flowline))
```

![](terrain-attributes_files/figure-gfm/vbw-example1-1.png)<!-- -->

``` r
valley_width_naive(valley_bottom(catchment, flowline, dem))
```

    ## [1] 666.9869

``` r
valley_polygon <- valley_bottom(catchment, flowline, dem)
valley_width_via_transects(valley_bottom(catchment, flowline, dem), flowline, max_width=1000)
```

    ## 148.7645 [m]

``` r
ggplot() + 
  geom_sf(data=catchment) +
  geom_sf(data=valley_polygon, fill="lightblue") +
  geom_sf(data=flowline, color="blue") +
  geom_sf(data=perpendicular_transects(flowline, length=1000) |> st_intersection(valley_polygon), color="red") 
```

![](terrain-attributes_files/figure-gfm/vbw-example1-2.png)<!-- -->

``` r
vb1 <- flowlines |> 
  filter(huc_8 %in% selected_huc_8) |> 
  filter(comid==8062593) |> # comment out the filter to process all
  mutate(flowline = st_set_crs(st_zm(geometry),project_crs)) |>
  st_drop_geometry() |>
  mutate(valley_bottom_geoms = st_as_sfc(bind_rows(map(comid, possibly(function(x) valley_bottom_by_comid(x), otherwise=NA))), crs=project_crs)) |>
  mutate(vb_width_simple = map_dbl(valley_bottom_geoms, possibly(function(x) valley_width_naive(x), otherwise=NA))) |> 
  mutate(vb_width_transect = map2_dbl(valley_bottom_geoms, flowline, function(x, y) valley_width_via_transects(x, y, 1000)),
         vb_width_transect = if_else(!is.nan(vb_width_transect),vb_width_transect,NA)) |> 
  select(comid, flowline, valley_bottom_geoms, vb_width_simple, vb_width_transect) 
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `vb_width_simple = map_dbl(...)`.
    ## Caused by warning:
    ## ! [perim] unknown CRS. Results can be wrong

``` r
ggplot() +  geom_sf(data=vb1$valley_bottom_geoms) + geom_sf(data=vb1$flowline)
```

![](terrain-attributes_files/figure-gfm/vbw-batch-1.png)<!-- -->

### Valley Bottom Width via VBET

Extract flowlines to use in vbet

``` r
flowlines_for_vbet <- flowlines |>
  filter(gnis_name %in% c("Yuba River")) |>
  st_transform(st_crs(dem)) |>
  st_zm() |> 
  st_crop(st_bbox(dem)) |>
  smoothr::densify(5)
```

    ## Warning: attribute variables are assumed to be spatially constant throughout
    ## all geometries