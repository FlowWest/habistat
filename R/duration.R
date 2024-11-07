#' Retrieve CDEC CSV
#'
#' @param sta 3-character alphanumeric code for CDEC station
#' @param sen Numeric code for CDEC sensor
#' @param dur 1-character code for duration type: `D`=daily, `H`=hourly, `E`=event (15 minute)
#' @param start_date Start date for query expressed as R date object e.g. lubridate::ymd("1900-10-01")
#' @param end_date End date for query expressed as R date object e.g. lubridate::ymd("2023-09-30")
#' @param dir Directory to store outputs in
#'
#' @returns The path to the downloaded CDEC file in `csv.gz` (compressed CSV) format
#' @md
#' @import httr
#' @importFrom readr read_csv
#'
#' @examples
#' @export
cdec_csv_retrieve <- function(sta=character(), sen=character(), dur=character(),
                              start_date=lubridate::ymd("1900-10-01"), end_date=lubridate::ymd("2023-09-30"),
                              dir="temp") {
  name <- str_to_upper(paste0(sta, "_", sen, "_", dur))
  dir.create(dir, recursive=T)
  filename <- file.path(dir, paste0(name,".csv.gz"))
  if(!file.exists(filename)){
    message(paste0("downloading to ",filename))
    dir.create(dir, recursive = TRUE)
    data_raw <-
      httr::GET(
        url="https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet",
        query=list(Stations=sta, SensorNums=sen, dur_code=dur,
                   Start=format(start_date,"%Y-%m-%d"), End=format(end_date,"%Y-%m-%d"))) |>
      httr::content("raw")
    gzf <- gzcon(file(filename, "wb"))
    data_raw |> writeBin(gzf)
    close(gzf)
  } else {
    message(paste0(filename, " exists, loading"))
  }
  return(filename)
}

#' Parse CDEC CSV
#'
#' @param filename Path to a downloaded CDEC file in the format output by the CDEC `CSVDataServlet` service, such as retrieved by the `retrieve_cdec_csv()` function. May be compressed (`.csv.gz`) or uncompressed (`.csv`).
#'
#' @returns A `tbl_df` data frame with cleaned tabular data
#' @md
#'
#' @examples
#' @export
cdec_csv_parse <- function(filename) {
  read_csv(filename,
           col_select = c(date_time = "DATE TIME", value = "VALUE"),
           col_types = list(col_datetime(), col_double())) |>
    janitor::clean_names() |>
    mutate(water_year = if_else(month(date_time)>=10, year(date_time)+1, year(date_time)))
}

#' Calculate Days Inundated
#'
#' @param q_gauge A vector of flow values, typically a time series
#' @param q_crit A numeric value against which to calculate the exceedance
#' @param stat Either "max" to calculate maximum days inundated (default) or "sum" to calculate total days inundated
#'
#' @returns A single numeric value for the specified statistic.
#' @md
#'
#' @examples
#' @export
duration_calc_days_inundated <- function(q_gauge, q_crit, stat="max") {

  n_days <- length(q_gauge)

  exceedence_intervals <-
    tibble(exceeds_q = q_gauge > q_crit,
           exceedence_event = with(rle(exceeds_q), rep(seq_along(lengths), lengths))) |>
    group_by(exceedence_event) |>
    mutate(cumulative_exceedence = if_else(exceeds_q, seq_along(exceedence_event), 0)) |>
    filter(exceeds_q) |>
    summarize(duration_days = max(cumulative_exceedence))

  if(stat=="sum"){
    return(sum(exceedence_intervals$duration_days))
  } else if(stat=="max") {
    return(pmax(max(exceedence_intervals$duration_days),0))
  }
}

#' Apply Duration HSI to Flow-to-Suitable-Area Curve
#'
#' @param fsa A `data.frame` or `tbl_df` containing a flow-to-suitable-area curve
#' @param drc A `data.frame` or `tbl_df` containing a flow-to-duration-HSI curve
#' @param fsa_q The unquoted name of the flow column in `fsa` (defaults to `q`)
#' @param fsa_wua The unquoted name of the suitable area column in `fsa` (defaults to `wua`)
#' @param drc_q The unquoted name of the flow column in `drc` (defaults to `q`)
#' @param drc_dhsi The unquoted name of the duration suitability factor column in `drc` (defaults to `dhsi`)
#'
#' @returns A `tbl_df` data frame with one row per flow (`q`) containing columns for the original suitable area (`wua`) and the new duration-weighted suitable area (`durwua`) along with interim calculation columns.
#' @md
#'
#' @export
duration_apply_dhsi_to_fsa_curve <- function(fsa, drc,
                                         fsa_q=q, fsa_wua=wua,
                                         drc_q=q, drc_dhsi=dhsi) {

  fun <- function(fsa, drc,
                  fsa_q={{fsa_q}}, fsa_wua={{fsa_wua}},
                  drc_q={{drc_q}}, drc_dhsi={{drc_dhsi}}) {

    full_join(fsa |> transmute(q={{fsa_q}}, wua={{fsa_wua}}),
              drc |> transmute(q={{drc_q}}, dhsi={{drc_dhsi}}),
              by = join_by(q)) |>
      # fill in any gaps with linear interpolation
      arrange(q) |>
      mutate(wua = zoo::na.approx(wua, q, na.rm=F),
             dhsi = zoo::na.approx(dhsi, q, na.rm=F)) |>
      drop_na() |>
      # start by computing the marginal WUA added at each step
      mutate(marginal_wua = if_else(row_number()==1, wua, wua - lag(wua))) |>
      # apply each flow's DHSI to whatever area is added or removed at that flow
      ### mutate(marginal_durwua = pmin(dhsi * marginal_wua, marginal_wua),
      ###       # pmin is a simplification to account for loss; this will end up overestimating the magnitude of the marginal drop in habitat
      ###       durwua = pmax(cumsum(marginal_durwua), 0))
      mutate(wua_added = if_else(marginal_wua > 0, marginal_wua, 0),
             marginal_durwua_added = dhsi * wua_added,
             wua_removed = (-1) * if_else(marginal_wua < 0, marginal_wua, 0),
             marginal_durwua_removed = cummean(dhsi) * wua_removed,
             ###        marginal_durwua_removed = cummean.na.rm(if_else(dhsi > 0, dhsi, NA)) * wua_removed,
             marginal_durwua = marginal_durwua_added - marginal_durwua_removed,
             durwua = pmin(pmax(cumsum(marginal_durwua), 0), wua))
  }

  fun <- possibly(fun, otherwise=tibble(q = numeric(0),
                                        dhsi = numeric(0),
                                        wua = numeric(0),
                                        durwua = numeric(0))) # apply error handling
  # note: this error handling method masks error details that occur within this function

  return(fun(fsa, drc,
             fsa_q={{fsa_q}}, fsa_wua={{fsa_wua}},
             drc_q={{drc_q}}, drc_dhsi={{drc_dhsi}}))
}

#' Remove Baseflow via Flow Threshold (Single)
#'
#' @param data A `data.frame` or `tbl_df` containing a flow-to-suitable-area curve
#' @param flow_threshold A numeric value for a flow threshold defining baseflow.
#' @param .flow_var Unquoted column name for the flow in cfs, defaults to `flow_cfs`
#' @param .habitat_var Unquoted column name for the suitable area, defaults to `wua_per_lf`
#'
#' @returns A `tbl_df` data frame with one row per flow containing two columns: the flow and the modified suitable habitat area with baseflow removed.
#' @md
#'
#' @export
flow_threshold_remove_baseflow <- function(data, flow_threshold,
                                           .flow_var = flow_cfs, .habitat_var = wua_per_lf) {
  # # reuses the code already written for the DSHI analysis to apply a baseflow threshold
  #
  # weights <- tibble(q = c(0, flow_threshold), w = c(0, 1))
  #
  # result <- duration_apply_dhsi_to_fsa_curve(drc = weights, drc_q = q, drc_dhsi = w,
  #                                            fsa = data, fsa_q = .flow_var, fsa_wua = .habitat_var)
  #
  # return(result)

  min_q <- 0
  max_q <- data |> pull({{.flow_var}}) |> max()

  baseflow_removal_drc <-
    data |>
    #complete({{.flow_var}} := flow_threshold) |>
    arrange({{.flow_var}}) |>
    mutate(dhsi = if_else({{.flow_var}} <= flow_threshold, 0, 1))

  result <-
    duration_apply_dhsi_to_fsa_curve(fsa=data, drc=baseflow_removal_drc,
                                     fsa_q={{.flow_var}}, fsa_wua={{.habitat_var}},
                                     drc_q={{.flow_var}}, drc_dhsi = dhsi)

  #if(!is.na(result[[1]])) {
     return(result |>
              transmute({{.flow_var}} := q,
                        {{.habitat_var}} := durwua) |>
              rename_with(function(x) paste0(x, "_nbfc"), everything()))
  #} else {
  #  return(tibble({{.flow_var}} := c(NA),
  #                {{.habitat_var}} := c(NA)))
  #}
}

#' Remove Baseflow via Flow Threshold (Multiple)
#'
#' @param data A `data.frame` or `tbl_df` containing flow-to-suitable-area curves for multiple groups. All variables except for the `.flow_var` and `.habitat_var` will be treated as groups.
#' @param .flow_threshold_var Unquoted column name for the flow threshold defining baseflow for each group.
#' @param .flow_var Unquoted column name for the flow in cfs, defaults to `flow_cfs`
#' @param .habitat_var Unquoted column name for the suitable area, defaults to `wua_per_lf`
#'
#' @returns A modified version of the original `data.frame` or `tbl_df` with additional flow and habitat columns identified by "_nbfc" ("no baseflow channel)
#' @md

#' @export
flow_threshold_remove_baseflow_apply <- function(data, .flow_threshold_var,
                                                 .flow_var = flow_cfs, .habitat_var = wua_per_lf) {
  data |>
    nest(fsa_raw = c({{.flow_var}}, {{.habitat_var}})) |>
    mutate(fsa_nbfc = pmap(list(fsa_raw, {{.flow_threshold_var}}),
                           function(x, q) flow_threshold_remove_baseflow(x, q, .flow_var = {{.flow_var}}, .habitat_var = {{.habitat_var}}))) |>
    unnest(c(fsa_raw, fsa_nbfc))

}

#' Split Flow-to-Suitable-Area Curve at Flow Threshold (Helper Function)
#'
#' @param suitable_area A numeric vector of suitable area values, e.g. in ft2 per lf
#' @param flow A numeric vector of flow values in cfs, same length as `suitable_area`
#' @param flow_threshold A numeric value for the instream-floodplain flow threshold in cfs
#'
#' @returns A named list including the original `suitable_area` and `flow` vectors, an identifier of instream versus floodplain (`above_threshold`), and suitable area values for instream and floodplain (`habitat_below` and `habitat_above`)
#' @md
#'
#' @export
flow_threshold_split <- function(suitable_area, flow, flow_threshold) {

  if (length(flow) > 0 &
      length(suitable_area) > 0 &
      all(!is.na(flow)) &
      all(!is.na(suitable_area)) &
      all(!is.na(flow_threshold))) {

    qt <- flow_threshold

    suitable_area_below_flow_threshold <- if (any(flow < qt)) suitable_area[which(flow < qt)]

    below_at_threshold <- if (any(flow < qt)) suitable_area_below_flow_threshold[[length(suitable_area_below_flow_threshold)]] else 0

    #suitable_area_below <- pmax(0, if_else(suitable_area <= below_at_threshold, suitable_area, below_at_threshold))

    #uitable_area_above <- pmax(0, suitable_area - suitable_area_below) #

    suitable_area_below <- if_else((suitable_area <= below_at_threshold) & (flow <= flow_threshold), suitable_area, below_at_threshold)

    suitable_area_above <- suitable_area - suitable_area_below #


    return(list(flow = flow,
                suitable_area = suitable_area,
                above_threshold = (flow >= qt),
                threshold = flow_threshold,
                habitat_below = suitable_area_below,
                habitat_above = suitable_area_above))

  } else {
    return(list(flow = flow,
                suitable_area = suitable_area,
                above_threshold = NA,
                threshold = flow_threshold,
                habitat_below = NA,
                habitat_above = NA))
  }
}

#' Split Flow-to-Suitable-Area Curve at Flow Threshold
#'
#' @param data A `data.frame` or `tbl_df` containing flow-to-suitable area curves (one row per flow per site, columns for flow and suitable area)
#' @param flow_threshold A numeric value for a flow threshold defining the instream-floodplain split
#' @param .flow_threshold_var, An unquoted column name for the flow threshold column in the data frame. Ignored if `flow_threshold` is provided.
#' @param .id_var Unquoted column name for the variable id, defaults to `comid`
#' @param .flow_var Unquoted column name for the flow in cfs, defaults to `flow_cfs`
#' @param .habitat_var Unquoted column name for the suitable area, defaults to `wua_per_lf`
#' @param .group_var Unquoted column name for an optional grouping variable
#'
#' @returns A new version of the `data.frame` or `tbl_df` containing additional columns for the instream and floodplain (below or above threshold) split
#' @md
#'
#' @export
flow_threshold_split_apply <- function(data,
                                       flow_threshold = NULL,
                                       .flow_threshold_var = NULL,
                                       .id_var = comid,
                                       .flow_var = flow_cfs,
                                       .habitat_var = wua_per_lf,
                                       .group_var = NULL) {

  if (!is.null(flow_threshold)) {

      summarized <-
        data |>
        group_by({{.id_var}}, {{.group_var}}) |>
        summarize(result = list(flow_threshold_split(suitable_area = {{.habitat_var}},
                                                     flow = {{.flow_var}},
                                                     flow_threshold = flow_threshold))) |>
        ungroup()

  } else {

    summarized <-
      data |>
      group_by({{.id_var}}, {{.group_var}}, {{.flow_threshold_var}}) |>
      summarize(result = list(flow_threshold_split(suitable_area = {{.habitat_var}},
                                                   flow = {{.flow_var}},
                                                   flow_threshold = {{.flow_threshold_var}}))) |>
      ungroup() |>
      select(-{{.flow_threshold_var}})

  }

  out <- summarized |>
           unnest_wider(result, names_sep="_") |>
           unnest(c(-{{.id_var}}, -{{.group_var}})) |>
           rename_with(function(x) str_replace(x, "result_", ""), starts_with("result_")) |>
           rename({{.flow_var}} := flow,
                  {{.habitat_var}} := suitable_area)

}

#' Rearing Duration HSI for Valley Lowland
#'
#' @param x A vector of values of number of days
#' @param baseflow_threshold_days An optional threshold number of days defining baseflow
#'
#' @returns A vector of values of duration HSI
#' @md
#'
#' @export
duration_hsi_rearing_vl <- function(x, baseflow_threshold_days = Inf) {
  case_when(x == 0  ~ 0,
            x < 18  ~ 0.66,
            x < 25 ~ 1.00,
            x < baseflow_threshold_days ~ 0.66,
            x >= baseflow_threshold_days ~ 0)
}

#' Rearing Duration HSI for Valley Foothill
#'
#' @param x A vector of values of number of days
#' @param baseflow_threshold_days An optional threshold number of days defining baseflow
#'
#' @returns A vector of values of duration HSI
#' @md
#'
#' @export
duration_hsi_rearing_vf <- function(x, baseflow_threshold_days = Inf) {
  case_when(x <= 0  ~ 0,
            x < 10 ~ 1.00,
            x < baseflow_threshold_days ~ 0.66,
            x >= baseflow_threshold_days ~ 0)
}

#' Spawning Duration HSI
#'
#' @param x A vector of values of number of days
#'
#' @returns A vector of values of duration HSI
#' @md
#'
#' @export
duration_hsi_spawning <- function(x) {
  case_when(x <= 0  ~ 0,
            x < 84 ~ 0,
            x >= 84 ~ 1.00)
}
