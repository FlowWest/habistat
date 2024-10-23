# Data not updating?
# (1) Check that the latest source hydraulic Rds outputs are present
# (2) Run model_cleaned.Rmd to re-export wua_hydraulic and wua_hydraulic_interp
# (3) Rebuild package and restart R session

function(input, output, session){

  # pre-load dynamic output items and make sure they do not suspend when hidden
  output$units_selector <- renderUI({})
  outputOptions(output, "units_selector", suspendWhenHidden = FALSE)
  output$out_spawning_toggle <- renderUI({})
  outputOptions(output, "out_spawning_toggle", suspendWhenHidden = FALSE)
  output$out_streamgage_selector <- renderUI({})
  outputOptions(output, "out_streamgage_selector", suspendWhenHidden = FALSE)
  output$out_flowscale_toggle <- renderUI({})
  outputOptions(output, "out_flowscale_toggle", suspendWhenHidden = FALSE)


  # WUA VARIABLE UNITS ---------------------------------------------------------
  wua_var <- reactive({
    switch(input$wua_var,
           "wua_per_lf_pred_SD" = paste0(input$wua_units,"_pred_SD"),
           "wua_per_lf_pred_SN" = paste0(input$wua_units,"_pred_SN"),
           "wua_per_lf_actual" = paste0(input$wua_units,"_actual"))
    })

  wua_lab <- reactive({
    switch(input$wua_units,
           "wua_per_lf" = "Suitable Habitat Area (ft2) per linear ft",
           "wua_acres" = "Suitable Habitat Area (acres)")
  })

  wua_suf <- reactive({
    switch(input$wua_units,
           "wua_per_lf" = "ft2 / ft",
           "wua_acres" = "ac")
  })

  output$units_selector <- renderUI({
    shinyWidgets::radioGroupButtons("wua_units", "Select Habitat Units",
                 choices=c("ft2 per linear ft" = "wua_per_lf",
                           "total acres" = "wua_acres"),
                 selected= if (input$flowline_scope=="comid") "wua_per_lf" else "wua_acres",
                 disabled = is.null(input$main_map_shape_click$id))
  })

  #input_parms <- reactiveValues(spawning_filter = FALSE)
  output$out_spawning_toggle <- renderUI({
    if (input$flowline_scope=="comid" & input$habitat_type=="spawning") {
      checkboxInput("spawning_toggle", "Use geomorphic spawning gravel likelihood", value = TRUE)
    }
  })
  #observe({
  #  input_parms$spawning_filter <- input$spawning_toggle
  #})

  # IDENTIFIERS FOR CLICKED ELEMENT --------------------------------------------

  most_recent_map_click <- reactiveValues(type = "none",
                                          lng = NULL,
                                          lat = NULL)

  clicked_item_label <- reactive({
    input$main_map_shape_click # so that this item updates on map click
    if (most_recent_map_click$type == "comid") {
      comid_name <- attr$gnis_name[[which(attr$comid==selected_point$comid)]]
      comid_watershed_name <- attr$watershed_level_3[[which(attr$comid==selected_point$comid)]]
      comid_name_suffix <- if_else(!is.na(comid_name),
                                   glue::glue(" ({comid_name})"),
                                   glue::glue(" (tributary of {comid_watershed_name})"))
      paste0("ComID ", selected_point$comid, comid_name_suffix)
    } else if (most_recent_map_click$type == "mainstem") {
      paste0(selected_mainstem$river_name, " Mainstem")
    } else if (most_recent_map_click$type == "watershed") {
      paste0(selected_watershed$watershed_name, " Watershed")
    }
  })

  output$clicked_item_heading <- renderUI({
    if (length(clicked_item_label()) > 0) {
      h4(clicked_item_label())
    } else {
      h4(paste0("no ",input$flowline_scope," selected"))
    }
  })

  # REACTIVE FLOW FILTER (COMID) -----------------------------------------------

  message(paste(names(predictions), collapse = ", "))

  active_predictions <- reactive({

    active_predictions <- predictions |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
      filter(habitat == input$habitat_type) |>
      inner_join(attr |> select(comid, spawning_gravel_either_method), by=join_by(comid))

    if(isTRUE(input$spawning_toggle)) {
      active_predictions <- active_predictions |>
        mutate(across(starts_with("wua_") & contains("_pred"), function(x) if_else((habitat == "spawning") & (!spawning_gravel_either_method), 0, x)))
    }

    return(active_predictions)
  })

  active_geom <- reactive({
    geom |>
      inner_join(active_predictions(),
                 by = join_by(comid)) |>
      mutate(wua_per_lf = !!sym(input$wua_var)) |>
      filter(!is.na(wua_per_lf)) |>
      mutate(object_label = paste0("<strong>", comid, " ", gnis_name, "</strong>", "<br />",
                                   "Scale-Dependent Model: ",
                                   round(wua_per_lf_pred_SD,2), " ft2/ft", "<br />",
                                   #"Scale-Dependent Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SD_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model: ",
                                   round(wua_per_lf_pred_SN,2), " ft2/ft", "<br />",
                                   #"Scale-Normalized Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SN_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Actual: ", round(wua_per_lf_actual,2), " ft2/ft"))
  })

  active_geom_mainstem <- reactive({
    mainstems |>
      inner_join(active_predictions_mainstem(),
                 by = join_by(river_cvpia)) |>
      mutate(wua_per_lf = !!sym(input$wua_var)) |>
      filter(!is.na(wua_per_lf)) |>
      mutate(object_label = paste0("<strong>", river_cvpia, "</strong>", "<br />",
                                   "Scale-Dependent Model: ",
                                   round(wua_per_lf_pred_SD,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model: ",
                                   round(wua_per_lf_pred_SN,2), " ft2/ft", "<br />"))
  })

  # SELECTED COMID OBSERVER ----------------------------------------------------

  selected_point <- reactiveValues(object_id = NULL,
                                   comid = NULL,
                                   lng = NULL,
                                   lat = NULL)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 6) == "comid_") {
        most_recent_map_click$type <- "comid"
        selected_point$object_id <- input$main_map_shape_click$id
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
        selected_point$comid <- as.numeric(str_replace(selected_point$object_id, "comid_", ""))
      }
    }
  })

  # PREDICTION TABLES, ATTRIBUTES, PLOTS ---------------------------------------

  active_pred_table <- reactive({
    if (most_recent_map_click$type == "comid") {
      active_predictions() |>
        select(any_of(names(var_names))) |>
        filter(comid == selected_point$comid) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        mutate(across(where(is.logical), function(x) if_else(x, "True", "False"))) |>
        mutate(across(-where(is.character), as.character)) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    } else if (most_recent_map_click$type == "watershed") {
      active_predictions_watershed() |>
        select(any_of(names(var_names))) |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        mutate(across(where(is.logical), function(x) if_else(x, "True", "False"))) |>
        mutate(across(-where(is.character), as.character)) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    } else if (most_recent_map_click$type == "mainstem") {
      active_predictions_mainstem() |>
        select(any_of(names(var_names))) |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        mutate(across(where(is.logical), function(x) if_else(x, "True", "False"))) |>
        mutate(across(-where(is.character), as.character)) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    }
  })

  output$pred_table <- DT::renderDT({
    if (most_recent_map_click$type != "none") {
      DT::datatable(active_pred_table(),
                    options = list(paging = F, searching = F),
                    selection = "none")
    }
  })

  active_attr_table <- reactive({ # for COMID only
   #if (!is.null(selected_point$object_id) & (length(selected_point$comid)>0)) {
     # if (!is.null(selected_point$object_id)) {
     if (most_recent_map_click$type == "comid") {
     result <- attr |>
       filter(comid == selected_point$comid) |>
       mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
       mutate(across(where(is.logical), function(x) if_else(x, "True", "False"))) |>
       mutate(across(-where(is.character), as.character)) |>
       pivot_longer(everything()) |>
       mutate(name = var_names[name])
     #gc()
     return(result)
     }
  })

  output$attr_table <- DT::renderDT({
     if (most_recent_map_click$type == "comid") {
        DT::datatable(active_attr_table(),
                      options = list(paging = F, searching = F),
                      selection = "none")
      }
    })

  output$fsa_plot <- renderPlot({
    palette_linetypes <- c("Unscaled" = "solid",
                           "Duration Scaled" = "dashed")
    palette_colors <- c("Scale-Dependent" = "#6388b4",
                        "Scale-Normalized" = "#8cc2ca",
                        "Actual" = "#ffae34")
    if (most_recent_map_click$type == "comid") { #& (length(selected_point$comid)>0)) {
    predictions |>
      filter(comid == selected_point$comid) |>
      filter(habitat == input$habitat_type) |>
      ggplot(aes(x = flow_cfs)) + #|> add_color_scale(type = input$habitat_type) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_actual")), color="Actual", linetype = "Unscaled")) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype = "Unscaled")) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype = "Unscaled")) +
      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
        input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
        input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
        input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
      #geom_line(aes(y = !!sym(var_not_selected), linetype = "Unscaled"), color = "gray") +
      #geom_line(aes(y = !!sym(input$wua_var), linetype = "Unscaled"), linewidth=1) +
      #geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = durwua)) +
      #geom_hline(aes(yintercept = chan_width_ft)) + #, linetype="Channel Width (ft)")) +
      #geom_vline(aes(xintercept = baseflow_cfs)) +
      #geom_text(aes(x = 1, y = chan_width_ft, label = chan_width_ft)) +
      scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
      scale_y_continuous(limits = c(0, NA)) +
      #scale_y_continuous(trans = ihs, labels = scales::label_comma(), limits = c(0, NA)) +
      theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
      xlab("Flow (cfs)") + ylab(wua_lab()) +
      scale_color_manual(name = "Model Type",
                         values = palette_colors) +
      scale_linetype_manual(name = paste("Duration Analysis", coalesce(paste0("(",str_to_upper(selected_gage()), ")"), "")),
                            values = palette_linetypes)
    } else if (most_recent_map_click$type == "watershed") {
      predictions_watershed |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        #geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Analysis")) +
        geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
          input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
          input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
          input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab(wua_lab()) +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
        scale_linetype_manual(name = paste("Duration Analysis", coalesce(paste0("(",str_to_upper(selected_gage()), ")"), "")),
                              values = palette_linetypes)
    } else if (most_recent_map_click$type == "mainstem") {
      predictions_mainstem |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
       # geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
          input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
          input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
          input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab(wua_lab()) +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
        scale_linetype_manual(name = paste("Duration Analysis", coalesce(paste0("(",str_to_upper(selected_gage()), ")"), "")),
                              values = palette_linetypes)
    } else {
      ggplot()
    }})

  # LEAFLET MAP FUNCTIONS ------------------------------------------------------

  layer_flowlines <- function(m, show = TRUE, type = "comid") {

    # first remove any existing flowlines
    m |> leaflet::removeShape(geom$object_id)
    m |> leaflet::removeShape(mainstems$object_id)

    if(type == "comid") {
    if(show) {
      # pal_limits <- c(min(active_geom()$wua_per_lf), max(active_geom()$wua_per_lf))
      message("plotting ", nrow(active_geom()), " flowlines")
      #m |> leafgl::addGlPolylines(data = active_geom(),layerId = ~object_id)
      m |> leaflet::addPolylines(data = active_geom(), # TODO switch to leafgl::addGlPolylines
                                 layerId = ~object_id,
                                 label = ~lapply(object_label, htmltools::HTML),
                                 color = ~pal(wua_per_lf, type=input$habitat_type), #pal(wua_per_lf),
                                 opacity = 1,
                                 weight = 2,
                                 options = leaflet::pathOptions(pane = "Flowlines"),
                                 group = "flowlines",
                                 highlightOptions = leaflet::highlightOptions(color = "#8B0000",
                                                                              weight = 3,
                                                                              bringToFront = TRUE)
      ) |>
        leaflet::addLegend(position = "bottomright",
                           colors = rev(flow_scale_colors[[input$habitat_type]]),
#                           colors = rev(pal(seq(pal_limits[[1]], pal_limits[[2]], (pal_limits[[2]]-pal_limits[[1]])/5))),
                           labels = lapply(paste("&ge;", rev(flow_scale_breaks[[input$habitat_type]])), htmltools::HTML),
#                           labels = c(round(pal_limits[[2]],2), rep("", 5-1), round(pal_limits[[1]],2)),
                           title = "Suitable Habitat Area (ft2) per linear ft",
                           layerId = "clegend")
      } else {
        m |> leaflet::removeShape(active_geom()$object_id) |> leaflet::removeControl("clegend")
        # TODO switch to leaflet::removeGlPolylines
      }
    } else if (type == "mainstem") {
    if(show) {
      pal_limits <- c(min(active_geom_mainstem()$wua_per_lf), max(active_geom_mainstem()$wua_per_lf))
      message("plotting ", nrow(active_geom_mainstem()), " flowlines")
      m |> leaflet::addPolylines(data = active_geom_mainstem(),
                                 layerId = ~mainstem_id,
                                 label = ~lapply(object_label, htmltools::HTML),
                                 color = ~pal(wua_per_lf, type=input$habitat_type), #pal(wua_per_lf),
                                 opacity = 1,
                                 weight = 2,
                                 options = leaflet::pathOptions(pane = "Flowlines"),
                                 group = "flowlines",
                                 highlightOptions = leaflet::highlightOptions(color = "#8B0000",
                                                                              weight = 3,
                                                                              bringToFront = TRUE)
      ) |>
        leaflet::addLegend(position = "bottomright",
                           colors = rev(flow_scale_colors[[input$habitat_type]]),
                           #                           colors = rev(pal(seq(pal_limits[[1]], pal_limits[[2]], (pal_limits[[2]]-pal_limits[[1]])/5))),
                           labels = lapply(paste("&ge;", rev(flow_scale_breaks[[input$habitat_type]])), htmltools::HTML),
                           #                           labels = c(round(pal_limits[[2]],2), rep("", 5-1), round(pal_limits[[1]],2)),
                           title = "Suitable Habitat Area (ft2) per linear ft",
                           layerId = "clegend")
    } else {
        m |> leaflet::removeShape(active_geom_mainstem()$object_id) |> leaflet::removeControl("clegend")
      }
    }
  }

  # LEAFLET RENDER -------------------------------------------------------------

  output$main_map <- renderLeaflet({

    shinyjs::showElement(id = 'loading_action')

    bbox <- c(xmin=-122.3, ymin=38.5, xmax=-121.3, ymax=39.7)

    leaflet::leaflet() |>
      leaflet::addMapPane("Basemap", zIndex = 400) |>
      leaflet::addMapPane("ValleyLowland", zIndex = 440) |>
      leaflet::addMapPane("Watersheds", zIndex = 445) |>
      leaflet::addMapPane("Flowlines", zIndex = 470) |>
      leaflet::addMapPane("Overlays", zIndex = 475) |>
      leaflet::addMapPane("AOI", zIndex = 480) |>
      leaflet::addMapPane("Reference", zIndex = 490) |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 1,
                                                       maxNativeZoom = 13,
                                                       maxZoom = 13,
                                                       pane = "Basemap"
                        ), group = "Terrain (default)") |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Reference_Overlay/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 1,
                                                       #minZoom = 10,
                                                       pane = "Reference"
                        ), group = "Terrain (default)") |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 1,
                                                       maxZoom = 13,
                                                       maxNativeZoom = 23,
                                                       pane = "Basemap"
                        ), group = "Aerial Imagery") |>
      leaflet::addTiles(urlTemplate = "", attribution = 'Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, IGN, IGP, UPR-EGP, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, Esri, and the GIS User Community') |>
      leaflet::fitBounds(lng1 = bbox[["xmin"]],
                         lat1 = bbox[["ymin"]],
                         lng2 = bbox[["xmax"]],
                         lat2 = bbox[["ymax"]]) |>
      leaflet::addPolygons(data = watersheds,
                           stroke = T,
                           weight = 1,
                           color = "#31a1b3",
                           opacity = 0.5,
                           fill = T,
                           fillColor = "#FFFFFF",
                           fillOpacity = 0,
                           layerId = ~watershed_id,
                           group = "watersheds",
                           label = ~lapply(watershed_label, htmltools::HTML),
                           highlightOptions = highlightOptions(stroke = T,
                                                               weight = 1,
                                                               color = "white",
                                                               opacity = 1,
                                                               fill = T,
                                                               fillColor = "#31a1b3",
                                                               fillOpacity = 0.5,
                                                               bringToFront = F)
      ) |>
      leaflet::addPolygons(data = hqt, group = "HQT - Valley Lowland",
                           popup = "Habitat Quantification Tool Boundary",
                           color = "darkgrey",
                           weight = 0,
                           fillColor = "grey",
                           fillOpacity = 0.33,
                           options = leaflet::pathOptions(pane = "ValleyLowland")) |>
      addLayersControl(
        baseGroups = c("Terrain (default)", "Aerial Imagery"),
        overlayGroups = c("flowlines", "watersheds", "streamgages", "HQT - Valley Lowland"),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      hideGroup("HQT - Valley Lowland") |>
      hideGroup("Aerial Imagery")
  })

  observe({
    proxy <- leaflet::leafletProxy("main_map", session = session)
    proxy |>
      layer_flowlines(type = input$flowline_scope)
  })

# SELECTED WATERSHED OBSERVER ----------------------------------------------------

selected_watershed <- reactiveValues(object_id = NA,
                                     watershed_name = NA,
                                     lng = NA,
                                     lat = NA)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 10) == "watershed_") {
        most_recent_map_click$type <- "watershed"
        selected_watershed$object_id <- input$main_map_shape_click$id
        selected_watershed$watershed_name <-watersheds$watershed_level_3[[which(watersheds$watershed_id==input$main_map_shape_click$id)]]
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
      }
    }
  })

  selected_watershed_geom <- reactive({
    watersheds |>
      filter(watershed_id == selected_watershed$object_id) |>
      mutate(object_id = "active_watershed")
  })

  selected_watershed_bbox <- reactive({
    st_bbox(selected_watershed_geom())
  })

  selected_watershed_name <- reactive({
    if(!is.na(selected_watershed$object_id)) {
      watersheds$watershed_level_3[[which(watersheds$watershed_id==selected_watershed$object_id)]]
    } else {
      NA
    }
  })

  # SELECTED MAINSTEM OBSERVER ----------------------------------------------------

  selected_mainstem <- reactiveValues(object_id = NA,
                                       river_name = NA,
                                       lng = NA,
                                       lat = NA)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 9) == "mainstem_") {
        most_recent_map_click$type <- "mainstem"
        selected_mainstem$object_id <- input$main_map_shape_click$id
        selected_mainstem$river_name <-mainstems$river_cvpia[[which(mainstems$mainstem_id==input$main_map_shape_click$id)]]
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
      }
    }
  })

  active_mainstem_geom <- reactive({
    watersheds |>
      filter(mainstem_id == selected_mainstem$object_id) |>
      mutate(object_id = "active_mainstem")
  })

  # ACTIVE MAP LOGIC ------------------------------------------------

  observe({

    proxy <- leaflet::leafletProxy("main_map", session = session)

    # Active watershed

    if (most_recent_map_click$type == "watershed") {
      # (nrow(selected_watershed_geom()) > 0) &
      proxy |>
        leaflet::removeShape("active_watershed") |>
        leaflet::addPolygons(data = selected_watershed_geom(),
                             stroke = T,
                             weight = 6,
                             color = "#f6b911",
                             opacity = 1,
                             fill = T,
                             fillColor = "#f6b911",
                             fillOpacity = 0.33,
                             group = "watersheds", # this will make the active watershed also show/hide with the layer toggle
                             layerId = "active_watershed",
                             label = ~lapply(watershed_label, htmltools::HTML)) |>
        leaflet::fitBounds(lng1 = selected_watershed_bbox()$xmin,
                           lat1 = selected_watershed_bbox()$ymin,
                           lng2 = selected_watershed_bbox()$xmax,
                           lat2 = selected_watershed_bbox()$ymax)

    } else {

      proxy |>
        leaflet::removeShape("active_watershed")

    }

    # Active flowline

    if (most_recent_map_click$type %in% c("comid", "mainstem")) {

      if (most_recent_map_click$type == "comid") {
        clicked_polyline <- active_geom() |>
          filter(comid == selected_point$comid)
      } else if (most_recent_map_click$type == "mainstem") {
        clicked_polyline <- active_geom_mainstem() |>
          filter(river_cvpia == selected_mainstem$river_name)
      }

      #bbox <- sf::st_bbox(clicked_polyline)

      proxy |>
        leaflet::removeShape("active_flowline") |>
        leaflet::addPolylines(data = clicked_polyline,
                              weight = 20,
                              color = "#f6b911",
                              opacity = 1,
                              layerId = "active_flowline",
                              label = ~lapply(object_label, htmltools::HTML)) #|>
      #fitBounds(bbox$ymin, bbox$xmin, bbox$ymax, bbox$xmax)

    } else {

      proxy |>
        leaflet::removeShape("active_flowline")

    }

  })

  # REACTIVE FLOW FILTER (WATERSHED AND MAINSTEM) ------------------------------

  active_flow_applied <- eventReactive(input$activeFlowApplyButton, {
    input$active_flow
  })

  active_map_params <- reactiveValues(flow = 1000) # make sure this matches the value specified in the "selected" attribute on sliderTextInput in ui.R
  observe({
    active_map_params$flow <- active_flow_applied()
  })

  active_predictions_watershed <- reactive({
    predictions_watershed |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
      filter(habitat == input$habitat_type)
  })

  active_predictions_mainstem <- reactive({
    predictions_mainstem |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
      filter(habitat == input$habitat_type)
  })

  # STREAMGAGE SELECTION -------------------------------------------------------

  # ACTIVE STREAMGAGE SELECTOR

  active_reach_info <- reactiveValues(is_mainstem = FALSE,
                                      watershed_name = NA,
                                      river_name = NA,
                                      multiplier = NA)

  # Update these attributes any time the map is clicked...
  observeEvent(input$main_map_shape_click, {
    message("--- Observe input$main_map_shape_click: Update active_reach_info values")
    if (most_recent_map_click$type == "comid") {
      active_reach_info$watershed_name <- (attr |>
                                             filter(comid == selected_point$comid) |>
                                             pull(watershed_level_3))[[1]]
      active_reach_info$river_name <- (attr |>
                                         filter(comid == selected_point$comid) |>
                                         pull(river_cvpia))[[1]]
      active_reach_info$is_mainstem <- !is.na((attr |>
                                                 filter(comid == selected_point$comid) |>
                                                 pull(river_cvpia))[[1]])
    } else if (most_recent_map_click$type == "watershed"){
      active_reach_info$watershed_name <- selected_watershed$watershed_name
      active_reach_info$is_mainstem <- FALSE
    } else if (most_recent_map_click$type == "mainstem"){
      active_reach_info$river_name <- selected_mainstem$river_name
      active_reach_info$is_mainstem <- TRUE
    }
    message("/// Observe input$main_map_shape_click: Update active_reach_info values")
  })

  # Whenever active river name or watershed name is updated...
  # ...Update list of possible streamgages...
  streamgage_options <- reactive({
    message("--- Reactive: Update streamgage_options")

    if(active_reach_info$is_mainstem) {
      message(paste("pulling streamgage options for", active_reach_info$river_name))
      streamgage_options <- streamgages[[active_reach_info$river_name]]
    } else if (most_recent_map_click$type == "watershed") {
      message(paste("pulling streamgage options for", active_reach_info$watershed_name))
      streamgage_options <- streamgages_by_watershed[[active_reach_info$watershed_name]]
    } else {
      streamgage_options <- list()
    }
    message(paste(length(streamgage_options), "streamgages found"))
    return(streamgage_options)
    message("/// Reactive: Update streamgage_options")
  })
  # ...Pull geometries for possible streamgages...
  streamgage_options_geom <- reactive({
    message("--- Reactive: Update streamgage_options_geom")
    if(isTRUE(length(streamgage_options()) > 0)){
      geom <- streamgage_pts |>
        filter(station_id %in% names(streamgage_options()))
    } else {
      geom <- st_sf(st_sfc())
    }
    message(paste(nrow(geom), "streamgage geometries filtered"))
    return(geom)
    message("/// Reactive: Update streamgage_options_geom")
  })
  # Pull out the ID of the nearest streamgage
  streamgage_nearest_id <- reactive({
    message("--- Reactive: Update streamgage_nearest_id")
    if(isTRUE(nrow(streamgage_options_geom()) > 0)) {
      geom <- streamgage_options_geom()
      streamgage_nearest_id <- geom$station_id[[
        st_nearest_feature(st_point(c(most_recent_map_click$lng,
                                      most_recent_map_click$lat)),
                           geom,
                           check_crs = FALSE)]]
    } else {
      streamgage_nearest_id <- NA
    }
    message(paste("streamgage_nearest_id =", streamgage_nearest_id))
    return(streamgage_nearest_id)
    message("/// Reactive: Update streamgage_nearest_id")
  })

  most_recent_gage_event <- reactiveValues(type = NA)

  # Identify the selected streamgage based on the selector menu and the available options
  selected_gage <- reactive({
    message("--- Reactive: Update selected_gage")
    if ((isTRUE(coalesce(input$streamgage_id, NA) %in% names(streamgage_options()))) & (most_recent_gage_event$type=="menu_select")) {
      selected_gage <- input$streamgage_id
    } else {
      selected_gage <- streamgage_nearest_id()
    }
    message(paste("selected_gage =", selected_gage))
    return(coalesce(selected_gage, NA))
    message("/// Reactive: Update selected_gage")
  })
  # Label the geometry with the selected streamgage
  streamgage_options_geom_labelled <- reactive({
    message("--- Reactive: Update streamgage_options_geom_labelled")
    geom <- streamgage_options_geom()
    message(paste(nrow(geom), "streamgage geometries retrieved"))
    if(isTRUE((nrow(geom) > 0))) {
      if(isTRUE(!is.na(selected_gage()))) {
        geom_lab <- geom |> mutate(selected = (station_id == coalesce(selected_gage(), NA)))
      } else {
        geom_lab <- geom |> mutate(selected = FALSE)
      }
    } else {
      geom_lab <- st_sf(st_sfc())
    }
    message(paste(nrow(geom_lab), "streamgage geometries returned"))
    return(geom_lab)
    message("/// Reactive: Update streamgage_options_geom_labelled")
  })
  # Update the map to include the latest streamgages including the active gage selection
  observe({
    message("--- Observe streamgage_options_geom_labelled: Update streamgage leaflet map layer")
    selected_gage() # This needs to be here for reactive to function
    proxy <- leaflet::leafletProxy("main_map")
    proxy |>
      leaflet::removeMarker(paste0("streamgage_", streamgage_pts$station_id))
    if (isTRUE(nrow(streamgage_options_geom_labelled()) > 0)) {
      proxy |>
        leaflet::addCircleMarkers(data = streamgage_options_geom_labelled(),
                                  layerId = ~paste0("streamgage_", station_id),
                                  group = "streamgages",
                                  label = ~station_label,
                                  color = ~if_else(selected, "#00A2E8", "#000000"),
                                  radius = 6,
                                  options = leaflet::markerOptions(pane = "Overlays"))
    }
    message("/// Observe streamgage_options_geom_labelled: Update streamgage leaflet map layer")
  })

  # On click of streamgage, update the lat/lon of most recent map click
  # which is used to update the active streamgage
  observeEvent(input$main_map_marker_click, {
    message("--- Observe input$main_map_marker_click: Update most_recent_map_click values")
    cat(input$main_map_marker_click$id)
    if (!is.null(input$main_map_marker_click$id)) {
      if(substr(input$main_map_marker_click$id, 1, 11) == "streamgage_") {
        message("clicked streamgage")
        most_recent_map_click$lng <- input$main_map_marker_click$lng
        most_recent_map_click$lat <- input$main_map_marker_click$lat
        most_recent_gage_event$type <- "map_click"
      }
    }
    message("/// Observe input$main_map_marker_click: Update most_recent_map_click values")
  })
  observeEvent(input$streamgage_id, {
    most_recent_gage_event$type <- "menu_select"
  })

  output$out_streamgage_selector <- renderUI({
    message("--- Render output$streamgage_selector")
    if(length(streamgage_options()) > 0) {
      message(paste(streamgage_options(), collapse=", "))
      el <- selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = setNames(names(streamgage_options()), streamgage_options()),
                  selected = selected_gage())
    } else {
      el <- selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = c())
      message("no stream gages")
    }
    return(el)
    message("/// Render output$streamgage_selector")
  })

  output$out_flowscale_toggle <- renderUI({
    message("--- Render output$flowscale_toggle")
    if (most_recent_map_click$type == "comid" & isTRUE(is.numeric(active_reach_info$multiplier))) {
      title <- paste0("Scale Flow by Drainage Area and Precipitation Ratio = ", round(active_reach_info$multiplier, 2))
      el <- checkboxInput("scale_flow", title, value=T)
      return(el)
    }
    message("/// Render output$flowscale_toggle")
  })

  # DURATION ANALYSIS ----------------------------------------------------------

  # STREAMGAGE RATING CURVE

  # for comid only
  streamgage_drc <- reactive({

    if (most_recent_map_click$type == "comid") {

      active_reach_attr <- attr |>
        filter(comid == selected_point$comid) |>
        as.list()

      active_streamgage_attr <-
        streamgage_attr |>
        filter(station_id == coalesce(selected_gage(), NA)) |>
        as.list()

      message(paste0("pulling streamgage_drc for ", selected_gage(),
                     " ", input$selected_run,
                     " ", input$habitat_type,
                     " ", input$selected_wyt))

      active_streamgage_data <-
        get_data(streamgage_duration_rating_curves, package = "habistat") |>
        filter((station_id == selected_gage()) &
                 (run == input$selected_run) &
                 (habitat == input$habitat_type) &
                 (wy_group == input$selected_wyt)) |>
        unnest(data)


      if (nrow(active_streamgage_data) > 0) {

        active_reach_info$multiplier <-
          (active_reach_attr$da_area_sq_km / active_streamgage_attr$da_gage) *
          (active_reach_attr$da_ppt_mean_mm / active_streamgage_attr$pc_gage)

        multiplier <- if (isTRUE(input$scale_flow)) active_reach_info$multiplier else 1

        message(paste("scaling streamgage flow data by", multiplier))

        active_streamgage_data |>
          mutate(model_q = model_q * multiplier) |>
          mutate(dhsi_selected = case_when(
            input$habitat_type == "spawning" ~ durhsi_spawning,
            active_reach_attr$hqt_gradient_class == "Valley Lowland" ~ durhsi_rearing_vl,
            TRUE ~ durhsi_rearing_vf))

      } else {

        tibble(q = list(), dhsi_selected = list(), avg_max_days_inundated = list())

      }
    }

  })

  # DURATION CURVE CALCULATION

  duration_curve <- reactive({

    if (length(selected_gage()) > 0) {

    if (most_recent_map_click$type == "comid") {

      fsa <-
        predictions |>
        filter(comid == selected_point$comid) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(wua_var()))

        if (nrow(streamgage_drc()) > 0) {

          message("duration hsi curve")

          duration_curve_result <-
            habistat::duration_apply_dhsi_to_fsa_curve(
            fsa = fsa,
            fsa_q = flow_cfs,
            fsa_wua = selected_wua,
            drc = streamgage_drc(),
            drc_q = model_q,
            drc_dhsi = dhsi_selected)

        } else {

          duration_curve_result <- tibble(q = list(), durwua = list())

        }

    } else if (most_recent_map_click$type %in% c("mainstem", "watershed")) {

      group_var <- switch(most_recent_map_click$type,
                          "mainstem" = "river_cvpia",
                          "watershed" = "watershed_level_3")

      predictions_filtered <- switch(most_recent_map_click$type,
                               "mainstem" = predictions |> filter(river_cvpia == selected_mainstem$river_name),
                               "watershed" = predictions |> filter(watershed_level_3 == selected_watershed$watershed_name))

      flow_xw <- switch(most_recent_map_click$type,
                        "mainstem" = cv_mainstems_flow_xw,
                        "watershed" = cv_watersheds_flow_xw)

      message("FSA")
      fsas <-
        predictions_filtered |>
        filter(habitat == input$habitat_type) |>
        transmute(!!sym(group_var), comid, flow_idx, flow_cfs, reach_length_ft, selected_wua = !!sym(input$wua_var)) |>
        nest(fsa_raw = c(flow_idx, flow_cfs, reach_length_ft, selected_wua), .by = c(!!sym(group_var), comid)) |>
        inner_join(flow_xw, by=join_by(!!sym(group_var), comid)) |>
        mutate(fsa_scaled = pmap(list(fsa_raw, multiplier), function(x, y) scale_fsa(x, y, .wua_var = selected_wua, .flow_var = flow_cfs))) |>
        rename(multiplier_fsa = multiplier) |>
        glimpse()

      message(paste("DRC", selected_gage()))
      drcs <-
        get_data(streamgage_duration_rating_curves, package = "habistat") |>
        filter((station_id == coalesce(selected_gage(), NA)) &
                 (run == input$selected_run) &
                 (habitat == input$habitat_type) &
                 (wy_group == input$selected_wyt)) |>
        rename(drc_raw = data) |>
        inner_join(streamgage_attr |> select(station_id, da_gage, pc_gage), by=join_by(station_id)) |>
        expand_grid(comid = fsas$comid) |>
        inner_join(flow_xw |>
                     select(comid, !!sym(group_var), da_reach, pc_reach),
                   by=join_by(comid)) |>
        mutate(multiplier = (da_reach / da_gage) * (pc_reach / pc_gage)) |>
        inner_join(attr |> select(comid, hqt_gradient_class), by=join_by(comid)) |>
        mutate(drc_scaled =
                 pmap(list(drc_raw, multiplier, hqt_gradient_class),
                      function(drc_raw, multiplier, hqt_gradient_class) {
                        drc_raw |>
                          mutate(model_q = model_q * multiplier) |>
                          mutate(dhsi_selected = case_when(input$habitat_type == "spawning" ~ durhsi_spawning,
                                                           hqt_gradient_class == "Valley Lowland" ~ durhsi_rearing_vl,
                                                           TRUE ~ durhsi_rearing_vf))
                      })) |>
        rename(multiplier_drc = multiplier) |>
        glimpse()

      if (nrow(drcs)>0) {

      message("FSA x DRC")
      joined <- inner_join(fsas, drcs, by=join_by(comid, !!sym(group_var), da_reach, pc_reach)) |>
        mutate(result = pmap(list(fsa_scaled, drc_scaled),
                             function(fsa, drc) {

          duration_applied <-
            habistat::duration_apply_dhsi_to_fsa_curve(
              fsa = fsa,
              fsa_q = flow_cfs,
              fsa_wua = selected_wua,
              drc = drc,
              drc_q = model_q,
              drc_dhsi = dhsi_selected)

          fsa |>
            left_join(duration_applied |> select(flow_cfs = q, wua, durwua), by=join_by(flow_cfs)) |>
            mutate(across(c(wua, durwua), function(y) zoo::na.approx(y, x = flow_cfs, na.rm=F, rule=2)))

        })) |>
        select(-fsa_raw, -drc_raw, -fsa_scaled, -drc_scaled) |>
        unnest(result) |>
        glimpse()

      message(paste("aggregated duration curve for", most_recent_map_click$type))

      duration_curve_result <-
        joined |>
        glimpse() |>
        group_by(habitat, run, wy_group,
                 !!sym(group_var), q = flow_idx) |> # better to use the float version
        summarize(across(c(wua, durwua), switch(input$wua_units,
                                                wua_per_lf = function(y) sum(y * reach_length_ft) / sum(reach_length_ft),
                                                wua_acres = function(y) sum(y * reach_length_ft) / 43560)),
                  .groups="drop")
      } else {

      duration_curve_result <- tibble(q = list(), wua = list(), durwua = list())

      }

    } else {

      fsa <- tibble(flow_cfs = c(),
                    flow_idx = c(),
                    selected_wua = c())
    }
    } else {

      duration_curve_result <- tibble(q = list(), wua = list(), durwua = list())

    }

  return(duration_curve_result)

  })

  output$dur_plot <- renderPlot({
    if(most_recent_map_click$type == "comid") {
      if((nrow(streamgage_drc())>0) & (nrow(duration_curve()) > 0)) {
      bind_rows("days" = streamgage_drc() |> transmute(x = model_q, y = avg_max_days_inundated),
                "dhsi" = streamgage_drc() |> transmute(x = model_q, y = dhsi_selected),
                "wua" = bind_rows("wua" = duration_curve() |> transmute(x = q, y = wua),
                                  "durwua" = duration_curve() |> transmute(x = q, y = durwua),
                                  .id = "grp"),
                .id = "varname") |>
        bind_rows(tribble(~varname, ~grp, ~x, ~y,
                          "dhsi", "ref", -Inf, 1,
                          "dhsi", "ref", Inf,  1)) |>
        mutate(varname = varname |> factor(levels = c("days", "dhsi", "wua"),
                                           labels = c("Max Length of Period Exceeding Flow per WY-Season (days)",
                                                      "Duration Suitability Factor",
                                                      paste0("Suitable Habitat Curve (", wua_suf(),")"))),
               grp = grp |> coalesce("None") |>
                            factor(levels = c("None", "ref", "wua", "durwua"),
                                   labels = c("None", "ref", "Original", "Duration-Weighted"))) |>
        ggplot() +
        geom_step(aes(x = x, y = y, linetype = grp)) +
        facet_wrap(~varname, ncol = 1, scales="free_y") +
        scale_x_log10(breaks = scales::breaks_log(8), labels = scales::label_comma()) +
        scale_y_continuous(breaks = scales::breaks_extended(8), labels = scales::label_comma(), limits=c(0, NA)) +
        annotation_logticks(sides = "b") +
        ylab("") + xlab("Flow (cfs)") +
        theme(legend.position = "none",
              panel.grid.minor = element_blank()) +
        scale_linetype_manual(name = "",
                              values = c("None" = "solid",
                                         "ref" = "dashed",
                                         "Original" = "solid",
                                         "Duration-Weighted" = "dashed"))
      } else {
        ggplot()
      }
    } else if (most_recent_map_click$type %in% c("mainstem", "watershed")) {
      duration_curve() |>
        ggplot() +
        geom_line(aes(x = q, y = wua, linetype="Original")) +
        geom_line(aes(x = q, y = durwua, linetype="Duration-Weighted")) +
        ylab(paste0("Suitable Habitat Area (", wua_suf(),")")) +
        xlab("Flow (cfs) at Outlet") +
        scale_x_log10(breaks = scales::breaks_log(8), labels = scales::label_comma()) +
        scale_y_continuous(breaks = scales::breaks_extended(8), labels = scales::label_comma(), limits=c(0, NA)) +
        annotation_logticks(sides = "b") +
        theme(legend.position = "top",
              panel.grid.minor = element_blank()) +
        scale_linetype_manual(name = "",
                              values = c("Original" = "solid",
                                         "Duration-Weighted" = "dashed")) +
        labs(caption = "Aggregated from duration-scaled habitat curves for outlet comid (at nominal cfs) and other comids (at cfs downscaled by drainage area and precipitation ratio).")
    }

#    plt_days <-
#      ggplot() +
#      geom_step(data=streamgage_drc(), aes(x = model_q, y = avg_max_days_inundated)) +
#      xlab("Flow (cfs)") +
#      ylab("Days per WY") + scale_y_continuous(limits = c(0, 365)) +
#      ggtitle("Max Length of Period Exceeding Flow per WY-Season")
#
#    plt_dhsi <-
#      ggplot() +
#      geom_step(data=streamgage_drc(), aes(x = model_q, y = dhsi_selected)) +
#      xlab("Flow (cfs)") +
#      ylab("Weight") + scale_y_continuous(limits = c(0, 1)) +
#      ggtitle("Duration Suitability Factor")
#
#    plt_dwua <- ggplot() +
#      geom_line(data=duration_curve(), aes(x = q, y = wua, linetype = "Original")) +
#      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype = "Duration-Weighted")) +
#      xlab("Flow (cfs)") +
#      ylab("WUA (ft2/ft)") +
#      ggtitle("Suitable Habitat Area") +
#      scale_linetype_manual(name = "",
#                            values = c("Original" = "solid",
#                                       "Duration-Weighted" = "dashed"))
#
#    # plt_dwua <- set_dim(plt_dwua, get_dim(plt_dhsi))
#
#    (plt_days + plt_dhsi + plt_dwua) +
#      plot_layout(heights = c(2, 1, 2), axes = "collect", ncol = 1) &
#      theme(legend.position = "bottom") &
#      scale_x_log10(labels = scales::label_comma()) &
#      annotation_logticks(sides = "b")

  })


#   # INTERFACE
#
#   hide(id = "pred_table")
#   hide(id = "attr_table")
#
#   observeEvent(input$showPredTable,{
#     toggle(id = "pred_table",anim = T)
#   })
#
#   observeEvent(input$showAttrTable,{
#     toggle(id = "attr_table",anim = T)
#   })

}
