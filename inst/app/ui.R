shinyUI(
  tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      #tags$style("body > .container-fluid {margin: 0px; padding: 0px;}")
    ),
    navbarPage(
      title = "HabiStat",
      id = "tabs",
      collapsible = TRUE,
      tabPanel("Interactive Map",
               sidebarPanel(
                 #style = "position:fixed; width:calc(inherit - 15px); overflow-y: scroll; height:100%; margin:0px; padding:0px; border: 0px; background: none;",
                 width = 6,
                 div(id = "mainControls",
                     div(style="display: inline-block;",
                         shinyWidgets::radioGroupButtons("habitat_type", "Select Habitat Type", choices=c("rearing", "spawning"), selected="rearing")),
                     div(style="display: inline-block;",
                         shinyWidgets::radioGroupButtons("flowline_scope", "Select Flowline Type", choices=c("comid", "mainstem"), selected="comid")),
                     div(style="display: inline-block;",
                     shinyWidgets::radioGroupButtons("wua_var", "Select Calculation Method", list("Scale-Dependent" = "wua_per_lf_pred_SD",
                                                                                                  "Scale-Normalized" = "wua_per_lf_pred_SN",
                                                                                                  "Actual" = "wua_per_lf_actual"))),
                 ),
                 div(id = "mapControls",
                     div(id = "control_active_flow_slider",
                         shinyWidgets::sliderTextInput("active_flow", "Select Flow (cfs) to show on map", choices=all_flows_idx, selected=1000, hide_min_max=T),
                         style="display:inline-block; width:85%"),
                     div(id = "control_active_flow_apply",
                         actionButton("activeFlowApplyButton" ,"Apply"),
                         style="display:inline-block; width:10%; vertical-align: bottom;"),
                     uiOutput("out_spawning_toggle"),
                 ),
                 h3("Results for Selected Item"),
                 uiOutput("clicked_item_heading"),
                 bslib::navset_tab(id = "tabset_sidebar",
                                   bslib::nav_panel("Suitable Habitat Area by Flow",
                                                    div(id = "controls_fsa",
                                                      uiOutput("units_selector"),
                                                    ),
                                                    shinycssloaders::withSpinner(plotOutput("fsa_plot"), hide.ui=F),
                                                    div(id = "predTable",
                                                        DT::DTOutput("pred_table")),
                                   ),
                                   bslib::nav_panel("Inundation Duration Analysis",
                                                    div(id = "controls_dur",
                                                      selectInput("selected_run", "Select Run", choices=c("fall", "late fall", "spring", "winter", "steelhead"), selected="fall"),
                                                      shinyWidgets::radioGroupButtons("selected_wyt", "Select Water Year Type", choices=c("Dry", "Wet"), selected="Dry"),
                                                    ),
                                                    uiOutput("out_streamgage_selector"),
                                                    uiOutput("out_flowscale_toggle"),
                                                    shinycssloaders::withSpinner(plotOutput("dur_plot"), hide.ui=F),
                                   ),
                                   bslib::nav_panel("Flowline Attributes",
                                                    div(id = "attrTable",
                                                        DT::DTOutput("attr_table")
                                                    )
                                   )
                 ),
               ),
               mainPanel(
                 width = 6, # main width plus sidebar width should add to 12
                 shinyjs::useShinyjs(),  # Initialize shinyjs
                 shinycssloaders::withSpinner(leafletOutput("main_map"), hide.ui=F),
                 #numericInput("active_flow", "Select Flow (cfs) to show on map", min=min(all_flows), max=max(all_flows), step=100, value=1000),
                 #shinycssloaders::withSpinner(leafgl::leafglOutput("main_map"))
               )
      )
    )
  )
)
