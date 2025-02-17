shinyUI(
  tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    navbarPage(
      title = tagList(
        tags$div(
          class = "title-container",
          tags$img(src = "habistat_wordmark.svg", class = "title-image", alt = "HabiStat")
        )
      ),
      id = "tabs",
      collapsible = TRUE,
      tabPanel("Interactive Map",
               sidebarPanel(
                 width = 6,
                 div(id = "mainControls",
                     div(style = "display: inline-block;",
                         shinyWidgets::radioGroupButtons("habitat_type",
                                                         label = "Select Habitat Type",
                                                         choices = c("rearing", "spawning"),
                                                         selected = "rearing"),
                         shinyBS::bsPopover(
                           id = "habitat_type",
                           title = "Habitat Info",
                           content = "Select rearing (in-channel + floodplain) or spawning habitat to flow relationships",
                           placement = "right",
                           trigger = "hover"
                         )
                     ),
                     div(style = "display: inline-block;",
                         shinyWidgets::radioGroupButtons("flowline_scope",
                                                         label = "Select Geographic Scope",
                                                         choices = c("comid", "mainstem", "watershed"),
                                                         selected = "comid"),
                         shinyBS::bsPopover(
                           id = "flowline_scope",
                           title = "Geographic Scope Info",
                           content = "comid represents sub-reach identifiers; mainstem refers to the entire mainstem reach; watershed aggregates to the watershed scale.",
                           placement = "right",
                           trigger = "hover"
                         )
                     ),
                     div(style = "display: inline-block;",
                         selectInput("wua_var",
                                     label = "Select Calculation Method",
                                     choices = list("Ensemble" = "wua_per_lf_pred",
                                                    "Scale-Dependent" = "wua_per_lf_pred_SD",
                                                    "Scale-Normalized" = "wua_per_lf_pred_SN",
                                                    "Actual" = "wua_per_lf_actual"),
                                     selected = "wua_per_lf_pred"),
                         shinyBS::bsPopover(
                           id = "wua_var",
                           title = "Calculation Info",
                           content = "Ensemble uses scale-dependent and -normalized averages; scale-dependent = Habitat Area/LF vs Flow; scale-normalized = Normalized Habitat Area/LF vs Normalized Flow; Actual uses observed modeled data.",
                           placement = "right",
                           trigger = "hover"
                         )
                     )
                 ),
                 div(id = "mapControls",
                     div(id = "control_active_flow_slider",
                         shinyWidgets::sliderTextInput("active_flow", "Select Flow (cfs) to show on map",
                                                       choices = seq(100, 10000, 100),
                                                       selected = 1000,
                                                       hide_min_max = TRUE),
                         style = "display:inline-block; width:85%"
                     ),
                     div(id = "control_active_flow_apply",
                         actionButton("activeFlowApplyButton", "Apply"),
                         style = "display:inline-block; width:10%; vertical-align: bottom;"
                     ),
                     uiOutput("out_spawning_toggle")
                 ),
                 h3("Results for Selected Item"),
                 uiOutput("clicked_item_heading"),
                 bslib::navset_tab(id = "tabset_sidebar",
                                   bslib::nav_panel("Suitable Habitat Area by Flow",
                                                    div(id = "controls_fsa",
                                                        uiOutput("units_selector")
                                                    ),
                                                    shinycssloaders::withSpinner(plotOutput("fsa_plot"), hide.ui = FALSE),
                                                    div(id = "predTable",
                                                        DT::DTOutput("pred_table")
                                                    )
                                   ),
                                   bslib::nav_panel("Inundation Duration Analysis",
                                                    div(id = "controls_dur",
                                                        selectInput("selected_run", "Select Run", choices = c("fall", "late fall", "spring", "winter", "steelhead"), selected = "fall"),
                                                        shinyWidgets::radioGroupButtons("selected_wyt", "Select Water Year Type", choices = c("Dry", "Wet"), selected = "Dry")
                                                    ),
                                                    uiOutput("out_streamgage_selector"),
                                                    uiOutput("out_flowscale_toggle"),
                                                    shinycssloaders::withSpinner(plotOutput("dur_plot"), hide.ui = FALSE)
                                   ),
                                   bslib::nav_panel("Flowline Attributes",
                                                    div(id = "attrTable",
                                                        DT::DTOutput("attr_table")
                                                    )
                                   )
                 )
               ),
               mainPanel(
                 width = 6,
                 shinyjs::useShinyjs(),
                 shinycssloaders::withSpinner(leafletOutput("main_map"), hide.ui = FALSE)
               )
      ),
      bslib::nav_item(a(href="https://flowwest.github.io/habistat/articles/interactive_map.html", "User Guide", target = "_blank")),
      bslib::nav_item(a(href="https://flowwest.github.io/habistat", "HabiStat Methods & Documentation", target = "_blank")),
    )
  )
)

