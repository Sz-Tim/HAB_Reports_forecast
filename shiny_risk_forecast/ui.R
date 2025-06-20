# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Shiny App: user interface

library(tidyverse)
library(shiny)
library(rprojroot)

dir(pattern="globals.R", recursive=T, full.names=T) |>
  source()

navbarPage(
  title="HAB risk forecasts", theme=shinythemes::shinytheme("flatly"),
  tabPanel(
    title="Weekly forecasts",
    tags$h1("Risk forecast for", textOutput("week.selected", inline=T)),
    sidebarLayout(
      sidebarPanel(
        tags$text("Select a monitoring site to view details."),
        tags$br(), tags$br(),
        selectInput(
          inputId="target",
          label="Select HAB or toxin",
          choices=c("All", targ_i$fig_long),
          selected="All",
          multiple=FALSE),
        selectInput(
          inputId="mapView",
          label="Reset map view",
          choices=c("All", "Shetland", "Outer Isles", "Inner Isles"),
          selected="All",
          multiple=F
        ),
        uiOutput("dates.map.weekSelectInput"),
      ),
      mainPanel(
        leafletOutput(outputId="dates.map.week.plot", 
                      height="700px")
      )
    )
  ),
  # tabPanel(
  #   title="Weekly forecasts (old)",
  #   tags$h1("Risk forecast for", textOutput("week.selected", inline=T)),
  #   sidebarLayout(
  #     sidebarPanel(
  #       tags$text("Select a monitoring site to view details."),
  #       tags$br(),
  #       tags$text("Click and drag to zoom, and double click to reset view."),
  #       tags$br(), tags$br(),
  #       selectInput(
  #         inputId="target",
  #         label="Select HAB or toxin",
  #         choices=targ_i$fig_long,
  #         selected=targ_i$fig_long[1],
  #         multiple=FALSE),
  #       selectInput(
  #         inputId="mapView",
  #         label="Reset map view",
  #         choices=c("All", "Shetland", "Outer Isles", "Inner Isles"),
  #         selected="All",
  #         multiple=F
  #       ),
  #       uiOutput("dates.map.weekSelectInput"),
  #       radioButtons(inputId="alertLag",
  #                    label="Separate forecasts by",
  #                    choices=list(
  #                      "Previous state"="previous",
  #                      "Actual state"="actual"
  #                    )),
  #       plotOutput("dates.map.week.pr_v_obs", 
  #                  height="300px"),
  #       tableOutput("dates.map.week.clicked"),
  #       width=5),
  #     mainPanel(
  #       tags$h3("Option 1:"),
  #       plotOutput(outputId="dates.map.week.plot.OLD",
  #                  brush=brushOpts(id="dates.map.week.brush",
  #                                  resetOnNew=T),
  #                  dblclick=clickOpts(id="dates.map.week.dblClick"),
  #                  click=clickOpts(id="dates.map.week.click"),
  #                  height="600px",
  #                  width="100%"),
  #       tags$h3("Option 2:"),
  #       plotlyOutput(outputId="dates.map.week.plotly",
  #                    height="500px",
  #                    width="100%"),
  #       width=7)
  #   )
  # ),
  tabPanel(
    title="Locations",
    tags$h1("Timeseries by monitoring location"),
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId="target",
          label="Select HAB or toxin",
          choices=targ_i$fig_long,
          selected=targ_i$fig_long[1],
          multiple=FALSE
        ),
        selectInput(
          inputId="mapView",
          label="Reset map view",
          choices=c("All", "Shetland", "Outer Isles", "Inner Isles"),
          selected="All",
          multiple=F
        ),
        uiOutput("sin.ts.sinSelectInput"),
        plotOutput(outputId="sin.map",
                   brush=brushOpts(id="sin.map.brush",
                                   resetOnNew=T),
                   dblclick=clickOpts(id="sin.map.dblClick"),
                   click=clickOpts(id="sin.map.click"),
                   height="500px",
                   width="100%")
      ),
      mainPanel(plotlyOutput(outputId="sin.timeseries",
                             height="600px",
                             width="100%"),
                width=8)
    )
  ),
  navbarMenu(
    title="Performance",
    tabPanel(
      title="By month",
      tags$h1("Past forecasting performance by month"),
      tags$text("A "), tags$em("skill score"), tags$text(" summarises the predictive ability of a model, where 1 is perfect and 0 is no added value compared to an overall mean."),
      tags$text("The model skill is assessed for two years of testing data using five performance metrics which each measure a different aspect of 'performance'. The values shown are the average among these metrics."),
      tags$text("Note that skill is only calculated for months with at least one 'hit' in the testing data."),
      # tags$h3("Option 1:"),
      # plotOutput(outputId="eval.month.plot.lines",
      #            height="300px", 
      #            width="100%"),
      # tags$h3("Option 2:"),
      plotlyOutput(outputId="eval.month.plot.lines.plotly",
                   height="400px", 
                   width="100%")
      # tags$h3("Option 3:"),
      # plotOutput(outputId="eval.month.plot",
      #            height="400px",
      #            width="100%")
    ),
    tabPanel(
      title="By location",
      tags$h1("Past forecasting performance by monitoring location"),
      tags$text("A "), tags$em("skill score"), tags$text(" summarises the predictive ability of a model, where 1 is perfect and 0 is no added value compared to an overall mean."),
      tags$text("The model skill is assessed for two years of testing data using five performance metrics which each measure a different aspect of 'performance'. The values shown are the average among these metrics."),
      tags$text("Note that skill is only calculated for sites with at least one 'hit' in the testing data."),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE
          ),
          selectInput(
            inputId="mapView",
            label="Reset map view",
            choices=c("All", "Shetland", "Outer Isles", "Inner Isles"),
            selected="All",
            multiple=F
          ),
          uiOutput("eval.sin.sinSelectInput"),
          plotOutput(outputId="eval.sin.plot",
                     height="250px",
                     width="100%"),
          width=5
        ),
        mainPanel(plotOutput(outputId="eval.sin.map",
                             brush=brushOpts(id="eval.sin.map.brush",
                                             resetOnNew=T),
                             dblclick=clickOpts(id="eval.sin.map.dblClick"),
                             click=clickOpts(id="eval.sin.map.click"),
                             height="600px",
                             width="100%"),
                  width=7)
      )
    )
  ),
  tabPanel(
    title="About",
    tags$br(),
    tags$h3("Funding"),
    tags$text("See current and previous"), tags$a("contributors", href="https://www.habreports.org/about.php"), tags$text("to HAB reports. Additional funding has been provided by"), tags$a("UHI Aquaculture Hub", href="https://www.uhi.ac.uk/en/research-enterprise/res-themes/mese/aquaculture/"), tags$text(""),
    tags$br(),
    tags$h3("Methods"),
    tags$text("See"), tags$a("Szewczyk et al. (2025)", href="https://doi.org/10.1016/j.hal.2024.102781"), tags$text("for full technical detail."),
    tags$br(),
    tags$h4("Data"),
    tags$text("Monitoring data is accessed through the"), tags$a("HAB Reports website", href="https://www.habreports.org"), 
    tags$text("which compiles observations performed by CEFAS, FSS, and SAMS."), 
    tags$text("Environmental data is accessed from the"), tags$a("Copernicus European North West Shelf biogeochemistry model", href="https://data.marine.copernicus.eu/product/NWSHELF_ANALYSISFORECAST_BGC_004_002/description"),
    tags$text("and the"), tags$a("Weather Research Forecast model", href="https://thredds.sams.ac.uk/thredds/catalog/scoats-wrf/catalog.html"), tags$text("run operationally by SAMS."),
    tags$h4("Models"),
    tags$text("The forecasts are from an ensemble of machine learning and Bayesian models. Forecasts are generated by each constituent model, and these forecasts are then used as inputs to the ensemble model."), 
    tags$h4("Performance"),
    tags$text("Model performance is assessed using five metrics which capture different definitions of 'performance'. These metrics are converted to a skill score, ranging from 0 (no information) to 1 (perfect prediction). The total model skill is the mean of these five scores.")
  )
)

