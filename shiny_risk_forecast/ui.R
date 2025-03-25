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
  navbarMenu(
    title="Latest",
    tabPanel(
      title="Map",
      tags$h3("Risk forecast for", textOutput("week.latest", inline=T)),
      sidebarLayout(
        sidebarPanel(
          tags$text("Select a monitoring site to view details."),
          tags$br(),
          tags$text("Click and drag to zoom, and double click to reset view."),
          tags$br(), tags$br(),
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE),
          plotOutput(outputId="all.map.current.pr_distr",
                     height="300px",
                     width="100%"),
          tableOutput("all.map.current.clicked"),
          width=5),
        mainPanel(
          plotOutput(outputId="all.map.current.plot",
                     brush=brushOpts(id="all.map.current.brush",
                                     resetOnNew=T),
                     dblclick=clickOpts(id="all.map.current.dblClick"),
                     click=clickOpts(id="all.map.current.click"),
                     height="600px",
                     width="100%"),
          width=7)
      )
    ),
    tabPanel(
      title="Total performance",
      tags$h3("Total forecasting performance"),
      tags$text("A "), tags$em("skill score"), tags$text(" summarises the predictive ability of a model, where 1 is perfect and 0 is no added value. The top row shows skill scores from 5 metrics across all past observations. The bottom row shows the original metric values."),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE),
          tags$br(), 
          tags$text("The ensemble forecasts (red) are constructed from individual constituent machine learning models (blue). While the constituent models are highly variable, the ensemble forecasts perform consistently. The grey diamond shows a null model based only on the day of year for comparison."),
          tags$br(), tags$br(),
          tags$text("Description of the metrics?")
        ),
        mainPanel(plotOutput(outputId="all.performance.plot",
                             height="600px",
                             width="100%"))
      )
    )
  ),
  navbarMenu(
    title="Past dates",
    tabPanel(
      title="Map",
      tags$h3("Risk forecast for", textOutput("week.selected", inline=T)),
      sidebarLayout(
        sidebarPanel(
          tags$text("Select a monitoring site to view details."),
          tags$br(),
          tags$text("Click and drag to zoom, and double click to reset view."),
          tags$br(), tags$br(),
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE),
          uiOutput("dates.map.weekSelectInput"),
          plotOutput("dates.map.week.pr_v_obs", 
                     height="300px"),
          tableOutput("dates.map.week.clicked"),
          width=5),
        mainPanel(
          plotOutput(outputId="dates.map.week.plot",
                     brush=brushOpts(id="dates.map.week.brush",
                                     resetOnNew=T),
                     dblclick=clickOpts(id="dates.map.week.dblClick"),
                     click=clickOpts(id="dates.map.week.click"),
                     height="600px",
                     width="100%"),
          width=7)
      )
    )
  ),
  navbarMenu(
    title="Monitoring locations",
    tabPanel(
      title="Forecast history",
      tags$h3("Timeseries by monitoring location"),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE
          ),
          uiOutput("sin.ts.sinSelectInput"),
          plotOutput(outputId="sin.map",
                     brush=brushOpts(id="sin.map.brush",
                                     resetOnNew=T),
                     dblclick=clickOpts(id="sin.map.dblClick"),
                     click=clickOpts(id="sin.map.click"),
                     height="600px",
                     width="100%")
        ),
        mainPanel(plotOutput(outputId="sin.timeseries",
                             height="600px",
                             width="100%"))
      )
    ),
    tabPanel(
      title="Past performance",
      tags$h3("Past forecasting performance by monitoring location"),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId="target",
            label="Select HAB or toxin",
            choices=targ_i$fig_long,
            selected=targ_i$fig_long[1],
            multiple=FALSE
          )
        ),
        mainPanel(plotOutput(outputId="sin.performance"))
      )
    )
  )
)


