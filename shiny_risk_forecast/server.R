# Project: HAB Reports Forecast
# www.habreports.org
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Shiny App: server logic


function(input, output, session) {


  # reactives: overview -----------------------------------------------------
  targ <- reactive({
    targ_i |> filter(fig_long==input$target)
  })
  
  targ_df <- reactive({
    all_df |>
      filter(y==targ()$abbr)
  })
  
  targ_available_weeks <- reactive({
    sort(unique(targ_df()$week), decreasing=T)
  })
  
  targ_available_sins <- reactive({
    sort(unique(targ_df()$sin), decreasing=T)
  })
  
  latest_targ <- reactive({
    targ_df() |>
      slice_max(week)
  })
  
  latest_sf <- reactive({
    inner_join(site_sf |> filter(type==targ()$type) |> select(sin, geometry), 
               latest_targ(), 
               by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")))
  })
  
  output$week.latest <- renderText({
    latest_targ()$week |> first() |> format("%Y-%b-%d")
  })
  
  

  # reactives: dates --------------------------------------------------------

  output$dates.map.weekSelectInput <- renderUI({
    selectInput(inputId="week", 
                label="Select week", 
                choices=targ_available_weeks(), 
                selected=first(targ_available_weeks()), 
                multiple=F
                )
  })
  
  selectedWeek_targ <- reactive({
    if(is.null(input$week)) {
      targ_df() |>
        filter(week==max(week))
    } else {
      targ_df() |>
        filter(week==input$week)
    }
  })
  
  selectedWeek_sf <- reactive({
    inner_join(site_sf |> filter(type==targ()$type) |> select(sin, geometry), 
               selectedWeek_targ(), 
               by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")),
             alert=factor(alert, 
                          levels=c("A0", "A1"),
                          labels=c("Below threshold", "Above threshold")))
  })
  
  output$week.selected <- renderText({
    selectedWeek_targ()$week |> first() |> format("%Y-%b-%d")
  })
  
  
  # reactives: SIN ----------------------------------------------------------

  sin_click <- reactiveValues(sin=NULL)
  observe({
    click <- input$sin.map.click
    if(!is.null(click)) {
      res <- nearPoints(site_means_sf |>
                          filter(y==targ()$abbr, 
                                 sin %in% targ_available_sins()) |> 
                          st_drop_geometry(), 
                        input$sin.map.click, "lon", "lat", 
                        maxpoints=1)
      if(nrow(res) > 0) {
        sin_click$sin <- res$sin[1]
      }
    }
  })
  
  output$sin.ts.sinSelectInput <- renderUI({
    if(is.null(sin_click$sin)) {
      selected_sin <- first(targ_available_sins())
    } else {
      selected_sin <- sin_click$sin
    }
    selectInput(inputId="sin", 
                label="Select site", 
                choices=targ_available_sins(), 
                selected=selected_sin, 
                multiple=F
    )
  })
  
  selectedSIN_targ <- reactive({
    if(is.null(input$sin)) {
      targ_df() |>
        filter(sin==first(sin))
    } else {
      targ_df() |> 
        filter(sin==input$sin)
    }
  })
  
  selectedSIN_sf <- reactive({
    inner_join(site_sf |> filter(type==targ()$type) |> select(sin, geometry), 
               selectedSIN_targ(), 
               by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")),
             alert=factor(alert, 
                          levels=c("A0", "A1"),
                          labels=c("Below threshold", "Above threshold")))
  })
  
  
  
  # reactives: validation ---------------------------------------------------
  
  val_targ <- reactive({
    validation_df |>
      filter(y==targ()$abbr) |>
      pivot_longer(c(".estimate", "skill"), names_to="estType") |>
      mutate(estType=factor(estType, 
                            levels=c("skill", ".estimate"),
                            labels=c("Skill score", "Metric value")))
  })
  


  # clicks and brushes ------------------------------------------------------

  map_ranges <- reactiveValues(x=st_bbox(scotland_sf)[c(1,3)],
                               y=st_bbox(scotland_sf)[c(2,4)])
    
  observe({
    brush <- input$all.map.current.brush
    if(!is.null(brush)) {
      map_ranges$x <- c(round(brush$xmin), round(brush$xmax))
      map_ranges$y <- c(round(brush$ymin), round(brush$ymax))
    }
    dblClick <- input$all.map.current.dblClick
    if(!is.null(dblClick)) {
      map_ranges$x=st_bbox(scotland_sf)[c(1,3)]
      map_ranges$y=st_bbox(scotland_sf)[c(2,4)]
    }
  })
  
  observe({
    brush <- input$dates.map.week.brush
    if(!is.null(brush)) {
      map_ranges$x <- c(round(brush$xmin), round(brush$xmax))
      map_ranges$y <- c(round(brush$ymin), round(brush$ymax))
    }
    dblClick <- input$dates.map.week.dblClick
    if(!is.null(dblClick)) {
      map_ranges$x=st_bbox(scotland_sf)[c(1,3)]
      map_ranges$y=st_bbox(scotland_sf)[c(2,4)]
    }
  })
  
  observe({
    brush <- input$sin.map.brush
    if(!is.null(brush)) {
      map_ranges$x <- c(round(brush$xmin), round(brush$xmax))
      map_ranges$y <- c(round(brush$ymin), round(brush$ymax))
    }
    dblClick <- input$sin.map.dblClick
    if(!is.null(dblClick)) {
      map_ranges$x=st_bbox(scotland_sf)[c(1,3)]
      map_ranges$y=st_bbox(scotland_sf)[c(2,4)]
    }
  })
  
  

  # Latest: Map -----------------------------------------------------------

  output$all.map.current.plot <- renderPlot({
    latest_sf() |>
      ggplot() + 
      geom_sf(data=scotland_sf) + 
      geom_sf(aes(colour=prA1, shape=prevAlert, size=prA1), stroke=1) +
      scale_colour_viridis_c("Forecasted\nrisk", 
                             option="inferno", end=0.9, limits=c(0, 1),
                             labels=label_percent()) +
      scale_shape_manual("Latest\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
      scale_size_continuous(limits=c(0, 1), range=c(0.7, 3), guide="none") +
      scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
      scale_y_continuous(limits=map_ranges$y, expand=c(0, 0))
  },
  res=100)
  
  output$all.map.current.pr_distr <- renderPlot({
    latest_sf() |>
      ggplot() + 
      geom_histogram(aes(prA1), fill="grey", colour="grey30", binwidth=0.1) + 
      scale_x_continuous("Forecasted risk", limits=c(-0.05, 1.05), 
                         labels=label_percent(), expand=c(0,0)) +
      ylab("Number of locations")
  }, 
  res=100)
  
  output$all.map.current.clicked <- renderTable({
    res <- nearPoints(latest_sf() |> st_drop_geometry(), 
                      input$all.map.current.click, "lon", "lat") |>
      mutate(Risk=paste0(round(prA1*100), "%"),
             Date=format(week, "%Y-%b-%d")) |>
      rename(SIN=sin, 
             `Last record`=prevAlert) |>
      select(SIN, Date, Risk, `Last record`)
    if (nrow(res) == 0)
      return(NULL)
    res
  })
  
  

  # Latest: Performance ---------------------------------------------------

  output$all.performance.plot <- renderPlot({
    val_targ() |>
      filter(estType=="Skill score") |>
      filter(modType != "Constituent") |>
      ggplot(aes(modType, value, colour=.metric, 
                 group=paste(.metric))) +
      geom_hline(yintercept=c(0, 0.5, 1), linewidth=0.2, colour="grey90") +
      geom_point() + geom_line() +
      scale_colour_brewer("Performance metric",
                          type="qual", palette=2,
                          labels=parse_format()) +
      scale_y_continuous("Skill score", limits=c(0, 1)) +
      scale_x_discrete(labels=parse_format()) +
      theme(axis.title.x=element_blank(),
            legend.position="inside",
            legend.background=element_blank(),
            legend.position.inside=c(0.175, 0.815))
  },
  res=100)
  
  

  # Past dates: Map ---------------------------------------------------------
  
  output$dates.map.week.plot <- renderPlot({
    ggplot(selectedWeek_sf()) + 
      geom_sf(data=scotland_sf) + 
      geom_sf(aes(colour=prA1, shape=alert, size=prA1), stroke=1) +
      scale_colour_viridis_c("Forecasted\nrisk", 
                             option="inferno", end=0.9, limits=c(0, 1),
                             labels=label_percent()) +
      scale_shape_manual("Actual\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
      scale_size_continuous(limits=c(0, 1), range=c(0.7, 3), guide="none") +
      scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
      scale_y_continuous(limits=map_ranges$y, expand=c(0, 0))
  },
  res=100)
  
  output$dates.map.week.pr_v_obs <- renderPlot({
    ggplot(selectedWeek_sf()) +
      stat_histinterval(aes(prA1, alert, fill=alert), 
                                scale=0.5, breaks=seq(0, 1, by=0.1)) +
      geom_dots(aes(prA1, alert, colour=prA1, shape=alert), 
                        side="bottom", scale=0.3, layout="swarm") +
      scale_colour_viridis_c("Forecasted risk", option="inferno", end=0.9) +
      scale_fill_viridis_d("Actual\nobservation", option="inferno", end=0.85) +
      scale_shape_manual("Actual\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
      scale_x_continuous("Forecasted risk", limits=c(-0.05, 1.05), 
                         labels=label_percent(), expand=c(0,0)) +
      scale_y_discrete(limits=paste(c("Below", "Above"), "threshold"),
                       breaks=paste(c("Below", "Above"), "threshold"),
                       labels=paste(c("Below", "Above"), "threshold", sep="\n")) +
      ggtitle("Actual vs. forecasted risk") +
      theme(legend.position="none",
            axis.title.y=element_blank())
  },
  res=100)
  
  output$dates.map.week.clicked <- renderTable({
    res <- nearPoints(selectedWeek_sf() |> st_drop_geometry(), 
                      input$dates.map.week.click, "lon", "lat") |>
      mutate(Risk=paste0(round(prA1*100), "%")) |>
      rename(SIN=sin, 
             Actual=alert) |>
      select(SIN, Risk, Actual)
    if (nrow(res) == 0)
      return(NULL)
    res
  })
  
  

  # Locations: Timeseries ---------------------------------------------------
  
  output$sin.map <- renderPlot({
    site_sf |> 
      filter(type==targ()$type, 
             sin %in% targ_available_sins()) |>
      mutate(selectedSIN=sin==selectedSIN_targ()$sin[1]) |>
      arrange(selectedSIN) |>
      ggplot() + 
      geom_sf(data=scotland_sf) + 
      geom_sf(aes(size=selectedSIN, shape=selectedSIN, colour=selectedSIN), stroke=1) +
      scale_size_manual(values=c(1.5, 3), guide="none") +
      scale_shape_manual(values=c(1, 19), guide="none") +
      scale_colour_manual(values=c("grey20", "blue"), guide="none") +
      scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
      scale_y_continuous(limits=map_ranges$y, expand=c(0, 0)) +
      theme(legend.position="bottom",
            legend.key.height=unit(2, "mm"),
            legend.key.width=unit(10, "mm"),
            legend.title=element_text(size=9),
            legend.title.position="top")
  },
  res=100)
  
  output$sin.timeseries <- renderPlot({
    selectedSIN_sf() |>
      ggplot() + 
      geom_ribbon(data=obs_df |> 
                    filter(y==targ()$abbr) |>
                    filter(sin==selectedSIN_targ()$sin[1]) |>
                    filter(between(week, min(selectedSIN_sf()$week), max(selectedSIN_sf()$week))),
                 aes(week, ymin=0, ymax=lnN_rel), fill="grey90", colour="grey90") +
      geom_point(aes(week, prA1), shape=1) +
      scale_y_continuous("Forecasted risk", limits=c(0, 1), 
                         labels=label_percent(), expand=c(0,0)) +
      scale_x_date(date_breaks="1 year", date_labels="%Y-%b") +
      theme(axis.title.x=element_blank(),
            panel.grid.major=element_line(colour="grey90", linewidth=0.15)) +
      ggtitle(paste0(targ()$fig_long, ": ", selectedSIN_sf()$sin[1]))
  },
  res=100)
  
  
}
