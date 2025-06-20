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
  
  targ_tl <- reactive({
    tl_i |> filter(abbr==targ()$abbr)
  })
  
  targ_df <- reactive({
    all_df |>
      filter(y==targ()$abbr)
  })
  
  targ_available_weeks <- reactive({
    sort(unique(targ_df()$week), decreasing=T)
  })
  
  targ_available_sins <- reactive({
    sort(unique(targ_df()$sin))
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

  sin_click <- reactiveValues(sin=NULL, eval.sin=NULL)
  
  observe({
    click <- input$sin.map.click
    if(!is.null(click)) {
      res <- nearPoints(site_sf |>
                          st_drop_geometry() |>
                          filter(sin %in% targ_available_sins()),
                        click, "lon", "lat",
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
  
  val_sin_targ <- reactive({
    validation_sin_df |>
      filter(y==targ()$abbr) |>
      pivot_longer(c(".estimate", "skill"), names_to="estType") |>
      mutate(estType=factor(estType, 
                            levels=c("skill", ".estimate"),
                            labels=c("Skill score", "Metric value")))
  })
  
  val_sin_mn_targ <- reactive({
    val_sin_targ() |> 
      filter(estType=="Skill score", 
             modType=="Ensemble") |> 
      group_by(y, sin) |> 
      summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) 
  })
  
  observe({
    click <- input$eval.sin.map.click
    if(!is.null(click)) {
      res <- nearPoints(site_sf |>
                          st_drop_geometry() |>
                          filter(sin %in% targ_available_sins()),
                        click, "lon", "lat",
                        maxpoints=1)
      if(nrow(res) > 0) {
        sin_click$eval.sin <- res$sin[1]
      }
    }
  })
  
  output$eval.sin.sinSelectInput <- renderUI({
    if(is.null(sin_click$eval.sin)) {
      selected_sin <- first(targ_available_sins())
    } else {
      selected_sin <- sin_click$eval.sin
    }
    selectInput(inputId="sin", 
                label="Select site", 
                choices=targ_available_sins(), 
                selected=selected_sin, 
                multiple=F
    )
  })
  
  selectedSIN_eval <- reactive({
    if(is.null(input$sin)) {
      val_sin_targ() |>
        filter(sin==first(sin))
    } else {
      val_sin_targ() |> 
        filter(sin==input$sin)
    }
  })
  


  # clicks and brushes ------------------------------------------------------

  map_ranges <- reactiveValues(x=st_bbox(scotland_sf)[c(1,3)],
                               y=st_bbox(scotland_sf)[c(2,4)])
  
  # all.map.current  
  observe({
    update_map_range(input$all.map.current.brush, map_ranges, "brush")
    update_map_range(input$all.map.current.dblClick, map_ranges, "dblClick", scot_bbox)
  })
  
  # dates.map.week
  observe({
    update_map_range(input$dates.map.week.brush, map_ranges, "brush")
    update_map_range(input$dates.map.week.dblClick, map_ranges, "dblClick", scot_bbox)
  })
  
  # sin.map
  observe({
    update_map_range(input$sin.map.brush, map_ranges, "brush")
    update_map_range(input$sin.map.dblClick, map_ranges, "dblClick", scot_bbox)
  })
  
  # eval.sin.map
  observe({
    update_map_range(input$eval.sin.map.brush, map_ranges, "brush")
    update_map_range(input$eval.sin.map.dblClick, map_ranges, "dblClick", scot_bbox)
  })
  
  

  # Weekly forecasts: Map ---------------------------------------------------

  output$dates.map.week.plot <- renderLeaflet({
    latest_df <- all_df |> 
      filter(year(date)==max(year(date))) |>
      group_by(y, sin) |>
      slice_max(date, n=1, with_ties=F) |>
      ungroup() |>
      left_join(targ_i |> select(abbr, fig_short, plotGroups, targ_ordered, col),
                by=join_by(y==abbr)) |>
      select(-lon, -lat) |>
      left_join(site_wgs |> select(sin, lon, lat)) |>
      inner_join(obs_df |> select(type, sin, y, week, lnN, tl)) |>
      left_join(tl_i |> group_by(targ_ordered) |> slice_head(n=1) |> select(targ_ordered, units)) |>
      mutate(units=str_replace_all(units, "\\\\", "/"),
             tl_col=case_when(tl=="TL0" ~ "#b2df8a",
                              tl=="TL1" ~ "#fed976",
                              tl=="TL2" ~ "#fd8d3c",
                              tl=="TL3" ~ "#e31a1c"),
             N=paste(format(round(expm1(lnN)), big.mark=",", trim=T), units)) |>
      filter(max(date) - date <= 14)
    
    
    fcst_cols <- viridis::inferno(102, end=0.85)
    tl_cols <- tibble(tl=paste0("TL", 0:3),
                      tl_col=c("#b2df8a", "#fed976", "#fd8d3c", "#e31a1c"))
    latest_leaflet_df <- latest_df |>
      mutate(fcst_pct=round(prA1*100, 0),
             fcst_chr=str_pad(fcst_pct, 2) |> str_replace(" ", "&#8194;")) |>
      mutate(y_lab=paste0("<td><span style='background-color: ", col, "'>&#8194;</span> ", fig_short, "</td>",
                          "<td><span style='color: ", fcst_cols[fcst_pct+1], "'>", fcst_chr, "%</span>&#8194;</td>", 
                          "<td><span style='color: ", tl_col, "'>■</span> ", N, "</td>",
                          "<td>", date, "</td>")) |> 
      arrange(targ_ordered) |>
      group_by(y, sin) |>
      slice_max(date) |>
      group_by(sin) |>
      summarise(lon=first(lon), 
                lat=first(lat),
                maxTL=max(tl),
                maxPr=max(prA1),
                lab=paste("<table><tr><th>Target</th><th>Risk</th><th>Latest</th><th>Date</th></tr><tr>", 
                          paste(y_lab, collapse="</tr><tr>"), 
                          "</tr></table>")) |>
      ungroup()
    
    tl_pal <- colorFactor(tl_cols$tl_col, domain=tl_cols$tl)
    pr_pal <- colorNumeric(viridis::inferno(100, end=0.9), c(0, 1))
    leaflet(data=latest_leaflet_df) |>
      addProviderTiles(providers$Esri.WorldTopoMap) |>
      fitBounds(-8.53, 54.3, 0, 61) |>
      addCircleMarkers(~lon, ~lat, radius=6,
                       popup=~paste0("<b>Site name</b><br>", 
                                     "<span style='font-size: 11px'>", sin, "</span>", 
                                     "<hr>", lab),
                       weight=2,
                       color="#252525",
                       fillColor=~tl_pal(maxTL), 
                       fillOpacity=0.8) |>
      addLegend("topleft", colors=tl_cols$tl_col, labels=c("Green", "Yellow", "Amber", "Red"),
                title="Highest<br>status", opacity=1)
  })
  output$dates.map.week.plot.OLD <- renderPlot({
    map_latest_sf <- inner_join(site_sf |> filter(type==targ()$type) |> select(sin, geometry), 
               targ_df() |>
                 select(-lon, -lat) |>
                 filter(week <= selectedWeek_targ()$week[1]) |>
                 group_by(sin) |>
                 slice_max(week) |>
                 filter(week >= (selectedWeek_targ()$week[1] - ddays(14))), 
               by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")),
             alert=factor(alert, 
                          levels=c("A0", "A1"),
                          labels=c("Below threshold", "Above threshold")))
    if(input$alertLag=="previous") {
      ggplot(map_latest_sf) +
      # ggplot(selectedWeek_sf()) + 
        geom_sf(data=scotland_sf) + 
        geom_sf(aes(colour=prA1, shape=prevAlert, size=prA1), stroke=1) +
        scale_colour_viridis_c("Forecasted\nrisk", 
                               option="inferno", end=0.9, limits=c(0, 1),
                               labels=label_percent()) +
        scale_shape_manual("Previous\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
        scale_size_continuous(limits=c(0, 1), range=c(0.7, 3), guide="none") +
        scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
        scale_y_continuous(limits=map_ranges$y, expand=c(0, 0))
    } else {
      ggplot(map_latest_sf) +
        # ggplot(selectedWeek_sf()) + 
        geom_sf(data=scotland_sf) + 
        geom_sf(aes(colour=prA1, shape=alert, size=prA1), stroke=1) +
        scale_colour_viridis_c("Forecasted\nrisk", 
                               option="inferno", end=0.9, limits=c(0, 1),
                               labels=label_percent()) +
        scale_shape_manual("Actual\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
        scale_size_continuous(limits=c(0, 1), range=c(0.7, 3), guide="none") +
        scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
        scale_y_continuous(limits=map_ranges$y, expand=c(0, 0))
    }
  },
  res=100)
  
  output$dates.map.week.pr_v_obs <- renderPlot({
    map_latest_sf <- inner_join(site_sf |> filter(type==targ()$type) |> select(sin, geometry), 
                                targ_df() |>
                                  select(-lon, -lat) |>
                                  filter(week <= selectedWeek_targ()$week[1]) |>
                                  group_by(sin) |>
                                  slice_max(week) |>
                                  filter(week >= (selectedWeek_targ()$week[1] - ddays(14))), 
                                by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")),
             alert=factor(alert, 
                          levels=c("A0", "A1"),
                          labels=c("Below threshold", "Above threshold")))
    if(input$alertLag=="previous") {
      ggplot(map_latest_sf) +
        # ggplot(selectedWeek_sf()) + 
        stat_histinterval(aes(prA1, prevAlert, fill=prevAlert), 
                          scale=0.5, breaks=seq(0, 1, by=0.1)) +
        geom_dots(aes(prA1, prevAlert, colour=prA1, shape=prevAlert), 
                  side="bottom", scale=0.3, layout="swarm") +
        scale_colour_viridis_c("Forecasted risk", option="inferno", end=0.9, limits=c(0, 1)) +
        scale_fill_viridis_d("Previous\nobservation", option="inferno", end=0.85) +
        scale_shape_manual("Previous\nobservation", values=c(1, 5), labels=label_wrap_gen(10)) +
        scale_x_continuous("Forecasted risk", limits=c(-0.05, 1.05), 
                           labels=label_percent(), expand=c(0,0)) +
        scale_y_discrete(limits=paste(c("Below", "Above"), "threshold"),
                         breaks=paste(c("Below", "Above"), "threshold"),
                         labels=paste(c("Below", "Above"), "threshold", sep="\n")) +
        ggtitle("Previous vs. forecasted risk") +
        theme(legend.position="none",
              axis.title.y=element_blank())
    } else {
      ggplot(map_latest_sf) +
        # ggplot(selectedWeek_sf()) + 
        stat_histinterval(aes(prA1, alert, fill=alert), 
                          scale=0.5, breaks=seq(0, 1, by=0.1)) +
        geom_dots(aes(prA1, alert, colour=prA1, shape=alert), 
                  side="bottom", scale=0.3, layout="swarm") +
        scale_colour_viridis_c("Forecasted risk", option="inferno", end=0.9, limits=c(0, 1)) +
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
    }
    
  },
  res=100)
  
  output$dates.map.week.clicked <- renderTable({
      res <- nearPoints(selectedWeek_sf() |> st_drop_geometry(), 
                        input$dates.map.week.click, "lon", "lat") |>
        mutate(Risk=paste0(round(prA1*100), "%")) |>
        rename(SIN=sin, 
               Previous=prevAlert,
               Actual=alert) |>
        select(SIN, Risk, Previous, Actual)
    if (nrow(res) == 0)
      return(NULL)
    res
  })
  
  output$dates.map.week.plotly <- renderPlotly({
    site_wgs <- site_sf |>
      st_transform(4326) |>
      sevcheck::add_lonlat(drop_geom=T)
    targ_latest_df <- targ_df() |>
      select(-lon, -lat) |>
      filter(week <= selectedWeek_targ()$week[1]) |>
      group_by(sin) |>
      slice_max(week) |>
      filter(week >= (selectedWeek_targ()$week[1] - ddays(14))) |>
      mutate(Forecast=paste0(round(prA1*100), "%"))
    selected_wk_plotly <- site_wgs |> 
      filter(type==targ_latest_df$type[1]) |> 
      select(sin, lon, lat) |>
      inner_join(targ_latest_df, by=join_by(sin)) |>
      arrange(prA1) |>
      mutate(prevAlert=factor(prevAlert, 
                              levels=c("A0", "A1"),
                              labels=c("Below threshold", "Above threshold")),
             alert=factor(alert, 
                          levels=c("A0", "A1"),
                          labels=c("Below threshold", "Above threshold"))) |>
      mutate(prevAlertNum=as.character(prevAlert))
    plot_ly(selected_wk_plotly, 
            type="scattermapbox", 
            mode="markers",
            lat=~lat, 
            lon=~lon, 
            marker=list(cmin=0, cmax=1, cauto=FALSE, color=~prA1, 
                        colorbar=list(title="Forecast", ticks="inside", dtick=0.25, tickformat=".0%"), 
                        colorscale="RdBu", allowoverlap=TRUE, reversescale=F),
            # color=~prA1,
            # colors=viridis::inferno(100),
            # symbol=~prevAlert,
            # symbols=c(18, 19),
            text=~sin,
            size=2,
            hovertemplate=paste("<b>%{text}</b><br>",
                                "Forecast: %{marker.color: .0%}")) |>
      layout(mapbox = list(
        style='open-street-map',
        zoom=4.5,
        center=list(lon = -5, lat = 58)))
  })
  
  

  # Locations: Timeseries ---------------------------------------------------
  
  output$sin.map <- renderPlot({
    site_sf |> 
      filter(type==targ()$type, 
             sin %in% targ_available_sins()) |>
      ggplot() + 
      geom_sf(data=scotland_sf) + 
      geom_sf(colour="grey20", size=1.5, shape=1, stroke=1) +
      geom_sf(data=site_sf |> filter(type==targ()$type, 
                                     sin==selectedSIN_targ()$sin[1]),
              colour="blue", size=3.5, shape=1, stroke=2) +
      scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
      scale_y_continuous(limits=map_ranges$y, expand=c(0, 0)) +
      theme(legend.position="bottom",
            legend.key.height=unit(2, "mm"),
            legend.key.width=unit(10, "mm"),
            legend.title=element_text(size=9),
            legend.title.position="top")
  },
  res=100)

  output$sin.timeseries <- renderPlotly({
    combined_df <- bind_rows(
      obs_df |> 
        filter(y==targ()$abbr,
               sin==selectedSIN_targ()$sin[1],
               between(week, min(selectedSIN_sf()$week), max(selectedSIN_sf()$week))) |>
        mutate(lab="Observed",
               yval=lnN_rel,
               lwidth=1),
      selectedSIN_sf() |> st_drop_geometry() |>
        mutate(lab="Forecast",
               yval=prA1,
               lwidth=2)
    ) |>
      arrange(lab, week)
    max_lnN <- combined_df |> slice_max(lnN_rel)
    tl_ticks <- targ_tl()$min_lnN / (max_lnN$lnN / max_lnN$lnN_rel)
    tl_labs <- targ_tl()$min_ge
    
    plot_ly(combined_df,
            x=~week, 
            y=~yval, 
            mode="lines+markers", 
            type="scatter", 
            color=~lab,
            colors=c("#084594", "#bdbdbd"),
            symbol=~lab, 
            symbols=c('circle', 'o'),
            # linetype=~lab,
            line=list(width=~lwidth),
            hovertemplate=paste("<b>%{x}</b><br>",
                                "%{y: .0%}<br>"),
            connectgaps=F) |>
      add_trace(x=~week, y=~yval, yaxis="y2", hoverinfo="none",
                colors=NULL, showlegend=F, inherit=F, mode="none") |>
      layout(xaxis=list(title="Date", tickmode="auto"),
             yaxis=list(title="Forecast", 
                        range=c(0, 1), 
                        dtick=0.25,
                        tickformat=".0%",
                        gridcolor="#FFFFFF"),
             yaxis2=list(title=ifelse(targ()$type=="hab", 
                                      "Observed density (cells/l)",
                                      "Observed concentration (ug/kg)"), 
                         range=c(0, 1), 
                         tickmode="array",
                         tickvals=tl_ticks,
                         ticktext=tl_labs,
                         gridcolor="#bdbdbd",
                         overlaying="y", side="right"))
    
  })
  
  
  
  
  
  # Performance: Total ------------------------------------------------------
  
  output$eval.total.plot <- renderPlot({
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
  
  
  

  # Performance: By month ---------------------------------------------------
  
  output$eval.month.plot <- renderPlot({
    p.month <- validation_month_df |>
      filter(modType=="Ensemble") |>
      group_by(plotGroups, y, modType, date_std) |>
      summarise(skill=mean(skill, na.rm=T)) |>
      ungroup() |>
      inner_join(targ_i |> select(abbr, targ_ordered), by=join_by(y==abbr)) |>
      ggplot(aes(date_std, targ_ordered, fill=skill)) +
      geom_raster() +
      scale_fill_viridis_b("Skill\nscore", limits=c(NA, 1), n.breaks=6) +
      scale_x_date("Month", date_breaks="month", date_labels="%b", 
                   expand=c(0, 0), oob=scales::oob_keep) +
      facet_grid(plotGroups~., scales="free_y", space="free_y") +
      ggtitle("Monthly") +
      theme(panel.background=element_rect(fill="grey97", colour="black"),
            strip.background=element_blank(),
            strip.text=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank())
    p.all <- validation_df |>
      filter(modType=="Ensemble") |>
      group_by(y) |>
      summarise(skill=mean(skill, na.rm=T)) |>
      ungroup() |>
      inner_join(targ_i |> select(abbr, targ_ordered, plotGroups), by=join_by(y==abbr)) |>
      mutate(yval=1) |>
      ggplot(aes(yval, targ_ordered, label=format(skill, digits=2), fill=skill)) +
      geom_raster(alpha=0.75) +
      geom_text() +
      scale_fill_viridis_c("Skill\nscore", limits=c(0, 1), n.breaks=5, guide="none") +
      scale_y_discrete(labels=label_wrap_gen(12)) +
      facet_grid(plotGroups~., scales="free_y", space="free_y") +
      ggtitle("Total") +
      theme(panel.background=element_rect(fill="grey97", colour="black"),
            strip.background=element_blank(),
            strip.text=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank())
    cowplot::plot_grid(p.all, p.month, align="h", axis="tb", nrow=1, rel_widths=c(0.2, 1))
  },
  res=100)
  
  output$eval.month.plot.lines <- renderPlot({
    monthly_skill_df <- validation_month_df |>
      filter(modType=="Ensemble") |>
      group_by(plotGroups, y, modType, date_std) |>
      summarise(skill=mean(skill, na.rm=T)) |>
      ungroup() |>
      inner_join(targ_i |> select(abbr, targ_ordered, fig_short), by=join_by(y==abbr)) |>
      mutate(month=month(date_std, label=F))
    monthly_skill_df |>
      ggplot(aes(month, skill, colour=targ_ordered)) +
      geom_point() +
      stat_smooth(se=F, method="loess", span=1) +
      geom_rug(data=monthly_skill_df |> group_by(plotGroups, targ_ordered, fig_short) |>
                 summarise(skill=median(skill),
                           month=1), sides="l",
               linewidth=1) +
      geom_text(data=monthly_skill_df |> group_by(plotGroups, targ_ordered, fig_short) |>
                  summarise(skill=median(skill)) |> 
                  mutate(month=1,
                         lab=fig_short),
                aes(label=lab), size=3, hjust=0, nudge_x=0.25, fontface="bold", colour="black") +
      scale_colour_brewer("", palette="Paired") +
      scale_x_continuous("Month", breaks=1:12, labels=str_sub(month.abb, 1, 1)) +
      scale_y_continuous("Skill score", limits=c(0, 1), breaks=c(0, 0.5, 1), 
                         labels=c("(no info)\n0\n", "0.5", "\n1\n(perfect)"), 
                         oob=scales::oob_keep) +
      facet_grid(.~plotGroups) +
      theme(legend.position="none",
            panel.grid.major.y=element_line(colour="grey90", linewidth=0.5),
            panel.background=element_rect(fill=NA, colour="black"),
            strip.background=element_blank(),
            strip.text=element_blank())
  },
  res=100)
  
  output$eval.month.plot.lines.plotly <- renderPlotly({
    total_skill_df <- validation_df |>
      filter(modType=="Ensemble") |>
      group_by(y) |>
      summarise(mnSkill=mean(skill, na.rm=T)) |>
      ungroup() 
    monthly_skill_df <- validation_month_df |>
      filter(modType=="Ensemble") |>
      group_by(plotGroups, y, modType, date_std) |>
      summarise(skill=mean(skill, na.rm=T)) |>
      ungroup() |>
      inner_join(targ_i |> select(abbr, targ_ordered, fig_short), by=join_by(y==abbr)) |>
      mutate(month=month(date_std)) |>
      left_join(total_skill_df)
    plot_ly(monthly_skill_df, x=~month) |>
      add_trace(y=~skill, 
                mode="markers+lines",
                color=~targ_ordered,
                colors="Paired",
                text=~mnSkill,
                type="scatter",
                hovertemplate=paste("Skill:<br>",
                                    "%{y: .2} (%{x})<br>",
                                    "%{text: .2} (Overall)")) |>
      layout(yaxis=list(title="Skill score",
                        range=c(0, 1), 
                        tickmode="array",
                        tickvals=c(0, 0.5, 1),
                        ticktext=c("(no info)\n0\n", "0.5", "\n1\n(perfect)")),
             xaxis=list(title="Month",
                        range=c(0.5, 12.5),
                        tickmode="array",
                        tickvals=1:12,
                        ticktext=month.abb))
                        # ticktext=str_sub(month.abb, 1, 1)))
  })
  

  
  # Performance: By SIN -----------------------------------------------------
  
  output$eval.sin.map <- renderPlot({
    val_sin_mn_targ() |>
      ggplot() + 
      geom_sf(data=scotland_sf) + 
      geom_point(aes(lon, lat, colour=value), size=2.5, shape=1, stroke=1.5) + 
      geom_point(data=val_sin_mn_targ() |> filter(sin==selectedSIN_eval()$sin[1]),
                 aes(lon, lat), colour="blue", size=6, shape=1, stroke=1) +
      scale_colour_viridis_c("Skill score", limits=c(0, 1), na.value="black",
                             breaks=c(0, 0.5, 1),
                             labels=c("(no info)\n0\n", "0.5", "\n1\n(perfect)")) +
      scale_x_continuous(limits=map_ranges$x, expand=c(0, 0)) +
      scale_y_continuous(limits=map_ranges$y, expand=c(0, 0)) +
      theme(axis.title=element_blank())
  },
  res=100)
  
  output$eval.sin.plot <- renderPlot({
    val_sin_mn_targ() |>
      ggplot(aes(value)) +
      geom_dots(side="bottom", layout="swarm") +
      stat_histinterval(breaks=seq(-1, 2, by=0.1)) +
      geom_vline(xintercept=filter(val_sin_mn_targ(), 
                                   sin==selectedSIN_eval()$sin[1])$value,
                 colour="blue") +
      scale_x_continuous("Skill score", limits=c(0, 1), oob=scales::oob_keep,
                         breaks=c(0, 0.5, 1),
                         labels=c("0\n(no info)", "0.5", "1\n(perfect)")) +
      theme(axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            legend.position="inside",
            legend.background=element_blank(),
            legend.position.inside=c(0.175, 0.815)) +
      ggtitle(paste0(targ()$fig_long, ": ", selectedSIN_eval()$sin[1]))
  },
  res=100)
  
  
}