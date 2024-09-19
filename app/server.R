server <- function(input, output, session) {
  
  # --- Set maximum upload size ---
  options(shiny.maxRequestSize = 500 * 1024^2)  # Set maximum upload size to 500 MB
  # End of Set maximum upload size
  
  # --- Define default NetCDF file path ---
  default_nc_path <- "./output_diversity.nc"
  # End of Define default NetCDF file path
  
  # --- Function to extract file name without extension ---
  extract_file_name <- function(file_path) {
    tools::file_path_sans_ext(basename(file_path))
  }
  # End of Function to extract file name without extension
  
  # --- Function to check file name for data type ---
  check_file_name <- function(file_name, data_type) {
    if (data_type == "Ensemble") {
      return(grepl("ensemble", tolower(file_name)) || grepl("ens", tolower(file_name)))
    } else if (data_type == "Diversity") {
      return(grepl("diversity", tolower(file_name)) || grepl("diverse", tolower(file_name)))
    }
    return(FALSE)
  }
  # End of Function to check file name for data type
  
  # --- Define custom color palette ---
  pal <- mako_pal(100)
  # End of Define custom color palette
  
  # --- Dynamically render UI for highlight type ---
  output$highlight_type_ui <- renderUI({
    req(input$highlight_high_sd)  # Show only if highlight_high_sd is checked
    if (input$highlight_high_sd) {
      selectInput("highlight_type", "Highlight Type", 
                  choices = c("Red pixel" = "red", "Area Outline" = "hatched"))
    }
  })
  # End of Dynamically render UI for highlight type
  
  # --- Reset file input when data_type changes ---
  observeEvent(input$data_type, {
    output$ncfile <- renderUI({
      fileInput("ncfile", "Choose the appropriate NetCDF File", accept = ".nc")
    })
  })
  # End of Reset file input when data_type changes
  
  # --- Reactive expression for NetCDF file ---
  nc_file <- reactive({
    if (is.null(input$ncfile)) {
      # Load default NetCDF file if no file is uploaded
      showModal(modalDialog(
        title ="Welcome in the visualization app!",
        paste("Default file loaded, first choose the appropriate data type and then load your own nc file."),
        footer = NULL,
        easyClose = TRUE
      ))
      print(default_nc_path)
      return(open.nc(default_nc_path))
    } else {
      return(open.nc(input$ncfile$datapath))
    }
  })
  # End of Reactive expression for NetCDF file
  
  # --- Function to read NetCDF file and update UI ---
  observeEvent(input$data_type, {
    req(input$data_type)
    nc <- nc_file()
    
    # Check if the file name contains appropriate keywords for the chosen data type
    file_name <- extract_file_name(if (is.null(input$ncfile)) default_nc_path else input$ncfile$name)
    data_type <- input$data_type
    if (!check_file_name(file_name, data_type)) {
      showModal(modalDialog(
        title = "Warning",
        paste("Be careful! You selected", data_type, "but it seems you loaded the wrong NetCDF file."),
        footer = NULL,
        easyClose = TRUE
      ))
      return()
    }
    
    output$feature_select <- renderUI({
      req(input$data_type)
      
      # Fetch choices based on feature_names attribute of mean_values
      feature_names <- att.get.nc(nc, variable = "mean_values", attribute = "feature_names")
      feature_names <- unlist(strsplit(feature_names, ","))
      
      if (input$data_type == "Diversity") {
        selectInput("feature", "3. Select Diversity Type", choices = setNames(1:length(feature_names), feature_names))
      } else if (input$data_type == "Ensemble") {
        selectInput("feature", "3. Select Species", choices = setNames(1:length(feature_names), feature_names))
      }
    })
  })
  # End of Function to read NetCDF file and update UI
  
  # --- Render UI for Ensemble Box ---
  output$ensemble_box <- renderUI({
    if (input$data_type == "Ensemble") {
      box(
        width = 12,
        p("All CEPHALOPOD outputs follow a strict quality check that includes", tags$br(),
          tags$b("- Environmental Feature Importance:"), "Assesses the potential of the set of environmental features in explaining the biological target", tags$br(),
          tags$b("- Predictive performance:"), "Evaluates the accuracy of each algorithm in predicting the observed target", tags$br(),
          tags$b("- Feature importance:"), "Measures how well key features emerge in explaining the observed pattern", tags$br(),
          tags$b("- Projection uncertainty:"), "Ensures low uncertainty in habitat suitability projections", tags$br(),
          "These quality checks are displayed in a comprehensive traffic-light system, that indicate their success, with red, yellow or green lights showing the number of quality checks passed. The white color indicate which quality checks did not pass."),
        tags$hr(style="border-color: black;"),
        plotOutput("ensemble_plot")
      )
    } else {
      box(
        width = 12,
        p("All CEPHALOPOD outputs follow a strict quality check that includes", tags$br(),
          tags$b("- Environmental Feature Importance:"), "Assesses the potential of the set of environmental features in explaining the biological target", tags$br(),
          tags$b("- Predictive performance:"), "Evaluates the accuracy of each algorithm in predicting the observed target", tags$br(),
          tags$b("- Feature importance:"), "Measures how well key features emerge in explaining the observed pattern", tags$br(),
          tags$b("- Projection uncertainty:"), "Ensures low uncertainty in habitat suitability projections", tags$br(),
          "These quality checks are only displayed at the taxa-level. The diversity estimate above only integrates taxa-level ensemble that successfully passed all quality checks."),
        tags$hr(style="border-color: black;")
      )
    }
  })
  # End of Render UI for Ensemble Box
  
  # --- Event handler for update button and highlight_high_sd checkbox ---
  observeEvent(c(input$update, input$highlight_high_sd), {
    req(input$data_type, input$feature, input$time, input$variable)
    
    nc <- nc_file()
    
    lon <- var.get.nc(nc, "longitude")
    lat <- var.get.nc(nc, "latitude")
    
    nlon <- length(lon)
    nlat <- length(lat)
    
    # Adjust the time variable to use the appropriate index for annual mean
    time_index <- ifelse(input$time == 13, 13, input$time)
    
    # Get the selected variable data for mean_values
    mean_values <- var.get.nc(nc, "mean_values", 
                              start = c(1, 1, time_index, as.numeric(input$feature)),
                              count = c(nlon, nlat, 1, 1))
    
    # Get the selected variable data for sd_values
    sd_values <- var.get.nc(nc, "sd_values", 
                            start = c(1, 1, time_index, as.numeric(input$feature)),
                            count = c(nlon, nlat, 1, 1))
    
    raster_map_mean <- create_raster_from_nc(mean_values, lon, lat)
    raster_map_sd <- create_raster_from_nc(sd_values, lon, lat)
    
    raster_map <- raster_map_mean  # Default to mean_values initially
    if (input$variable == 'sd_values') {
      raster_map <- raster_map_sd  # Switch to sd_values if selected
    }
    
    raster_map[raster_map == 0] <- NA
    
    if (input$contour_checkbox) {
      # Generate contours
      contours <- tryCatch({
        raster::rasterToContour(raster_map, nlevels = 4)
      }, error = function(e) {
        NULL
      })
    } else {
      contours <- NULL
    }
    
    # Create land mask where raster_map is NA or 0
    land_mask <- is.na(raster_map) | raster_map == 0
    # Create land raster: Non-zero at land, NA elsewhere
    land_raster <- raster_map
    land_raster[land_mask] <- -9999
    land_raster[!land_mask] <- NA  # Set non-land values to NA
    message(extent(land_raster))
    
    output$map <- renderLeaflet({
      
      base_map <- leaflet(options = leafletOptions(
        crs = leafletCRS(crsClass = "L.CRS.EPSG4326"),
        minZoom = 1,  # Set minimum zoom level
        maxZoom = 5  # Set maximum zoom level
      )) %>%
        setMaxBounds(lng1 = -180, lat1 = -90, lng2 = 180, lat2 = 90) %>%
        fitBounds(lng1 = -180, lat1 = -90, lng2 = 180, lat2 = 90)
      
      if (input$highlight_high_sd) {
        # Calculate the threshold for high SD values
        sd_threshold <- 0.5
        high_sd_raster <- raster_map_sd
        high_sd_raster[values(raster_map_sd) > sd_threshold] <- 999
        
        if (!is.null(input$highlight_type)) {
          if (input$highlight_type == "red") {
            high_sd_raster[values(raster_map_sd) <= sd_threshold] <- NA
            
            base_map %>%
              addRasterImage(raster_map, colors = pal, opacity = 0.7) %>%
              # addRasterImage(land_raster, colors = "black") %>%
              {
                if (!is.null(contours)) {
                  addPolylines(., data = contours, color = "yellow", weight = 2)
                } else {
                  .
                }
              } %>%
              addRasterImage(high_sd_raster, colors = 'red', opacity = 0.8)
            
          } else if (input$highlight_type == "hatched") {
            # Add a box around the raster at the minimum value to avoid polygons outside the raster limit
            tmp <- as.matrix(high_sd_raster)
            tmp[1,] <- min(tmp, na.rm = TRUE)
            tmp[nrow(tmp),] <- min(tmp, na.rm = TRUE)
            tmp[,1] <- min(tmp, na.rm = TRUE)
            tmp[,ncol(tmp)] <- min(tmp, na.rm = TRUE)
            high_sd_raster <- setValues(high_sd_raster, tmp)
            
            # Create a contour line at the threshold
            rc0 <- raster::rasterToContour(high_sd_raster > 1, nlevels = 1)
            
            all_lines <- rc0@lines[[1]]@Lines
            
            base_map %>%
              addRasterImage(raster_map, colors = pal, opacity = 0.7) %>%
              # addRasterImage(land_raster, colors = 'black') %>%
              {
                if (!is.null(contours)) {
                  addPolylines(., data = contours, color = "yellow", weight = 2)
                } else {
                  .
                }
              } %>%
              addPolylines(data = rc0, color = "red", weight = 2)  # Add hatched lines for high SD values
            
          }
        }
      } else {
        base_map %>%
          addRasterImage(raster_map, colors = pal, opacity = 0.7, project = FALSE) %>%
          # addRasterImage(land_raster, colors = 'black') %>%
          {
            if (!is.null(contours)) {
              addPolylines(., data = contours, color = "yellow", weight = 2)
            } else {
              .
            }
          }
      }
    })
    
    # Output the legend plot
    output$legend_plot <- renderPlot({
      # Extract min and max values from raster_map, excluding NA values
      valid_values <- values(raster_map)
      valid_values <- valid_values[!is.na(valid_values)]
      min_val <- min(valid_values)
      max_val <- max(valid_values)
      
      par(mar = c(3, 1, 1, 1))  # Adjust margins to provide space for axis labels
      
      # Generate the legend plot using image.plot
      fields::image.plot(
        zlim = c(min_val, max_val),
        legend.only = TRUE,
        col = pal,  # Use custom palette 'pal' here
        horizontal = TRUE,  # Make the legend horizontal
        legend.line = 2.5,  # Adjust legend label position
        legend.cex = 0.8  # Adjust legend text size
      )
      title_str <- paste("Legend:", input$data_type, " values")
      title(main = title_str, cex.main = 1.2, line = -0.5)
    })
  })
  # End of Event handler for update button and highlight_high_sd checkbox
  
  # --- Update the custom time label based on slider input ---
  output$time_label <- renderText({
    month_names <- c("January", "February", "March", "April", "May", "June", 
                     "July", "August", "September", "October", "November", "December")
    if (input$time == 13) {
      "Selected month: Annual Mean"
    } else {
      paste("Selected month:", month_names[input$time])
    }
  })
  # End of Update the custom time label based on slider input
  
  # --- Event handler for ensemble data type and feature selection ---
  observeEvent(c(input$data_type, input$feature), {
    req(input$data_type, input$feature, input$ncfile)  # Ensure necessary inputs are available
    
    nc <- nc_file()
    file_name <- extract_file_name(if (is.null(input$ncfile)) default_nc_path else input$ncfile$name)
    
    if (input$data_type == "Ensemble" & (grepl("ensemble", tolower(file_name)) || grepl("ens", tolower(file_name)))) {
      
      # Fetch the 'reco' attribute for the selected species
      reco <- att.get.nc(nc, variable = "mean_values", attribute = "reco")
      a <- unlist(strsplit(reco, ";"))
      lines_ <- a[as.numeric(input$feature)]
      lines <- strsplit(lines_, "\n")[[1]]
      
      # From the string -> reconstructed dataframe
      df <- data.frame(
        PRE_VIP = integer(length(lines)),
        FIT = integer(length(lines)),
        CUM_VIP = integer(length(lines)),
        DEV = integer(length(lines)),
        Recommandation = character(length(lines)),
        COL = character(length(lines)),
        stringsAsFactors = FALSE
      )
      
      row_names <- character(length(lines))
      
      for (i in seq_along(lines)) {
        parts <- strsplit(lines[i], ": ")[[1]] #isolate model name, present in the beginning of each line
        row_names[i] <- parts[1]  # Set row name
        values <- unlist(strsplit(parts[2], ", "))
        
        df[i, "PRE_VIP"] <- as.integer(values[1])
        df[i, "FIT"] <- as.integer(values[2])
        df[i, "CUM_VIP"] <- as.integer(values[3])
        df[i, "DEV"] <- as.integer(values[4])
        df[i, "Recommandation"] <- values[5]
        df[i, "COL"] <- values[6]
      }
      
      
      rownames(df) <- row_names
      
      # Create traffic_col based on traffic_val
      traffic_col <- rep(df$COL, each = 4)
      traffic_val <- df[, 1:4] %>% as.matrix() %>% t() %>% c()
      traffic_col[which(traffic_val == 0)] <- "white"
      
      # Plot the traffic lights and recommendations
      output$ensemble_plot <- renderPlot({
        plot.new()
        par(mar = c(1, 3, 7, 1), xpd = NA)
        plot(x = rep(1:4, nrow(df)), y = rep(nrow(df):1, each = 4), axes = FALSE, cex = 3,
             xlim = c(0, 5), ylim = c(0, nrow(df) + 1), ylab = "", xlab = "",
             pch = 21, col = "black", bg = traffic_col)
        axis(side = 3, at = 1:4, labels = c("A priori \n var. imp.", "Predictive \n performance", "Cumulative \n var. imp.", "Projection \n uncertainty"),
             tick = FALSE, line = NA, cex.axis = 1, las = 2)
        axis(side = 2, at = nrow(df):1, labels = rownames(df), tick = FALSE, line = -0.5, las = 2, cex.axis = 1)
        axis(side = 4, at = nrow(df):1, labels = df$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
      })
    } else {
      output$ensemble_plot <- renderPlot({
        plot.new()  # Render nothing when data_type is not "Ensemble"
      })
    }
  })
  # End of Event handler for ensemble data type and feature selection
  
  # --- Function to create raster from NetCDF data ---
  create_raster_from_nc <- function(nc_data, lon, lat) {
    nlon <- length(lon)
    nlat <- length(lat)
    raster_data <- matrix(nc_data, nrow = nlat, ncol = nlon, byrow = TRUE)
    r <- raster(raster_data)
    crs(r) <- CRS("+proj=longlat +datum=WGS84")
    extent(r) <- c(min(lon), max(lon), min(lat), max(lat))
    return(r)
  }
  # End of Function to create raster from NetCDF data
  
}
