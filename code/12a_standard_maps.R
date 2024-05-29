#' =============================================================================
#' @name standard_maps
#' @description Simple function for building standard mean and SD maps from
#' projection data
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param MONTH list corresponding to the groups of month to computes maps from.
#' @return plots mean and uncertainty maps per model or ensemble

standard_maps <- function(FOLDER_NAME = NULL,
                          SUBFOLDER_NAME = NULL,
                          MONTH = list(c(10,11,12,1,2,3),
                                       4:9)){
  if(!exists("SEED")){
    SEED <- 123
  }
  set.seed(SEED)
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
                       message(paste(Sys.time(), "******************** START : standard_maps ********************"))

  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))

  # --- 1.3. Check for projections
  if((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE){
    message("No validated algorithms to display projections from")
    # Stop logs
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  }

  # --- 1.3. Initialize loop over algorithm or targets for later
  if(CALL$DATA_TYPE == "proportions"){loop_over <- 1:ncol(QUERY$Y)
  }else if(CALL$FAST == FALSE){loop_over <- CALL$HP$MODEL_LIST
  }else {loop_over <- MODEL$MODEL_LIST}

  # --- 1.4. Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/05_standard_maps.pdf"))

  # --- 1.5. Set initial plot layout & requirements
  par(mfrow = c(4,3), mar = c(2,2,7,1))
  r0 <- CALL$ENV_DATA[[1]][[1]]

  # --- 1.6 Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA

  # --- 2. Plot the quality checks
  # --- 2.1. Compute the recommendation table
  rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
  traffic_col <- rep(rec$COL, each = 4)
  traffic_val <- rec[,1:4] %>% as.matrix() %>% t() %>% c()
  traffic_col[which(traffic_val == 0)] <- "white"
  # --- 2.2. Plot the traffic lights and recommendations
  plot.new()
  if(CALL$DATA_TYPE != "proportions"){mtext(paste("QUALITY CHECK \n", QUERY$annotations$scientificname, "\n ID:", QUERY$annotations$worms_id))}
  if(CALL$DATA_TYPE == "proportions"){mtext("QUALITY CHECK \n for proportions")}
  par(mar = c(1,3,7,1), xpd = NA)
  plot(x = rep(1:4, nrow(rec)), y = rep(nrow(rec):1, each = 4), axes = FALSE, cex = 3,
       xlim = c(0,5), ylim = c(0,nrow(rec)+1), ylab = "", xlab = "",
       pch = 21, col = "black", bg = traffic_col)
  axis(side = 3, at = 1:4, labels = c("A priori \n var. imp.","Perdictive \n performance","Cumulative \n var. imp.","Projection \n uncertainty"), tick = FALSE, line = NA, cex.axis = 1, las = 2)
  axis(side = 2, at = nrow(rec):1, labels = rownames(rec), tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
  axis(side = 4, at = nrow(rec):1, labels = rec$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
  abline(h = 0) # PDF wide separator
  plot.new()
  par(mar = c(2,2,4,1))

  # --- 3. Plot the legends
  # --- 3.1. biomass or habitat suitability legend
  hsi_pal <- inferno_pal(100)
  plot.new()
  # --- 3.1.1. Extract the plot true scale (average max across bootstrap)
  if(CALL$DATA_TYPE == "continuous"){
    plot_scale <- lapply(loop_over, FUN = function(z){
        z = MODEL[[z]]$proj$y_hat %>% apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    }) %>%
      unlist() %>%
      quantile(0.95, na.rm = TRUE)
  } else {
    plot_scale <- 1
  }

  # --- 3.1.2. Plot the colorbar
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.7,
                col = hsi_pal, border = "black")
  # --- 3.1.3. Apply the scale to the axis caption
  # This is only informative and the raster will be rescaled by the maximum
  axis(side = 1, at = seq(0, 1, length.out = 5), labels = round(seq(0, 1 * plot_scale, length.out = 5), 2))
  text(x = 0.5, y = 0.3, "Habitat Suitability Index", adj = 0.5)

  # --- 3.2. Observation vs 75% quartile
  plot.new()
  points(x = 0.1, y = 0.4, pch = 22, col = "black", bg = "gray80", cex = 5)
  text(x = 0.2, y = 0.4, "Q75 Habitat Suitability Index", pos = 4)
  points(x = 0.1, y = 0.1, pch = 20, col = "black", cex = 2)
  text(x = 0.2, y = 0.1, "Observation", pos = 4)
  # --- 3.3. MESS x SD 2D color scale
  par(mar = c(5,5,3,2), xpd = FALSE)
  bivar_pal <- colmat(nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Standard deviation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = round(seq(0, 1 * plot_scale*0.5, length.out = 6), 2))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(4,2.5,3,1))

  # --- 4. Build model-level outputs
  for(i in loop_over){
    # --- 4.1. Compute the different layers
    # --- 4.1.1. Extract projection data for the iteration
    if(CALL$DATA_TYPE == "proportions"){val_raw <- MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
    }else{val_raw <- MODEL[[i]][["proj"]][["y_hat"]]}

    # --- 4.2. Loop over month for maps
    for(j in 1:length(MONTH)){
      m <- MONTH[[j]]
      # --- 4.2.1. Mean value
      # Rescaled by the maximum to match the colorbar
      val <- apply(val_raw[,,m], 1, function(x)(x = mean(x, na.rm = TRUE)))
      r_m <- r0 %>% setValues(val / plot_scale)
      r_m[r_m>1] <- 1 # set the maximum at Q95

      # --- 4.2.2. Coefficient of variation
      # Computes mean SD across bootstrap and than average across month
      if(length(m) > 1){
        val <- apply(val_raw[,,m], c(1,3), function(x)(x = sd(x, na.rm = TRUE))) %>%
          apply(1, function(x)(x = mean(x, na.rm = TRUE)))
      } else {
        val <- apply(val_raw[,,m], 1, function(x)(x = sd(x, na.rm = TRUE)))
      }

      r_sd <- r0 %>% setValues(val)
      r_sd[r_sd > plot_scale*25] <- plot_scale*25
      r_sd[r_sd <= 0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included

      # --- 4.2.3. MESS
      r_mess <- QUERY$MESS[[m]]*-1
      if(nlayers(r_mess) > 1){r_mess <- calc(r_mess, mean, na.rm = TRUE)}
      r_mess[r_mess<0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included
      r_mess[r_mess>100] <- 100

      # --- 4.3. Plot the corresponding maps
      # --- 4.3.1. Plot the habitat suitability map
      if(CALL$DATA_TYPE == "proportions"){
        # Proportions
        tmp <- which(QUERY$annotations$worms_id == colnames(QUERY$Y)[i])
        plot(r_m, col = hsi_pal,
             legend=FALSE, cex.main = 1,
             main = paste("Projection for", QUERY$annotations$scientificname[tmp], "\n Month:", paste(m, collapse = ",")))
        mtext(text = paste("Predictive performance (", names(MODEL[["MBTR"]][["eval"]])[1], ") =", MODEL[["MBTR"]][["eval"]][[1]],
                           "\n Colorbar scale:", format(round(max(getValues(r_m), na.rm = TRUE), 5), scientific = TRUE)),
              side = 1, line = 3, cex = 0.6)
      } else {
        # biomass or habitat suitability values
        plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))],
             legend=FALSE, cex.main = 1,
             main = paste("Projection (", i, ") \n Month:", paste(m, collapse = ",")))
        mtext(text = paste("Predictive performance (", names(MODEL[[i]][["eval"]])[1], ") =", MODEL[[i]][["eval"]][[1]]),
              side = 1, line = 2, cex = 0.7)
      }

      # Land mask
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      box("figure", col="black", lwd = 1)

      # --- 4.3.2. Plot the observations
      # --- 4.3.2.1. Top biomass quartile as contour
      plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
      # --- 4.3.2.2. Land mask
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      box("figure", col="black", lwd = 1)
      # --- 4.3.2.3. Observations location
      if(CALL$DATA_TYPE == "proportions"){tmp <- QUERY$S
      } else {tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]}
      # --- 4.3.2.4. Observation colors
      if(CALL$DATA_TYPE == "continuous"){
        points(tmp$decimallongitude, tmp$decimallatitude,
               col = col_numeric("inferno", domain = range(0:plot_scale, na.rm = TRUE), alpha = 0.2, na.color = hsi_pal[100])(QUERY$Y$measurementvalue), pch = 20, cex = 0.6)
      } else if(CALL$DATA_TYPE == "proportions") {
        points(tmp$decimallongitude, tmp$decimallatitude,
               col = col_numeric("inferno", domain = range(QUERY$Y[,i], na.rm = TRUE), alpha = 0.2)(QUERY$Y[,i]), pch = 20, cex = 0.6)
      } else {
        points(tmp$decimallongitude, tmp$decimallatitude,
               col = "black", pch = 20, cex = 0.6)
      }

      # --- 4.3.3. Plot the uncertainties
      # Land mask
      plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
      box("figure", col="black", lwd = 1)
      # SD x MESS plot
      r <- bivar_map(rasterx = r_sd, rastery = r_mess, colormatrix = bivar_pal,
                     cutx = seq(0,max(getValues(r_sd), na.rm = TRUE), length.out = 101), cuty = 0:100)
      plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
      # Subtitle display of NSD
      if(CALL$DATA_TYPE == "proportions"){
        mtext(text = paste("Projection uncertainty (", names(MODEL[["MBTR"]][["eval"]])[4], ") =", round(MODEL[["MBTR"]][["eval"]][[4]],2)),
              side = 1, line = 2, cex = 0.7)
      } else {
        mtext(text = paste("Projection uncertainty (", names(MODEL[[i]][["eval"]])[4], ") =", round(MODEL[[i]][["eval"]][[4]],2)),
              side = 1, line = 2, cex = 0.7)
      }


    } # End m month loop
  } # End i model loop

  # --- 5. Build ensemble quality checks
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    # --- 5.1. Compute the recommendation table
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE, ENSEMBLE = TRUE)
    traffic_col <- rep(rec$COL, each = 4)
    traffic_val <- rec[,1:4] %>% as.matrix() %>% t() %>% c()
    traffic_col[which(traffic_val == 0)] <- "white"
    # --- 2.2. Plot the traffic lights and recommendations
    plot.new()
    par(mar = c(1,3,7,1), xpd = NA)
    plot(x = rep(1:4, nrow(rec)), y = rep(nrow(rec):1, each = 4), axes = FALSE, cex = 3,
         xlim = c(0,5), ylim = c(0,nrow(rec)+1), ylab = "", xlab = "",
         pch = 21, col = "black", bg = traffic_col)
    axis(side = 3, at = 1:4, labels = c("A priori \n var. imp.","Perdictive \n performance","Cumulative \n var. imp.","Projection \n uncertainty"), tick = FALSE, line = NA, cex.axis = 1, las = 2)
    axis(side = 2, at = nrow(rec):1, labels = rownames(rec), tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
    axis(side = 4, at = nrow(rec):1, labels = rec$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
    plot.new()
    par(mar = c(4,2.5,3,1))
  } # End ensemble QC

  # --- 6. Build ensemble outputs
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    # --- 6.1. Build ensemble array, weighted by evaluation values, re-scale max=1
    y_ens <- NULL
    for(i in MODEL$MODEL_LIST){
      y_ens <- abind(y_ens, MODEL[[i]][["proj"]][["y_hat"]],  along = 2)
    }
    # --- 6.2. Compute the different layers
    # --- 6.2.1. Mean value
    # Rescaled by the maximum to match the colorbar
    val <- apply(y_ens, 1, function(x)(x = mean(x, na.rm = TRUE)))
    r_m <- r0 %>% setValues(val / plot_scale)
    r_m[r_m>1] <- 1 # set maximum at Q95

    # --- 6.2.2. Coefficient of variation
    val <- apply(y_ens, 1, function(x)(x = sd(x, na.rm = TRUE)))
    r_sd <- r0 %>%
      setValues(val)
    r_sd[r_sd > plot_scale*0.25] <- plot_scale*0.25
    r_sd[r_sd <= 0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included

    # --- 6.3. Plot the corresponding maps
    # --- 6.3.1. Plot the biomass
    # biomass or habitat suitability values
    plot(r_m, col = hsi_pal[max(1, floor(r_m@data@min*100)):min(100, ceiling(r_m@data@max*100))], legend=FALSE,
         main = "Projection ( Ensemble )")
    mtext(text = paste("Predictive performance (", names(MODEL[["ENSEMBLE"]][["eval"]])[1], ") =", round(MODEL[["ENSEMBLE"]][["eval"]][[1]],2)),
          side = 1, line = 2, cex = 0.7)
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    box("figure", col="black", lwd = 1)

    # --- 6.3.2. Plot the observations
    # --- 6.3.2.1. Top biomass quartile as contour
    plot(r_m > quantile(r_m, 0.75), col = c("white","gray80"), legend=FALSE, main = "Observations")
    # --- 6.3.2.2. Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    box("figure", col="black", lwd = 1)
    # --- 6.3.2.3. Observations
    tmp <- QUERY$S[which(QUERY$Y$measurementvalue > 0),]
    if(CALL$DATA_TYPE == "continuous"){
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = col_numeric("inferno", domain = range(0:plot_scale, na.rm = TRUE), alpha = 0.2, na.color = hsi_pal[100])(QUERY$Y$measurementvalue), pch = 20, cex = 0.6)
    } else {
      points(tmp$decimallongitude, tmp$decimallatitude,
             col = "black", pch = 20)
    }

    # --- 6.3.3. Plot the uncertainties
    # Land mask
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    box("figure", col="black", lwd = 1)
    # SD x MESS plot
    r <- bivar_map(rasterx = r_sd, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = seq(0,max(getValues(r_sd), na.rm = TRUE), length.out = 101), cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
    mtext(text = paste("Projection uncertainty (", names(MODEL[["ENSEMBLE"]][["eval"]])[3], ") =", round(MODEL[["ENSEMBLE"]][["eval"]][[3]],2)),
          side = 1, line = 2, cex = 0.7)
  } # End if ENSEMBLE = TRUE

  # --- 7. Wrap up and save
  # --- 7.1. Stop PDF saving
  dev.off()
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 7.3. Save model ensemble
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"),
       compress = "gzip", compression_level = 6)
  # --- 7.4. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION

