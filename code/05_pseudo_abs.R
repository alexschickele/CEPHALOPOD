#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

pseudo_abs <- function(FOLDER_NAME = NULL,
                       SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  set.seed(123)
  
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : pseudo_absences ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  if(is.null(CALL$NB_PA)){CALL$NB_PA = nrow(QUERY$S)}
  
  # --- 1.3. Base raster and land
  r <- CALL$ENV_DATA[[1]][[1]]
  r[!is.na(r)] <- 0
  
  land <- r
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Double check data type - plot observation locations anyway
  if(CALL$DATA_TYPE != "presence_only"){
    message("No Pseudo-absence generation necessary for this data type")
    pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_observations.pdf"))
    par(mfrow = c(2,2), mar = c(7,3,7,3))
    
    # --- 2.1.1. Geographical plot for continuous data
    # We plot an artificial land with the scale of the observation first, for the legend
    if(CALL$DATA_TYPE == "continuous"){
      plot_scale <- quantile(QUERY$Y$measurementvalue, 0.95)
      tmp <- (land-9998)*plot_scale
      tmp[1] <- 0
      plot(tmp, col = inferno_pal(100), main = paste("OBSERVATIONS \n samples for ID:", SUBFOLDER_NAME), 
           sub = paste("Nb. of observations after binning:",  nrow(QUERY$S)))
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      # Now we scale the observation values
      tmp <- QUERY$Y$measurementvalue
      tmp[tmp>plot_scale] <- plot_scale
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = col_numeric("inferno", domain = c(0, plot_scale))(tmp),
             pch = 20, cex = 0.6)
      box("figure", col="black", lwd = 1)
    }

    # --- 2.1.2. Geographical plot for proportion type
    if(CALL$DATA_TYPE == "proportions"){
      plot(land, col = "antiquewhite4", legend=FALSE, main = paste("OBSERVATIONS \n samples for ID:", SUBFOLDER_NAME), 
           sub = paste("Nb. of observations after binning:",  nrow(QUERY$S)))
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = "black", pch = 20, cex = 0.6)
      box("figure", col="black", lwd = 1)
    }
    
    # --- 2.2. Histogram of values
    hist(unlist(QUERY$Y), breaks = 25, col = scales::alpha("black", 0.5), main = "OBSERVATIONS \n Histogram of raw values", 
         xlab = paste(CALL$DATA_TYPE, "observed values"))
    box("figure", col="black", lwd = 1)
    
    # --- 2.3. Longitude and latitude profile
    hist(QUERY$S$decimallongitude, breaks = seq(-200,200, 20), col = scales::alpha("black", 0.5), 
         main = "OBSERVATIONS \n Longitudinal spectrum", xlab = "Longitude")
    box("figure", col="black", lwd = 1)
    
    hist(QUERY$S$decimallatitude, breaks = seq(-100,100, 10), col = scales::alpha("black", 0.5), 
         main = "OBSERVATIONS \n Latitudinal spectrum", xlab = "Latitude")
    box("figure", col="black", lwd = 1)
    
    dev.off()
    
    log_sink(FILE = sinkfile, START = FALSE)
    return(SUBFOLDER_NAME)
  } # if datatype !presence_only
  
  # --- 3. Initialize background definition
  # --- 3.1. Extract the monthly bias
  # To know how many pseudo-absences to draw from each month background
  month_freq <- summary(as.factor(QUERY$S$month))
  month_freq <- trunc(month_freq/(sum(month_freq)/CALL$NB_PA))
  
  # --- 3.2. Initialize lon x lat x time df
  xym <- NULL
  
  # --- 4. Extract pseudo-absence
  for(m in names(month_freq)){
    # --- 4.1. Extract presence points
    presence <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude, month) %>% 
      dplyr::filter(month == m) %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 4.2. Background definition
    # Need to redefine base raster here in case of different NA between months
    r <- CALL$ENV_DATA[[as.numeric(m)]][[1]]
    r[!is.na(r)] <- 0
    
    # --- 4.2.1. Based on a geographical distance
    # Defining the background as cells distant from more than n-km from presence
    if(CALL$METHOD_PA == "mindist"){
      # Calculate distance to presences in raster object
      background <- raster::rasterize(presence, r)
      background[background < 1] <- NA
      background <- raster::distance(background)
      
      # Filter by minimum distance to presence
      background <- synchroniseNA(stack(background, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame() %>% 
        dplyr::filter(layer < CALL$DIST_PA & !is.na(layer)) %>% 
        dplyr::select(x, y)
    } # End if mindist
    
    # --- 4.2.2. Biased by cumulative-distance to presence
    if(CALL$METHOD_PA == "cumdist"){
      
      # Compute cumulative distance
      background <- r %>% rasterToPoints() %>% .[,1:2]
      cumdist <- pointDistance(presence, background, lonlat = TRUE) %>% 
        apply(2, sum)
      
      # Define weighted background raster
      background <- r
      background[!is.na(background)] <- cumdist
      background <- (background/max(getValues(background), na.rm = TRUE)-1)*-1
      
      background <- synchroniseNA(stack(background, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame()
    } # End if cumdist
    
    # --- 4.2.3. Density of presence within a buffer
    if(CALL$METHOD_PA == "density"){
      # Prepare base raster
      presence <- raster::rasterize(presence, r)
      presence[!is.na(presence)] <- 1
      
      # Compute density
      focal_w <- focalWeight(r, 20, type = "Gauss")
      dens <- raster::focal(presence, focal_w, fun = function(x){sum(x, na.rm = TRUE)}, pad = TRUE)
      # dens <- raster::focal(presence, focal_w, fun = function(x){mean(x, na.rm = TRUE)}, pad = TRUE)
      dens <- (dens/max(getValues(dens), na.rm = TRUE))
      dens[!is.na(presence)] <- NA
      
      # Define the weighted background raster
      background <- synchroniseNA(stack(dens, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame()
    } # End if density
    
    # --- 4.2.4. Custom background filter
    if(!is.null(CALL$BACKGROUND_FILTER)){
      if(is.data.frame(CALL$BACKGROUND_FILTER)){
        background <- CALL$BACKGROUND_FILTER
      } else {
        background <- raster(CALL$BACKGROUND_FILTER) %>% 
          rasterToPoints() %>% 
          as.data.frame()
      }
    } # END if background_filter TRUE
    
    # --- 4.2.5. Additional Environmental distance (MESS)
    if(CALL$PA_ENV_STRATA == TRUE){
      env_strata <- dismo::mess(x = CALL$ENV_DATA[[as.numeric(m)]], v = QUERY$X)
      env_strata <- synchroniseNA(stack(env_strata, CALL$ENV_DATA[[1]][[1]]))[[1]] # mess has values on all pixels, so we remove land
      env_strata[env_strata > 0] <- 0 # 0 proba better than NA to avoid errors
      env_strata[env_strata == -Inf] <- min(env_strata[!is.na(env_strata) & env_strata != -Inf])
      
      env_strata <- (env_strata*(-1)) %>% 
        rasterToPoints() %>% 
        as.data.frame()
      
      env_strata$mess <- env_strata$mess/max(env_strata$mess, na.rm = TRUE)

      background <- merge(background, env_strata, all.x = TRUE) # use merge instead of join to avoid parallel error
      if(ncol(background) == 4){
        background[,3] <- apply(background[,3:4], 1, function(x)(x = mean(x, na.rm = TRUE)))
        background <- background[,-4]
      }
    } # end if env strata
    
    # --- 4.2.6. Additional Percentage Random
    if(!is.null(CALL$PER_RANDOM) == TRUE & ncol(background) == 3){
      background[,3] <- background[,3]*(1-CALL$PER_RANDOM)+CALL$PER_RANDOM # try to get a fix random PA generation
      background <- background %>% 
        dplyr::filter(!is.na(3))
    } # if random
    
    # --- 4.3. Sample within the background data
    # Add a resample option if there is not enough background available
    if(nrow(background) > 0){
      if(ncol(background) == 3){
        if(nrow(background) < month_freq[m]){
          message(" PSEUDO-ABS : background too small, selection with replacement !")
          tmp <- sample(x = 1:nrow(background), size = month_freq[m], replace = TRUE, prob = background[,3]) # sqr the distance to bias more
        } else {
          tmp <- sample(x = 1:nrow(background), size = month_freq[m], prob = background[,3]) # sqr the distance to bias more
        } 
      } else {
        if(nrow(background) < month_freq[m]){
          message(" PSEUDO-ABS : background too small, selection with replacement !")
          tmp <- sample(x = 1:nrow(background), size = month_freq[m], replace = TRUE)
        } else {
          tmp <- sample(x = 1:nrow(background), size = month_freq[m])
        } 
      } # end if col background
      
      # --- 4.4. Concatenate over month
      xym0 <- background[tmp,1:2] %>% cbind(rep(as.numeric(m), month_freq[m]))
      xym <- rbind(xym, xym0)
      
    } # end if nrow background

  } # m month loop
  
  # --- 5.Fast PDF to check the absences location
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_pseudo_abs.pdf"))
  par(mfrow = c(2,2), mar = c(7,3,7,3))
  
  # --- 5.1. Geographical distribution
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("OBSERVATIONS \n presence samples for ID:", SUBFOLDER_NAME), 
       sub = paste("Nb. of observations after binning:",  nrow(QUERY$S)))
  points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, col = scales::alpha("black", 0.2), pch = 20)
  box("figure", col="black", lwd = 1)
  
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("PSEUDO_ABSENCES \n location for ID:", SUBFOLDER_NAME), 
       sub = paste("NB_PA :", CALL$NB_PA, "// METHOD_PA :", CALL$METHOD_PA))
  points(xym[,1:2], col = scales::alpha("red", 0.2), pch = 20)
  box("figure", col="black", lwd = 1)
  
  # --- 5.2. Longitude and latitude profile
  hist(QUERY$S$decimallongitude, breaks = seq(-200,200, 20), col = scales::alpha("black", 0.5), 
       main = "TARGET \n Longitudinal spectrum", xlab = "Longitude")
  hist(xym$x, breaks = seq(-200,200, 20), col = scales::alpha("red", 0.5), add = TRUE)
  box("figure", col="black", lwd = 1)
  
  hist(QUERY$S$decimallatitude, breaks = seq(-100,100, 10), col = scales::alpha("black", 0.5), 
       main = "TARGET \n Latitudinal spectrum", xlab = "Latitude")
  hist(xym$y, breaks = seq(-100,100, 10), col = scales::alpha("red", 0.5), add = TRUE)
  box("figure", col="black", lwd = 1)
  dev.off()
  
  # --- 6. Append the query
  # --- 6.1. Feature table
  # Extract from the corresponding monthly raster
  X <- NULL
  for(i in 1:nrow(xym)){
    tmp <- raster::extract(CALL$ENV_DATA[[xym[i,3]]], xym[i,1:2]) %>% 
      as.data.frame()
    X <- rbind(X, tmp)
  }
  QUERY$X <- rbind(QUERY$X, X)
  
  # --- 6.2. Target table - replace by 0 and 1's
  QUERY$Y <- data.frame(measurementvalue = c(rep(1, nrow(QUERY$Y)), 
                                             rep(0, nrow(xym))))
  
  # --- 6.3. Sample table
  S <- data.frame(decimallongitude = xym$x,
                  decimallatitude = xym$y,
                  month = xym[,3],
                  measurementtype = "Pseudo-absence") %>% 
    mutate(ID = row_number()+nrow(QUERY$S))
  QUERY$S <- QUERY$S %>% 
    bind_rows(S)

  # --- 7. Wrap up and save
  # --- 7.1. Save QUERY object
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 7.3. Pretty return
  return(SUBFOLDER_NAME)
  
} # END FUNCTION
