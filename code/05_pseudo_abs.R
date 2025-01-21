#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' 
#' @param FOLDER_NAME A character string indicating the folder containing the output.
#' @param SUBFOLDER_NAME A list of subfolders to parallelize the process on.
#' 
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' @return Updates the output in a QUERY.RData and CALL.Rdata files
#' 

pseudo_abs <- function(FOLDER_NAME = NULL,
                       SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Start logging process
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : pseudo-absence ********************"))
  
  # --- 1.2. Load metadata and previous query
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))
  # If no pseudo-absence number is set, use the number of presence points
  if (is.null(CALL$NB_PA)) CALL$NB_PA <- nrow(QUERY$S)
  
  # --- 1.3. Base raster and land
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r <- CALL$ENV_DATA[[1]][[1]]
  r[!is.na(r)] <- 0
  
  land <- r
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Double check data type - plot observation locations anyway
  if(CALL$DATA_TYPE != "presence_only"){
    message("No Pseudo-absence generation necessary for this data type")
    
    # --- 2.1. Set up plot file and layout
    pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_observations.pdf"))
    par(mfrow = c(2,2), mar = c(7,3,7,3))
    
    # --- 2.2. Geographical plot for continuous data
    # We plot an artificial land with the scale of the observation first, for the legend
    if(CALL$DATA_TYPE == "continuous"){
      plot_scale <- quantile(QUERY$Y$measurementvalue, 0.95)
      tmp <- (land-9998)*plot_scale
      tmp[1] <- 0
      plot(tmp, col = inferno_pal(100), main = paste("OBSERVATIONS \n samples for ID:", SUBFOLDER_NAME), 
           sub = paste("Nb. of observations after binning:",  nrow(QUERY$S)), type = "continuous")
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      
      # Plot observations with color scale based on measurement value
      tmp <- QUERY$Y$measurementvalue
      tmp[tmp>plot_scale] <- plot_scale
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = col_numeric("inferno", domain = c(0, plot_scale))(tmp),
             pch = 20, cex = 0.6)
      box("figure", col="black", lwd = 1)
    } # end if continuous

    # --- 2.3. Geographical plot for proportion type
    if(CALL$DATA_TYPE == "proportions"){
      plot(land, col = "antiquewhite4", legend=FALSE, main = paste("OBSERVATIONS \n samples for ID:", SUBFOLDER_NAME), 
           sub = paste("Nb. of observations after binning:",  nrow(QUERY$S)), type = "classes")
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = "black", pch = 20, cex = 0.6)
      box("figure", col="black", lwd = 1)
    } # end if proportions
    
    # --- 2.4. Histogram of values
    hist(unlist(QUERY$Y), breaks = 25, col = scales::alpha("black", 0.5), main = "OBSERVATIONS \n Histogram of raw values", 
         xlab = paste(CALL$DATA_TYPE, "observed values"))
    box("figure", col="black", lwd = 1)
    
    # --- 2.5. Longitude and latitude profile
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
    if(is.null(CALL$BACKGROUND_FILTER)){
      # --- 4.2.1. By default as the density of presence within a buffer
      presence <- terra::rasterize(as.matrix(presence), r)
      presence[!is.na(presence)] <- 1
      
      # Compute density
      background <- terra::focal(presence, terra::focalMat(r, 20, type = "Gauss"),
                                 fun = sum, na.rm = TRUE, expand = TRUE)
      
      background <- background / max(terra::values(background), na.rm = TRUE) # Rescale to 1
      background[!is.na(terra::values(presence))] <- NA # Probability at presence location is NA
      
      # Define the weighted background raster
      background <- terra::as.data.frame(terra::mask(background, r), xy = TRUE, na.rm = TRUE)
      colnames(background) <- c("x","y","layer")
    } else {
      # --- 4.2.2. Using a custom background filter if provided
      if(is.data.frame(CALL$BACKGROUND_FILTER)){
        background <- CALL$BACKGROUND_FILTER
      } else {
        background <- terra::rast(CALL$BACKGROUND_FILTER) 
        background <- terra::as.data.frame(terra::mask(background, r), xy = TRUE, na.rm = TRUE)
        colnames(background) <- c("x","y","layer")
      } # if provided a raster or data.frame
    } # if background filter

    # --- 4.3. Additional Environmental distance (MESS)
    if(CALL$PA_ENV_STRATA == TRUE){
      # Compute the mess analysis
      env_strata <- predicts::mess(x = CALL$ENV_DATA[[as.numeric(m)]], v = QUERY$X)
      env_strata <- terra::mask(env_strata, r) # mess has values on all pixels, so we remove land
      env_strata[env_strata > 0] <- 0 # 0 proba better than NA to avoid errors
      env_strata[env_strata == -Inf] <- min(env_strata[!is.na(env_strata) & env_strata != -Inf])
      
      # Convert to a dataframe
      env_strata <- terra::as.data.frame(env_strata*(-1), xy = TRUE, na.rm = TRUE) # As dataframe
      env_strata$mess <- env_strata$mess/max(env_strata$mess, na.rm = TRUE) # Rescale to 1

      # Merge with the current background
      background <- merge(background, env_strata, all.x = TRUE) # use merge instead of join to avoid parallel error
      background[,3] <- apply(background[,3:4], 1, function(x)(x = mean(x, na.rm = TRUE)))
      background <- background[,-4]
      
      rm(env_strata)
      gc()
    } # end if env strata
    
    # --- 4.4. Additional Percentage Random
    if(!is.null(CALL$PER_RANDOM) == TRUE){
      background[,3] <- background[,3]*(1-CALL$PER_RANDOM)+CALL$PER_RANDOM # try to get a fix random PA generation
      background <- background %>% dplyr::filter(!is.na(3)) # NA safe check
    } # if random
    
    # --- 4.5. Sample within the background data
    # Add a resample option if there is not enough background available
    if(nrow(background) < month_freq[m]){
      message(" PSEUDO-ABS : background too small, selection with replacement !")
      tmp <- sample(x = 1:nrow(background), size = month_freq[m], replace = TRUE, prob = background[,3]) # sqr the probability to bias more
    } else {
      tmp <- sample(x = 1:nrow(background), size = month_freq[m], prob = background[,3]) # sqr the probability to bias more
    } # if resample or not
      
    # --- 4.4. Concatenate over month
    xym0 <- background[tmp,1:2] %>% cbind(rep(as.numeric(m), month_freq[m]))
    xym <- rbind(xym, xym0)
    
    rm(presence, background, xym0)
    gc()
  } # m month loop
  
  # --- 5. Fast PDF to check the absences location
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
    tmp <- terra::extract(CALL$ENV_DATA[[xym[i,3]]], xym[i,1:2])[-1] %>% 
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
