#' =============================================================================
#' @name diversity_maps
#' @description Function extracting all projections by species and bootstrap.
#' Then computes a set of diversity metrics and plots
#' @param FOLDER_NAME name of the corresponding folder
#' @param HILL a value of Hill number to use for diversity computation
#' @param MONTH a list of month to concatenate together
#' @return plots mean and uncertainty maps per diversity metric
#' @return a diversity and mess object per projection, species and bootstrap.
#' @return a .nc matching the EMODnet standards with the same informations as 
#' above. - Saved in FOLDERNAME.

diversity_maps <- function(FOLDER_NAME = NULL,
                           HILL = 1,
                           MONTH = list(c(10,11,12,1,2,3),
                                        4:9)){

  # --- 1. Global parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # Baseline raster
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r0 <- CALL$ENV_DATA[[1]][[1]]
  
  # Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Predict first - assemble later estimation
  message(paste0(Sys.time(), "--- DIVERSITY : build the ensembles - START"))
  
  # --- 2.1. Which subfolder list and which model in the ensemble
  all_files <- list.files(paste0(project_wd, "/output/", FOLDER_NAME), recursive = TRUE)
  model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]
  ensemble_files <- mclapply(model_files, function(x){
    memory_cleanup() # low memory use
    
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/MODEL.RData"))
    if(length(MODEL$MODEL_LIST) >= 1){
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/QUERY.RData"))
      return(list(SUBFOLDER_NAME = x, MODEL_LIST = MODEL$MODEL_LIST, Y = QUERY$Y, MESS = QUERY$MESS, REC = MODEL$recommandations))
    } else {
      return(NULL)
    } # if model list
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.preschedule = FALSE) %>% 
    .[lengths(.) != 0] %>% 
    .[grep("Error", ., invert = TRUE)] # to exclude any API error or else
  
  
  # --- 2.2. Loop over the files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - loop over files"))
  tmp <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    
    # --- 2.2.1. Load MODEL files
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
    
    # --- 2.2.2. Extract projections in a matrix
    # If there is more than 1 algorithm, we extract and average across algorithm
    # Output matrix is cell x bootstrap x month
    if(length(x$MODEL_LIST) > 1){
      m <- lapply(x$MODEL_LIST, function(y){
        MODEL[[y]][["proj"]]$y_hat
      }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
    } else {
      m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
    } # end if
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.cleanup = TRUE)
  
  # --- 2.3. Re-arrange the array
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - format to array"))
  data <- tmp %>% 
    abind(along = 4) %>% 
    aperm(c(1,4,2,3))
  
  # --- 3. Initialize diversity computing
  full_cell_id <- which(!is.na(data[,1,1,1]))
  cell_vector <- terra::values(r0) %>% as.numeric()
  
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/06_hill_diversity_",HILL,".pdf"))

  for(m in seq_along(MONTH)){
    # --- 3.1. Format the data
    tmp <- data[,MONTH[[m]],,]
    tmp <- apply(tmp, c(1,4), function(x)(x = mean(x, na.rm = TRUE)))
    df0 <- tmp[full_cell_id,] %>% as.data.frame() # subset full cells
    
    # --- 3.2. Compute the diversity
    div0 <- hill_taxa(comm = df0, q = HILL) # compute diversity
    div <- cell_vector # base
    div[full_cell_id] <- div0 # fill non empty cells
    
    # --- 3.3. Plot
    r <- terra::setValues(r0, div)
    plot(r, col = parula_pal(100), main = paste("Hill diversity of", HILL, "(nb. of species)", "\n Month:", paste(MONTH[[m]], collapse = ",")))
    plot(land, col = "black", add = TRUE, legend = FALSE)
    
  } # end month loop
  
  dev.off()
  
  
  
} # END FUNCTION


