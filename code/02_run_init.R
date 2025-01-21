#' =============================================================================
#' @name run_init
#' @description Initializes the output directory and parameters for a run. 
#' 1. Creates an output folder corresponding to the current run (based on species 
#' and selection criteria) where results will be stored.
#' 2. Initializes all parameters required for the run and saves them in a CALL object.
#' 
#' @param FOLDER_NAME A string representing the name of the run folder to be used.
#' @param SP_SELECT A vector of IDs for the species targeted for parallel processing.
#' @param FAST Boolean - Projections and plots are not computed for algorithms that 
#' do not pass the Quality Checks.
#' @param LOAD_FROM A string indicating the path to a previous list_bio object from 
#' another folder. This avoids redoing the lengthy initial list_bio step.
#' 
#' @param SP_SELECT A vector of species identifiers (e.g., Aphia ID for traditional 
#' data - updated if WORMS_CHECK is set to TRUE.
#' @param WORMS_CHECK Boolean - matches SP_SELECT against the WoRMS taxonomy.
#' 
#' @param DATA_TYPE Type of output data, influencing the sub-folder structure and function
#' triggered. Either "presence_only", "continuous" or "proportions".
#' 
#' @param ENV_VAR A list of .nc filenames to extract the environmental features from.
#' @param ENV_PATH A string or vector indicating the path to the root directory where 
#' the .nc files are located.
#' 
#' @param METHOD_PA Method for generating pseudo-absences: "mindist", "cumdist", 
#' or "density".
#' @param NB_PA The number of pseudo-absences to generate.
#' @param DIST_PA If METHOD_PA is "mindist", the distance (in meters) from presences 
#' used to define background data (for expert use only).
#' @param BACKGROUND_FILTER A data frame (2 columns: longitude and latitude) or a path 
#' to a raster object to sample pseudo-absences from environmentally distinct cells. 
#' This filter allows finer tuning for pseudo-absence selection.
#' @param PER_RANDOM The proportion of pseudo-absences sampled randomly in addition 
#' to the specified background.
#' @param PA_ENV_STRATA Pseudo-absences are sampled from environmentally distinct 
#' cells relative to presences using a mess analysis, in addition to the background choice.
#' 
#' @param OUTLIER Boolean - Biological outliers will be removed in the target.
#' @param RFE Boolean - Performs recursive feature elimination for predictor pre-selection.
#' @param ENV_COR A numeric threshold for removing correlated environmental features from 
#' QUERY objects and CALL. If NULL, no removal occurs.
#' 
#' @param NFOLD An integer specifying the number of folds to train each algorithm.
#' @param FOLD_METHOD A method for creating folds - either random split "kfold" 
#' or longitudinal split "lon".
#' 
#' @param MODEL_LIST A list of algorithms to train.
#' @param LEVELS The hyperparameter grid size for each parameter.
#' @param TARGET_TRANSFORMATION A path to a function(x, REVERSE = T/F) for transforming 
#' the target variable. Useful for left-skewed distribution in the target.
#' 
#' @param ENSEMBLE Boolean - Computes an ensemble during evaluation and projection steps.
#' @param N_BOOTSTRAP The number of bootstrap iterations for projections and partial 
#' dependency plots.
#' 
#' @param CUT A numeric value or NULL; if numeric, specifies the quantile (0 to 1) 
#' at which projections are considered zero. Projections without observations are excluded.
#' It avoids projecting hotspots in geographical areas non-connected to observed distribution.
#' 
#' @details Various data transformations between DATA_SOURCE and DATA_TYPE include:
#' - "occurrence" to "presence_only" (default; not recommended for > 50 targets)
#' - "biomass" or "abundance" to "continuous" (default; not recommended for > 50 targets)
#' - "biomass" or "abundance" to "proportions" (not recommended if sampling stations differ)
#' - "MAG" to "proportions" (default; not recommended for > 50 targets)
#' - "MAG" to "continuous" represented as Shannon-diversity
#' - "MAG" to "presence_only" represented as presence-only
#' 
#' @return Creates a subfolder for each species in SP_SELECT and returns a vector of 
#' sub-directory names for subsequent parallel computing.
#' @return Saves all global parameters in a CALL.RData object.


run_init <- function(FOLDER_NAME = "test_run",
                     SP_SELECT = NULL,
                     WORMS_CHECK = TRUE,
                     FAST = FALSE,
                     LOAD_FROM = NULL,
                     DATA_TYPE = NULL,
                     ENV_VAR = NULL,
                     ENV_PATH = "/net/meso/work/nknecht/Masterarbeit/General_Pipeline/Data/environmental_climatologies",
                     METHOD_PA = "density",
                     NB_PA = NULL,
                     PER_RANDOM = 0.25,
                     DIST_PA = NULL,
                     BACKGROUND_FILTER = NULL,
                     PA_ENV_STRATA = FALSE,
                     OUTLIER = TRUE,
                     RFE = TRUE,
                     ENV_COR = 0.8,
                     NFOLD = 3,
                     FOLD_METHOD = "lon",
                     MODEL_LIST = c("GLM","GAM","RF","MLP","SVM","BRT"),
                     LEVELS = 3,
                     TARGET_TRANSFORMATION = NULL,
                     ENSEMBLE = TRUE,
                     N_BOOTSTRAP = 10,
                     CUT = NULL){
  
  # --- 1. Initialization and setup
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Start logs - create file
  # Create a log file to save all text output during the run
  if(is.null(LOAD_FROM)){
    sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                         START = TRUE)
    message(paste(Sys.time(), "******************** START : run_init ********************"))
  }
  
  # --- 1.2. Define the different folder path
  # --- 1.2.1. Work in an existing directory
  in_path <- out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  
  # --- 1.2.2. Duplicate to a new directory
  if(!is.null(LOAD_FROM)){
    in_path <- paste0(project_wd, "/output/", LOAD_FROM)
    out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  }
  
  # --- 1.3. Load parameters
  load(paste0(in_path, "/CALL.RData"))
  
  # --- 2. Create new directory if needed
  # If we want to use an old list_bio without changing the directory
  if(!is.null(LOAD_FROM)){
    if(file.exists(out_path)==TRUE){
      stop("--- This new foldername is already used")
    } else {
      dir.create(out_path)
      sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                           START = TRUE)
      message(paste(Sys.time(), "******************** START : run_init ********************"))
    } # if exists
  } # if LOAD_FROM
  
  # --- 3. Check the taxonomic assignation - WORMS_CHECK
  # Only if TRUE and data_source != proportions // you may want to turn it off if the data are not adapted (e.g., functional, MAGs)
  
  # --- 3.1. Initialize object
  SP_SELECT_INFO <- NULL # Initialize as NULL; it will store info on species if WORMS_CHECK is TRUE
  SP_SELECT <- SP_SELECT[!is.na(SP_SELECT)] # Remove NA values from SP_SELECT
  
  # Check if WORMS_CHECK is enabled and the data source is not "MAG"
  if (WORMS_CHECK && CALL$DATA_SOURCE != "MAG") {
    message("WORMS_CHECK: TRUE - checking for taxonomically unaccepted names and synonyms against the WoRMS taxonomic backbone\n")
    
    # --- 3.2. Extract WoRMS information from SP_SELECT
    SP_SELECT <- as.numeric(SP_SELECT) # Convert SP_SELECT to numeric for worms ID compatibility
    SP_SELECT_INFO <- lapply(SP_SELECT, function(x) worms_check(ID = x, MARINE_ONLY = TRUE))
    
    # Filter out unaccepted names (keep only non-empty results)
    SP_SELECT_INFO <- SP_SELECT_INFO[lengths(SP_SELECT_INFO) > 0] 
    
    # --- 3.3. Identify duplicates
    valid_ids <- sapply(SP_SELECT_INFO, function(x) x$VALID) # Extract valid names
    duplicates <- duplicated(valid_ids) # Identify duplicates
    duplicates_indices <- which(duplicates) # Indices of duplicate entries
    
    # --- 3.4. Merge the duplicates
    if (length(duplicates_indices) > 0) {
      while (length(duplicates_indices) > 0) {
        # --- 3.4.1. Find indices of duplicated VALID vector
        indices <- which(valid_ids == valid_ids[duplicates_indices[1]]) # Find all indices of the first duplicate
        message(paste("WORMS: subfolder", paste(valid_ids[indices], collapse = " & "), "are duplicates - merging them"))
        
        # --- 3.4.2. Merge SYNONYM into the first element of indices
        SP_SELECT_INFO[[indices[1]]]$SYNONYM <- unique(c(SP_SELECT_INFO[[indices[1]]]$SYNONYM, unlist(lapply(indices, function(idx) SP_SELECT_INFO[[idx]]$SYNONYM))))
        
        # --- 3.4.3. Remove merged entries
        SP_SELECT_INFO[indices[-1]] <- NULL # Use NULL instead of NA to free up memory
        
        # --- 3.4.4. Update duplicates_indices
        duplicates_indices <- duplicates_indices[duplicates_indices != indices[-1]] # Remove processed duplicates
      } # end while
    } # end if duplicates
    
    # --- 3.5. Add nice names and update SP_SELECT
    names(SP_SELECT_INFO) <- sapply(SP_SELECT_INFO, function(x) x[["VALID"]]) # Set names based on VALID names
    SP_SELECT <- names(SP_SELECT_INFO) # Update SP_SELECT to contain valid names
  } # end if
  
  # --- 4. Create species sub-directories
  # Named by the content of the "worms_id" column OR "proportions" OR "shannon" if relevant. 
  # Contains a QUERY object with the corresponding species selection to consider in each sub folder.
  
  if(DATA_TYPE == "proportions"){
    dir.create(paste0(out_path, "/proportions"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/proportions/QUERY.RData"))
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "MAG"){
    dir.create(paste0(out_path, "/shannon"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/shannon/QUERY.RData"))
  } else {
    for(i in SP_SELECT){
      dir.create(paste0(out_path, "/", i))
      QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = i))
      save(QUERY, file = paste0(out_path, "/", i, "/QUERY.RData"))
    }
  } # if DATA_TYPE and SOURCE
  
  # --- 5. Extract environmental raster from NCDF
  # --- 5.1. Get the list of files from the ENV_PATH
  list_nc <- list.files(ENV_PATH, pattern = "\\.nc$", full.names = TRUE)
  
  # --- 5.2. Get the list of variables to include and exclude
  var_out <- grep("!", ENV_VAR, value = TRUE) %>% gsub("!", "", .)
  var_in <- grep("!", ENV_VAR, value = TRUE, invert = TRUE)
  
  # Determine which variables to keep
  ENV_VAR <- sub("\\.nc$", "", basename(list_nc))
  if (length(var_in) > 0) ENV_VAR <- ENV_VAR[ENV_VAR %in% var_in]
  if (length(var_out) > 0) ENV_VAR <- ENV_VAR[!ENV_VAR %in% var_out]
  
  # --- 5.3. Provide a list per month, containing raster stack of all variables
  # Also plot the native predictor distribution for each variable
  
  ENV_DATA <- vector("list", 12)
  pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/Available_predictors.pdf"))
  par(mfrow = c(4,2), oma = c(0,0,2,0))
  
  for (m in 1:12) {
    # --- 5.3.1. Extract data from the .nc files using nc_to_raster (now with terra)
    stack_month <- lapply(list_nc[str_detect(list_nc,paste(ENV_VAR,collapse="|"))], 
                          function(x) {
                            result <- nc_to_raster(MONTH = m,
                                                   NC = x,
                                                   MIN_DEPTH = CALL$SAMPLE_SELECT$FEATURE_MIN_DEPTH,
                                                   MAX_DEPTH = CALL$SAMPLE_SELECT$FEATURE_MAX_DEPTH)
                            return(result)
                          })
    
    # Extract dimension information
    dim_info <- lapply(stack_month, function(x) {
      paste0(names(x[[2]]), "(", x[[2]], ")") %>% paste(collapse = " ; ")
    })
    
    # Extract the actual raster from the nc_to_raster output
    stack_month <- lapply(stack_month, function(x) x[[1]])
    
    # Convert list of SpatRasters to a SpatRaster collection (equivalent to stack in raster)
    stack_month <- rast(stack_month)
    
    # --- 5.3.2. Assign proper names to the raster layers
    names(stack_month) <- ENV_VAR
    
    # --- 5.3.3. Plot the native range in a PDF
    # Enables the user to see if one variable has a very different range than others
    # Loops over each variable to ensure all are plotted by the terra library
    for(f in 1:nlyr(stack_month)){
      plot(stack_month[[f]], main = paste(names(stack_month)[f], "\n month n°:", m, "-", unlist(dim_info)),
           col = viridis::viridis(100), cex.main = 0.6) 
    } # end f feature layer loop
    
    # --- 5.3.4. Synchronize NAs across predictors and assign to list
    ENV_DATA[[m]] <- synchroniseNA(stack_month)
    message(paste(Sys.time(), "--- ENV. STACK : month", m, "done \t"))
  }  # End month loop
  
  dev.off()
  
  # --- 5.4. Synchronize NAs across months
  # --- 5.4.1. Extract the base layer (first raster of each month)
  base_r <- lapply(ENV_DATA, function(x) x[[1]]) %>% 
    rast() %>% 
    synchroniseNA() %>% 
    .[[1]]
  
  # Normalize base raster to be used for NA alignment
  base_r <- base_r / base_r
  
  # --- 5.4.2. Multiply each monthly stack by the base raster to align NAs
  ENV_DATA <- lapply(ENV_DATA, function(x) {
    x <- x * base_r  # Apply the normalization
    names(x) <- ENV_VAR  # Ensure layer names are preserved
    x <- terra::wrap(x) # Ensure that we can save it in a RData file
    return(x)
  })
  
  # --- 6. Update CALL object
  # --- 6.1. Append CALL with ENV_DATA
  CALL[["ENV_DATA"]] <- ENV_DATA
  
  # --- 6.2. Append CALL with DATA_TYPE
  # Define it the same as DATA_SOURCE if left blank
  if(is.null(DATA_TYPE)){DATA_TYPE <- CALL$DATA_SOURCE}
  CALL[["DATA_TYPE"]] <- DATA_TYPE
  
  # --- 6.3. Append CALL with all other objects
  CALL[["SP_SELECT"]] <- SP_SELECT
  CALL[["WORMS_CHECK"]] <- WORMS_CHECK
  CALL[["SP_SELECT_INFO"]] <- SP_SELECT_INFO
  CALL[["FAST"]] <- FAST
  CALL[["ENV_VAR"]] <- ENV_VAR
  CALL[["ENV_PATH"]] <- ENV_PATH
  CALL[["NB_PA"]] <- NB_PA
  CALL[["PER_RANDOM"]] <- PER_RANDOM
  CALL[["BACKGROUND_FILTER"]] <- BACKGROUND_FILTER
  CALL[["PA_ENV_STRATA"]] <- PA_ENV_STRATA
  CALL[["OUTLIER"]] <- OUTLIER
  CALL[["RFE"]] <- RFE
  CALL[["ENV_COR"]] <- ENV_COR
  CALL[["NFOLD"]] <- NFOLD
  CALL[["FOLD_METHOD"]] <- FOLD_METHOD
  CALL[["MODEL_LIST"]] <- MODEL_LIST
  CALL[["LEVELS"]] <- LEVELS
  CALL[["TARGET_TRANSFORMATION"]] <- TARGET_TRANSFORMATION
  CALL[["ENSEMBLE"]] <- ENSEMBLE
  CALL[["N_BOOTSTRAP"]] <- N_BOOTSTRAP
  CALL[["CUT"]] <- CUT
  
  # --- 7. Wrap up and save
  # --- 7.1. Save file(s)
  save(CALL, file = paste0(out_path, "/CALL.RData"))
  
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
  # --- 7.3. List sub-directories to return
  # To be returned as a vector, further used to define parallel runs
  if(DATA_TYPE == "proportions"){
    parallel <- "proportions"
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "MAG"){
    parallel <- "shannon"
  } else {
    parallel <- CALL$SP_SELECT
  }
  return(parallel)
  
} # END FUNCTOIN

