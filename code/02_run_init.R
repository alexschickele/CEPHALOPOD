#' =============================================================================
#' @name run_init
#' @description (1) initialize an output folder corresponding to the run call (i.e., species
#' and selection criteria) in which the output will be saved.
#' @description (2) initialize all parameters of the run and stores them in
#' a CALL object

#' @param FOLDER_NAME name of the run folder we want to work in
#' @param SP_SELECT vector of IDs, corresponding to the species to parallelize on
#' @param FAST TRUE or FALSE; if TRUE, does not compute projections and plot for algorithms
#' that did not pass the Quality Checks
#' @param LOAD_FROM load a previous list_bio object from another folder to be
#' duplicated in the new FOLDER_NAME. It avoids re-doing the initial list_bio step
#' that can be long for omics data

#' @param SP_SELECT species to run the analysis for, in form of Aphia ID for traditional data
#' and OTU ID for omics - updated if WORMS_CHECK = TRUE
#' @param WORMS_CHECK check for the WoRMS taxonomic backbone and updates SP_SELECT accordingly
#' @param DATA_TYPE the output type of the data, which can influence the sub-folder
#' architecture. See details.
#'
#' @param ENV_VAR a list of .nc files to extract the main variable from, located in ENV_PATH
#' @param ENV_PATH string or vector of path to the root where the .nc are.
#' @param METHOD_PA method of pseudo-absence, either "mindist" or "cumdist" or "density"
#' @param NB_PA number of pseudo-absences to generate
#' @param DIST_PA if METHOD_PA = "mindist", distance from presences (in meters),
#'  from which to define the background data. Expert use only.
#' @param BACKGROUND_FILTER alternative background filter for finer tuning, such
#' as selecting pseudo-absences within the sampled background of a given campaign
#' or instrument deployment. Passed by the user in the form of a 2 column
#' data frame, x = longitude and y = latitude where the pseudo-absences
#' can be sampled. Or a path to a raster object where pseudo-absences are sampled in
#' non NA cells, weighted by the cell values.
#' @param PER_RANDOM ratio of pseudo-absences that are sampled randomly in the background
#'  - in addition to the background choice
#' @param PA_ENV_STRATA pseudo absences are also sampled in environmentally distinct cells from
#' the presences using a mess analysis - in addition to the background choice


#' @param OUTLIER if TRUE, remove outliers
#' @param RFE if true, performs a recursive feature exclusion predictor pre-selection
#' @param ENV_COR numeric, removes the correlated environmental values from the
#' query objects and CALL according to the defined threshold. Else NULL.

#' @param NFOLD number of folds, used defined integer
#' @param FOLD_METHOD method used to create the folds, integer between "kfold"
#' and "lon"

#' @param MODEL_LIST list of algorithms from which to compute hyperparameter
#' selection
#' @param LEVELS maximum number of parameter values to test in each of the grids
#' @param TARGET_TRANSFORMATION path to a function(x, REVERSE = T/F) to transform the target variable

#' @param ENSEMBLE TRUE or FALSE; if TRUE, computes an ensemble at the evaluation and projection steps
#' @param N_BOOTSTRAP number of bootstrap to do for the projections and partial dependency plots
#' @param CUT numeric or NULL; if numeric, quantile (between 0 and 1) at which the projections are considered to be 0
#' Projection patches without observation are then removed.
#' @param PROJ_PATH (optional) path to a environmental raster, potentially
#' different than the one given in the QUERY object. This is the case for
#' supplementary projections in other time and climate scenarios for example.
#' To your own risk and only for expert users !

#' @details Different data transformation between DATA_SOURCE and DATA_TYPE are implemented, including:
#' - "occurrence" to "binary" (default)
#' - "biomass" to "continuous" (default, not recommended for more than 50 targets)
#' - "biomass" to "proportions" (not recommended if the sampling stations are not the same)
#' - "omic" to "proportions" (default, not recommended for more than 50 targets)
#' - "omic" to "continuous" in form of richness
#' - "omic" to "binary" in form of presence-only
#' @return creates one subfolder per SP_SELECT and a vector of sub-directories names
#' used as an argument for parallel computing later
#' @return all global parameters in a CALL.RData object

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
                     TARGET_TRANSFORMATION = "/net/meso/work/aschickele/Bluecloud_WB_local/function/target_transformation_yj_auto.R",
                     ENSEMBLE = TRUE,
                     N_BOOTSTRAP = 10,
                     CUT = 0.1,
                     PROJ_PATH = NULL){

  # --- 1. Initialize function
  # --- 1.1. Start logs - create file
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

  # --- 3. Check the taxonomic assignation
  # Only if TRUE and data_source != proportions // you may want to turn it off if the data are not adapted (e.g., functional, omics)
  # --- 3.1. Initialize object
  SP_SELECT_INFO <- NULL # if FALSE stays like this
  SP_SELECT <- SP_SELECT[!is.na(SP_SELECT)] # remove NA from the list

  if(WORMS_CHECK == TRUE & CALL$DATA_SOURCE != "omic"){
    message("WORMS_CHECK: TRUE - checking for taxonomic unnaccepted names and synonyms against the WoRMS taxonomic backbone \n")
    # --- 3.2. Extract worms information from the SP_SELECT
    SP_SELECT <- as.numeric(SP_SELECT) # numeric for worms ID
    SP_SELECT_INFO <- lapply(SP_SELECT, function(x){
      out = worms_check(ID = x, MARINE_ONLY = TRUE)
    })
    SP_SELECT_INFO <- SP_SELECT_INFO[lengths(SP_SELECT_INFO) > 0] # filter out unnaccepted names

    # --- 3.3. Identify duplicates
    duplicates <- duplicated(sapply(SP_SELECT_INFO, "[[", "VALID")) %>% which(. == TRUE)
    valid_id <- lapply(SP_SELECT_INFO, function(x){x <- x$VALID}) %>% unlist()
    duplicates_id <- which(valid_id %in% valid_id[duplicates])

    # --- 3.4. Merge the duplicates
    if(length(duplicates_id) > 0){
      while(length(duplicates_id) > 0) {
        # --- 3.4.1. Find indices of duplicated VALID vector
        # indices <- which(sapply(SP_SELECT_INFO, function(x) all(x$VALID == SP_SELECT_INFO[[i]]$VALID)) == TRUE)

        indices <- which(valid_id == valid_id[duplicates_id])
        message(paste("WORMS: subfolder", paste(valid_id[indices], collapse = " & "), "are a dupplicates - merging it"))

        # --- 3.4.2. Merge SYNONYM into the first element of indices
        SP_SELECT_INFO[[indices[1]]]$SYNONYM <- unlist(sapply(indices, function(idx) SP_SELECT_INFO[[idx]]$SYNONYM)) %>% as.vector() %>% unique()
        SP_SELECT_INFO[indices[-1]] <- NA

        # --- 3.4.3. Update duplicates_id
        duplicates_id <- duplicates_id[-indices]

      } # end while
    } # if duplicates

    SP_SELECT_INFO <- SP_SELECT_INFO[!is.na(SP_SELECT_INFO)] # filter out old duplicates names

    # --- 3.5. Add nice names and update SP_SELECT
    names(SP_SELECT_INFO) <- lapply(SP_SELECT_INFO, function(x)(x = x[["VALID"]]))
    SP_SELECT <- names(SP_SELECT_INFO)
  } # end if

  # --- 4. Create species sub-directories
  # Named by their Worms ID, OTU ID, or named proportion or richness depending on
  # the data source and type. Contains a QUERY object with the corresponding
  # species selection to consider in each sub folder.
  if(DATA_TYPE == "proportions"){
    dir.create(paste0(out_path, "/proportions"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/proportions/QUERY.RData"))
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
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
  # --- 5.1. Get the list of files
  list_nc <- list.files(ENV_PATH) %>%
    .[grep(pattern = ".nc", x = .)]

  # --- 5.2. Get the list of variables
  var_out <- ENV_VAR %>% .[grep("!", . , invert = FALSE)] %>% gsub("!", "", .)
  var_in <- ENV_VAR %>% .[grep("!", . , invert = TRUE)]

  var_all <- str_sub(list_nc, 1, -4)
  if(length(var_in) != 0){var_all <- var_all[var_all %in% var_in]}
  if(length(var_out) != 0){var_all <- var_all[!c(var_all %in% var_out)]}
  ENV_VAR <- var_all

  # --- 5.3. Provide a list per month, containing raster stack of all variables
  # Add plot of the native predictor distribution
  ENV_DATA <- list()
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/Available_predictors.pdf"))

  for(m in 1:12){
    # --- 5.3.1. Extract the data from the .nc files
    stack_month <- lapply(paste0(ENV_PATH, "/", ENV_VAR, ".nc"),
                          FUN = function(x){x = nc_to_raster(MONTH = m,
                                                             NC = x,
                                                             MIN_DEPTH = CALL$SAMPLE_SELECT$FEATURE_MIN_DEPTH,
                                                             MAX_DEPTH = CALL$SAMPLE_SELECT$FEATURE_MAX_DEPTH)
                          })
    dim_info <- lapply(stack_month, FUN = function(x)(x = x[[2]])) %>%
      lapply(FUN = function(x)(x = paste0(names(x), "(", x, ")")) %>% paste(collapse = " ; "))
    stack_month <- lapply(stack_month, FUN = function(x)(x = x[[1]])) %>%
      raster::stack()

    # --- 5.3.2. Pretty names
    names(stack_month) <- ENV_VAR

    # --- 5.3.3. Plot native range in a PDF
    # Enables user to see if one variables has a very different range than others
    raster::plot(stack_month, main = paste(names(stack_month), "\n month n°:", m, "-", unlist(dim_info)),
         cex.main = 0.8, col = viridis_pal(100), nr = 4, nc = 2)

    # --- 5.3.4. SynchroniseNA across predictors - and assign to list
    ENV_DATA[[m]] <- synchroniseNA(stack_month)
    message(paste(Sys.time(), "--- ENV. STACK : month", m, "done \t"))
  } # end month loop
  dev.off()

  # --- 5.4. SynchroniseNA across month
  # --- 5.4.1. Extract the base layer of each month
  base_r <- lapply(ENV_DATA, FUN = function(x)(x = x[[1]])) %>%
    raster::stack() %>%
    synchroniseNA() %>%
    .[[1]]
  base_r <- base_r/base_r

  # --- 5.4.2. Multiply each month stack by the base
  ENV_DATA <- lapply(ENV_DATA, FUN = function(x){x = x*base_r
                                                 names(x) = ENV_VAR
                                                 return(x)})

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
  CALL[["METHOD_PA"]] <- METHOD_PA
  CALL[["NB_PA"]] <- NB_PA
  CALL[["PER_RANDOM"]] <- PER_RANDOM
  CALL[["DIST_PA"]] <- DIST_PA
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
  CALL[["PROJ_PATH"]] <- PROJ_PATH

  # --- 7. Wrap up and save
  # --- 7.1. Save file(s)
  save(CALL, file = paste0(out_path, "/CALL.RData"))

  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)

  # --- 7.3. List sub-directories to return
  # To be returned as a vector, further used to define parallel runs
  if(DATA_TYPE == "proportions"){
    parallel <- "proportions"
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
    parallel <- "shannon"
  } else {
    parallel <- CALL$SP_SELECT
  }
  return(parallel)

} # END FUNCTOIN

