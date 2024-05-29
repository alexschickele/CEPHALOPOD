#' =============================================================================
#' @name query_bio_wrapper
#' @description wrapper around functions to extract biological data according
#' to a user defined set of parameters, among the available species in list_bio.
#' The extracted data is formatted to be directly usable by the models available
#' in this workbench.
#' @param FOLDER_NAME name of the corresponding work folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic annotation
#' for each target species.
#' @return Saves the output in a QUERY.RData file

query_bio_wrapper <- function(FOLDER_NAME = NULL,
                              SUBFOLDER_NAME = NULL){
  set.seed(123)
  # --- 1. Initialize function
  # --- 1.1. Start logs - new file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "wt"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_bio_wrapper ********************"))
  # --- 1.3. Load the run metadata in the CALL.RData object
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  # --- 1.4. Load the parallel node metadata in the QUERY.RData object
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 2. Redirection to the OBIS query
  # For occurrence source data
  if(CALL$DATA_SOURCE == "occurrence"){

    # --- 2.1. Run function
    QUERY <- query_occurrence(FOLDER_NAME = FOLDER_NAME,
                              QUERY = QUERY)
  } # End ATLANTECO redirection

  # --- 3. Redirection to the ATLANTECO query
  # For continuous  source data
  if(CALL$DATA_SOURCE == "biomass"){

    # --- 3.1. Run function
    QUERY <- query_biomass(FOLDER_NAME = FOLDER_NAME,
                             QUERY = QUERY)
  } # End ATLANTECO redirection

  # --- 4. Redirection to the MGNIFY query
  if(CALL$DATA_SOURCE == "omic"){

    # --- 4.2. Run function
    QUERY <- query_omic(FOLDER_NAME = FOLDER_NAME,
                          QUERY = QUERY)
  } # End MGNIFY redirection

  # --- 5. Redirection to the CUSTOM query
  if(CALL$DATA_SOURCE != "omic" & CALL$DATA_SOURCE != "biomass" & CALL$DATA_SOURCE != "occurrence"){
    # --- 5.1. Load function
    source(file = paste0(project_wd, "/code/03d_query_custom.R"))

    # --- 5.2. Run function
    QUERY <- query_custom(FOLDER_NAME = FOLDER_NAME,
                          QUERY = QUERY)
  } # End MGNIFY redirection

  # --- 6. Wrap up and save
  # --- 6.1 Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/QUERY.RData"))
  # --- 6.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 6.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
