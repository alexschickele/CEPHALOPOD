#' =============================================================================
#' @name query_bio_wrapper
#' @description Wrapper around functions to extract biological data according 
#' to a user-defined set of parameters. Extracted data is formatted for direct
#' use in models. Handles occurrence, biomass, abundance, and metagenomics (MAG) data.
#' 
#' @param FOLDER_NAME Name of the main working folder.
#' @param SUBFOLDER_NAME Subfolder name to parallelize on.
#' 
#' @return SUBFOLDER_NAME The function returns the subfolder name.
#' @details The function saves a 'QUERY.RData' file with the processed data.
#' 

query_bio_wrapper <- function(FOLDER_NAME = NULL, SUBFOLDER_NAME = NULL) {
  
  # --- 1. Initialize function and environment
  set.seed(123)  # Ensure reproducibility
  
  # --- 1.1. Start logging (open a new log file for this subfolder)
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/log.txt"), open = "wt"),
    START = TRUE)
  
  message(Sys.time(), "******************** START : query_bio_wrapper ********************")
  
  # --- 1.2. Load essential metadata
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))  # Load metadata for the current run
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))  # Load query-specific metadata
  
  # --- 2. Data Source Redirection
  # Use conditional blocks to call the relevant query function based on the DATA_SOURCE
  
  if (CALL$DATA_SOURCE == "occurrence") {
    # --- 2.1. Occurrence data query (e.g., from OBIS)
    QUERY <- query_occurrence(FOLDER_NAME = FOLDER_NAME, QUERY = QUERY)
    
  } else if (CALL$DATA_SOURCE %in% c("biomass", "abundance")) {
    # --- 3.1. Biomass or abundance data query
    QUERY <- query_abundance_biomass(FOLDER_NAME = FOLDER_NAME, QUERY = QUERY)
    
  } else if (CALL$DATA_SOURCE == "MAG") {
    # --- 4.1. Metagenomics (MAG) query (from MGNIFY)
    QUERY <- query_MAG(FOLDER_NAME = FOLDER_NAME, QUERY = QUERY)
    
  } else {
    # --- 5.1. Custom data query (if the data source is not one of the predefined types)
    source(paste0(project_wd, "/code/03d_query_custom.R"))  # Source the custom query script
    
    # --- 5.2. Run the custom query
    QUERY <- query_custom(FOLDER_NAME = FOLDER_NAME, QUERY = QUERY)
  }
  
  # --- 6. Wrap-up and save the results
  # Save the QUERY object with updated data
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 6.1 Stop logging and close the log file
  log_sink(FILE = sinkfile, START = FALSE)
  
  # --- 6.2 Return the subfolder name for downstream processes
  return(SUBFOLDER_NAME)
  
}  # END FUNCTION

