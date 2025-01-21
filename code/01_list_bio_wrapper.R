#' =============================================================================
#' @name list_bio_wrapper
#' @description Wrapper function that handles different data types and services to list
#' available species or taxonomic data, applying user-specified sample selection criteria.
#' The function redirects to appropriate helper functions depending on the data source.
#'
#' @param FOLDER_NAME Name of the folder to create, corresponding to the run.
#' @param DATA_SOURCE The type of data to query, which can be:
#' "occurrence" (OBIS data), "biomass" or "abundance" (AtlantECO data), "MAG" (MATOU/MAG data),
#' or a file path (e.g., "custom") for a user-provided dataset.
#' @param SAMPLE_SELECT A list of sample selection criteria, including:
#' - MIN_SAMPLE: Minimum number of geographical points required (i.e., non-zero records).
#' - TARGET_MIN_DEPTH, TARGET_MAX_DEPTH: Minimum and maximum target depth (in meters).
#' - FEATURE_MIN_DEPTH, FEATURE_MAX_DEPTH: Minimum and maximum feature depth (in meters).
#' - START_YEAR, STOP_YEAR: Start and stop years for the sampling period.
#'
#' @details For custom tables, the data must contain specific columns (e.g., scientificname, worms_id, etc.).
#' For other data sources, this function calls appropriate helper functions to retrieve and process data.
#'
#' @return Depending on the data source, this function returns:
#' - For "occurrence", "abundance", or "biomass" data: A data frame with species occurrences.
#' - For "MAG" data: A full list of metadata, samples, and taxonomic annotations.
#' - For custom data: A formatted data frame based on the input file.
#' The result is also saved in the specified folder to avoid re-running queries.

list_bio_wrapper <- function(FOLDER_NAME = "test_run",
                             DATA_SOURCE = "biomass",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){
  
  # --- 1. Initialization and setup
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Assign default feature depth range if not provided
  if(is.null(SAMPLE_SELECT$FEATURE_MIN_DEPTH)) SAMPLE_SELECT$FEATURE_MIN_DEPTH <- SAMPLE_SELECT$TARGET_MIN_DEPTH
  if(is.null(SAMPLE_SELECT$FEATURE_MAX_DEPTH)) SAMPLE_SELECT$FEATURE_MAX_DEPTH <- SAMPLE_SELECT$TARGET_MAX_DEPTH
  
  # --- 1.2. Parameter validation: Check that DATA_SOURCE is a character string
  if(!is.character(DATA_SOURCE)){
    stop("The data source should be one of 'biomass', 'abundance', 'occurrence', 'MAG', or a file path for custom data (e.g., .csv, .txt, .xlsx).")
  }
  
  # --- 1.3. Folder creation for the run
  # Construct the folder path based on the provided FOLDER_NAME
  folderpath <- paste0(project_wd, "/output/", FOLDER_NAME)
  
  # If the folder already exists, stop the function to prevent overwriting
  if(file.exists(folderpath)){
    stop("--- OVERWRITE ALARM: This folder name is already used. \nPlease choose another FOLDER_NAME or delete the existing folder.")
  } else {
    dir.create(folderpath)  # Create the folder if it does not exist
  }
  
  # --- 2. Data source redirection
  # Depending on the data source type, call the appropriate function
  
  # --- 2.1. OBIS/GBIF (Occurrence) Data
  if(DATA_SOURCE == "occurrence"){
    LIST_BIO <- list_occurrence(DATA_SOURCE = DATA_SOURCE, SAMPLE_SELECT = SAMPLE_SELECT)
  }
  
  # --- 3. ATLANTECO (Biomass or Abundance) Data
  if(DATA_SOURCE == "biomass" | DATA_SOURCE == "abundance"){
    LIST_BIO <- list_abundance_biomass(DATA_SOURCE = DATA_SOURCE, SAMPLE_SELECT = SAMPLE_SELECT)
  }
  
  # --- 4. MATOU (MAG) Data
  if(DATA_SOURCE == "MAG"){
    LIST_BIO <- list_MAG(SAMPLE_SELECT = SAMPLE_SELECT)
  }
  
  # --- 5. Custom Data Access
  # If the data source is not one of the predefined types, assume it is a file path for custom data.
  if(DATA_SOURCE != "MAG" & DATA_SOURCE != "biomass" & DATA_SOURCE != "abundance" & DATA_SOURCE != "occurrence"){
    LIST_BIO <- list_custom(DATA_SOURCE = DATA_SOURCE, SAMPLE_SELECT = SAMPLE_SELECT, FOLDER_NAME = FOLDER_NAME)
  }
  
  # --- 6. Wrap up and save the result
  # Create a list object to store the data source, sample selection criteria, and result (LIST_BIO)
  CALL <- list(DATA_SOURCE = DATA_SOURCE, SAMPLE_SELECT = SAMPLE_SELECT, LIST_BIO = LIST_BIO)
  
  # Save the result in the specified folder to avoid re-running the query in future runs
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  
  # Return the list of biological data (LIST_BIO)
  return(CALL$LIST_BIO)
  
} # END FUNCTION
