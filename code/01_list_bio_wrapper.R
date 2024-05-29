#' =============================================================================
#' @name list_bio_wrapper
#' @description wrapper around functions to list the available species within the
#' different data types and access service sources.
#' @param FOLDER_NAME name of the folder to create, corresponding to the run
#' @param DATA_SOURCE type of data to query from : string among "occurrence" (i.e. OBIS data),
#' "biomass" (i.e. Ecotaxa data) or "omic" (i.e. MGnify) data;
#' this parameter can also be a string path to a custom input table (i.e. "custom")
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection (i.e. non-zero records)
#' TARGET_MIN_DEPTH and TARGET_MAX_DEPTH: minimum and maximum depth levels in m (for target)
#' FEATURE_MIN_DEPTH and FEATURE_MAX_DEPTH: minimum and maximum depth levels in m (for feature)
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @details For custom tables, the data table must contain the following columns and one row per sample :
#' - scientificname : name of the taxa
#' - worms_id : aphiaID corresponding to the taxa name. Fill with another identified if none is available or not relevant
#' - decimallatitude : decimal, latitude of the sample (-90 to +90)
#' - decimallongitude : decimal, longitude of the sample (-180 to +180)
#' - depth : integer, depth of the sample
#' - year : integer, year of sampling
#' - measurementvalue : numeric or string (e.g. present)
#' - measurementunit : information on the unit of the measurement
#' - taxonrank : taxonomic ranking corresponding to the scientific name (e.g. species, gender, order...)
#' @return for "occurrence" and "biomass" data, returns a data frame with the number
#' of occurrences per worms_ID available
#' @return for "omic" data, returns a complete list of metadata, samples, taxonomic
#' annotations available.
#' @return for custom data from file path, returns the formatted file in memory
#' @return the returned object is saved in the run file to avoid re-running the query

list_bio_wrapper <- function(FOLDER_NAME = "test_run",
                             DATA_SOURCE = "biomass",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){

  # --- 1. Initialize
  # --- 1.1. Add default predictor depth range if not specified
  if(is.null(SAMPLE_SELECT$FEATURE_MIN_DEPTH)){SAMPLE_SELECT$FEATURE_MIN_DEPTH <- SAMPLE_SELECT$TARGET_MIN_DEPTH}
  if(is.null(SAMPLE_SELECT$FEATURE_MAX_DEPTH)){SAMPLE_SELECT$FEATURE_MAX_DEPTH <- SAMPLE_SELECT$TARGET_MAX_DEPTH}

  # --- 1.1. Parameter checking
  if(!is.character(DATA_SOURCE)){
    stop("The specified data source should be 'biomass', 'occurrence', 'omic' or a path to file (.csv, .txt, .xlsx) for custom data")
  }

  # --- 1.2. Folder creation
  # If the directory exist, stops the creation and inform the user to manually delete the
  # previous runs, or name the directory differently
  folderpath <- paste0(project_wd,"/output/",FOLDER_NAME)
  if(file.exists(folderpath)==TRUE){
    stop("--- OVERWRITE ALARM: This foldername is already used. \n Please choose another FOLDER_NAME or delete the concerned folder")
  } else {
    dir.create(folderpath)
  }

  # --- 2. Redirection to OBIS data access
  # For occurrence source data
  if(DATA_SOURCE == "occurrence"){

    # --- 2.1. Run function
    LIST_BIO <- list_occurrence(DATA_SOURCE = DATA_SOURCE,
                                SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection

  # --- 3. Redirection to ATLANTECO data access
  # For biomass source data
  if(DATA_SOURCE == "biomass"){

    # --- 3.1. Run function
    LIST_BIO <- list_biomass(DATA_SOURCE = DATA_SOURCE,
                               SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection

  # --- 4. Redirection to MGNIFY data access
  # For omic source data
  if(DATA_SOURCE == "omic"){

    # --- 4.1. Run function
    LIST_BIO <- list_omic(SAMPLE_SELECT = SAMPLE_SELECT)
  } # End MGNIFY redirection

  # --- 5. Redirection to CUSTOM data access
  # For any type of data
  if(DATA_SOURCE != "omic" & DATA_SOURCE != "biomass" & DATA_SOURCE != "occurrence"){

    # --- 5.1. Run function
    LIST_BIO <- list_custom(DATA_SOURCE = DATA_SOURCE,
                            SAMPLE_SELECT = SAMPLE_SELECT,
                            FOLDER_NAME = FOLDER_NAME)
  } # End MGNIFY redirection

  # --- 6. Wrap up and save
  CALL <- list(DATA_SOURCE = DATA_SOURCE,
               SAMPLE_SELECT = SAMPLE_SELECT,
               LIST_BIO = LIST_BIO)
  save(CALL, file = paste0(project_wd,"/output/", FOLDER_NAME, "/CALL.RData"))
  return(CALL$LIST_BIO)

} # END FUNCTION
