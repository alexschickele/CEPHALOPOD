#' =============================================================================
#' @name query_occurrence
#' @description Extracts biological data from the OBIS database based on a
#' user-defined list of species, time, and depth range. The extracted data is 
#' formatted to be directly usable by models in this workbench.
#' @param FOLDER_NAME Name of the folder where output is stored.
#' @param QUERY The QUERY.RData object containing the list of species.
#' 
#' @return QUERY Updated QUERY object with data frames for Y (target values),
#' S (sample stations), and annotations (taxonomic information).
#' 

query_occurrence <- function(FOLDER_NAME = NULL, QUERY = NULL) {
  
  # --- 1. Load run metadata and species list
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))  # Load metadata for the current run
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT  # Extract species list from the QUERY object
  
  # --- 1.1. Handle species synonyms and children if WORMS_CHECK is enabled
  # Replace the list of species by the one matched against the taxonomy
  if (CALL$WORMS_CHECK == TRUE) {
    SP_SELECT <- CALL$SP_SELECT_INFO[[SP_SELECT]] %>% 
      unlist() %>% as.numeric() %>% .[!is.na(.)]  # Get all related species IDs
  }
  
  # --- 2. Query OBIS occurrences for the selected species
  target <- occurrence(taxonid = SP_SELECT) %>% 
    dplyr::filter(aphiaID %in% !!SP_SELECT) %>%  # Filter by species IDs
    dplyr::filter(basisOfRecord %in% c("Occurrence", "HumanObservation", "LivingSpecimen", "MachineObservation")) %>%
    dplyr::filter(decimalLatitude != 0 & decimalLatitude >= -90 & decimalLatitude <= 90) %>% 
    dplyr::filter(decimalLongitude != 0 & decimalLongitude >= -180 & decimalLongitude <= 180) %>%
    dplyr::filter(occurrenceStatus == "present") %>% 
    dplyr::mutate(date_month = format(as.Date(eventDate, format = "%Y-%m-%d"), "%m")) %>%  # Extract month from date
    dplyr::select(any_of(c("scientificName", "aphiaID", "decimalLatitude", "decimalLongitude", "depth", "date_year", "date_month", "occurrenceStatus", "basisOfRecord", "taxonRank"))) %>% 
    dplyr::filter(date_year >= CALL$SAMPLE_SELECT$START_YEAR & date_year <= CALL$SAMPLE_SELECT$STOP_YEAR) %>% 
    dplyr::filter(depth >= CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH & depth <= CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH) %>%
    distinct()
  
  # --- 2.1. Ensure taxonRank column exists
  if (!"taxonRank" %in% names(target)) {
    target <- target %>% mutate(taxonRank = "undefined")
  }
  
  # --- 2.2. Standardize column names
  colnames(target) <- c("scientificname", "worms_id", "decimallatitude", "decimallongitude", "depth", "year", "month", "measurementvalue", "measurementunit", "taxonrank")
  
  # --- 3. Query GBIF occurrences for the selected species
  YEAR_0 <- seq(CALL$SAMPLE_SELECT$START_YEAR, CALL$SAMPLE_SELECT$STOP_YEAR)
  SNAME_0 <- unique(target$scientificname)
  
  target_gbif <- mapply(function(YEAR, SNAME) {
    occ_data(scientificName = SNAME,
             year = YEAR,
             depth = paste0(CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH, ",", CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH),
             occurrenceStatus = 'PRESENT',
             limit = 99000)$data
  }, 
  YEAR = rep(YEAR_0, length(SNAME_0)), 
  SNAME = rep(SNAME_0, length(YEAR_0))) %>%
    bind_rows() %>% 
    distinct()
  
  colnames(target_gbif) <- tolower(colnames(target_gbif)) # make sure they are lower case
  
  # --- 3.1. Reformat GBIF data if available
  if (nrow(target_gbif) > 0) {
    target_gbif <- target_gbif %>%
      dplyr::select(any_of(c("decimallatitude", "decimallongitude", "depth", "year", "month"))) %>%
      dplyr::filter(!is.na(decimallatitude) & !is.na(decimallongitude)) %>%
      distinct() %>%
      mutate(scientificname = target$scientificname[1],  # Use the first species name
             worms_id = QUERY$SUBFOLDER_INFO$SP_SELECT,  # Use original species ID
             measurementvalue = "present",
             measurementunit = "Occurrence",
             taxonrank = target$taxonrank[1]) %>%
      dplyr::select(any_of(c("scientificname", "worms_id", "decimallatitude", "decimallongitude", "depth", "year", "month", "measurementvalue", "measurementunit", "taxonrank")))
  }
  
  # --- 4. Combine OBIS and GBIF data
  if (nrow(target_gbif) > 0) {
    target <- rbind(target, target_gbif) %>%
      distinct() %>%
      mutate(nb_occ = n())
  } else {
    target <- target %>% mutate(nb_occ = n())
  }
  
  # --- 5. Create Y (target) table
  Y <- target %>% dplyr::select(measurementvalue)
  
  # --- 6. Create S (sample stations) table
  S <- target %>%
    dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>%
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           month = as.numeric(month),
           ID = row_number())  # Add unique row IDs
  
  # --- 7. Create an Annotation table for species information
  annotations <- target %>%
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>%
    distinct()
  
  # --- 8. Save results in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)  # Return the updated QUERY object
  
}  # END FUNCTION
