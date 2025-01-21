#' =============================================================================
#' @name list_occurrence
#' @description Extracts available Aphia_IDs (Worms IDs) based on user-defined criteria
#' from the OBIS and GBIF databases. Filters species by year, depth, and minimum number of occurrences.
#' 
#' @param DATA_SOURCE Passed from the wrapper function (typically "occurrence").
#' @param SAMPLE_SELECT List of criteria passed from the wrapper function, including:
#' - START_YEAR, STOP_YEAR: The time period for filtering the data (year range).
#' - TARGET_MIN_DEPTH, TARGET_MAX_DEPTH: Depth range for the query (in meters).
#' - MIN_SAMPLE: Minimum number of occurrences for a species to be included.
#' 
#' @return A data frame of available Worms IDs (Aphia_IDs) with the number of occurrences
#' from OBIS and GBIF combined, and taxonomic information. Each row represents a unique species.
#' 
list_occurrence <- function(DATA_SOURCE,
                            SAMPLE_SELECT){
  
  # --- 1. Retrieve taxa from OBIS
  # This section queries the OBIS database for species based on the user-defined date and depth range.
  # Filters out species with fewer than 10 records to speed up processing (pre-filtering).
  message("--- LIST BIO: Retrieving species occurrences from OBIS and GBIF")
  
  data_list <- checklist(
    startdate = as.Date(as.character(SAMPLE_SELECT$START_YEAR), format = "%Y"),
    enddate = as.Date(as.character(SAMPLE_SELECT$STOP_YEAR), format = "%Y"),
    startdepth = SAMPLE_SELECT$TARGET_MIN_DEPTH,
    enddepth = SAMPLE_SELECT$TARGET_MAX_DEPTH
  ) %>%
    dplyr::filter(records >= 10) %>%  # Only keep species with more than 2 records
    dplyr::select(acceptedNameUsageID, taxonRank, acceptedNameUsage, records) %>%
    distinct()  # Ensure uniqueness
  
  # --- 1.1. Rename columns for clarity
  colnames(data_list) <- c("worms_id", "taxonrank", "scientificname", "obis_occ")
  
  # --- 1.2. Clean OBIS data
  # Remove species names with spurious characters, such as brackets (e.g., updated OBIS records might include brackets)
  data_list <- data_list %>% 
    dplyr::filter(!grepl("\\[", scientificname))  # Filter out names with brackets
  
  # --- 2. Retrieve data from GBIF
  # Queries GBIF for the number of occurrences based on species name. Filtering by year and depth will happen later.
  message(paste(Sys.time(), "--- Retrieving GBIF data: START"))
  
  gbif_occ <- mclapply(
    data_list$scientificname,
    FUN = function(NAME){
      occ_count(
        scientificName = NAME,
        year = paste0(SAMPLE_SELECT$START_YEAR, ",", SAMPLE_SELECT$STOP_YEAR),  # Query over the same year range
        depth = paste0(SAMPLE_SELECT$TARGET_MIN_DEPTH, ",", SAMPLE_SELECT$TARGET_MAX_DEPTH),  # Query over the same depth range
        occurrenceStatus = 'PRESENT'  # Only consider species with 'PRESENT' status
      )
    },
    mc.cores = MAX_CLUSTERS  # Run in parallel to speed up the process
  ) %>% unlist()  # Unlist the results to get a vector
  
  message(paste(Sys.time(), "--- Retrieving GBIF data: DONE"))
  
  # --- 3. Combine OBIS and GBIF data
  # Add the GBIF occurrences to the OBIS occurrences and calculate total occurrences (nb_occ).
  data_list <- data.frame(
    data_list,
    gbif_occ = gbif_occ,
    nb_occ = data_list$obis_occ + gbif_occ  # Sum of OBIS and GBIF occurrences
  ) %>%
    dplyr::filter(nb_occ >= SAMPLE_SELECT$MIN_SAMPLE)  # Filter species with fewer than the minimum sample size
  
  # --- 4. Return the final list of species with occurrence data
  return(data_list)
  
} # END FUNCTION
