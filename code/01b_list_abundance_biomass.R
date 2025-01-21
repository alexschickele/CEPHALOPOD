#' =============================================================================
#' @name list_abundance_biomass
#' @description Extracts available Aphia_IDs corresponding to user-defined criteria
#' from the available data in the AtlantECO database.
#' @param DATA_SOURCE The data type to query: either 'abundance' or 'biomass'.
#' @param SAMPLE_SELECT List of sample selection criteria, including:
#' TARGET_MIN_DEPTH, TARGET_MAX_DEPTH, START_YEAR, STOP_YEAR, and MIN_SAMPLE.
#' @return A list of available Worms ID (Aphia ID) and the number of occurrences
#' within the selected data type and sample criteria.

list_abundance_biomass <- function(DATA_SOURCE,
                                   SAMPLE_SELECT) {
  
  # --- 1. Establish connection to the database
  # Connect to the remote AtlantECO database. Ensure credentials are secure.
  # This function uses PostgreSQL to access the data. Details for the database
  # and connection credentials are available in the source documentation.
  db <- dbConnect(
    drv = PostgreSQL(),
    host = "postgresql-srv.d4science.org",
    dbname = "bc2026_wb3_db",
    user = "bluecloud_wb3_reader",
    password = "1a6414f89d8a265c8bdd",
    port = 5432
  )
  
  # --- 2. Query the database for abundance or biomass data
  # Retrieve species (Aphia ID) and number of occurrences within the specified
  # sample criteria. Only non-zero measurements are kept.
  message("--- LIST BIO: retrieving species data from ATLANTECO (", DATA_SOURCE, ")")
  
  data_list <- tbl(db, paste0(DATA_SOURCE, "_data")) %>%
    dplyr::filter(
      depth >= !!SAMPLE_SELECT$TARGET_MIN_DEPTH,
      depth <= !!SAMPLE_SELECT$TARGET_MAX_DEPTH,
      year >= !!SAMPLE_SELECT$START_YEAR,
      year <= !!SAMPLE_SELECT$STOP_YEAR,
      measurementvalue != 0,               # Exclude zero measurement values
      !is.na(measurementvalue),            # Exclude missing values
      worms_id != "Not found"              # Exclude invalid Worms IDs
    ) %>%
    group_by(worms_id) %>% # Group by Worms ID (Aphia ID) and calculate occurrences
    summarize(nb_occ = n(), .groups = 'drop') %>%
    dplyr::filter(nb_occ >= !!SAMPLE_SELECT$MIN_SAMPLE) %>% # Only select relevant rows (skip unneeded columns)
    collect()  # Bring the data into R memory only after filtering
  
  # --- 3. Disconnect from the database
  # Disconnecting from the database to free up the connection after query execution.
  dbDisconnect(db)
  
  # --- 4. Wrap up and return the result
  # Return the filtered and processed species data list
  return(data_list)
  
} # END FUNCTION
