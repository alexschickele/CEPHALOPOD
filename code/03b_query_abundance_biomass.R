#' =============================================================================
#' @name query_abundance_biomass
#' @description extracts biological from the Atlanteco database based on a
#' user-defined list of species, time, and depth range. The extracted data is 
#' formatted to be directly usable by models in this workbench.
#' @param FOLDER_NAME Name of the folder where output is stored.
#' @param QUERY The QUERY.RData object containing the list of species.
#' 
#' @return QUERY Updated QUERY object with data frames for Y (target values),
#' S (sample stations), and annotations (taxonomic information).
#' 

query_abundance_biomass <- function(FOLDER_NAME = NULL,
                                    QUERY = NULL){
  
  # --- 1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 1.1. Single query if no WORMS_CHECK
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  # --- 1.2. Multiple query if WORMS_CHECK = TRUE; i.e., also requesting children and synonyms
  if(CALL$WORMS_CHECK == TRUE){
    SP_SELECT <- CALL$SP_SELECT_INFO[[SP_SELECT]] %>% unlist() %>% as.numeric() %>% .[!is.na(.)]
  }
  
  # --- 2. Connect to database (Make sure to close connection after data extraction)
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bc2026_wb3_db",
    user="bluecloud_wb3_reader",
    password="1a6414f89d8a265c8bdd",
    port=5432
  )
  
  # --- 3. Extract data efficiently, only collect necessary columns
  target <- tbl(db, paste0(CALL$DATA_SOURCE, "_data")) %>%
    dplyr::filter(worms_id %in% !!SP_SELECT) %>%
    dplyr::select(decimallatitude, decimallongitude, depth, year, month, measurementvalue, worms_id) %>%
    collect() %>%
    mutate(month = str_pad(month, 2, pad = "0")) # Efficiently pad month
  
  # --- 4. Apply filtering and conditions in a single filter operation
  target <- target %>%
    dplyr::filter(
      depth >= CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH &
        depth <= CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH &
        year >= CALL$SAMPLE_SELECT$START_YEAR &
        year <= CALL$SAMPLE_SELECT$STOP_YEAR &
        !is.na(measurementvalue) & measurementvalue != 0) %>%
    group_by(worms_id) %>%
    distinct() %>%
    ungroup() %>%
    mutate(nb_occ = n())  # count occurrences
  
  # --- 4.2. Optimize proportions calculation using matrix operations
  if(CALL$DATA_TYPE == "proportions"){
    target_proportions <- target %>%
      pivot_wider(names_from = "worms_id", values_from = "measurementvalue", values_fn = mean) %>%
      mutate(across(everything(), ~ replace_na(.x, 0))) %>%
      rowwise() %>%
      mutate(row_sum = sum(c_across(where(is.numeric)))) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), ~ .x / row_sum)) %>%
      dplyr::select(-row_sum)
  }
  
  # --- 5. Create Y target table
  Y <- if (CALL$DATA_TYPE == "proportions") {
    target_proportions
  } else {
    target %>%
      dplyr::select(measurementvalue)
  }
  
  # --- 6. Create S sample table with column selection
  S <- target %>%
    dplyr::select(decimallatitude, decimallongitude, depth, year, month) %>%
    mutate(
      decimallatitude = as.numeric(decimallatitude),
      decimallongitude = as.numeric(decimallongitude),
      month = as.numeric(month),
      ID = row_number()
    ) %>%
    distinct()
  
  # --- 7. Create Annotation table
  annotations <- target %>%
    dplyr::select(worms_id, nb_occ) %>%
    distinct()
  
  # --- 8. Disconnect from database after data is fully extracted
  dbDisconnect(db)
  
  # --- 9. Save results in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
}
