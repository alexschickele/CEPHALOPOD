#' =============================================================================
#' @name query_custom
#' @description extracts biological from the provided custom input, according to 
#' a user provided list of species, time and depth range. The extracted data is 
#' formatted to be directly usable by the models available in this workbench.
#' 
#' @param FOLDER_NAME Name of the folder where output is stored.
#' @param QUERY The QUERY.RData object containing the list of species.
#' 
#' @return QUERY Updated QUERY object with data frames for Y (target values),
#' S (sample stations), and annotations (taxonomic information).
#' 

query_custom <- function(FOLDER_NAME = NULL, QUERY = NULL){
  
  # --- 1. Load CALL.RData, which contains metadata and biological data
  load(file.path(project_wd, "output", FOLDER_NAME, "CALL.RData"))
  
  # --- 1.1. Extract species list (SP_SELECT) from the QUERY object
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  
  # --- 1.2. If WORMS_CHECK is TRUE, retrieve species, including synonyms and child taxa
  if (CALL$WORMS_CHECK) {
    SP_SELECT <- CALL$SP_SELECT_INFO[[SP_SELECT]] %>% 
      unlist() %>% 
      as.numeric() %>% 
      na.omit()  # Remove NA values to clean up the species list
  }
  
  # --- 2. Create the Y target table based on the data type in CALL
  if (CALL$DATA_TYPE == "presence_only") {
    # --- 2.1. For presence data, filter and normalize measurement values
    Y <- CALL$LIST_BIO %>% 
      dplyr::filter(worms_id %in% !!SP_SELECT) %>% 
      dplyr::select(measurementvalue) %>%
      mutate(measurementvalue = ifelse(is.numeric(measurementvalue), 1, measurementvalue))  # Normalize to 1 for presence
  } else if (CALL$DATA_TYPE == "continuous") {
    # --- 2.2. For continuous biomass data, filter measurement values
    Y <- CALL$LIST_BIO %>% 
      dplyr::filter(worms_id %in% !!SP_SELECT) %>% 
      dplyr::select(measurementvalue)
  } else if (CALL$DATA_TYPE == "proportions") {
    # --- 2.3. For proportions data, pivot the table to wide format and normalize
    target_proportions <- CALL$LIST_BIO %>% 
      dplyr::filter(worms_id %in% !!SP_SELECT) %>% 
      dplyr::select("decimallatitude", "decimallongitude", "depth", "year", "month", "measurementvalue", "worms_id") %>%
      pivot_wider(names_from = "worms_id", values_from = "measurementvalue", values_fn = mean, values_fill = 0)
    
    # Prepare Y by removing irrelevant columns and normalizing across rows
    Y <- target_proportions %>%
      dplyr::select(-c(decimallatitude, decimallongitude, depth, year, month)) %>%
      apply(1, function(x) x / sum(x, na.rm = TRUE)) %>%
      t() %>%
      as.data.frame()
  }
  
  # --- 3. Create the S sample table for spatial and temporal information
  if (CALL$DATA_TYPE %in% c("presence_only", "continuous")) {
    # --- 3.1. Filter sample data, exclude non-relevant columns, and convert to numeric
    S <- CALL$LIST_BIO %>% 
      dplyr::filter(worms_id %in% !!SP_SELECT) %>% 
      dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>%
      mutate(across(c(decimallatitude, decimallongitude, month), as.numeric),
             ID = row_number())  # Add a unique ID for each sample
  } else if (CALL$DATA_TYPE == "proportions") {
    # --- 3.2. Handle proportions data and exclude Y columns
    S <- target_proportions %>%
      dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ", names(Y)))) %>%
      mutate(across(c(decimallatitude, decimallongitude, month), as.numeric),
             ID = row_number())
  }
  
  # --- 4. Create the annotation table for species metadata (e.g., taxonomic information)
  annotations <- CALL$LIST_BIO %>%
    dplyr::filter(worms_id %in% !!SP_SELECT) %>%
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>%
    distinct()  # Ensure unique rows for species information
  
  # --- 5. Update and return the QUERY object with Y, S, and annotations tables
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
}
