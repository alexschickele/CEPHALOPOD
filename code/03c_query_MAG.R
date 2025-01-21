#' =============================================================================
#' @name query_MAG
#' @description extracts biological from the MATOU data downloaded in list_bio, 
#' according to a user provided list of species, time and depth range. 
#' The extracted data is formatted to be directly usable by the models available
#' in this workbench.
#' 
#' @param FOLDER_NAME Name of the folder where output is stored.
#' @param QUERY The QUERY.RData object containing the list of species.
#' 
#' @return QUERY Updated QUERY object with data frames for Y (target values),
#' S (sample stations), and annotations (taxonomic information).
#' 

query_MAG <- function(FOLDER_NAME = NULL, QUERY = NULL){
  
  # --- 1. Initialize Database Connection
  db <- dbConnect(
    drv = PostgreSQL(),
    host = "postgresql-srv.d4science.org",
    dbname = "bluecloud_demo2",
    user = "bluecloud_demo2_u",
    password = "6a26c54a05ec5dede958a370ca744a",
    port = 5432
  )
  
  # --- 2. Load CALL Parameters
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  
  # --- 3. Extract species correspondences from LIST_BIO
  target_input <- CALL$LIST_BIO %>%
    dplyr::filter(scientificname %in% SP_SELECT)
  
  # --- 3.1. Find appropriate taxonomic rank that matches SP_SELECT
  taxo_rank_nb <- table(target_input$taxonrank)
  trank <- names(taxo_rank_nb)[taxo_rank_nb == length(SP_SELECT)]
  
  if(length(trank) == 0){
    message("QUERY_MAG: No unique taxonomic rank found. Please select valid species names.")
    return(NULL)  # Early exit if no matching rank
  }
  
  # --- 4. Raw Query to Fetch Data from the Database
  message(paste(Sys.time(), "QUERY_MAG: Querying the raw data..."))
  
  # Query and filter directly within the database to reduce data transferred to R
  target <- tbl(db, "data") %>%
    dplyr::select(readCount, Station, Latitude, Longitude, Phylum, Class, Order, Family, Genus, Genes) %>%
    mutate(MAG = str_sub(Genes, 1, 22)) %>% 
    dplyr::select(-Genes) %>% 
    mutate(MAG = str_replace(paste0(MAG, "tmp"), "_tmp|tmp", "")) %>%
    filter_at(vars(trank), any_vars(. %in% SP_SELECT)) %>%
    collect() %>%  # Collect the data after filtering
    group_by(scientificname = !!sym(trank), taxonrank = trank, Station) %>%
    summarize(measurementvalue = sum(readCount), .groups = "drop")
  
  message(paste(Sys.time(), "DONE"))
  
  # --- 5. Normalize by Total Reads per Station
  sum_station <- tbl(db, "sum_station") %>%
    collect()
  
  target_norm <- target %>%
    left_join(sum_station, by = "Station") %>%
    mutate(measurementvalue = measurementvalue / sum_reads) %>%
    group_by(scientificname, taxonrank, Station) %>%
    summarize(measurementvalue = mean(measurementvalue), .groups = "drop")
  
  # --- 6. Join with Locations Data and Fix Column Names
  locs_w_time <- tbl(db, "locs_w_time") %>%
    collect()
  
  target_pretty <- target_norm %>%
    left_join(locs_w_time, by = "Station") %>%
    mutate(
      measurementunit = "Relative proportion of metagenomic reads; size class 0.8 to 5 micrometer",
      depth = 1,  # Default depth of 1
      worms_id = scientificname # No WoRMS taxonomy for MAGs, so we use the scientific name as an ID
    ) %>%
    dplyr::select(scientificname, worms_id, Latitude, Longitude, depth, year, month, measurementvalue, measurementunit, taxonrank, Station) %>%
    dplyr::rename(decimallatitude = Latitude, decimallongitude = Longitude)
  
  # --- 7. Remove rows with missing coordinates
  target_pretty <- target_pretty %>%
    dplyr::filter(!is.na(decimallatitude) & !is.na(decimallongitude))
  
  # --- 8. Build Sample Table (S)
  S <- target_pretty %>%
    dplyr::select(-measurementvalue, -worms_id, -taxonrank, -scientificname) %>%
    distinct() %>%
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           month = as.numeric(month),
           ID = row_number())
  
  # --- 9. Build Target Table (Y)
  # Use a list to extract measurement values for each species and bind columns efficiently
  Y <- lapply(SP_SELECT, function(species){
    target_pretty %>%
      dplyr::filter(scientificname == species) %>%
      dplyr::select(measurementvalue)
  }) %>%
    bind_cols()
  
  # Rename columns in Y with species names
  colnames(Y) <- SP_SELECT
  
  # --- 10. Remove Rows that Sum to Zero
  to_remove <- which(rowSums(Y) == 0)
  if(length(to_remove) > 0){
    Y <- Y[-to_remove, ]
    S <- S[-to_remove, ]
  }
  
  # --- 11. Transform Data Based on CALL$DATA_TYPE
  if(CALL$DATA_TYPE == "proportions"){
    Y <- apply(Y, 1, function(x) x / sum(x)) %>%
      t() %>%
      as.data.frame()
  }
  
  if(CALL$DATA_TYPE == "presence_only"){
    Y <- (Y > 0) %>% as.data.frame()  # Convert to presence data
    names(Y) <- SP_SELECT
  }
  
  if(CALL$DATA_TYPE == "continuous"){
    Y <- apply(Y, 1, function(x) vegan::diversity(x, "shannon")) %>%
      as.data.frame()
    names(Y) <- "measurementvalue"
  }
  
  # --- 12. Build Annotations Table
  annotations <- target_pretty %>%
    dplyr::select(worms_id, taxonrank, scientificname) %>%
    distinct()
  
  # --- 13. Save Data to the QUERY Object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  # --- 14. Disconnect from Database
  dbDisconnect(db)
  
  return(QUERY)
  
}  # END FUNCTION
