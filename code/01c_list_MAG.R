#' =============================================================================
#' @name list_MAG
#' @description Extracts available metadata, taxonomy, and read counts per MAG 
#' based on user-defined criteria. Data is retrieved from the MATOU service 
#' available through Bluecloud Plankton Genomics VLAB.
#' @param SAMPLE_SELECT List of user-defined selection criteria, including
#' MIN_SAMPLE, START_YEAR, STOP_YEAR, TARGET_MIN_DEPTH, and TARGET_MAX_DEPTH.
#' @return A data frame of metadata, samples, and taxonomic annotations.
#' The output contains the number of occurrences for each taxonomic rank.
list_MAG <- function(SAMPLE_SELECT){
  
  # --- 1. Establish Database Connection
  # Connect to the MATOU PostgreSQL database hosted on Bluecloud. 
  # Ensure credentials are handled securely.
  db <- dbConnect(
    drv = PostgreSQL(),
    host = "postgresql-srv.d4science.org",
    dbname = "bluecloud_demo2",
    user = "bluecloud_demo2_u",
    password = "6a26c54a05ec5dede958a370ca744a",
    port = 5432
  )
  
  # --- 2. Raw Query for MAG Data
  # Retrieve MAG data with the relevant taxonomic ranks and read counts.
  message(paste(Sys.time(), "LIST_MAG: Starting raw query"))
  
  data_w_taxo <- tbl(db, "data") %>%
    mutate(MAG = str_sub(Genes, 1, 22)) %>% # Extract MAG ID from Genes and filter by positive read counts
    dplyr::filter(readCount > 0) %>%
    dplyr::select(Station, Phylum, Class, Order, Family, Genus, MAG) %>% # Select relevant columns
    distinct() %>%
    collect()
  
  message(paste(Sys.time(), "LIST_MAG: Query complete"))
  
  # Clean MAG names by removing extra underscores
  data_w_taxo$MAG <- str_replace(paste0(data_w_taxo$MAG, "tmp"), "_tmp|tmp", "")
  
  # --- 3. Filter Samples Based on User Input
  # Retrieve and filter station data from MATOU based on user-defined year and depth ranges.
  message(paste(Sys.time(), "LIST_MAG: Filtering stations based on criteria"))
  
  locs_w_time <- tbl(db, "locs_w_time") %>%
    collect() %>%
    mutate(depth = 1) %>% # MATOU data are surface samples (depth = 1)
    dplyr::filter(
      year >= SAMPLE_SELECT$START_YEAR,
      year <= SAMPLE_SELECT$STOP_YEAR,
      depth >= SAMPLE_SELECT$TARGET_MIN_DEPTH,
      depth <= SAMPLE_SELECT$TARGET_MAX_DEPTH
    ) %>%
    dplyr::select(Station)
  
  # Join the filtered stations with the MAG data
  data_w_taxo <- data_w_taxo %>%
    inner_join(locs_w_time, by = "Station")
  
  # --- 4. Format and Summarize the Data by Taxonomic Rank
  # Create a list of unique scientific names and their taxonomic rank.
  message(paste(Sys.time(), "LIST_MAG: Formatting taxonomy and calculating occurrences"))
  
  list_raw <- lapply(c("Phylum", "Class", "Order", "Family", "Genus", "MAG"), function(rank) {
    id <- which(colnames(data_w_taxo) == rank)
    data_w_taxo %>%
      dplyr::select(all_of(1:id)) %>%
      mutate(taxonrank = rank, scientificname = .data[[rank]])
  }) %>%
    bind_rows() %>%
    distinct()
  
  # --- 5. Count the Number of Observations
  # Group by taxonomic rank and scientific name, and count the number of occurrences.
  message(paste(Sys.time(), "LIST_MAG: Counting number of observations"))
  
  list_bio <- list_raw %>%
    group_by(taxonrank, scientificname) %>%
    summarize(nb_occ = n(), .groups = 'drop') %>%
    dplyr::filter(nb_occ >= SAMPLE_SELECT$MIN_SAMPLE) %>%
    distinct() %>% 
    mutate(worms_id = scientificname) # worms_id do not exist for MAGs, thus adding scientificname instead
  
  # add the full taxonomy as supplementary columns
  taxo_info <- list_raw %>% 
    dplyr::select(taxonrank, scientificname, Class, Order, Family, Genus, MAG) %>% 
    distinct()
  
  # get the full taxo in list_bio
  list_bio <- list_bio %>% 
    left_join(taxo_info) %>% 
    dplyr::filter(!is.na(scientificname)) # remove the NA scientificnames
  
  # --- 6. Disconnect from the Database
  # Ensure that the database connection is closed.
  dbDisconnect(db)
  
  # --- 7. Return the Final Data
  # Return the filtered and summarized list of taxonomic ranks and occurrences.
  return(list_bio)
  
} # END FUNCTION
