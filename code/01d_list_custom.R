#' =============================================================================
#' @name list_custom
#' @description Extracts available Aphia_ID (Worms ID) corresponding to user-defined criteria
#' from the custom database provided by the user. The function filters the data based on
#' the user's sample selection criteria (e.g., depth, year) and removes unnecessary metadata.
#' 
#' @param DATA_SOURCE File path to the data source (.xlsx, .txt, or .csv format)
#' @param SAMPLE_SELECT A list containing user-defined filtering criteria, including depth and time ranges
#' @param FOLDER_NAME Folder name where output files should be saved (if needed)
#' @return A filtered list of available Worms/Aphia IDs and the number of occurrences within the selected samples

list_custom <- function(DATA_SOURCE,
                        SAMPLE_SELECT,
                        FOLDER_NAME){
  
  # --- 1. Open the file based on its extension
  # --- 1.1. Extract the file extension
  ext <- tools::file_ext(DATA_SOURCE)
  
  # --- 1.2. Conditional data loading based on file type
  if (ext == "xlsx") {
    # For .xlsx files, use readxl (no vroom support for Excel)
    df <- readxl::read_xlsx(DATA_SOURCE)
    colnames(df) <- tolower(colnames(df)) # Standardize column names to lowercase
  } else if (ext == "txt" | ext == "csv") {
    # For .txt or .csv files, use vroom (more efficient for large datasets)
    df <- vroom::vroom(DATA_SOURCE)
    colnames(df) <- tolower(colnames(df)) # Standardize column names to lowercase
  } else {
    # Return a message and exit if file format is not supported
    message("The file extension is not recognized. Please use .xlsx, .txt, or .csv.")
    return(NULL)
  }
  
  # --- 2. Verify presence of mandatory columns
  # --- 2.1. Define the required columns for analysis
  names_qc <- c("scientificname", "worms_id", "decimallatitude", "decimallongitude", 
                "depth", "year", "month", "measurementvalue", "measurementunit", "taxonrank")
  
  # --- 2.2. Check which required columns are present in the data
  names_df <- df %>% dplyr::select(any_of(names_qc)) %>% names()
  
  # --- 2.3. Inform the user if mandatory columns are missing and stop the function
  if (length(names_qc) != length(names_df)) {
    message("The data table must contain the following columns, with one row per sample: \n
        - scientificname : taxonomic name of the species/taxon \n
        - worms_id : AphiaID or other identifier for the species/taxon \n
        - decimallatitude : latitude of the sample in decimal degrees (-90 to +90) \n
        - decimallongitude : longitude of the sample in decimal degrees (-180 to +180) \n
        - depth : sample depth in meters \n
        - year : year of sampling (integer) \n
        - month : month of sampling (integer) \n
        - measurementvalue : recorded value for the sample (numeric or string, e.g., 'present') \n
        - measurementunit : units of the measurement value \n
        - taxonrank : taxonomic rank (e.g., species, genus, order...)")
    
    return(NULL) # Exit the function if columns are missing
  }
  
  # --- 3. Reduce metadata if the number of columns exceeds the memory limit (e.g., 25 columns)
  # This is necessary for large datasets to prevent memory issues during further processing.
  if (ncol(df) > 25) {
    # Save the full dataset for traceback purposes
    save(df, file = file.path(project_wd, "output", FOLDER_NAME, "LIST_BIO_RAW.RData"))
    
    # Reduce the dataset to only the required columns
    df <- df %>% dplyr::select(any_of(names_qc))
    
    # Inform the user that metadata has been reduced
    message("LIST_CUSTOM: reduced metadata to the strict minimum due to memory constraints. 
            The full dataset has been saved separately for reference.")
  }
  
  # --- 4. Filter the data based on user-defined sample selection criteria
  # Only samples that match the depth and time ranges, and have a valid measurement value, are retained.
  data_list <- df %>%
    dplyr::filter(
      depth >= SAMPLE_SELECT$TARGET_MIN_DEPTH & depth <= SAMPLE_SELECT$TARGET_MAX_DEPTH, # Filter by depth range
      year >= SAMPLE_SELECT$START_YEAR & year <= SAMPLE_SELECT$STOP_YEAR, # Filter by year range
      !is.na(measurementvalue), # Remove rows with missing measurement values
      !grepl("abs|Abs", measurementvalue) # Remove rows where the measurement value indicates absence
    ) %>%
    distinct() %>% # Ensure unique rows (if duplicates are present)
    group_by(scientificname) %>% # Group by taxon name to count occurrences
    mutate(nb_occ = n()) %>% # Count occurrences for each taxon
    ungroup()
  
  # --- 5. Return the filtered data list with occurrence counts
  return(data_list)
  
} # end function
