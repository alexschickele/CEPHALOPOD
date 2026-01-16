#' =============================================================================
#' @name query_env
#' @description This function appends the query output with environmental values 
#' extracted from rasters at the sampling stations. It also updates the query object 
#' with cleaned data (removing duplicate stations, handling missing values) and 
#' filters samples based on a minimum sample size criterion.
#'
#' The function utilizes the terra library for raster extraction and spatial operations.
#' It ensures environmental data is correctly matched to the biological observations.
#' If necessary, the function retrieves environmental data from the nearest valid cell.
#'
#' @param FOLDER_NAME A character string indicating the folder containing the output.
#' @param SUBFOLDER_NAME A list of subfolders to parallelize the process on.
#'
#' @return A data frame (X) of environmental values at the sampling stations.
#' @return Updated data frames (Y and S) with duplicate stations removed.
#' @return The updated QUERY object saved to a QUERY.RData file.
#' @return An updated list of subfolders based on the minimum sample size criterion.
#'

query_env <- function(FOLDER_NAME = NULL, SUBFOLDER_NAME = NULL) {
  
  # --- 1. Initialize function
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Start logging process
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_env ********************"))
  
  # --- 1.2. Load metadata and previous query
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 1.3. Get environmental feature names
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  features_name <- names(CALL$ENV_DATA[[1]])

  # --- 2. Re-grid sample to raster resolution and filter invalid rows
  res <- terra::res(CALL$ENV_DATA[[1]])[1]  # Resolution of raster
  digit <- nchar(sub('^0+','',sub('\\.','', res)))-1  # Calculate number of decimal places in the resolution
  sample <- QUERY$S %>%
    cbind(QUERY$Y) %>%
    mutate(
      decimallatitude = round(decimallatitude + 0.5 * res, digits = digit) - 0.5 * res,
      decimallongitude = round(decimallongitude + 0.5 * res, digits = digit) - 0.5 * res
    ) %>%
    dplyr::filter(!is.na(decimallatitude) & !is.na(decimallongitude) & !is.na(month))  # Remove rows with NA in lat, long, or month
  
  # --- 2.3. Early return if insufficient data
  if (nrow(sample) <= CALL$SAMPLE_SELECT$MIN_SAMPLE) {
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  }
  
  # --- 3. Group by coordinates and month, aggregate duplicate samples
  S <- sample %>%
    dplyr::select(-names(QUERY$Y)) %>%
    group_by(decimallongitude, decimallatitude, month) %>%
    reframe(across(everything(), ~ str_flatten(unique(.x), collapse = ";"))) %>% 
    mutate(ID = row_number())  # Generate new row IDs
  
  # --- 4. Average biological measurement per group of identical coordinates and month
  Y <- sample %>%
    dplyr::select(decimallongitude, decimallatitude, month, names(QUERY$Y)) %>%
    group_by(decimallongitude, decimallatitude, month) %>%
    reframe(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::select(names(QUERY$Y))
  
  # --- 5. Extract environmental data for the sampling stations
  X <- matrix(NA, nrow = nrow(S), ncol = length(features_name))
  to_remove <- NULL
  
  for (j in seq_len(nrow(S))) {
    month <- as.numeric(S$month[j])
    features <- CALL$ENV_DATA[[month]]
    
    xy <- S[j, c("decimallongitude", "decimallatitude")] %>% as.data.frame()
    
    # Extract environmental data using terra::extract
    env_values <- terra::extract(features, xy)[-1]
    
    # If any NA values, find the nearest non-NA cell
    if (is.na(sum(env_values))) {
      r_dist <- terra::distance(features[[1]], 
                                terra::vect(xy, geom = c("decimallongitude", "decimallatitude"), crs="+proj=longlat +datum=WGS84"))  # Distance to points
      non_na_mask <- !is.na(terra::values(features[[1]]))
      r_dist[!non_na_mask] <- NA  # Exclude areas with NA values
      r_dist[r_dist > 500e3] <- NA
      
      if (any(!is.na(terra::values(r_dist)))) {
        nearest_idx <- where.min(r_dist) # Find the nearest valid cell
        env_values <- terra::extract(features, nearest_idx[1,2]) # Extract from the first min cell ID
      } else {
        to_remove <- c(to_remove, j)  # Mark row for removal if too far inland
      }
    } # End if NA
    
    X[j, ] <- unlist(env_values)
  } # End j nrow(S) loop
  
  X <- as.data.frame(X)
  colnames(X) <- features_name
  
  # --- 6. Remove rows without valid environmental data
  if (length(to_remove) > 0) {
    S <- S[-to_remove, ]
    Y <- Y[-to_remove, ]
    X <- X[-to_remove, ]
    message(paste("Removed", length(to_remove), "rows due to missing environmental data."))
  }
  
  # --- 7. Environmental thinning
  # We keep one observation per environmental space bin
  # Recommended for large datasets only, alleviates bias in the environmental space
  
  # --- 7.1. Testing if inputs are right
  input_test <- lapply(CALL$ENV_THINNING, function(x)(x %in% CALL$ENV_VAR))
  
  # --- 7.2. Perform the thinning
  if(length(grep(FALSE, input_test)) == 0 & !is.null(CALL$ENV_THINNING)){
    # --- 7.2.1. Build the environmental space
    env_space <- terra::values(features) %>% as.data.frame() %>% 
      dplyr::select(CALL$ENV_THINNING) %>% 
      apply(., 2, function(x){quantile(x, seq(0,1,0.01), na.rm = T)}) %>% 
      as.data.frame()
    
    # --- 7.2.2. Identify unique bins
    X_env_space_id <- X %>% 
      dplyr::select(CALL$ENV_THINNING) %>% 
      mutate_all(., function(x){
        bin = quantile(x, seq(0,1, length.out = CALL$ENV_THINNING_BINS), na.rm = T)
        cut(x, breaks = bin, include.lowest = T, labels = F)
      })
    
    # --- 7.2.3. Subset the occurrence per bin
    to_group <- X_env_space_id %>% 
      mutate(id = row_number()) %>% 
      dplyr::group_by_at(CALL$ENV_THINNING) %>% 
      mutate(group_id = cur_group_id()) %>% .$group_id
    
    # --- 7.2.4. Apply the thinning on the input data
    # Biological target and features value is averaged per bin
    # Location and metadata is drawn from the first sample as they have the same feature associated anyway
    if (length(to_group) > 0) {
      S_thinned <- cbind(S, to_group) %>% 
        group_by(to_group) %>% 
        slice_sample(n = 1) %>% 
        ungroup() %>% 
        dplyr::select(-to_group)
      S <- S_thinned
      
      Y_thinned <- cbind(Y, to_group) %>% 
        group_by(to_group) %>% 
        summarize_all(median) %>% 
        ungroup() %>% 
        dplyr::select(-to_group)
      Y <- Y_thinned
      
      X_thinned <- cbind(X, to_group) %>% 
        group_by(to_group) %>% 
        summarize_all(median) %>% 
        ungroup() %>% 
        dplyr::select(-to_group)
      X <- X_thinned
      
      message(paste("Keeping only", nrow(S), "rows in unique environmental space bins."))
    }
  } else {
    message("ENVIRONMENTAL THINNING: could not perform the environmental thinning as the specified features are not part of the inputs \n")
  }
  
  # --- 8. Filter out rare species (only for proportion data)
  if (CALL$DATA_TYPE == "proportions") {
    target_filter <- apply(Y, 2, function(x) sum(!is.na(x)))
    valid_targets <- which(target_filter >= CALL$SAMPLE_SELECT$MIN_SAMPLE)
    
    sample_filter <- rowSums(!is.na(Y[, valid_targets]))
    valid_samples <- which(sample_filter > 0)
    
    Y <- Y[valid_samples, valid_targets]
    X <- X[valid_samples, ]
    S <- S[valid_samples, ]
    
    # Normalize proportions data
    Y <- t(apply(Y, 1, function(row) row / sum(row, na.rm = TRUE)))
    
    # Update species selection in CALL object
    CALL$SP_SELECT <- colnames(Y)
    CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::wrap(x)) # Pack the rasters before saving
    
    save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  }
  
  # --- 9. Save updated query data
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["X"]] <- X
  
  if (nrow(Y) >= CALL$SAMPLE_SELECT$MIN_SAMPLE && (nrow(Y) / ncol(Y)) > 1) {
    QUERY[["eval"]][["SAMPLE_SIZE"]] <- TRUE
  } else {
    QUERY[["eval"]][["SAMPLE_SIZE"]] <- FALSE
  }
  
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 10. End logging and return
  log_sink(FILE = sinkfile, START = FALSE)
  
  if (QUERY[["eval"]][["SAMPLE_SIZE"]]) {
    return(SUBFOLDER_NAME)
  } else {
    message(paste("Insufficient sample size or bad row/column ratio for proportions data in", SUBFOLDER_NAME))
    return(NA)
  }
  
} # END FUNCTION

