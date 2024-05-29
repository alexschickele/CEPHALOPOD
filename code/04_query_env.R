#' =============================================================================
#' @name query_env
#' @description appends the query_bio output with a list of environmental
#' values at the sampling stations and a path to the environmental raster that
#' will be used for projections
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return X: a data frame of environmental values at the sampling stations
#' @return Y and S updated with duplicate stations removed
#' @return Updates the output in a QUERY.RData file
#' @return an updated list of subfolders according to the minimum number of occurrence criteria

query_env <- function(FOLDER_NAME = NULL,
                      SUBFOLDER_NAME = NULL){
  set.seed(123)
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_env ********************"))
  # --- 1.2. Load the run metadata and query
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 1.3. Get feature names
  features_name <- CALL$ENV_DATA[[1]] %>% names()

  # --- 2. Re-grid sample on the raster resolution and filter
  # --- 2.1. Do the aggregation
  # The cell centers are at .5, thus it is re-gridded to the nearest .5 value
  res <- res(CALL$ENV_DATA[[1]])[[1]]
  digit <- nchar(sub('^0+','',sub('\\.','',res)))-1
  sample <- QUERY$S %>%
    cbind(QUERY$Y) %>%
    mutate(decimallatitude = round(decimallatitude+0.5*res, digits = digit)-0.5*res) %>%
    mutate(decimallongitude = round(decimallongitude+0.5*res, digits = digit)-0.5*res)

  # --- 2.2. Remove NA in coordinates x month
  sample <- sample %>%
    dplyr::filter(!is.na(decimallatitude) & !is.na(decimallongitude) & !is.na(month))

  # --- 2.3. Early return in case of no biological data
  if(nrow(sample) <= CALL$SAMPLE_SELECT$MIN_SAMPLE){
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  }

  # --- 3. Select one sample per group of identical coordinates x month
  # Among each group of identical lat and long, concatenates description
  # Updates the ID as the previous one is overwritten (we do not keep the raw data)
  S <- sample %>%
    dplyr::select(-names(QUERY$Y)) %>%
    group_by(decimallongitude, decimallatitude, month) %>%
    reframe(across(everything(), ~ str_flatten(unique(.x), collapse = ";"))) %>%
    mutate(ID = row_number())

  # --- 4. Average measurement value per group of identical coordinates x month
  # The corresponding biological value is averaged across all row of S
  Y0 <- sample %>%
    dplyr::select(decimallongitude, decimallatitude, month, names(QUERY$Y))

  Y <- NULL
  for(n in 1:nrow(S)){
    tmp <- Y0 %>%
      inner_join(S[n,], by = c("decimallongitude", "decimallatitude", "month")) %>%
      dplyr::select(names(QUERY$Y))
    tmp <- apply(tmp, 2, mean)
    Y <- rbind(Y, tmp)
  }
  Y <- as.data.frame(Y)

  # --- 5. Extract the environmental data in the data frame
  # If there is an NA, extract from nearest non-NA cells
  # within a certain radius - two grid cells
  X <- NULL
  to_remove <- NULL
  for(j in 1:nrow(S)){
    # --- 5.1. Point toward the right monthly raster
    month <- as.numeric(S$month[j])
    features <- CALL$ENV_DATA[[month]]

    # --- 5.2. First try to extract environmental data
    xy <- S[j,] %>% dplyr::select(x = decimallongitude, y = decimallatitude)
    tmp <- raster::extract(features, xy) %>%
      as.data.frame()

    # --- 5.3. Use neighbor cell if NA is not far inland
    if(is.na(sum(tmp))){
      r_dist <- distanceFromPoints(features, xy) # Compute distance to NA point
      r_dist <- synchroniseNA(stack(r_dist, features[[1]]))[[1]] # Synchronize NA
      r_dist[r_dist > 200e3] <- NA
      min_dist <- which.min(getValues(r_dist)) # Get closest non-NA point ID
      tmp <- features[min_dist] %>%
        as.data.frame()

      # --- 5.4. Remove if NA is too far inland
      if(sum(tmp)==0){
        to_remove <- c(to_remove, j)
      }
    }

    X <- rbind(X, tmp)
  } # End for j

  # --- 6. Early return if no environmental data matching
  if(length(X) == 0){
    message("No environmental data could be extracted for these observation locations \n
            Please check the coverage of the environmental predictors and the location of the observations")
    return(NA)
    }
  colnames(X) <- features_name

  # --- 7. Remove rows that are still NA - i.e. on land
  # Already done for X in the previous step
  if(length(to_remove) != 0){
    Y <- dplyr::slice(Y, -to_remove)
    S <- dplyr::slice(S, -to_remove)
    message(paste("--- ENV EXTRACT : Removed row number", to_remove, "more than 2 grid cells on land \n"))
  }

  # --- 8. Wrap up and save
  # --- 8.1. Remove rare targets and update sample list
  # Designed to clean the input table of proportion data (i.e. avoid sum lines = 0 in test of train sets)
  if(CALL$DATA_TYPE == "proportions"){
    target_filter <- apply(Y, 2, function(x)(x = sum(x/x, na.rm = TRUE)))
    sample_filter <- Y %>% dplyr::select(which(target_filter >= CALL$SAMPLE_SELECT$MIN_SAMPLE)) %>%
      apply(1, sum)
    Y <- Y[which(sample_filter > 0), which(target_filter >= CALL$SAMPLE_SELECT$MIN_SAMPLE)] %>%
      apply(1, function(x)(x = x/sum(x))) %>%
      aperm(c(2,1)) %>%
      as.data.frame()
    X <- X[which(sample_filter > 0),]
    S <- S[which(sample_filter > 0),]
  # Necessary to update SP_SELECT in the CALL object for later...
  CALL$SP_SELECT <- names(Y)
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  } # if proportions

  # --- 8.2. Append QUERY with the environmental values and save
  # And updated Y and S tables with duplicate coordinate removed
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["X"]] <- X
  # And a QC on the sample size or col/row ratio after env. binning
  if(nrow(Y) >= CALL$SAMPLE_SELECT$MIN_SAMPLE & (nrow(Y)/ncol(Y)) > 1){
    QUERY[["eval"]][["SAMPLE_SIZE"]] <- TRUE
  } else {
    QUERY[["eval"]][["SAMPLE_SIZE"]] <- FALSE
  } # end if

  # Save
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 8.3. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)

  # --- 8.4. Update list of SUBFOLDER_NAME
  if(nrow(Y) >= CALL$SAMPLE_SELECT$MIN_SAMPLE & (nrow(Y)/ncol(Y)) > 1){
    return(SUBFOLDER_NAME)
  } else {
    message(paste("QUERY ENV:", SUBFOLDER_NAME, "The sample does not match the minimum sample size - or has a row/col ratio under 1:1 (for proportions data) \n
            Please work on the data to increase the sample size"))
    return(NA)
  }

} # END FUNCTION
