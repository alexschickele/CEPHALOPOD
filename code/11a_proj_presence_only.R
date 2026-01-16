#' =============================================================================
#' @name proj_presence_only
#' @description computes spatial projections for the presence_only data sub-pipeline
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_presence_only <- function(QUERY, MODEL, CALL){

  # --- 1. Initialize function
  # --- 1.1. Open base raster and values
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r0 <- CALL$ENV_DATA[[1]][[1]] # Base raster 

  # --- 1.2. Define the projections to compute
  # Only the algorithm that passed QC if the argument FAST equals TRUE
  loop_over <- if(CALL$FAST) MODEL$MODEL_LIST else CALL$HP$MODEL_LIST

  # --- 2. Define bootstraps with tidymodels
  boot_split <- bootstraps(cbind(QUERY$Y, QUERY$X), 
                           times = CALL$N_BOOTSTRAP, strata = "measurementvalue")

  # --- 3. Start the loop over algorithms
  for(i in loop_over){
    
    # --- 3.1. Fit model on bootstrap
    # fit_resamples() does not save models by default. Thus the control_resamples()
    boot_fit <- MODEL[[i]][["final_wf"]] %>%
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x),
                    verbose = FALSE)) %>%
      unnest(.extracts)


    # --- 4. Loop over month for predictions
    y_hat <- NULL
    for(m in seq_along(CALL$ENV_DATA)){

      # --- 4.1. Load the right features
      features <- terra::subset(CALL$ENV_DATA[[m]], QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
        terra::as.data.frame()
      
      # --- 4.2. Compute one prediction per bootstrap
      # As we extracted the model information in a supplementary column, we can
      # directly compute the bootstrap within the synthetic resample object.
      boot_proj <- boot_fit %>%
        mutate(proj = purrr::map(.extracts, function(x)(x = predict(x, features))))

      # --- 4.3. First transform the object into a cell x bootstrap matrix
      # /!\ Need to create a unique row identifier for pivot_wider to work...
      boot_proj <- boot_proj %>%
        dplyr::select(id, proj) %>%
        unnest(c(id, proj)) %>%
        as.data.frame() %>%
        group_by(id) %>%
        mutate(row = row_number()) %>%
        pivot_wider(names_from = id, values_from = .pred) %>%
        dplyr::select(-row) 

      # --- 4.4. Assign the desired values to the non-NA cells in the base raster
      boot_proj <- apply(boot_proj, 2, function(x){
        r <- terra::values(r0) # Base raster values
        r[!is.na(r)] <- x
        x <- r
      })

      # --- 4.5. Concatenate with previous month
      y_hat <- abind(y_hat, boot_proj, along = 3)
      message(paste("--- PROJ :", i,"- month", m, "done \t"))
      
      rm(boot_proj, features)
      gc()
      
    } # for m month

    # --- 5. Cut spatial discontinuities
    if (!is.null(CALL$CUT)) {
      y_hat <- apply(y_hat, c(2, 3), function(x) {
        xy <- QUERY$S %>% dplyr::select(decimallongitude, decimallatitude) %>%
          .[QUERY$Y == 1, ] # Extract values
        x2 <- x
        x[x < (CALL$CUT * 1)] <- 0 # Set threshold
        r_patch <- terra::patches(setValues(r0, x), zeroAsNA = TRUE) # Compute patches
        id_patch <- terra::extract(r_patch, xy)$patches %>% unique() %>% .[!is.na(.)] # Verify if each patch has observations
        
        if(length(id_patch) > 0) { # Make sure there is patches
          r_patch[!(r_patch %in% id_patch)] <- 0 # Remove the ones that dont
          r_patch[r_patch > 0] <- 1 # scale
          x <- x2 * terra::values(r_patch)
        } else {
          x <- x2
          message(paste("--- CUT : could not been applied to this projection of", i, 
                        "because it has a maximum suitability index inferior to the cutting threshold \n"))
        }
        
        return(x)
      })
    } # if cut

    # --- 6. Compute the average CV across bootstrap runs as a QC
    if(dim(y_hat)[[2]] == CALL$N_BOOTSTRAP) {
      NSD <- apply(y_hat, c(1,3), function(x)(x = sd(x, na.rm = TRUE))) %>%
        mean(na.rm = TRUE)
      NSD <- NSD/mean(y_hat, na.rm = TRUE)
    } else {
      NSD <- NA
      message(paste("--- PROJ: Model", i, " discarded, bootstrap did not complete"))
    }
    
    # --- 7. Append the MODEL object
    # --- 7.1. Save the evaluation metric and projections
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    MODEL[[i]][["eval"]][["NSD"]] <- NSD

    # --- 7.2. Discard low quality models according to NSD
    # Fixed at 50% average across the projection
    if(MODEL[[i]][["eval"]][["NSD"]] > 0.5 | is.na(MODEL[[i]][["eval"]][["NSD"]]) | is.null(MODEL[[i]][["eval"]][["NSD"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to NSD =", MODEL[[i]][["eval"]][["NSD"]], "> 0.5 \n"))
    }

    # --- 7.3. Discard if the PRE_VIP QC in the query is too low
    # In case of FAST = TRUE; the species was kept until here but we do not want
    # to compute an ensemble for it.
    if(QUERY$eval$PRE_VIP < 0.05){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to PRE_VIP =", QUERY$eval$PRE_VIP, "< 0.05 \n"))
    }
  } # for i model loop

  # --- 8. Build the ensemble NSD if there is an ensemble
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    MODEL[["ENSEMBLE"]][["eval"]][["NSD"]] <- lapply(MODEL$MODEL_LIST,
                                                     FUN = function(x){
                                                       x <- MODEL[[x]]$eval$NSD
                                                     }) %>%
      unlist() %>% mean()
  } # End ENSEMBLE QC

  return(MODEL)

} # END FUNCTION


