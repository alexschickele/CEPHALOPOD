#' =============================================================================
#' @name proj_proportions
#' @description computes spatial projections for the proportions data sub-pipeline
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the model object from the master pipeline
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_proportions <- function(QUERY,
                             MODEL,
                             CALL){

  # --- 1. Initialize function
  # --- 1.1. Open base raster and values
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r0 <- CALL$ENV_DATA[[1]][[1]] # Base raster 

  # --- 1.2. Early stop function if model did not pass QC and fast = TRUE
  if(CALL$FAST == TRUE & (length(MODEL$MODEL_LIST) == 0)){
    return(MODEL)
  }

  # --- 1.3. Source the MBTR functions
  library(reticulate)
  source_python(paste0(project_wd,"/function/mbtr_function.py"))

  # --- 2. Define bootstraps
  # Run the bootstrap generation from tidy models
  boot_split <- bootstraps(cbind(QUERY$Y, QUERY$X), times = CALL$N_BOOTSTRAP)

  # --- 3. Fit model on bootstrap
  # --- 3.1. Extract bootstrap input
  for(b in 1:CALL$N_BOOTSTRAP){
    X_tr <- boot_split$splits[[b]] %>% analysis() %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_tr, paste0(project_wd, "/data/MBTR_cache/", b, "_X_tr.feather"))

    Y_tr <- boot_split$splits[[b]] %>% analysis() %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_tr, paste0(project_wd, "/data/MBTR_cache/", b, "_Y_tr.feather"))

    X_val <- boot_split$splits[[b]] %>% assessment() %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_val, paste0(project_wd, "/data/MBTR_cache/", b, "_X_val.feather"))

    Y_val <- boot_split$splits[[b]] %>% assessment() %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_val, paste0(project_wd, "/data/MBTR_cache/", b, "_Y_val.feather"))
  } # for b
  
  rm(X_tr, Y_tr, X_val, Y_val) # Intermediate memory cleanup
  gc()

  # --- 3.2. Fit model
  # --- 3.2.1. Load .py functions & instance
  # NOTE: I know this is ugly but reticulate loses the connection at every parallel
  # worker AND every for loop turn... therefore it is re-sourced many times in
  # the script

  library(reticulate)
  source_python(paste0(project_wd,"/function/mbtr_function.py"))

  # --- 3.2.2. Computing the model objects
  message(paste0(Sys.time(), "--- MBTR: re-fit model for bootstrap"))
  boot_fit <- mcmapply(FUN = mbtr_fit,
                       path = paste0(project_wd, "/data/MBTR_cache/",1:CALL$N_BOOTSTRAP,"_"),
                       hp_id = "0",
                       loss_type='mse',
                       n_boosts = as.integer(MODEL$MBTR$final_wf$NBOOST),
                       min_leaf= MODEL$MBTR$final_wf$MEAN_LEAF,
                       learning_rate=MODEL$MBTR$final_wf$LEARNING_RATE,
                       lambda_weights=MODEL$MBTR$final_wf$LEARNING_RATE/100,
                       lambda_leaves=0,
                       n_q= as.integer(MODEL$MBTR$final_wf$N_Q),
                       early_stopping_rounds = 10,
                       SIMPLIFY = FALSE,
                       USE.NAMES = FALSE,
                       mc.cores = 20)

  # --- 3.2.3. Reload them in R because of "previous session invalidity"
  for(b in 1:CALL$N_BOOTSTRAP){
    boot_fit[[b]] <- py_load_object(paste0(project_wd, "/data/MBTR_cache/",b,"_0_m"), pickle = "pickle")[[1]]
  }
  message("Bootstrap fitting - DONE")

  # --- 4. Loop over month for predictions
  y_hat <- NULL
  for(m in 1:length(CALL$ENV_DATA)){

    # --- 4.1. Load the right features
    features <- terra::subset(CALL$ENV_DATA[[m]], QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
      terra::as.data.frame()

    # --- 4.2. Computing one prediction per bootstrap
    library(reticulate)
    source_python(paste0(project_wd,"/function/mbtr_function.py")) # needed for some reason...

    boot_proj <- mclapply(boot_fit,
                          function(a_boot){
                            library(reticulate)
                            source_python(paste0(project_wd,"/function/mbtr_function.py"))
                            out <- mbtr_predict(model = a_boot, X_pred = features, n_boosts = as.integer(MODEL$MBTR$final_wf$NBOOST))
                            return(out)},
                          mc.cores = 20) %>%
      abind(along = 3)

    # --- 4.3. Assign the desired values to the non-NA cells in the list
    boot_proj <- apply(boot_proj, c(2,3), function(x){
      r <- terra::values(r0) # Base raster values
      r[!is.na(r)] <- x
      x <- r
    }) %>%
      aperm(c(1,3,2))

    # --- 4.4. Concatenate with previous month
    y_hat <- abind(y_hat, boot_proj, along = 4)
    message(paste("--- PROJ : month", m, "done \t"))
    
    rm(boot_proj, features)
    gc()
  } # for m month

  # --- 6. Compute the average CV across bootstrap runs as a QC
  NSD <- apply(y_hat, c(1,3,4), function(x)(x = sd(x, na.rm = TRUE))) %>%
    mean(na.rm = TRUE)
  NSD <- NSD/mean(y_hat, na.rm = TRUE)

  # --- 7. Append the MODEL object
  # --- 7.1. Save the evaluation metric and projections
  MODEL[["MBTR"]][["proj"]][["y_hat"]] <- y_hat
  MODEL[["MBTR"]][["eval"]][["NSD"]] <- NSD

  # --- 7.2. Discard low quality models according to NSD
  # Fixed at 50% average across the projection
  if(MODEL[["MBTR"]][["eval"]][["NSD"]] > 0.5 | is.na(MODEL[["MBTR"]][["eval"]][["NSD"]])){
    MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != "MBTR"]
    message(paste("--- EVAL : discarded", "MBTR", "due to NSD =", MODEL[["MBTR"]][["eval"]][["NSD"]], "> 0.5 \n"))
  }

  return(MODEL)

} # END FUNCTION

