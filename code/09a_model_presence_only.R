#' =============================================================================
#' =============================================================================
#' @name model_presence_only
#' @description Sub-pipeline for model fitting on presence-only data.
#' It is invoked via the model_wrapper.R function, hence no default parameters.
#' 
#' @param CALL the call object from the master pipeline. 
#' @param QUERY the query object from the master pipeline.
#' 
#' @return A list containing the best model, hyperparameters, and predictions
#' for each resampling fold.
#' 

model_presence_only <- function(CALL, QUERY){
  
  # --- 1. Initialize
  # --- 1.1. Store hyperparameter set in MODEL object
  MODEL <- CALL$HP
  
  # --- 1.2. Get model list
  model_name <- CALL$MODEL_LIST
  rm(CALL) # Remove CALL for memory use
  gc()

  # --- 2. Loop with all selected models
  for(i in seq_along(model_name)){
    
    message(paste(Sys.time(), "--- Start hyper parameter tuning for", model_name[i], "---"))
    
    # --- 2.1. Define the formula
    # General formula for most models; special handling for GAM
    if(model_name[i] != "GAM"){
      formula_vars <- QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = " + ")
    } else {
      # GAM requires spline terms
      formula_vars <- paste0("s(", QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = ", k = 3) + s("), ", k = 3)")
    }
    formula <- paste0("measurementvalue ~ ", formula_vars) %>% as.formula()
    
    # --- 2.2. Define workflow adapted to hyper parameter tuning
    model_wf <- workflow() %>% 
      add_variables(outcomes = "measurementvalue", predictors = QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
      add_model(MODEL[[model_name[i]]][["model_spec"]], formula = formula)

    # --- 2.3. Run the model for each fold x (hyper parameter grid rows)
    # Runs hyper parameter tuning if a grid is present in the HP (e.g. no GLM tune)
    if(!is.null(MODEL[[model_name[i]]][["model_grid"]])){
      
      model_res <- model_wf %>% 
        tune_grid(resamples = QUERY$FOLDS$resample_split,
                  grid = MODEL[[model_name[i]]][["model_grid"]],
                  metrics = yardstick::metric_set(rmse), 
                  control = control_grid(verbose = TRUE, allow_par = FALSE))
    } else {
      model_res <- model_wf %>% 
        fit_resamples(resamples = QUERY$FOLDS$resample_split,
                      control = control_resamples(verbose = TRUE, allow_par = FALSE))
    }

    # --- 2.4. Select best hyper parameter set
    # Based on RMSE values per model run (rsq does not work with 0's)
    model_best <- model_res %>% select_best(metric = "rmse")
    # Retrieve the corresponding RMSE as well
    MODEL[[model_name[i]]][["best_fit"]] <- model_res %>% show_best(metric = "rmse") %>% .[1,]
    
    # --- 2.5. Define final workflow
    final_wf <- model_wf %>% finalize_workflow(model_best)
    
    # --- 2.6. Run the model on same cross validation splits
    # We have one fit per cross validation saved in a list, to be passed to further steps
    final_fit <- lapply(seq_along(QUERY$FOLDS$resample_folds), function(x){
      out <- final_wf %>% last_fit(QUERY$FOLDS$resample_split$splits[[x]])
      return(out)
    }) # loop over the same cross validation folds
    
    MODEL[[model_name[i]]][["final_wf"]] <- final_wf
    MODEL[[model_name[i]]][["final_fit"]] <- final_fit
    
    # --- 2.7. Display information
    message(paste(Sys.time(), "--- DONE ---"))
  }

  return(MODEL)
  
} # END FUNCTION



