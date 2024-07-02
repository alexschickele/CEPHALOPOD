#' =============================================================================
#' @name model_binary
#' @description sub-pipeline corresponding to the model fitting procedure for
#' binary data. Called via the model_wrapper.R function, thus no default
#' parameter values are defined.
#' @param CALL the call object from the master pipeline.
#' @param QUERY the query object from the master pipeline
#' @return returns a list object containing the best model, associated hyper
#' parameters and predicted values per re sampling folds

model_binary <- function(CALL,
                         QUERY){

  # --- 1. Initialize
  # --- 1.1. Storage in MODEL object
  MODEL <- CALL$HP

  # --- 2. Loop with all selected models
  for(i in 1:length(CALL$MODEL_LIST)){
    # --- 2.1. Display information
    message(paste(Sys.time(), "--- Start hyper parameter tuning for", CALL$MODEL_LIST[i], "---"))

    # --- 2.2. Define the formula
    # General formula for all algorithms vs Specific for GAM (including the spline)
    if(CALL$MODEL_LIST[i] != "GAM"){tmp <- QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = " + ")
    } else {tmp <- paste0("s(", QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = ", k = 3) + s("), ", k = 3)")}
    formula <- paste0("measurementvalue ~ ", tmp) %>% as.formula()

    # --- 2.3. Define workflow adapted to hyper parameter tuning
    model_wf <- workflow() %>%
      add_variables(outcomes = "measurementvalue", predictors = QUERY$SUBFOLDER_INFO$ENV_VAR) %>%
      add_model(MODEL[[CALL$MODEL_LIST[i]]][["model_spec"]], formula = formula)

    # --- 2.4. Run the model for each fold x (hyper parameter grid rows)
    # Runs hyper parameter tuning if a grid is present in the HP (e.g. no GLM tune)
    if(!is.null(MODEL[[CALL$MODEL_LIST[i]]][["model_grid"]])){
      model_res <- model_wf %>%
        tune_grid(resamples = QUERY$FOLDS$resample_split,
                  grid = MODEL[[CALL$MODEL_LIST[i]]][["model_grid"]],
                  metrics = yardstick::metric_set(rmse),
                  control = control_grid(verbose = TRUE, allow_par = FALSE))
    } else {
      model_res <- model_wf %>%
        fit_resamples(resamples = QUERY$FOLDS$resample_split,
                      control = control_resamples(verbose = TRUE, allow_par = FALSE))
    }

    # --- 2.5. Select best hyper parameter set
    # Based on RMSE values per model run (rsq does not work with 0's)
    model_best <- model_res %>%
      select_best(metric="rmse")
    # Retrieve the corresponding RMSE as well
    rmse_best <- model_res %>% show_best(metric="rmse") %>% .[1,]
    MODEL[[CALL$MODEL_LIST[i]]][["best_fit"]] <- rmse_best

    # --- 2.6. Define final workflow
    final_wf <- model_wf %>%
      finalize_workflow(model_best)

    # --- 2.7. Run the model on the initial split
    # Save the workflow and model object in a list to be passed to further steps
    final_fit <- final_wf %>%
      last_fit(QUERY$FOLDS$init_split)

    MODEL[[CALL$MODEL_LIST[i]]][["final_wf"]] <- final_wf
    MODEL[[CALL$MODEL_LIST[i]]][["final_fit"]] <- final_fit

    # --- 2.8. Display information
    message(paste(Sys.time(), "--- DONE ---"))
  }

  return(MODEL)

} # END FUNCTION



