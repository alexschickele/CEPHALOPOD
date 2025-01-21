#' =============================================================================
#' @name eval_presence_only
#' @description sub-pipeline for model evaluation corresponding to presence_only data
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @param PDF_PATH preset path for the pdf save
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return variable importance plots as PDF file

eval_presence_only <- function(CALL, QUERY, MODEL, PDF_PATH){

  # --- 1. Model performance assessment
  for(i in MODEL$MODEL_LIST){
    
    # --- 1.1. Load final model data
    # Loop over the cross validation runs
    
    model_data <- lapply(1:(length(MODEL[[i]][["final_fit"]])), function(x){
      # Extract and collect final predictions
      final_fit <- MODEL[[i]][["final_fit"]][[x]] %>%
        collect_predictions()
      
      # Create a data frame for y (obs) and y_hat (predicted)
      df <- data.frame(y = final_fit$measurementvalue, 
                       y_hat = final_fit$.pred)
      return(df)
    }) %>% bind_rows()
    
    # --- 1.2. Compute Continuous Boyce Index into MODELS object
    MODEL[[i]][["eval"]][["CBI"]] <- ecospat.boyce(fit = model_data$y_hat,
                                                   obs = model_data$y_hat[which(model_data$y == 1)], # predictions at observed locations
                                                   PEplot = FALSE,
                                                   method = "pearson",
                                                   rm.duplicate = FALSE,
                                                   res = 100) %>%
      .$cor
    
    # --- 1.3. Memory cleanup after each iteration
    rm(model_data)
    gc()

  } # for each model loop

  # --- 2. Variable importance - algorithm level
  # --- 2.1. Initialize function
  # --- 2.1.1. Storage
  var_imp <- list()
  
  # --- 2.1.2. General features - used later for plots
  features <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[1]]$data %>%
    dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  
  # --- 2.1.3. Loop over models
  for(i in MODEL$MODEL_LIST){
    
    # --- 2.2. Loop over the cross-validations
    var_imp[[i]][["Raw"]] <- lapply(1:CALL$NFOLD, function(x){
      
      # --- 2.2.1. Extract related target and feature
      # Model and cross-validation specific
      id <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]][["in_id"]]
      
      features_x <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]]$data[id,] %>%
        dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
      
      target <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]]$data[id,] %>%
        dplyr::select(measurementvalue)
      
      # --- 2.2.2. Extract final model fit
      m <- extract_fit_parsnip(MODEL[[i]][["final_fit"]][[x]])
      
      # --- 2.2.3. Build model explainer
      explainer <- explain_tidymodels(model = m,
                                      data = features_x,
                                      y = target)
      
      # --- 2.2.4. Intermediate cleanup
      rm(features_x, target, m)
      gc()
      
      # --- 2.2.5. First in terms of RMSE, i.e., raw var importance for later ensemble computing
      message(paste("--- VAR IMPORTANCE : compute for", i))
      out <- model_parts(explainer = explainer,
                         loss_function = loss_root_mean_square) %>%
        dplyr::filter(permutation != 0) %>%
        dplyr::filter(variable != "_baseline_") %>%
        group_by(permutation) %>%
        mutate(value = dropout_loss - dropout_loss[variable == "_full_model_"]) %>%
        dplyr::filter(variable != "_full_model_") %>% 
        ungroup() %>% 
        mutate(cv = x) # add the cv information for percentage
      
      return(out)
      
    }) %>% bind_rows() # end cross validation loop
    
    # --- 2.3. Further compute it as percentage for model-level plot
    # Security if a model did not fit, hence var_imp = 0 for all predictors
    # Avoids an error leading to a function stop, while other models could be OK
    if(sum(var_imp[[i]][["Raw"]][["value"]]) > 0){
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]] %>%
        group_by(cv, permutation) %>% 
        mutate(value = value / sum(value) * 100)  %>%
        ungroup() %>%
        dplyr::select(variable, value) %>%
        dplyr::filter(!is.na(value)) %>%
        mutate(variable = fct_reorder(variable, value, .desc = TRUE))
    } else {
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]]
    }
    
    # --- 2.4. Compute cumulative variable importance
    MODEL[[i]][["eval"]][["CUM_VIP"]] <- var_imp[[i]][["Percent"]] %>%
      group_by(variable) %>%
      summarise(average = mean(value)) %>%
      dplyr::slice(1:3) %>%
      dplyr::select(average) %>%
      sum()
    
    # --- 2.5. Save row VIP
    MODEL[[i]][["vip"]] <- var_imp[[i]][["Percent"]]
    
  } # for each model loop
  
  # --- 3. Removing low quality algorithms
  for(i in MODEL$MODEL_LIST){
    # --- 3.1. Based on model performance
    # Fixed at 0.3 for CBI value or NA (in case of a 0 & 1 presence_only model prediction)
    if(MODEL[[i]][["eval"]][["CBI"]] < 0.5 | is.na(MODEL[[i]][["eval"]][["CBI"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CBI =", MODEL[[i]][["eval"]][["CBI"]], "< 0.5 \n"))
    }

    # --- 3.2. Based on cumulative variable importance
    # Fixed at 30% cumulative importance for the top three predictors
    if(MODEL[[i]][["eval"]][["CUM_VIP"]] < 50 | is.na(MODEL[[i]][["eval"]][["CUM_VIP"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CUM_VIP =", MODEL[[i]][["eval"]][["CUM_VIP"]], "< 50% \n"))
    }
  } # for each model loop

  # --- 4. Variable importance - Ensemble level
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    # --- 4.1. Aggregate and weight raw data
    ens_imp <- NULL
    message("--- VAR IMPORTANCE : compute ensemble")

    for(i in MODEL$MODEL_LIST){
      # Concatenate raw variable importance
      ens_imp <- rbind(ens_imp, var_imp[[i]][["Raw"]])
    } # End i model loop

    # --- 4.2. Further compute it as percentage
    var_imp[["ENSEMBLE"]][["Percent"]] <- ens_imp %>%
      group_by(cv, permutation) %>%
      mutate(value = value / sum(value) * 100 * length(MODEL$MODEL_LIST))  %>%
      ungroup() %>%
      dplyr::select(variable, value) %>%
      mutate(variable = fct_reorder(variable, value, .desc = TRUE))
  } # End if ensemble TRUE

  # --- 5. Build the ensemble QC if there is an ensemble
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    # --- 5.1. Ensemble CBI
    MODEL[["ENSEMBLE"]][["eval"]][["CBI"]] <- lapply(MODEL$MODEL_LIST,
                                                     FUN = function(x){
                                                       x <- MODEL[[x]]$eval$CBI
                                                     }) %>%
      unlist() %>% mean()

    # --- 5.2. Ensemble cumulative VIP
    MODEL[["ENSEMBLE"]][["eval"]][["CUM_VIP"]] <- var_imp[["ENSEMBLE"]][["Percent"]]  %>%
      group_by(variable) %>%
      summarise(average = mean(value)) %>%
      dplyr::slice(1:3) %>%
      dplyr::select(average) %>%
      sum()

    # --- 5.3. Ensemble raw VIP
    MODEL[["ENSEMBLE"]][["vip"]] <- var_imp[["ENSEMBLE"]][["Percent"]]
  } # End ENSEMBLE QC

  # --- 6. Variable importance - Plot

  # --- 6.2. Define plots to display
  # All if FAST == FALSE; those that passed QC if there is more than 1
  if(CALL$FAST == FALSE){
    plot_display <- CALL$HP$MODEL_LIST
  } else if(length(MODEL$MODEL_LIST) >= 1){
    plot_display <- MODEL$MODEL_LIST
  } else {
    plot_display <- NULL
  }

  # --- 6.3. Plot algorithm level and ensemble variable importance
  if(is.null(plot_display) == FALSE){
    # --- 6.1. Graphical specification
    pdf(PDF_PATH)
    par(mfrow = c(3,1), mar = c(6,2,3,20))
    # --- 6.3.1. Algorithm level plot
    for(i in plot_display){
      # Define the color (green = QC passed; red = no)
      if(i %in% MODEL$MODEL_LIST == TRUE){pal <- "#1F867B"
      } else {pal <- "#B64A60"}

      # Do the plot
      tmp <- var_imp[[i]][["Percent"]]
      boxplot(tmp$value ~ tmp$variable, axes = FALSE, horizontal = TRUE,
              main = paste("PREDICTOR IMPORTANCE (", i, ")"),
              sub = paste("Predictive performance (CBI) =", round(MODEL[[i]]$eval$CBI, 2), "; Cumulated var. importance (%; top 3) =", round(MODEL[[i]]$eval$CUM_VIP, 0)),
              col = pal, ylab = "", xlab = "Variable importance (%)")
      axis(side = 4, at = 1:ncol(features), labels = levels(tmp$variable), las = 2, cex.axis = 0.6)
      axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10))
      abline(v = seq(0, 100, 10), lty = "dotted")
      box()
      box("figure", col = "black", lwd = 1)
    } # End i model loop

    # --- 6.3.2. Ensemble level plot
    if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
      tmp <- var_imp[["ENSEMBLE"]][["Percent"]]
      boxplot(tmp$value ~ tmp$variable, axes = FALSE, horizontal = TRUE,
              main = "PREDICTOR IMPORTANCE ( Ensemble )",
              sub = paste("Predictive performance (CBI) =", round(MODEL[["ENSEMBLE"]]$eval$CBI, 2), "; Cumulated var. importance (%; top 3) =", round(MODEL[["ENSEMBLE"]]$eval$CUM_VIP, 0)),
              col = "gray50",
              ylab = "", xlab = "Variable importance (%)")
      axis(side = 4, at = 1:ncol(features), labels = levels(tmp$variable), las = 2, cex.axis = 0.6)
      axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10))
      abline(v = seq(0, 100, 10), lty = "dotted")
      box()
      box("figure", col = "black", lwd = 1)


    } # End if ensemble TRUE

  dev.off()
  } # End if model list > 1 or FAST == FALSE

  # --- 7. Wrap up and save

  return(MODEL)

} # END FUNCTION
