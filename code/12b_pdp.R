#' =============================================================================
#' @name pdp
#' @description Compute and plot partial dependency plots at the model and
#' ensemble level. The proportion data have their own PDP function embedded in
#' MBTR. (FOR NOW)
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.

pdp <- function(FOLDER_NAME = NULL,
                SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  set.seed(123)
  
  # --- 1.1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.2. Early return - if no model have fitted
  if(length(MODEL$MODEL_LIST) == 0){
    message("No PDP are computed for this species - no algorithms passed the QC or fitted")
    return(NULL)
  }
  
  # --- 1.3. Source the MBTR functions
  if(CALL$DATA_TYPE == "proportions"){
    library(reticulate)
    source_python(paste0(project_wd,"/function/mbtr_function.py"))
  } # if proportions
  
  # --- 1.4. Multivariate PDP function - for proportions
  mpartial <- function(object, train, pred_var, nboosts, grid_resolution=10, cores=10) {
    #' @param object a fitted model for which a predict() method exists
    #' @param train training set on which the model was fitted
    #' @param pred_var name (or index) of the variable of the training set for which the pdp is to be computed
    #' @param grid_resolution number of points along the values of `pred_var` at which the pdp is to be computed
    #' @param nboosts number of boosting rounds for the validated model to predict
    #' @return A data.frame with the values of `pred_var` and the predicted value of all response variables
    # define the grid of pred.var values at which to compute the pdp
    
    library(reticulate)
    source_python(paste0(project_wd,"/function/mbtr_function.py"))
    
    x <- dplyr::select(train, all_of(pred_var))
    grid <- seq(min(x, na.rm=T), max(x, na.rm=T), length.out=grid_resolution)
    
    # compute the pdp for each variable
    yhat <- parallel::mclapply(grid, function(val) {
      X <- train
      X[,pred_var] <- val %>% as.data.frame()
      this_yhat <- mbtr_predict(model = object, X_pred = X, n_boosts = nboosts)
      this_yhat <- apply(this_yhat, 2, mean)
    }, mc.cores=cores)
    
    # combine the results for all grid values
    yhat <- data.frame(x=grid, do.call(rbind, yhat))
    
    return(yhat)
  }
  
  # --- 1.5. For loop parameters
  if(CALL$DATA_TYPE == "proportions"){loop_over <- 1:ncol(QUERY$Y)
  }else{loop_over <- MODEL$MODEL_LIST}
  
  # --- 1.6. Initialize global PDP storage
  pdp_all <- NULL
  
  # --- 2. Define bootstraps
  # --- 2.1. Target transformation
  if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
    message("PROJ --- Transforming the target variable according to the provided function")
    source(CALL$TARGET_TRANSFORMATION)
    tmp <- target_transformation(QUERY$Y$measurementvalue, REVERSE = FALSE)
    Y <- data.frame(tmp$out)
    colnames(Y) <- "measurementvalue"
    QUERY[["target_transformation"]][["yj_obj"]] <- tmp$yj_obj
  } else {
    Y <- QUERY$Y
  }
  # --- 2.2. Re-assemble all query tables
  tmp <- cbind(Y, QUERY$X, QUERY$S)
  
  # --- 2.3. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = CALL$N_BOOTSTRAP)
  
  # --- 2.4. Save bootstrap on disk for proportions
  if(CALL$DATA_TYPE == "proportions"){
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
  } # if proportions
  
  # ----------------------------------------------------------------------------
  for(i in MODEL$MODEL_LIST){
    # --- 3. Fit the bootstrap
    if(CALL$DATA_TYPE == "proportions"){
      # --- 3.1. Fit MBTR on bootstrap - for proportions only
      # --- 3.1.1. Source the MBTR functions
      library(reticulate)
      source_python(paste0(project_wd,"/function/mbtr_function.py"))
      # --- 3.1.2. Compute the fit
      message(paste0(Sys.time(), "--- MBTR: re-fit model for bootstrap"))
      boot_fit <- mcmapply(FUN = mbtr_fit,
                           path = paste0(project_wd, "/data/MBTR_cache/",1:CALL$N_BOOTSTRAP,"_"),
                           hp_id = "0",
                           loss_type='mse',
                           n_boosts = as.integer(1000),
                           min_leaf= MODEL$MBTR$final_wf$MEAN_LEAF,
                           learning_rate=MODEL$MBTR$final_wf$LEARNING_RATE,
                           lambda_weights=MODEL$MBTR$final_wf$LEARNING_RATE/100,
                           lambda_leaves=0,
                           n_q= as.integer(MODEL$MBTR$final_wf$N_Q),
                           early_stopping_rounds = 10,
                           SIMPLIFY = FALSE,
                           USE.NAMES = FALSE,
                           mc.cores = CALL$N_BOOTSTRAP)
      
      # --- 3.1.3. Reload them in R because of "previous session invalidity"
      for(b in 1:CALL$N_BOOTSTRAP){
        boot_fit[[b]] <- py_load_object(paste0(project_wd, "/data/MBTR_cache/",b,"_0_m"), pickle = "pickle")[[1]]
      } # boot_fit reload
      
    } else {
      # --- 3.2. Fit tidymodel on bootstrap
      # We use control_resample > extract_fit_parsnip to extract the fit object
      boot_fit <- MODEL[[i]][["final_wf"]] %>% 
        fit_resamples(resamples = boot_split,
                      control = control_resamples(extract = function (x) extract_fit_parsnip(x))) %>% 
        unnest(.extracts)
    } # if proportions or not
    
    # --- 4. Compute the raw partial dependency 
    # --- 4.1. Initialize the pdp object
    pdp_m <- NULL
    message(paste(Sys.time(), "--- PDP : computing for", i))
    
    for(b in 1:CALL$N_BOOTSTRAP){
      if(CALL$DATA_TYPE == "proportions"){
        
        # --- 4.2. Compute PDP for MBTR - proportions data
        # --- 4.2.1 Loop over environmental variables
        pdp_b <- NULL
        for(j in QUERY$SUBFOLDER_INFO$ENV_VAR){
          tmp <- mpartial(object = boot_fit[[b]],
                          pred_var = j,
                          grid_resolution = 20,
                          train = QUERY$X[QUERY$SUBFOLDER_INFO$ENV_VAR],
                          nboosts = as.integer(MODEL$MBTR$final_wf$NBOOST),
                          cores = 20)
          colnames(tmp) <- c(j, CALL$SP_SELECT)
          pdp_b <- abind(pdp_b, tmp, along = 3)
        } # j env variable
        
        # --- 4.2.2. Aggregation at the model level
        pdp_m <- abind(pdp_m, pdp_b, along = 4)
        
      } else {
        # --- 4.3. Compute PDP for tidymodels
        # --- 4.3.1. Build model explainer
        explainer <- explain_tidymodels(model = boot_fit$.extracts[[b]],
                                        data = QUERY$X[QUERY$SUBFOLDER_INFO$ENV_VAR],
                                        y = QUERY$Y,
                                        verbose = FALSE)
        
        # --- 4.3.2. Do the computation
        tmp <- model_profile(explainer = explainer,
                             variables = QUERY$SUBFOLDER_INFO$ENV_VAR,
                             N = 100,
                             type = "partial") 
        tmp <- tmp[["agr_profiles"]] %>% 
          dplyr::select(-"_label_", -"_ids_")
        
        # --- 4.3.3. Aggregation at the model level
        names(tmp) <- c("var","x","yhat") # Clean naming
        if(b == 1){pdp_m <- tmp
        }else {pdp_m <- cbind(pdp_m, tmp$yhat)}
        
      } # if proportions or not
    } # End j bootstrap loop
    
    # --- 4.4. Reverse transformation
    if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
      message("--- PROJ : reverse transformation of the target")
      pdp_m[,-c(1,2)] <- apply(pdp_m[,-c(1,2)], 2, function(x){
        x <- target_transformation(x, REVERSE = TRUE, PARAM = QUERY$target_transformation)
      }) %>% as.data.frame()
    } # end if transformation
    
    # --- 5. Compute mean and coefficient of variation
    if(CALL$DATA_TYPE == "proportions"){
      # --- 5.1. Compute for MBTR - proportion data
      pdp_all <- apply(pdp_m, 2, list)[-1] %>% 
        lapply(function(x)(x = data.frame(var = rep(QUERY$SUBFOLDER_INFO$ENV_VAR, each = dim(x[[1]])[1]),
                                          y_hat_m = apply(x[[1]], c(1,2), mean) %>% as.vector(),
                                          y_hat_cv = apply(x[[1]], c(1,2), function(z)(z = sd(z)/mean(z)*100)) %>% as.vector())))
      
      
      x <- apply(pdp_m, 2, list)[[1]] %>% 
        .[[1]] %>% 
        .[,,1] %>% 
        as.vector()
      
      pdp_all <- lapply(pdp_all, function(z)(z = cbind(z, x)))
      MODEL$MBTR$pdp_all <- pdp_all
      
    } else {
      # --- 5.2. Compute for tidymodels
      # --- 5.2.1. At the model level
      tmp_m <- apply(pdp_m[,-c(1,2)], 1, mean)
      tmp_cv <- apply(pdp_m[,-c(1,2)], 1, cv)
      pdp_m <- cbind(pdp_m[,c(1,2)], tmp_m, tmp_cv)
      colnames(pdp_m) <- c("var","x","y_hat_m","y_hat_cv")
      
      # --- 5.2.2. Write in model object
      MODEL[[i]][["pdp"]] <- pdp_m
      
      # --- 5.2.3. Aggregate in an all_pdp object
      pdp_all[[i]] <- pdp_m
    } # if proportions or not
  } # End model loop
  
  # --- 5.3. Save model object
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"),
       compress = "gzip", compression_level = 6)
  
  # ----------------------------------------------------------------------------
  
  # --- 6. Plotting PDPs
  # --- 6.1 Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/06_pdp.pdf"))
  par(mfrow = c(4,3), mar = c(4,4,3,4))
  
  # --- 6.2. Initialize 
  # --- 6.2.1. Color palette
  if(CALL$DATA_TYPE == "proportions"){m_names <- names(pdp_all)
  } else {m_names <- MODEL$MODEL_LIST}
  pal <- brewer.pal(length(m_names), "Paired")
  if(CALL$DATA_TYPE != "proportions"){pal <- pal[which(m_names %in% MODEL$MODEL_LIST)]}
  
  # --- 6.2.2. Plot scaling (average Q95 across bootstrap)
  # Because continuous data are not between 0 and 1
  if(CALL$DATA_TYPE == "continuous"){
    plot_scale <- lapply(MODEL$MODEL_LIST, FUN = function(z){
      z = MODEL[[z]]$proj$y_hat %>% apply(1, function(x)(x = mean(x, na.rm = TRUE)))
    }) %>% 
      unlist() %>% 
      quantile(0.95, na.rm = TRUE)
    
  } else {
    plot_scale <- 1
  }
  
  
  # --- 6.3. Iteratively compute the plots
  for(i in QUERY$SUBFOLDER_INFO$ENV_VAR){
    for(j in 1:length(pal)){
      # --- 6.3.1. Prepare data table
      tmp <- pdp_all[[j]] %>% 
        dplyr::filter(var == i)
      
      # --- 6.3.2. Plot
      if(j == 1){
        plot(tmp$x, tmp$y_hat_m, type = 'l', lwd = 1, col = pal[j],
             ylim = c(0, plot_scale),
             xlab = "Environmental value", ylab = "Norm. response", main = paste("PARTIAL DEPENDENCE \n (", i, ")"), cex.main = 0.8)
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = scales::alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = m_names[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
        grid(col = "gray20")
        box("figure", col="black", lwd = 1)
      } else {
        lines(tmp$x, tmp$y_hat_m, type = 'l', ylim = c(0,1), lwd = 1, col = pal[j],
              xlab = "", ylab = "", main = "")
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = scales::alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = m_names[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
      } # End if
    } # End j model/target loop
  } # End i plot loop
  
  # --- 7. Wrap up and save
  # --- 7.1. Stop PDF saving
  dev.off()
  # --- 7.2. Pretty return
  return(SUBFOLDER_NAME)
  
} # END FUNCTION

