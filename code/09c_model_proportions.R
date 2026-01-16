#' =============================================================================
#' @name model_proportions
#' @description sub-pipeline corresponding to the model fitting procedure for
#' proportions data. Called via the model_wrapper.R function, thus no default
#' parameter values are defined.
#' @param CALL the call object from the master pipeline. 
#' @param QUERY the query object from the master pipeline
#' @return returns a path to the python model file

model_proportions <- function(CALL,
                              QUERY){
  
  # --- 1. Initialize function
  # --- 1.1. Storage in MODEL object
  MODEL <- CALL$HP
  
  # --- 1.2. Source the MBTR functions
  library(reticulate)
  source_python(paste0(project_wd,"/function/mbtr_function.py"))
  
  # --- 1.3. Create a clean cache folder in data
  # i.e. to write train and test files to pass to python (necessary for parallel)
  dir.create(paste0(project_wd, "/data/MBTR_cache"))
  to_clean <- list.files(paste0(project_wd, "/data/MBTR_cache"), full.names = TRUE)
  file.remove(to_clean)
  
  # --- 2. Run the model for each fold x hyper parameter
  # --- 2.1. Initialize loss & hp object
  loss <- list()
  hp <- 1:nrow(MODEL$MBTR$model_grid)
  
  # --- 2.2. Initialize the hyperparameter fit function
  # i.e. to use mcmapply across hyperparameters for faster computing
  mbtr_hp_fit <- function(hp, path){
    source_python(paste0(project_wd,"/function/mbtr_function.py"))
    m <- mbtr_fit(path,
                  hp_id = as.character(hp),
                  loss_type='mse',
                  n_boosts = as.integer(10000),
                  min_leaf= MODEL$MBTR$model_grid$MEAN_LEAF[hp],
                  learning_rate=MODEL$MBTR$model_grid$LEARNING_RATE[hp],
                  lambda_weights=MODEL$MBTR$model_grid$LEARNING_RATE[hp]/100,
                  lambda_leaves= 0,
                  n_q= as.integer(MODEL$MBTR$model_grid$N_Q[hp]),
                  early_stopping_rounds = 10)
    # tmp <- m[[2]] %>% unlist()
    m[[2]] <- m[[2]] %>% unlist()
    return(m)
  } # end function
  
  # --- 2.3. Do the resample fit
  for(cv in 1:CALL$NFOLD){
    # --- 2.3.1. Load .py functions & instance
    # NOTE: I know this is ugly but reticulate loses the connection at every parallel
    # worker AND every for loop turn... therefore it is re-sourced many times in
    # the script
    
    library(reticulate)
    source_python(paste0(project_wd,"/function/mbtr_function.py"))
    
    # --- 2.3.2. Prepare the inputs
    message(paste0(Sys.time(), "--- MBTR: writing data to cache for cv = ", cv))
    X_tr <- QUERY$FOLDS$resample_folds[[cv]]$analysis %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_tr, paste0(project_wd, "/data/MBTR_cache/", cv, "_X_tr.feather"))
    
    Y_tr <- QUERY$FOLDS$resample_folds[[cv]]$analysis %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_tr, paste0(project_wd, "/data/MBTR_cache/", cv, "_Y_tr.feather"))
    
    X_val <- QUERY$FOLDS$resample_folds[[cv]]$assessment %>%
      dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
    write_feather(X_val, paste0(project_wd, "/data/MBTR_cache/", cv, "_X_val.feather"))
    
    Y_val <- QUERY$FOLDS$resample_folds[[cv]]$assessment %>%
      dplyr::select(as.character(CALL$SP_SELECT))
    write_feather(Y_val, paste0(project_wd, "/data/MBTR_cache/", cv, "_Y_val.feather"))
    
    # --- 2.3.3. Fitting hyperparameters
    message(paste0(Sys.time(), "--- MBTR: fit hyperparameters for cv = ", cv))
    loss[[cv]] <- mcmapply(FUN = mbtr_hp_fit,
                           path = paste0(project_wd, "/data/MBTR_cache/",cv,"_"),
                           hp = hp,
                           mc.cores = 20)
    message("Algorithm fitting - DONE")
  } # cv loop
  
  # --- 2.4. Compute the minimum loss and corresponding nb. of boost rounds
  min_loss <- nboost <- list()
  for(i in 1:nrow(MODEL$MBTR$model_grid)){
    max_boost <- lapply(loss, function(x)(x = length(x[2,][[i]]))) %>% unlist() %>% min()
    min_loss[[i]] <- lapply(loss, function(x)(x = x[2,][[i]][1:max_boost])) %>% as.data.frame() %>% apply(1, mean) %>% min()
    nboost[[i]] <- which(lapply(loss, function(x)(x = x[2,][[i]][1:max_boost])) %>% as.data.frame() %>% apply(1, mean) == min_loss[[i]])
  }
  best_hp <- which(unlist(min_loss) == min(unlist(min_loss)))[1] # Takes the first in case there is two equal
  nboost <- nboost[[best_hp]]
  
  # --- 2.5. Retrieve the corresponding RMSE as well
  MODEL[["MBTR"]][["best_fit"]] <- sqrt(min_loss[[best_hp]])

  # --- 3. Write best model in a file
  # --- 3.1. Best workflow
  MODEL[["MBTR"]][["final_wf"]] <- MODEL$MBTR$model_grid[best_hp,] %>% mutate(NBOOST = nboost)
  
  # --- 3.2. Best model
  # It is already written by default in the training process
  # We just write the path to the data cache corresponding to the trained models
  MODEL[["MBTR"]][["final_fit"]] <- paste0(project_wd, "/data/MBTR_cache/",1:CALL$NFOLD, "_", best_hp)
  return(MODEL)
  
} # END FUNCTION
