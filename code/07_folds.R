#' =============================================================================
#' @name folds
#' @description this functions adds an FOLD list to the query. It corresponds to
#' the different folds between train and test splits and associated parameters
#' 
#' @param FOLDER_NAME A character string indicating the folder containing the output.
#' @param SUBFOLDER_NAME A list of subfolders to parallelize the process on.
#' 
#' @return Updates the FOLDS in a QUERY.RData file

folds <- function(FOLDER_NAME = NULL,
                  SUBFOLDER_NAME = NULL){

  # --- 1. Initialize function
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Start logging process
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : folds ********************"))
  
  # --- 1.2. Load metadata and previous query
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 1.3. Target transformation
  if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
    message("--- FOLDS : Transforming the target variable according to the provided function")
    source(CALL$TARGET_TRANSFORMATION)
    tmp <- target_transformation(x = QUERY$Y$measurementvalue, REVERSE = FALSE)
    Y <- data.frame(tmp$out)
    colnames(Y) <- "measurementvalue"
    QUERY[["target_transformation"]][["yj_obj"]] <- tmp$yj_obj
    
    rm(tmp)
    gc()
    
  } else {
    Y <- QUERY$Y
  }

  # --- 2. Define cross-validation split
  # We use the same for both hyper parameter tuning and model evaluation
  
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(Y, QUERY$X, QUERY$S)
  
  # --- 2.2. Parameter check
  if(CALL$FOLD_METHOD != "kfold" & CALL$FOLD_METHOD != "lon"){
    stop("FOLD_METHOD not implemented or incorrect. It should be 'kfold' or 'lon'")
  }

  # --- 2.3. Normal k-fold re sampling
  # --- 2.3.1. For univariate data - strata is possible so we do it
  if(CALL$FOLD_METHOD == "kfold" & CALL$DATA_TYPE != "proportions"){
    folds <- vfold_cv(data = tmp,
                      strata = measurementvalue,
                      v = CALL$NFOLD)
  }
  
  # --- 2.3.2. For multivariate data - strata is not possible
  if(CALL$FOLD_METHOD == "kfold" & CALL$DATA_TYPE == "proportions"){
    folds <- vfold_cv(data = tmp,
                      v = CALL$NFOLD)
  }

  # --- 2.4. Longitudinal block re sampling
  if(CALL$FOLD_METHOD == "lon"){
    folds <- group_vfold_cv(data = tmp,
                            group = c(decimallongitude),
                            v = CALL$NFOLD)
  }

  # --- 3. Append QUERY and CALL objects
  # S$id keeps track of the initial row numbers within each split/re sample
  
  # --- 3.1. Append QUERY
  QUERY$FOLDS$resample_split <- folds
  for(i in 1:nrow(folds)){
    fold_name <- folds$id[i]
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["assessment"]] <- folds$splits[[i]] %>% assessment()
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["analysis"]] <- folds$splits[[i]] %>% analysis()
  }

  # --- 4. Wrap up and save
  # --- 4.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 4.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 4.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
