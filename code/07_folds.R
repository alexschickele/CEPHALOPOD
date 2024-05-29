#' =============================================================================
#' @name folds
#' @description this functions adds an FOLD list to the query. It corresponds to
#' the different folds between train and test splits and associated parameters
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return Updates the FOLDS in a QUERY.RData file

folds <- function(FOLDER_NAME = NULL,
                  SUBFOLDER_NAME = NULL){

  # --- 1. Initialize function
  if(!exists("SEED")){
    SEED <- 123
  }
  set.seed(SEED)
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : folds ********************"))

  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 1.3. Target transformation
  if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
    message("--- FOLDS : Transforming the target variable according to the provided function")
    source(CALL$TARGET_TRANSFORMATION)
    tmp <- target_transformation(x = QUERY$Y$measurementvalue, REVERSE = FALSE)
    Y <- data.frame(tmp$out)
    colnames(Y) <- "measurementvalue"
    QUERY[["target_transformation"]][["yj_obj"]] <- tmp$yj_obj
  } else {
    Y <- QUERY$Y
  }

  # --- 2. Initial split
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(Y, QUERY$X, QUERY$S)

  # --- 2.2. Do the initial split
  # The initial split is done on 30% of the data, stratified,
  # to have a well representative evaluation sample
  # --- 2.2.1. For univariate data - strata is possible so we do it
  if(CALL$DATA_TYPE != "proportions"){
    init_split <- tmp %>%
      initial_split(prop = 0.7,
                    strata = measurementvalue)
  }
  # --- 2.2.2. For multivariate data - strata is not possible
  if(CALL$DATA_TYPE == "proportions"){
    init_split <- tmp %>%
      initial_split(prop = 0.7)
  }

  # --- 2.2.3. Append FOLD object
  QUERY$FOLDS$init_split <- init_split
  QUERY$FOLDS$test <- testing(init_split)
  QUERY$FOLDS$train <- training(init_split)

  # --- 3. Train set resampling
  # --- 3.1. Parameter check
  if(CALL$FOLD_METHOD != "kfold" & CALL$FOLD_METHOD != "lon"){
    stop("FOLD_METHOD not implemented or incorrect. It should be 'kfold' or 'lon'")
  }

  # --- 3.2. Normal k-fold re sampling
  # --- 3.2.1. For univariate data - strata is possible so we do it
  if(CALL$FOLD_METHOD == "kfold" & CALL$DATA_TYPE != "proportions"){
    folds <- vfold_cv(data = QUERY$FOLDS$train,
                      strata = measurementvalue,
                      v = CALL$NFOLD)
  }
  # --- 3.2.2. For multivariate data - strata is not possible
  if(CALL$FOLD_METHOD == "kfold" & CALL$DATA_TYPE == "proportions"){
    folds <- vfold_cv(data = QUERY$FOLDS$train,
                      v = CALL$NFOLD)
  }

  # --- 3.3. Longitudinal block re sampling
  if(CALL$FOLD_METHOD == "lon"){
    folds <- group_vfold_cv(data = QUERY$FOLDS$train,
                            group = c(decimallongitude),
                            v = CALL$NFOLD)
  }

  # --- 4. Append QUERY and CALL objects
  # S$id keeps track of the initial row numbers within each split/re sample
  # --- 4.1. Append QUERY
  QUERY$FOLDS$resample_split <- folds
  for(i in 1:nrow(folds)){
    fold_name <- folds$id[i]
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["assessment"]] <- folds$splits[[i]] %>% assessment()
    QUERY$FOLDS[["resample_folds"]][[fold_name]][["analysis"]] <- folds$splits[[i]] %>% analysis()
  }

  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 5.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
