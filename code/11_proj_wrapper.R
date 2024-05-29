#' =============================================================================
#' @name proj_wrapper
#' @description wrapper function redirecting towards the projection sub-pipeline
#' corresponding to the type of data.
#' @description In case of proportion data, an input converting section is run,
#' to properly pass the inputs to python library MBTR
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.
#' @return outputs are saved in the MODEL.RData object

proj_wrapper <- function(FOLDER_NAME = NULL,
                         SUBFOLDER_NAME = NULL){
  set.seed(123)
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : proj_wrapper ********************"))

  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))

  # --- 2. Redirection to the PRESENCE model projections
  if(CALL$DATA_TYPE == "binary"){

    # --- 2.1. Run function
    MODEL <- proj_binary(QUERY = QUERY,
                         MODEL = MODEL,
                         CALL = CALL)
  } # END if binary

  # --- 3. Redirection to the CONTINUOUS model projections
  if(CALL$DATA_TYPE == "continuous"){

    # --- 3.1. Run function
    MODEL <- proj_continuous(QUERY = QUERY,
                             MODEL = MODEL,
                             CALL = CALL)
  } # END if continuous

  # --- 4. Redirection to the PROPORTIONS model projections
  if(CALL$DATA_TYPE == "proportions"){

    # --- 4.1. Run function
    MODEL <- proj_proportions(QUERY = QUERY,
                              MODEL = MODEL,
                              CALL = CALL)
  } # END if proportions

  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/MODEL.RData"),
       compress = "gzip", compression_level = 6)
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 5.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
