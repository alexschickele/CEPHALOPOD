#' =============================================================================
#' @name model_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @description In case of proportion data, an input converting section is run,
#' to properly pass the inputs to python library MBTR
#' 
#' @param FOLDER_NAME Name of the main working folder.
#' @param SUBFOLDER_NAME Subfolder name to parallelize on.
#' 
#' @return a model list object containing the different model objects
#' @return in case of proportion data, the model list object contains the path
#' to the model files as it cannot be passed as an object in memory
#' @return outputs are saved in a MODEL.RData object

# TO DO : implement the input converter for MBTR

model_wrapper <- function(FOLDER_NAME = NULL,
                          SUBFOLDER_NAME = NULL){

  # --- 1. Initialize function and environment
  set.seed(123)  # Ensure reproducibility
  
  # --- 1.1. Start logging (open a new log file for this subfolder)
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  
  message(Sys.time(), "******************** START : model_wrapper ********************")
  
  # --- 1.2. Load essential metadata
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))  # Load metadata for the current run
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))  # Load query-specific metadata

  # --- 2. Redirection to presence_only model
  if(CALL$DATA_TYPE == "presence_only"){

    # --- 2.1. Run function
    MODEL <- model_presence_only(CALL, QUERY = QUERY)
  } # END if presence_only

  # --- 3. Redirection to CONTINUOUS model
  if(CALL$DATA_TYPE == "continuous"){

    # --- 3.1. Run function
    MODEL <- model_continuous(CALL, QUERY = QUERY)
  } # END if continuous

  # --- 4. Redirection to PROPORTIONS model
  if(CALL$DATA_TYPE == "proportions"){

    # --- 4.1.. Run function
    MODEL <- model_proportions(CALL, QUERY = QUERY)
  } # END if proportions

  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 5.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
