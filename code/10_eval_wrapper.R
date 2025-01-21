#' =============================================================================
#' @name eval_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return outputs are saved in the MODEL.RData object

eval_wrapper <- function(FOLDER_NAME = NULL,
                         SUBFOLDER_NAME = NULL){

  # --- 1. Initialize function and environment
  set.seed(123)  # Ensure reproducibility
  
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : eval_wrapper ********************"))

  # --- 1.2. Load essential metadata
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))  # Load metadata for the current run
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))  # Load query-specific metadata
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData")) # Load model-specific metadata

  # --- 1.3. Setting up a path for the PDF file
  pdf_path <- paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/04_variable_importance.pdf")

  # --- 2. Redirection to PRESENCE model evaluation
  if(CALL$DATA_TYPE == "presence_only"){

    # --- 2.1. Run function
    MODEL<- eval_presence_only(CALL = CALL, QUERY = QUERY, MODEL = MODEL, PDF_PATH = pdf_path)
  } # END if presence_only

  # --- 3. Redirection to CONTINUOUS model evaluation
  if(CALL$DATA_TYPE == "continuous"){

    # --- 3.1. Run function
    MODEL <- eval_continuous(CALL = CALL, QUERY = QUERY, MODEL = MODEL, PDF_PATH = pdf_path)
  } # END if continuous

  # --- 4. Redirection to PROPORTIONS model evaluation
  if(CALL$DATA_TYPE == "proportions"){

    # --- 4.1. Run function
    MODEL <- eval_proportions(CALL = CALL, QUERY = QUERY, MODEL = MODEL, PDF_PATH = pdf_path)
  } # END if proportions

  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 5.3. Pretty return
  return(SUBFOLDER_NAME)

} # END FUNCTION
