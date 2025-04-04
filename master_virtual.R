#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' A. Schickele 2024
#' =============================================================================

# --- 0. Start up and load functions
# All will be called in the config file later
rm(list=ls())
closeAllConnections()
setwd("/net/meso/work/aschickele/CEPHALOPOD")
source(file = "./code/00_config.R")
run_name <- "VIRTUALSPECIES_continuous_v2"
MAX_CLUSTERS <- 100

# --- 0. Generate the virtual species
source(file = "./preparation/generate_virtual.R")
# VIRTUAL <- generate_virtual(ENV_VAR = NULL,
#                             ENV_PATH = "/net/meso/work/clercc/Predictors/PIPELINE_SET/TEST_SET",
#                             MONTH = 6,
#                             NOISE_SD = c(1,2,3),
#                             DATA_TYPE = "proportions")

# --- 1. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "/net/meso/work/aschickele/CEPHALOPOD/data/virtualspecies/virtual_binary.csv",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 20, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 50, START_YEAR = 1950, STOP_YEAR = 2020))

# Define the list of species to consider
sp_list <- list_bio$worms_id %>% unique()
# save(VIRTUAL, file = paste0(project_wd, "/output/", run_name,"/VIRTUAL.RData"))

# --- 2. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the 
# global parameters in an object, (3) Builds a local list of monthly raster
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           WORMS_CHECK = FALSE,
                           FAST = FALSE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "presence_only",
                           ENV_VAR = NULL,
                           ENV_PATH = "/net/meso/work/clercc/Predictors/PIPELINE_SET/TEST_SET",
                           METHOD_PA = "density",
                           PER_RANDOM = 0,
                           PA_ENV_STRATA = TRUE,
                           OUTLIER = TRUE,
                           RFE = TRUE,
                           ENV_COR = 0.8,
                           NFOLD = 5,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","MLP","BRT","GAM","SVM","RF"),
                           LEVELS = 3,
                           TARGET_TRANSFORMATION = NULL,
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0)

# --- 3. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 4. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
subfolder_list <- mcmapply(FUN = query_env,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  .[grep("Error", ., invert = TRUE)] %>% # to exclude any API error or else
  as.vector()

# --- 5. Generate pseudo-absences if necessary
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 6. Outliers, Environmental predictor and MESS check 
# This functions returns an updated subfolder_list with meaningful feature set
subfolder_list <- mcmapply(FUN = query_check,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  as.vector()

# --- 7. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 8. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name)

# --- 9. Model fit
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 10. Model evaluation
# Performance metric and variable importance
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# ---11. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12. Output plots
# --- 12.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12.4 User synthesis
user_synthesis(FOLDER_NAME = run_name)

# --- 0. Evaluate virtual species
load("/net/meso/work/aschickele/CEPHALOPOD/output/00_METHOD_PAPER/runs/VIRTUALSPECIES_binary/VIRTUAL.RData") # load the same as before
save(VIRTUAL, file = paste0(project_wd, "/output/", run_name,"/VIRTUAL.RData")) # save for copy in the current folder

evaluate_virtual(FOLDER_NAME = run_name) # evaluate

save(VIRTUAL, file = paste0(project_wd, "/output/", run_name,"/VIRTUAL.RData")) # save again with the evaluation

# --- END --- 