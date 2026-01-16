#' =============================================================================
#' @name hyperparameter
#' @description this function creates an hyperparameter grid for each available
#' algorithm in the pipeline and returns those selected by the user
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @return a hyperparameter list of the chosen models, including spec and grid
#' per model
#' @return the general HP list is saved in the CALL.RData object

hyperparameter <- function(FOLDER_NAME = NULL){

  # --- 1. Initialize the  object
  set.seed(123)
  
  # --- 1.1. Output object
  HP <- list()
  # --- 1.2. Load CALL object
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 2. Fill up the object with a set of possible grids
  # --- 2.3. GENERALIZED LINEAR MODELS
  # --- 2.3.1. Define the model specifications
  if(CALL$DATA_TYPE == "presence_only"){
    HP$GLM$model_spec <- linear_reg(mode = "regression") %>% 
      set_engine("glm", 
                 family = stats::binomial(link = "logit")) %>% 
      step_normalize() %>% 
      translate()
  } else {
    HP$GLM$model_spec <- linear_reg(mode = "regression",
                                    engine = "glm") %>% 
      step_normalize() %>% 
      translate()
  }
  
  # --- 2.3.2. Define the grid according to built in functions
  # No parameters to tune for engine "glm"
  
  # --- 2.2. GENERALIZED ADDITIVE MODELS 
  # --- 2.2.1. Define the model specifications
  if(CALL$DATA_TYPE == "presence_only"){
    HP$GAM$model_spec <- gen_additive_mod(mode = "regression",
                                          adjust_deg_free = tune(),
                                          select_features = FALSE) %>% 
      set_engine("mgcv",
                 family = stats::binomial(link = "logit")) %>% 
      step_normalize() %>% 
      translate()
  } else {
    HP$GAM$model_spec <- gen_additive_mod(mode = "regression",
                                          engine = "mgcv",
                                          adjust_deg_free = tune(),
                                          select_features = FALSE) %>% 
      step_normalize() %>% 
      translate()
  }
  
  # --- 2.2.2. Define the grid according to built in functions
  HP$GAM$model_grid <- grid_regular(adjust_deg_free(range = c(0.5, 2.5)),
                                    levels = CALL$LEVELS)
  
  
  # --- 2.1. RANDOM FOREST 
  # --- 2.1.1. Define the model specifications
  HP$RF$model_spec <- rand_forest(mode = "regression",
                                  engine = "randomForest",
                                  trees = tune(),
                                  min_n = tune())
  
  # --- 2.1.2. Define the grid according to built in functions
  HP$RF$model_grid <- grid_regular(trees(range = c(200, 1000)),
                                   min_n(range = c(2, 20)),
                                   levels = CALL$LEVELS)
  
  # --- 2.4. SINGLE LAYER NEURAL NETWORK
  # --- 2.4.1. Define the model specifications
  HP$MLP$model_spec <- mlp(mode = "regression",
                           engine = "nnet",
                           hidden_units = tune(),
                           penalty = 0.01,
                           dropout = 0,
                           epochs = tune(),
                           activation = NULL,
                           learn_rate = 1e-2) %>% 
    step_normalize() %>% 
    translate()
  
  
  # --- 2.4.2. Define the grid according to built in functions
  HP$MLP$model_grid <- grid_regular(hidden_units(range = c(1, 5)),
                                    epochs(range = c(10, 50)),
                                    levels = CALL$LEVELS)
  
  # --- 2.5. BOOSTED REGRESSION TREES
  # --- 2.5.1. Define the model specifications
  HP$BRT$model_spec <- boost_tree(mode = "regression",
                                  engine = "xgboost",
                                  min_n = tune(),
                                  tree_depth = tune(),
                                  learn_rate = tune(),
                                  stop_iter = 10)
  
  # --- 2.5.2. Define the grid according to built in functions
  HP$BRT$model_grid <- grid_regular(min_n(range = c(2, 20)),
                                    tree_depth(range = c(3, 10)),
                                    learn_rate(range = c(1e-1, 1e-3)),
                                    levels = CALL$LEVELS)
  
  # --- 2.6. SUPPORT VECTOR MACHINE
  # --- 2.6.1. Define the model specifications
  HP$SVM$model_spec <- svm_rbf(mode = "regression",
                               engine = "kernlab",
                               cost = tune(),
                               rbf_sigma = tune(),
                               margin = tune()) %>% 
    step_normalize() %>% 
    translate()

  # --- 2.6.2. Define the grid according to built in functions
  HP$SVM$model_grid <- grid_regular(cost(range = c(-3, 3)),       # 10^-3 to 10^3
                                    rbf_sigma(range = c(-6, -1)), # 10^-6 to 10^-1
                                    svm_margin(range = c(-3, -1)),
                                    levels = CALL$LEVELS)
  
  # --- 2.7. MULTIVARIATE BOOSTED TREE REGRESSOR
  # Specific to proportions data
  HP$MBTR$model_grid <- data.frame(LEARNING_RATE = seq(1e-2, 1e-3, length.out = CALL$LEVELS),
                                   N_Q = seq(5,25, length.out = CALL$LEVELS),
                                   MEAN_LEAF = seq(2, 100, length.out = CALL$LEVELS)) %>% 
    expand.grid() %>% 
    unique()
  
  # --- 3. Hyper parameter selection
  # --- 3.1. For univariate data : According to the specified model list
  if(CALL$DATA_TYPE != "proportions"){
    HP <- HP[CALL$MODEL_LIST]
    HP$MODEL_LIST <- CALL$MODEL_LIST
  } else {
    HP <- HP["MBTR"]
    HP$MODEL_LIST <- "MBTR"
  }
  
  # --- 4. Wrap up and save
  # --- 4.1. Append CALL and save
  CALL[["HP"]] <- HP
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
} # END FUNCTION
