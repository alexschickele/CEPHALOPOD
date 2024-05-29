#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_check <- function(FOLDER_NAME = NULL,
                        SUBFOLDER_NAME = NULL){

  # --- 1. Initialize function
  if(!exists("SEED")){
    SEED <- 123
  }
  set.seed(SEED)
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_check ********************"))
  # --- 1.2. Load the run metadata and query
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))

  # --- 1.3. Moving average function
  # Short and only used here, thats why it is not in the function folder
  ma <- function(x, n = 10){
                 if(length(x) > n){ma = stats::filter(x, rep(1 / n, n), sides = 2)
                 } else {ma = NA}
                 return(ma)
    } # end function

  # --- 1.4. Target transformation - if continuous
  # Same transformation as at the model stage; used for outliers, spearman, mutual information and target vs. feature plot
  # As the RFE is a random forest, it works on quantiles, so we do not care about target transformation there.
  if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
    source(CALL$TARGET_TRANSFORMATION)
    Y <- target_transformation(x = QUERY$Y$measurementvalue, REVERSE = FALSE)$out %>% as.data.frame()
  } else {
    Y <- QUERY$Y %>% as.data.frame()
  } # end if

  # --- 2. Target outlier analysis
  # Outlier check on the query based on z-score (from Nielja code)
  if(CALL$OUTLIER == TRUE){
    if(CALL$DATA_TYPE == "binary"){
      message("--- Cannot perform outlier analysis on presence - pseudo absence data ---")
    } else {
      to_remove <- outlier_iqr_col(Y, n = 2.5) %>% as.vector()
      if(length(to_remove > 0)){
        Y <- dplyr::slice(Y, -to_remove) # we slice the transformed dataset as well for the plots
        QUERY$Y <- dplyr::slice(QUERY$Y, -to_remove)
        QUERY$X <- dplyr::slice(QUERY$X, -to_remove)
        QUERY$S <- dplyr::slice(QUERY$S, -to_remove)
      } # if to remove !NULL

      message(paste("--- OUTLIERS : Removed row number", to_remove, "\n"))
    }
  } # END if outlier TRUE

  # --- 2. Univariate feature check
  # Remove features that are not explaining the target better than a random one
  # We use both Mutual Information and Spearman as they are slightly different processes

  # --- 2.1. Prepare data for the plot // by feature x metric x random
  # Loops over columns if data are proportions
  univ_data <- lapply(1:ncol(QUERY$X), function(x){
    x <- lapply(1:ncol(Y), function(y){
      y <- lapply(1:100, function(b){
        b = feature_selection(V1 = Y[,y], V2 = QUERY$X[,x], METHOD = "emp")
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% abind(along = 3) %>% aperm(c(3,2,1))
  dimnames(univ_data)[[2]] <- c("mutual_information", "spearman")

  # Security for NA
  univ_data[is.na(univ_data)] <- 0

  # --- 2.2. Compute the average and SD per bins // for the error bars
  univ_data_med <- apply(univ_data, c(1,2), median)
  univ_data_sd <- apply(univ_data, c(1,2), sd)

  # --- 2.3. Nice output table
  univ_feature_check <- univ_data_med %>%
    as.data.frame() %>%
    mutate(ID = row_number(),
           varname = colnames(QUERY$X))

  # --- 2.4. List features to keep
  univ_to_keep <- univ_feature_check %>%
    dplyr::filter(mutual_information > 0 | spearman > 0)
  univ_to_remove <- colnames(QUERY$X)[-univ_to_keep$ID] # useful for the message to user

  # --- 2.5. Update QUERY with feature names to keep
  # --- 2.5.1. Inform user and log file
  if(length(univ_to_keep) == 0){
    message(paste("FEATURE UNIVARIATE CHECK:", SUBFOLDER_NAME, "no features explain the target distribution better than a random one.
            Please check the chosen feature and your target distribution."))
    return(NA) # early return if no features (avoid crash later)
  } else if (length(univ_to_remove) > 1){
    message(paste("FEATURE UNIVARIATE CHECK: discard", univ_to_remove,
                  "- no more variance explained than a NULL \n"))
  }

  # --- 2.5.2. Do the update
  QUERY$SUBFOLDER_INFO$ENV_VAR <- colnames(QUERY$X)[univ_to_keep$ID]

  # --- 3. Environmental variable correlation check
  # Removing correlated environmental variables to avoid correlated model features
  if(is.numeric(CALL$ENV_COR) == TRUE){
    # --- 3.1. Opening environmental value at presence points
    features <- QUERY$X[,QUERY$SUBFOLDER_INFO$ENV_VAR] # accounting for univariate check
    if(CALL$DATA_TYPE == "binary"){
      features <- features[which(QUERY$Y$measurementvalue != 0),]
    }

    # --- 3.2. Filter out unique features
    features_unique <- apply(features, 2, function(x)(x = length(unique(x))))
    features <- features[, which(features_unique > 1)]

    # --- 3.2. Check correlation/distance between variables
    features_dist <- cor(features, method = "pearson")
    features_dist <- as.dist(1-abs(features_dist))

    # --- 3.3. Do a clustering and cut
    features_clust <- hclust(features_dist) %>% as.dendrogram()
    features_group <- cutree(features_clust, h = 1-CALL$ENV_COR)

    # --- 3.4. Choose one variable within each inter-correlated clusters
    # Keeping the one with the highest Spearman and Mutual Information rel. to NULL
    cor_to_keep <- NULL
    for(i in 1:max(features_group)){
      tmp <- univ_feature_check %>%
        dplyr::filter(varname %in% names(which(features_group == i))) %>%
        mutate(var_expl = max(c(mutual_information, spearman))) %>%
        dplyr::arrange(var_expl) %>%
        dplyr::select(varname) %>% .[,1] %>% rev()
      message(paste("--- ENV_COR : Cluster", i, ": Keeping", tmp[1], "\n"))
      cor_to_keep <- c(cor_to_keep, tmp[1])
      if(length(tmp) > 1){
        out <- tmp[-1]
        message(paste("--- ENV_COR : Cluster", i, ": Removing", out, "\n"))
      }
    }

    # --- 3.5. Plot the corresponding dendrogram
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/02_env_cor.pdf"))
    par(mfrow = c(1,1), mar = c(10,2,4,15), cex = 0.6)
    pal <- rep("darkorange", length(features))
    pal[get_leaves_attr(features_clust, "label") %in% cor_to_keep] <- "#1F867B"
    labels_colors(features_clust) = pal
    plot(features_clust, horiz = TRUE, axes = FALSE,
         main = "ENVIRONMENTAL PREDICTORS \n Pearson's correlation (r) at the sampling stations", cex = 0.5,
         xlab = "Pearson correlation")
    mtext(side = 1, line = 7, "    Dendrogram of feature correlation. The orange line represents the threshold at which two features
    are considered as correlated. The one with the higher Mutual Information and Spearman correlation
    with the target is kept (green). The others are discarded for the rest of the analysis (orange)", cex = 0.6, adj = 0)
    axis(side = 1, at = seq(0,1,0.2), labels = seq(1,0,-0.2), las = 1, cex.axis = 1)
    abline(v = seq(0,1,0.2), col = "gray50", lty = "dashed")
    abline(v = 1-CALL$ENV_COR, col = "darkorange")
    dev.off()

    # --- 3.6. Update ENV_VAR
    QUERY$SUBFOLDER_INFO$ENV_VAR <- cor_to_keep
  } # END if env_cor TRUE

  # # --- 4. RFE (recursive feature elimination) importance analysis
  # Done with a Random forest using the method developed in the "Caret" library
  if(CALL$RFE == TRUE){
    if(CALL$DATA_TYPE == "proportions"){
      message("A RFE predictor selection is not possible for proportion data,
              please select carefully your predictors")
    } else {
      # --- 4.1. Initialize data and control parameters
      features <- QUERY$X[,QUERY$SUBFOLDER_INFO$ENV_VAR]
      rfe_df <- cbind(QUERY$Y, features)
      # Modify rfFuncs : we modify the selectSize code to use pickSizeBest instead of a 1.5 tolerance around the best
      rfFuncs$selectSize <- function(...) pickSizeBest(...) # weird functioning :-)
      rfe_control <- rfeControl(functions=rfFuncs, method="cv", number=5, rerank = FALSE)
      # --- 4.2. Run the RFE algorithm
      message(paste(Sys.time(), "--- RFE : Fitting the Recursive Feature Exclusion algorithm \n"))
      rfe_fit <- rfe(rfe_df[,-1], rfe_df[,1], sizes = c(1:ncol(rfe_df[,-1])), rfeControl = rfe_control)

      # --- 4.3. Extract the relevant predictors
      # --- 4.3.1. Compute the moving average loss
      loss_ma <- ma(rfe_fit$results$RMSE, n = 10)
      # --- 4.3.2. Compute the percentage loss
      loss_ma_pct <- (loss_ma[-1] - loss_ma[-length(loss_ma)])/loss_ma[-length(loss_ma)]*100
      # --- 4.3.3. Find the first minimum or <1% loss percentage
      # Consider all variables if "id" is NA, due to moving average not working for low variable number
      id <- which(loss_ma_pct > -1)[1]
      if(is.na(id) == TRUE){
        message("--- The moving window could not find a minimum loss, please considering
                removing this option from your run due to insufficient number of environmental variables")
        id <- ncol(features)
      }

      # --- 4.4. Compute variable importance
      rfe_vip <- rfe_fit$variables %>%
        dplyr::select(var, Overall) %>%
        mutate(var = fct_reorder(as.factor(var), Overall, .desc = TRUE))

      # --- 4.5. Extract the environmental variables
      rfe_to_keep <- rfe_vip$var[1:id] %>% as.character()
      message(paste("--- RFE : Selecting", rfe_to_keep, "\n"))

      # --- 4.6. Produce an information plot
      pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/03_feature_pre_selection.pdf"))
      par(mfrow = c(2,1), mar = c(4,3,2,15))
      # --- Variable importance
      boxplot(rfe_vip$Overall ~ rfe_vip$var,
              main = paste("ENVIRONMENTAL PREDICTORS \n A-priori importance for ID:", SUBFOLDER_NAME),
              xlab = "", ylab = "", axes = FALSE, outline = FALSE, horizontal = TRUE,
              col = c(rep("#1F867B", id), rep("gold", ncol(features)-id)), cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7)
      axis(side = 4, at = 1:ncol(features), labels = levels(rfe_vip$var), las = 2, cex.axis = 0.6)
      axis(side = 1, at = c(seq(0, 15, 5), seq(0, 100, 20)), labels = c(seq(0, 15, 5), seq(0, 100, 20)), cex.axis = 0.7)
      title(xlab = "Estimated importance (%)", line = 2, cex.lab = 0.7)
      abline(v = c(seq(0, 15, 5), seq(0, 100, 20)), lty = "longdash", col = "gray50")
      box()
      box("figure", col="black", lwd = 1)
      # --- Number of variables
      par(mar = c(6,3,2,15))
      plot(rfe_fit$results$RMSE, rfe_fit$results$Variables, pch = 21, cex = 2, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7,
           bg = c(rep("#1F867B", id), rep("gold", ncol(features)-id)), col = "black",
           main = paste("ENVIRONMENTAL PREDICTORS \n Optimal number for ID:", SUBFOLDER_NAME),
           ylab = "", xlab = "")
      axis(side = 4, at = 1:ncol(features), labels = levels(rfe_vip$var), las = 2, cex.axis = 0.6)
      title(xlab = "Loss metric (RMSE)", ylab = "Nb. of considered predictors", line = 2, cex.lab = 0.7)
      grid(col = "gray50")
      mtext(side = 1, line = 5, "Feature selection by recursive feature exclusion procedure (Random Forest algorithm). The upper panel presents the ranked feature
importance (%) in explaining the target. The lower panel represents the loss (RMSE) in function of the number of features considered.
The optimal number of features considered conresponds to the number after which adding a feature does not decrease the loss (moving
average of 5) by 1 %. The considered features are in green. The ones that do not bring additional information are discarded and
displayed in yellow.", cex = 0.6, adj = 0)
      dev.off()

      # --- 4.7. Update ENV_VAR
      QUERY$SUBFOLDER_INFO$ENV_VAR <- rfe_to_keep

    }
  } # END if RFE TRUE

  # --- 5. MESS analysis
  r_mess <- NULL
  for(m in 1:length(CALL$ENV_DATA)){
    # --- 5.1. Load necessary data
    features <- CALL$ENV_DATA[[m]] %>%
      raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR)
    names(features) <- QUERY$SUBFOLDER_INFO$ENV_VAR

    # --- 5.2. Compute the mess analysis
    # --- 5.2.1. Load environmental samples data
    tmp <- QUERY$X %>% dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))

    # --- 5.2.2. Remove the pseudo-absences from the analysis
    # We should not consider them as true input data as they are user defined
    if(CALL$DATA_TYPE == "binary"){
      tmp <- tmp[which(QUERY$Y != 0),]
    } # if pres remove pseudo abs

    # --- 5.2.3. Analysis
    r_mess[[m]] <- dismo::mess(x = features, v = as.data.frame(tmp), full = FALSE)

  } # for m month

  # --- 5.3. Append to query
  QUERY$MESS <- stack(r_mess)

  # --- 6. Verification of feature pre-selection
  # --- 6.1. Initialize QC and necessary data
  post_check <- univ_feature_check %>%
    dplyr::filter(varname %in% !!QUERY$SUBFOLDER_INFO$ENV_VAR)

  # --- 6.2. Compute the quality check
  MI_diff <- mean(post_check$mutual_information)
  S_diff <- mean(post_check$spearman)
  PRE_VIP <- max(c(MI_diff, S_diff))

  # --- 6.3. Append QUERY
  QUERY[["eval"]][["PRE_VIP"]] <- PRE_VIP

  # --- 7. Verification plot
  # --- 7.1. Initialize plot
  pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/04_feature_verification.pdf"))

  # --- 7.2. Initialize colors
  pal <- rep("#B64A60", ncol(QUERY$X))
  pal[which(names(QUERY$X) %in% univ_to_keep$varname)] <- "darkorange" # passed univariate selection
  pal[which(names(QUERY$X) %in% cor_to_keep)] <- "gold" # passed univariate selection

  if(CALL$RFE == TRUE){
    pal[which(names(QUERY$X) %in% QUERY$SUBFOLDER_INFO$ENV_VAR)] <- "#1F867B" # final selection OR after RFE
  } # end if

  # --- 7.3. Mutual information and Spearman plot
  # --- 7.3.1. Initialize plot
  par(mfrow = c(1,1), mar = c(10,6,4,4))

  # --- 7.3.2. Plot points
  plot(univ_data_med[,1], univ_data_med[,2],
       xlim = c(min(univ_data_med[,1])-0.1, max(univ_data_med[,1]))+0.1,
       ylim = c(min(univ_data_med[,2])-0.1, max(univ_data_med[,2]))+0.1,
       xlab = "Mutual Information rel. to NULL feature", ylab = "Spearman correlation rel. to NULL feature",
       main = "Target variance explained by each feature")
  mtext(side = 1, line = 7, "Each dot represents feature (numbered) that was not containing more information than random (red), correlated with another feature
(orange), not bringing additional information relative to other selected features (yellow), considered in the final feature set (green). The
lines associated to each feature represent a 1 standard deviation error bar. The red shading represent the threshold at which a feature
explains less information than random. The yellow shading represent the threshold over which, the average variance explained in the
feature set is considered as sufficient", cex = 0.5, adj = 0)
  polygon(x = c(-10,10,10,-10), y = c(0.05,0.05,-10,-10), col = scales::alpha("orange", 0.5), border = NA, density = 20)
  polygon(x = c(0.05,0.05,-10,-10), y = c(-10,10,10,-10), col = scales::alpha("orange", 0.5), border = NA, density = 20)
  polygon(x = c(-10,10,10,-10), y = c(0,0,-10,-10), col = scales::alpha("darkred", 0.5), border = NA, density = 20)
  polygon(x = c(0,0,-10,-10), y = c(-10,10,10,-10), col = scales::alpha("darkred", 0.5), border = NA, density = 20)

  grid(col = "black")
  abline(h = c(0, 0.05), v = c(0, 0.05) , lwd = 1, col = c("darkred", "orange"))

  # --- 7.3.3. Plot error bars
  for(i in 1:ncol(QUERY$X)){
    lines(x = c(univ_data_med[i,1]-univ_data_sd[i,1], univ_data_med[i,1]+univ_data_sd[i,1]),
          y = rep(univ_data_med[i,2], 2),
          lwd = 2)
    lines(x = rep(univ_data_med[i,1], 2),
          y = c(univ_data_med[i,2]-univ_data_sd[i,2], univ_data_med[i,2]+univ_data_sd[i,2]),
          lwd = 2)
  } # end for

  # --- 7.3.4. Add labels
  points(univ_data_med[,1], univ_data_med[,2], pch = 20, cex = 6, col = pal) # background
  points(univ_data_med[,1], univ_data_med[,2], cex = 4, lwd = 2) # circle
  text(univ_data_med[,1], univ_data_med[,2], labels = seq(1:ncol(QUERY$X)), offset = 0, col = "black", cex = 1, lwd = 10)

  # --- 7.4. Univariate target vs. feature plot
  # --- 7.4.1. initialize plot
  par(mfrow = c(4,3), mar = c(6,4,3,4))

  # --- 7.4.2. Loop over the features (and target if proportions)
  for(x in 1:ncol(QUERY$X)){
    for(y in 1:ncol(Y)){
      # --- 7.4.2.1. Prepare the data and bin the feature
      feature_bin <- discretize(QUERY$X[,x], nbins = 10, disc = "equalwidth")
      bin_mid <- seq(min(QUERY$X[,x]),max(QUERY$X[,x]), length.out = 10) %>% signif(2)

      if(ncol(QUERY$Y) > 1){yname <- paste("Target:", colnames(QUERY$Y)[y])} else{yname <- paste("Target:", SUBFOLDER_NAME)}

      # --- 7.4.2.2. Do the plot
      if(CALL$DATA_TYPE != "binary"){
        boxplot(Y[,y]~feature_bin[,1], axes = FALSE, col = pal[x],
                main  = paste("(", x, ")", names(QUERY$X[x])), cex.main = 0.8,
                sub = paste("Mutual Information:", signif(univ_feature_check$mutual_information[x],2), "| Spearman cor.:", signif(univ_feature_check$spearman[x],2)), cex.sub = 0.8,
                xlab = "Feature values", cex.lab = 0.8,
                ylab = yname)
        box()
        axis(side = 1, at = 1:10, labels = bin_mid, cex.axis = 0.7, las = 2)
        axis(side = 2, las = 2, cex.axis = 0.7)
        box("figure", col="black", lwd = 1)
      } else {
        if(length(unique(QUERY$X[,x][which(Y[,y] == 1)])) > 1){
          message("QUERY CHECK (information): features with unique values are not plotted")
          hist(QUERY$X[,x][which(Y[,y] == 1)], breaks = seq(min(QUERY$X[,x]), max(QUERY$X[,x]), length.out = 25), col = scales::alpha(pal[x], 0.3),
               main  = paste("(", x, ")", names(QUERY$X[x])), cex.main = 0.8, axes = FALSE,
               sub = paste("Mutual Information:", signif(univ_feature_check$mutual_information[x],2), "| Spearman cor.:", signif(univ_feature_check$spearman[x],2)), cex.sub = 0.8,
               xlab = "Feature values", cex.lab = 0.8,
               ylab = paste(yname, "(Frequency) \n Pseudo-abs. (gray) ; Presence (colored)"))
          hist(QUERY$X[,x][which(Y[,y] == 0)], breaks = seq(min(QUERY$X[,x]), max(QUERY$X[,x]), length.out = 25), col = scales::alpha("black", 0.3), add = TRUE)

          box()
          axis(side = 1, las = 2, cex.axis = 0.7)
          axis(side = 2, las = 2, cex.axis = 0.7)
          box("figure", col="black", lwd = 1)
          grid(col = "gray20")
        } # if unique value
      } # if binary data
    } # for x features
  } # for y target

  # --- 7.5. Wrap up plot
  dev.off()

  # --- 8. Wrap up and save
  # --- 8.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 8.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 8.3. Pretty return
  # Updates the list of species to model depending on the PRE_VIP (> 0.05) to avoid non meaningful feature ensembles
  if(QUERY$eval$PRE_VIP >= 0.05 | CALL$FAST == FALSE){
    return(SUBFOLDER_NAME)
  } else {
    message(paste("QUERY CHECK:", SUBFOLDER_NAME, "The selected features do not present significant trends to the observations; please work on the predictors and data. \n
            The species is discarded from further analysis"))
    return(NA)
  }

} # END FUNCTION
