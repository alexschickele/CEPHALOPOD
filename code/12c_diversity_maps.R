#' =============================================================================
#' @name diversity_maps
#' @description Function extracting all projections by species and bootstrap.
#' Then computes a set of diversity metrics and plots
#' @param FOLDER_NAME name of the corresponding folder
#' @param MONTH a list of month to concatenate together
#' @return plots mean and uncertainty maps per diversity metric
#' @return a diversity and mess object per projection, species and bootstrap.
#' @return a .nc matching the EMODnet standards with the same informations as 
#' above. - Saved in FOLDERNAME.

diversity_maps <- function(FOLDER_NAME = NULL,
                           MONTH = list(c(10,11,12,1,2,3),
                                        4:9)){

  # --- 1. Initialize function
  # --- 1.1. Global parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))

  # --- 1.2. Baseline raster
  r0 <- CALL$ENV_DATA[[1]][[1]]

  # --- 1.3. Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA

  # --- 2. Build the ensemble(s)
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - START"))
  
  # --- 2.1. Extract file information
  # Which subfolder list and which model in the ensemble
  all_files <- list.files(paste0(project_wd, "/output/", FOLDER_NAME), recursive = TRUE)
  model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]
  ensemble_files <- mclapply(model_files, function(x){
    memory_cleanup() # low memory use
    
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/MODEL.RData"))
    if(length(MODEL$MODEL_LIST) >= 1){
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/QUERY.RData"))
      return(list(SUBFOLDER_NAME = x, MODEL_LIST = MODEL$MODEL_LIST, Y = QUERY$Y, MESS = QUERY$MESS, REC = MODEL$recommandations))
    } else {
      return(NULL)
    } # if model list
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.preschedule = FALSE) %>% 
    .[lengths(.) != 0] %>% 
    .[grep("Error", ., invert = TRUE)] # to exclude any API error or else
  
  # --- 2.2. Loop over the files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - loop over files"))
  tmp <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    
    # --- 2.2.1. Load MODEL files
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
    
    # --- 2.2.2. Extract projections in a matrix
    # If there is more than 1 algorithm, we extract and average across algorithm
    # Output matrix is cell x bootstrap x month
    if(length(x$MODEL_LIST) > 1){
      m <- lapply(x$MODEL_LIST, function(y){
        MODEL[[y]][["proj"]]$y_hat
      }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
    } else {
      m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
    } # end if
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.cleanup = TRUE)
  
  # --- 2.3. Stack in a cell x species x bootstrap x month matrix
  # --- 2.3.1. Re-arrange the array
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - format to array"))
  all_ens <- tmp %>% 
    abind(along = 4) %>% 
    aperm(c(1,4,2,3))
  # --- 2.3.2. Pretty dimensions
  dimnames(all_ens)[[2]] <- lapply(ensemble_files, function(x){out <- x$SUBFOLDER_NAME}) %>% unlist() %>% as.character()
  dimnames(all_ens)[[4]] <- 1:12 %>% as.character()
  
  # --- 2.4. Memory cleanup
  rm(tmp)
  gc() # clean garbage and temporary files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - DONE"))
  
  # --- 2.5. Return if not enough ensembles
  if(dim(all_ens)[[2]] < 5){
    message("--- DIVERSITY : not enough ensembles to compute diversity from")
    return(NA)
  }
  
  # --- 2.6. Convert recommendations to list
  # --- 2.6.1. Extract recommandations as a list
  tmp <- lapply(ensemble_files, function(x){
    out <- setNames(as.list(apply(x$REC, 1, function(row) paste(row, collapse = ", "))), rownames(x$REC))
    return(out)
  }) # end lapply
  
  # --- 2.6.2. Flatten for each ensemble
  if(CALL$DATA_TYPE == "proportions"){
    flat_list <- paste(names(tmp), ": ", tmp[["MBTR"]])
  } else {
    flat_list <- lapply(tmp, function(df) {
      lapply(names(df), function(algo) {
        paste(algo, ": ", df[[algo]], sep = "")
      }) # end lapply
    }) # end lapply
  } # end if
  
  # --- 2.6.3. Flatten again
  flat_list <- lapply(flat_list, function(x) paste(unlist(x), collapse = "\n"))
  recommandation_list <- paste(flat_list, collapse = "; ")
  
  # --- 3. Compute habitat suitability - for the .nc in the shiny app.
  # --- 3.1. Across bootstrap for each month
  a_m <- mclapply(1:dim(all_ens)[1], function(i) {
    apply(all_ens[i,,,], c(1,3), function(x) mean(x, na.rm = TRUE))
  }, mc.cores = MAX_CLUSTERS) %>% abind(along = 3) %>% aperm(c(3,1,2))
  a_sd <- mclapply(1:dim(all_ens)[1], function(i) {
    apply(all_ens[i,,,], c(1,3), function(x) sd(x, na.rm = TRUE))
  }, mc.cores = MAX_CLUSTERS) %>% abind(along = 3) %>% aperm(c(3,1,2))
  
  # --- 3.2. Add a 13th month: mean across the year
  mean_across_12_a_m_expanded <- array(apply(a_m, c(1, 2), function(x) mean(x, na.rm = TRUE)),
                                       dim = c(dim(a_m)[1], dim(a_m)[2], 1))
  sd_across_12_a_sd_expanded <- array(apply(a_sd, c(1,2), function(x) sd(x, na.rm = TRUE)),
                                      dim = c(dim(a_m)[1],  dim(a_m)[2],1))
  
  a_m <- abind::abind(a_m, mean_across_12_a_m_expanded, along = 3) %>% aperm(., c(1, 3, 2))
  a_sd <- abind::abind(a_sd, sd_across_12_a_sd_expanded, along = 3) %>% aperm(., c(1, 3, 2))
  
  # --- 4. Construct diversity
  # --- 4.1. Compute alpha diversity indices
  a_shannon <- apply(all_ens, c(1,3,4), function(x)(x = vegan::diversity(x, "shannon")))
  a_richness <- apply(all_ens, c(1,3,4), function(x)(x = sum(x, na.rm = TRUE)))
  a_evenness <- a_shannon/(a_richness)
  a_invsimpson <- apply(all_ens, c(1,3,4), function(x)(x = vegan::diversity(x, "invsimpson")))

  # --- 4.2. Stack diversities
  div_all <- abind(a_shannon, a_richness, a_evenness, a_invsimpson, along = 4)
  div_m <- apply(div_all, c(1,3,4), function(x)(x = mean(x, na.rm = TRUE)))
  div_sd <- apply(div_all, c(1,3,4), function(x)(x = sd(x, na.rm = TRUE)))
  div_names <- c("a_shannon", "a_richness", "a_evenness", "a_invsimpson")

  # --- 4.3. Add a 13th month: mean across the year
  mean_across_12_div_m_expanded <- array(apply(div_m, c(1, 3), function(x) mean(x, na.rm = TRUE)),
                                         dim = c(dim(div_m)[1], 1, dim(div_m)[3]))
  sd_across_12_div_sd_expanded <- array(apply(div_sd, c(1,3), function(x) sd(x, na.rm = TRUE)),
                                        dim = c(dim(div_m)[1], 1, dim(div_m)[3]))
  
  div_m <- abind::abind(div_m, mean_across_12_div_m_expanded, along = 2)
  div_sd <- abind::abind(div_sd, sd_across_12_div_sd_expanded, along = 2)
  
  # --- 5. Extract MESS analysis
  message(paste(Sys.time(), "--- Extract MESS"))
  mess_all <- mclapply(model_files, FUN = function(s){
    # --- 6.1. Load corresponding objects
    # We have a file exist security here in case the diversity is run later
    model_file <- paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData")
    load(model_file)

    if(length(MODEL$MODEL_LIST) != 0){
      query_file <- gsub('MODEL.RData', 'QUERY.RData', model_file)
      load(query_file)
      mess_s <- getValues(QUERY$MESS*-1)
      }
  },
  mc.cores = min(length(model_files), MAX_CLUSTERS), mc.preschedule = FALSE) %>%
  abind(along = 3) %>%
  apply(c(1,2), mean)

  # --- 5.3. Final adjustments
  mess_all[mess_all > 100] <- 100

  # --- 5.4. Intermediate save
  save(div_all, file = paste0(project_wd, "/output/", FOLDER_NAME,"/DIVERSITY.RData"))

  # --- 6. Diversity plots
  # --- 6.1. Initialize pdf and plot
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/diversity_maps.pdf"))
  par(mfrow = c(3,2), mar = c(2,3,3,3))

  # --- 6.2. Plot the legends
  # --- 6.2.1. Diversity legend
  plot.new()
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.1,
                col = inferno_pal(100), border = "black")
  axis(side = 1)
  text(x = 0.5, y = 0.3, "Diversity value [0 - 1]", adj = 0.5)

  # --- 6.2.2. MESS x CV 2D color scale
  par(mar = c(5,5,1,2))
  bivar_pal <- colmat(nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Standard deviation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = round(seq(0, 0.25, length.out = 6), 2))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(5,3,3,1))

  # --- 6.3. Plot diversity maps
  # Loop over diversity maps and projections
  for(i in 1:dim(div_m)[[3]]){
    for(m in MONTH){

      # --- 6.3.1. Prepare the data
      tmp_m <- div_m[,m,] %>% apply(c(1,3), mean) %>% .[,i]
      tmp_sd <- div_sd[,m,] %>% apply(c(1,3), mean) %>% .[,i]
      plot_scale <- quantile(tmp_m, 0.95, na.rm = TRUE)

      # --- 6.3.2. Assign to raster
      # Average raster is cut at plot scale (Q95)
      r_m <- r0 %>% setValues(tmp_m)
      r_m[r_m > plot_scale] <- plot_scale

      r_sd <- r0 %>% setValues(tmp_sd)

      r_mess <- r0 %>% setValues(apply(mess_all[, m], 1, mean))
      r_mess[r_mess<0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included
      r_mess[r_mess>100] <- 100

      # --- 6.3.3. Average projection
      plot(r_m, col = inferno_pal(100), legend=FALSE,
           main = paste("DIVERSITY (",div_names[i], ") \n Month:", paste(m, collapse = ",")),
            cex.main = 1)
      mtext(text = paste("Projection rescaling factor: ", round(plot_scale,2),
                         "\n Number of species ensembles:", dim(all_ens)[[2]]),
            side = 1, line = 3, cex = 0.7)
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)

      # try to create contours
      contours <- tryCatch({
        raster::rasterToContour(r_m, nlevels = 4)
      }, error = function(e) {
          NULL
      }) # end tryCatch
      if (!is.null(contours)){plot(contours, add = TRUE)}

      box("figure", col="black", lwd = 1)

      # --- 6.3.4. Uncertainties projection
      # First we draw a land mask
      plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties", cex.main = 1)
      # Then compute and plot the 2-dimensional color scale for SD x MESS
      r <- bivar_map(rasterx = r_sd, rastery = r_mess, colormatrix = bivar_pal,
                     cutx = seq(0,max(getValues(r_sd), na.rm = TRUE), length.out = 101), cuty = 0:100)
      plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)

      # try to create contours
      contours <- tryCatch({
        raster::rasterToContour(r_sd, nlevels = 4)
      }, error = function(e) {
          NULL
      }) # end tryCatch
      if (!is.null(contours)){plot(contours, add = TRUE)}

      box("figure", col="black", lwd = 1)

    } # m month
  } # i diversity

  # --- 6.4. Stop PDF
  dev.off()

  # --- 7. Save diversity maps
  save(div_all, file = paste0(project_wd, "/output/", FOLDER_NAME,"/DIVERSITY.RData"))

  # --- 8. Generate the corresponding .nc files
  # --- 8.1. Species ensemble
  file = paste0(project_wd, "/output/", FOLDER_NAME,"/","output_ensemble.nc")
  CEPHALOPOD_to_netcdf(r0, a_m, a_sd, dimnames(a_m)[[3]], recommandation_list, file)
  
  # --- 8.2. Diversity level
  file = paste0(project_wd, "/output/", FOLDER_NAME,"/","output_diversity.nc")
  CEPHALOPOD_to_netcdf(r0, div_m, div_sd, div_names, recommandation_list, file)
  
} # END FUNCTION


