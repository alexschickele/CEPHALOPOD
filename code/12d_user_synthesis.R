#' =============================================================================
#' @name user_synthesis
#' @description This function merges all PDF produced during the analysis into
#' a synthetic document summarizing the input, quality checks and results.
#' @param FOLDER_NAME name of the corresponding folder
#' @return a PDF document in the FOLDER_NAME directory

user_synthesis <- function(FOLDER_NAME){


  # --- 1. Initialize function
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  #library(pdftools) # déjà installée comme dépendence d'un autre truc

  # --- 2. General information PDF
  # --- 2.1. Title and initial description
  pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/00_general_informations.pdf"))
  par(mar = c(2,1,2,1))
  plot.new()
  text(x = 0.5, y = 1, cex = 2, font = 2, "Detailed summary for user")
  text(x = 0.5, y = 0.95, cex = 0.7, font = 3, paste("For the following run:", FOLDER_NAME))
  text(x = 0, y = 0.75, cex = 0.7, font = 2, adj = c(0,1), "Document description:")
  text(x = 0, y = 0.7, cex = 0.7, font = 1, adj = c(0,1),
       "The present .pdf file provides detailed information related to the analysis in- and outputs. The input description
displayed hereafter correspond to the user-defined parameters related to the biological sample selection
criteria, the type of data used in the analysis and the choices for algorithm training. They are common to all
considered species. The outputs are displayed at the species level, in form of graphical outputs for each pipeline
steps, and only including those passing the embedded quality checks. The graphical outputs provide information
concerning:
     - The biological and environmental data location, selection, and importance in the observations.
     - The environmental variable importance in the model training.
     - The spatial projections and uncertainties
     - The biological response to each environmental variable in form of partial dependency plots
Depending on the user selection, a selection of diversity projections across species can be displayed at the end
of this document. It only considers species passing the quality checks.

Please note that the Bluecloud2026 Ecosystem Workbench development team did its best to provide a habitat
modelling pipeline, including the latest advances and quality checks in this field of research, and quality
assessment guidelines to a large audience. The interpretation of the results and final quality assessment is up to
the user, however.")

  # --- 2.2. Input description
  plot.new()
  text(x = 0, y = 1, cex = 0.7, font = 2, adj = c(0,1), "Input description:")
  text(x = 0, y = 0.95, cex = 0.7, font = 1, adj = c(0,1), paste("Only contains output that passed the quality checks:", CALL$FAST))

  text(x = 0, y = 0.9, cex = 0.7, font = 1, adj = c(0,1), "The following parameters are related to the biological sample selection:")
  text(x = 0, y = 0.85, cex = 0.7, font = 1, adj = c(0,1), paste("     - Minimum number of sample per species:", CALL$SAMPLE_SELECT$MIN_SAMPLE))
  text(x = 0, y = 0.8, cex = 0.7, font = 1, adj = c(0,1), paste("     - Depth range of the samples:", CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH, "to", CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH, "m"))
  text(x = 0, y = 0.75, cex = 0.7, font = 1, adj = c(0,1), paste("     - Time range of the samples:", CALL$SAMPLE_SELECT$START_YEAR, "to", CALL$SAMPLE_SELECT$STOP_YEAR))

  text(x = 0, y = 0.65, cex = 0.7, font = 1, adj = c(0,1), "The following parameters are related to the type of data used:")
  text(x = 0, y = 0.6, cex = 0.7, font = 1, adj = c(0,1), paste("     - Type of raw data query from the data access service:", CALL$DATA_SOURCE))
  text(x = 0, y = 0.55, cex = 0.7, font = 1, adj = c(0,1), paste("     - Type of raw data considered in the analysis (after eventual transformation of the raw data):", CALL$DATA_TYPE))

  text(x = 0, y = 0.45, cex = 0.7, font = 1, adj = c(0,1), "The following parameters are related to the algorithm training:")
  text(x = 0, y = 0.4, cex = 0.7, font = 1, adj = c(0,1), paste("     - Number of splits for the cross-validation folds:", CALL$NFOLD))
  text(x = 0, y = 0.35, cex = 0.7, font = 1, adj = c(0,1), paste("     - Type of cross-validation performed:", CALL$FOLD_METHOD))
  text(x = 0, y = 0.3, cex = 0.7, font = 1, adj = c(0,1), paste("     - Hyperparameter grid size:", CALL$LEVELS))
  text(x = 0, y = 0.25, cex = 0.7, font = 1, adj = c(0,1), paste("     - Compute an ensemble model:", CALL$ENSEMBLE))
  text(x = 0, y = 0.2, cex = 0.7, font = 1, adj = c(0,1), paste("     - Number of bootstrap in projections and partial dependency plots:", CALL$N_BOOTSTRAP))
  text(x = 0, y = 0.15, cex = 0.7, font = 1, adj = c(0,1), paste("     - Projection discontinuity cutoff:", CALL$CUT))

  dev.off()

  # --- 3. List files
  # --- 3.1. Listing all pdfs within the directory
  output_wd <- paste0(project_wd, "/output/", FOLDER_NAME)
  all_files <- list.files(output_wd, pattern = ".pdf", recursive = TRUE)
  # --- 3.2. Remove empty pdf
  to_remove <- file.size(paste0(output_wd, "/", all_files)) != 0
  all_files <- all_files[to_remove]

  # --- 4. Combine and save
  batch_size <- 100

  file_batches <- split(all_files, ceiling(seq_along(all_files) / batch_size))

  batched_pdfs <- list()

  # # Loop over the batches
  for (i in seq_along(file_batches)) {
  # Construct the input file paths for the current batch
  input_files <- paste0(output_wd, "/", file_batches[[i]])

  # Construct the output file path for the current batch
  output_file <- paste0(output_wd, "/BATCH_", i, ".pdf")

  # Combine the PDFs in the current batch
  pdf_combine(input = input_files, output = output_file)

  # Add the path of the batched PDF to the list
  batched_pdfs <- c(batched_pdfs, output_file)
  }

  pdf_combine(input = batched_pdfs, output = paste0(output_wd, "/USER_SUMMARY.pdf"))

  # Loop over the batched_pdfs
  for (file in batched_pdfs) {
    # Remove the file
    file.remove(file)
    }

} # END FUNCTION


