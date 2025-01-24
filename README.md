## **CEPHALOPOD - Comprehensive Ensemble Pipeline for Habitat modelling Across Large-scale Ocean Pelagic Observation Datasets**

This repository contains all code and technical documentation to understand the basis of CEPHALOPOD and install it on your local machine or institution computing cluster. CEPHALOPOD is a standardized and flexible framework for marine habitat modeling across various types and sources of input data. **It is targeted at - but not limited to - plankton observation datasets.**

CEPHALOPOD is also designed as the main tool for the Ecosystem Workbench of the Bluecloud2026 E.U. project. This initiative aims to enhance the accessibility, accuracy and inter-comparability of extensive plankton *in-situ* observation datasets, derived from traditional counting methods, quantitative imaging, and cutting-edge genomic (Metagenome Assembled Genome; MAG) techniques available from repositories such as OBIS, GBIF, AtlantECO or MATOU. It offers a sustainable working for generating ecosystem-level Essential Ocean and Biodiverisity Variables (EOVs & EBVs), adhering to rigorous QA/QC procedures and best practices in the habitat modeling field.

## **How to install it ?**

You can download and install CEPHALOPOD on your local machine or computing cluster by following the steps below. We recommend a minimum of 8 cores and 32 Gb of memory for a general usage across most case studies.

1.  Initiate an R session (version 4.2. or above) on your local machine or computing cluster

2.  Access a terminal at the designated installation path for CEPHALOPOD.

3.  Execute the following command: `git -clone https://github.com/alexschickele/CEPHALOPOD.git` or download CEPHALOPOD as an archive and decompress it in the designated installation path.

4.  Navigate to the file located at: ./code/00_config.R

5.  On the first use, execute the R library installation commands line by line, responding affirmatively to any interactive prompts. This ensures the installation of all necessary libraries and corresponding versions.

6.  CEPHALOPOD should now be ready for utilization !

## How is it structured ?

Now that you have CEPHALOPOD at your designated installation path, you should see several sub-folders. Here we provide an overview of their content and use.

-   The `./master.R` script contains an example script demonstrating the utilization of CEPHALOPOD. This file serves as a template to be customized or replicated for each usage instance of CEPHALOPOD (e.g., by modifying the data source or model parameters; we will come to that later).

-   The `./code` folder contains the collection of R scripts, each corresponding to the function associated with a distinct step within the CEPHALOPOD workflow. All R scripts are numbered in their order of use, matching the numbering in the `./master.R` script.

-   The `./function` folder contains custom R functions of general use across the different steps (e.g., definition of colorbars, operations on spatial data).

-   The `./data` folder corresponds to the directory containing input data and model cache. Please note that some input data are local to the ETH Zurich cluster where CEPHALOPOD was developed, hence not present in this GitHub repository. Refer to the "**Current status**" section below.

-   The `./output` designates the directory where CEPHALOPOD automatically stores its output for each instance or case study.

## How to use it ?

Now that you understand how to install CEPHALOPOD and what it contains, let's dive into an example from a master script. Please note that the following example does not correspond to a real case study and only aims to illustrate how to run CEPHALOPOD in a practical manner. For further information, please refer to associated publications and the technical documentation present in the R scripts in `./code`.

The `./master.R` script starts by a cleanup of your current R session and defines your CEPHALOPOD instance in the `run_name` object. The call to `./code/00_config.R` will load all necessary libraries and functions in your R environment.

```{r}
rm(list=ls())
closeAllConnections()
setwd("/your_installation_path")
source(file = "./code/00_config.R")
run_name <- "test"
```

We kickstart the modeling process by identifying the species available for analysis, from a defined source. We select these species based on specific criteria, such as a minimum sample size, the depth range for sampling, and the temporal range of the data. This selection ensures that the subsequent modeling efforts focus on species that meet the desired criteria. Parameters are the following :

`DATA_SOURCE`: either "occurrence", "biomass", "abundance", "MAG", "custom". This parameter will redirect the query to the appropriate database to source the data from. The latter corresponds to an already formatted .csv, .txt or .xlsx file.

`SAMPLE_SELECT`: a list containing the minimum sample size, target depth range and sample temporal range.

```{r}
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "occurrence",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 100, START_YEAR = 1950, STOP_YEAR = 2020))

```

Now we can look at the `LIST_BIO` object and build a vector of species identifier to extrapolate the distribution from. Here, we select all species in the *Thalassiosira* taxa. This code block can be customed depending on the taxa selection needed for your case study.

```{r}
sp_list <- list_bio %>%
  dplyr::filter(grepl("Thalassiosira ", scientificname)) %>%
  dplyr::select(worms_id) %>%
  unique() %>% pull()
```

In the second step of the modeling process, several critical actions are taken to prepare for subsequent stages. First, a dedicated output folder is created to systematically store the results of each species-level run. Global parameters for modeling are defined, including environmental variables, data types (e.g., proportions), and the selection of specific modeling methods. These parameters are stored into a `CALL.Rdata` object, which is referenced by all subsequent steps. This approach minimizes memory consumption and simplifies the execution of each step. Finally, a local repository of monthly environmental climatologies is established, utilizing `.nc` data. This step is crucial to ensure a uniform set and format of environmental predictors while also optimizing memory usage. Parameters are the following:

-   `FOLDER_NAME` name of the run folder we want to work in

-   `SP_SELECT` vector of identifiers (i.e., usually WoRMS identifiers), corresponding to the species to parallelize on

-   `FAST` `TRUE` or `FALSE`; if `TRUE`, CEPHALOPOD stops the processing of taxa that did not success at a quality check. It avoids unnecessary CPU and memory usage.

-   `LOAD_FROM` load a previous `LIST_BIO` object from another folder to be duplicated in the new `FOLDER_NAME`. It avoids repeating the initial steps of CEPHALOPOD.

-   `DATA_TYPE` the output type of the data, which can influence the sub-folder architecture. See details in the corresponding function.

-   `ENV_VAR` a list of .`nc` files to extract the main variable from, located in `ENV_PATH` . The names correspond to the variable name of the corresponding `.nc` file. If a `!` is specified at the start of the variable name, the variable is excluded from the list.

-   `ENV_PATH` string or vector of path to the root where the .`nc` are.

-   `METHOD_PA` method of pseudo-absence, either "`mindist`" or "`cumdist`" or "`density`" (recommended)

-   `NB_PA` number of pseudo-absences to generate

-   `PER_RANDOM` ratio of pseudo-absences that are sampled randomly in the background

-   `PA_ENV_STRATA` if `TRUE`, increases the pseudo-absence sampling probability in geographical cells with environmental conditions distant from the conditions of the observations.

-   `DIST_PA` if `METHOD_PA = "mindist"`, distance from presences (in meters), from which to define the background data. Expert use only

-   `BACKGROUND_FILTER` additional background filter for finer tuning, such as selecting pseudo-absences within the sampled background of a given campaign or instrument deployment. Passed by the user in the form of a 2 column data frame, x = longitude and y = latitude where the pseudo-absences can be sampled. Or a path to a raster object where pseudo-absences are sampled in non NA cells, weighted by the cell values.

-   `OUTLIER` if `TRUE`, remove outliers further than 2.5 standard deviation from the mean of the observations

-   `RFE` if `TRUE`, performs a Recursive Feature Elimination procedure that discards predictors that do not present significant importance in explaining the patterns in the observations, relative to the ones already selected.

-   `ENV_COR` numeric, removes the correlated environmental predictors from the query objects and `CALL` according to the defined threshold. Else `NULL`.

-   `NFOLD` number of folds, used defined integer

-   `FOLD_METHOD` method used to create the folds, integer between "`kfold`" and "`lon`"; respectively for normal k-fold or longitudinal block of observations (recommended to deal with spatial autocorrelation issues).

-   `MODEL_LIST` list of algorithms to use for calibration and hyperparameter selection

-   `LEVELS` maximum number of parameter values to test in each of the hyperparameter grids

-   `TARGET_TRANSFORMATION` path to a `function(x, REVERSE = T/F)` to transform the target variable prior to the model calibration step.

-   `ENSEMBLE` `TRUE` or `FALSE`; if `TRUE`, computes an ensemble at the evaluation and projection steps

-   `N_BOOTSTRAP` number of bootstrap to do for the projections and partial dependency plots

-   `CUT` numeric or `NULL`; if numeric, quantile (between 0 and 1) at which the projections are considered to be 0. Projection patches without observation are then removed.

The function is assigned to `subfolder_list` as it also returns the list of species to parallelize on.

```{r}
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           FAST = FALSE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "presence_only",
                           ENV_VAR = c("!climatology_s_0_50","!climatology_s_200_300"),
                           ENV_PATH = NULL, # replace by local path to environmental predictors,
                           METHOD_PA = "density",
                           PER_RANDOM = 0.05,
                           OUTLIER = TRUE,
                           RFE = TRUE,
                           ENV_COR = 0.8,
                           NFOLD = 3,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","GAM","RF","MLP","BRT","SVM"),
                           LEVELS = 3,
                           TARGET_TRANSFORMATION = NULL,
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0)
```

Now we can retrieve the biological data for the selected taxa to extrapolate the distribution from. These data form the foundation for training and validating species distribution models. The availability and quality of this data are crucial for the success of the modeling process.

The function is built in parallel over each taxa considered. It does not provide output in the console, as all the retrieve data are saved in a `QUERY.RData` object. As for all subsequent functions, the parameters are following:

`FOLDER_NAME`: the name of the folder in which the run is saved

`SUBFOLDER_NAME`: the name of all sub folders corresponding to each species to parallelize over

```{r}
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

We can now extract the environmental data corresponding to each biological sample at the time and location of the sample. In addition, the data is regridded and binned on a 1 x 1° monthly resolution. Taxa that do not meet the minimum number of occurrence after binning are discarded. Therefore, we assign the ouput of this function to an updated `subfolder_list` object.

```{r}
subfolder_list <- mcmapply(FUN = query_env,
                  FOLDER_NAME = run_name,
                  SUBFOLDER_NAME = subfolder_list,
                  mc.cores = min(length(subfolder_list), MAX_CLUSTERS)) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  .[grep("Error", ., invert = TRUE)] %>%
  as.vector()
```

For occurrence data, we need to artificially generate pseudo-absences to have a balanced dataset between 0 and 1's. This balance is crucial for accurate model training and predictive performance and best practices recommend pseudo-absences following the same biases as presences. They are by default generated following the nearby presence density. The function provides a .`PDF` file with the presence and pseudo-absence distribution in the geographical space (or the observations for continuous or proportion data types).

```{r}
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

In step 6, we rigorously check the biological and environmental dataset quality. First, we identify and handle biological outliers in the data to avoid introducing bias in the model training and extrapolated maps. Then, we assess the quality and relevance of environmental predictor variables used by discarding correlated environmental predictors and discarding predictors that do not explain a significant portion of the observed data. This step is crucial for a parsimonious feature set. The function updates `subfolder_list` with the species passing the quality check on the selected predictors. Finally, we perform a Multivariate Environmental Similarity Surface (MESS) analysis to identify the geographical areas outside the range of environmental values in which the model has been trained (i.e., environmental extrapolation).

Several .`PDF` files are provided in this steps, including a dendrogram of environmental predictor correlation, an a priori predictor importance and ranking of the predictors, as well as the corresponding loss after which, adding a predictor would not significantly increase the model performance.

```{r}
subfolder_list <- mcmapply(FUN = query_check,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS)) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  as.vector()
```

Now that the input data have been thoroughly quality checked, we can create splits and resampling folds to facilitate the training and assessment of species distribution models. Proper partitioning of data is essential for robust modeling results, and assessing model performance in reproducing the observed patterns. The model training design is oriented around a *n-time* cross validation between train and test set to find the best hyperparameters for each algorithm. Then, each algorithm is tested against a final evaluation set to assess its performance against an independent dataset.

```{r}
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

Each algorithm is associated to a set of hyperparameter to best reproduce the observed biological target. This step focuses on configuring hyperparameters for training the species distribution models algorithms. These hyperparameters dictate the model's behavior during training and significantly influence model performance.

```{r}
hyperparameter(FOLDER_NAME = run_name)
```

At this stage, we defined all necessary parameters and quality checked all input data. Therefore, we can start to fit species distribution models to the data. We employ various modeling algorithms, such as Generalized Linear Models (`GLM`), Generalized Additive Models (`GAM`), Random Forest (`RF`), Support Vector Machines (`SVM`), Boosted Regression Trees (`BRT`) and Multilayer Perceptrons (`MLP`), to capture the relationships between environmental features and the biological target(s).

From this steps on, all outputs are saved in a `MODEL.Rdata` object in each species corresponding sub folder. The input parameters remain the same as in the previous steps.

```{r}
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

To assess each of the algorithm's performance, we perform a serie of quality checks on the evaluation dataset.

-   `FIT`: the metric quantifying the fitting performance of the algorithm. How well does it reproduces the patterns in the observed data?

-   `VIP`: the variable importance (PDF output), that quantifies the contribution of each environmental predictor.

```{r}
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

All algorithm that successfully passed the evaluation step (if `FAST = TRUE`) are considered of sufficient quality to built spatial projections. We retrain a full model with the parameters and hyperparameters previously selected and perform spatial projections for *n-bootstrap* resamples to estimate the associated uncertainty. A quality check is performed on this uncertainty estimation.

-   `NSD`: the normalized standard deviation, quantifying the uncertainty between bootstraps at each grid-cell, averaged over all grid-cells of the projection domain.

```{r}
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

Following the projections, CEPHALOPOD provides a summarized `PDF` output of the each projections and quality checks for each taxa, algorithm and ensembles. The projections can be provided in a monthly resolution or any lower temporal resolution.

```{r}
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

In addition, we are also performing partial dependency plots for each species and environmental predictors. It shows the marginal response of the habitat suitability to each environmental feature. The corresponding outputs are also saved in a `PDF` file.

```{r}
mcmapply(FUN = pdp,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))
```

Finally, all PDF files present in each species folder can be summarized in a unique summary for users, concatenating all quality checks and outputs of the present workbench ecosystem run.

```{r}
user_synthesis(FOLDER_NAME = run_name)
```

Supplementary analysis such as diversity estimates can be performed by any user based on the information, data and output stored in the `QUERY` and `MODEL.Rdata` objects.

## Current status (Jan, 24th, 2025)

Here you will find updated informations concerning the code status from a development perspective.

-   The predictor set can be found online at: <https://data.d4science.net/m9WC>
