#' =============================================================================
#' @name CEPHALOPOD_to_netcdf
#' @title Write mean and standard deviation values to NetCDF
#' @description
#'   This function saves geographical and time-series data in NetCDF format,
#'   including mean and standard deviation values for multiple features.
#' @param r0   Baseline raster
#' @param a_m  A three-dimensional array (geo_coordinates x time x feature) containing mean values.
#' @param a_sd A three-dimensional array (geo_coordinates x time x feature) containing standard deviation values.
#' @param output_file Name of the NetCDF file to be created.
#' @return
#'   A NetCDF file containing mean and standard deviation values of features at geographical and time scales,
#'   saved in the specified output folder.
#' Saved in FOLDERNAME.

CEPHALOPOD_to_netcdf <- function(r0, a_m, a_sd, names, reco_list, output_file = "output_abundance_test.nc") {

  MONTH <- as.list(1:13)
  nlon <- ncol(r0)  # Assuming r0 is defined elsewhere
  nlat <- nrow(r0)  # Assuming r0 is defined elsewhere
  ntime <- length(MONTH)
  nfeatures <- dim(a_m)[3]  # Number of variables (formerly species)
  names <- as.character(names)
  print(names)
  lon_values <- seq(-180, 180, length.out = nlon)
  lat_values <- seq(-90, 90, length.out = nlat)
  
  # Define dimensions
  lon_dim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_values)
  lat_dim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_values)
  time_dim <- ncdim_def(name = "time", units = "months since 1970-01-01", vals = 1:ntime)
  feature_dim <- ncdim_def(name = "feature", units = "", vals = 1:nfeatures)
  
  # Define variables
  mean_var <- ncvar_def(name = "mean_values", units = "NA", dim = list(lon_dim, lat_dim, time_dim, feature_dim), missval = -9999, longname = "Mean Values")
  sd_var <- ncvar_def(name = "sd_values", units = "NA", dim = list(lon_dim, lat_dim, time_dim, feature_dim), missval = -9999, longname = "Standard Deviation Values")
  
  # Create NetCDF file
  nc_file <- nc_create(output_file, vars = list(mean_var, sd_var))
  
  # Add attributes
  ncatt_put(nc_file, "mean_values", "reco", reco_list)
  ncatt_put(nc_file, "mean_values", "feature_names", paste(names, collapse = ","))
  
  # Write coordinate data
  ncvar_put(nc_file, "longitude", lon_values)
  ncvar_put(nc_file, "latitude", lat_values)
  ncvar_put(nc_file, "time", 1:ntime)
  
  # Write data to NetCDF
  for (i in 1:nfeatures) {
    for (m in 1:ntime) {
      tmp_m <- a_m[, m, i]
      tmp_sd <- a_sd[, m, i]
      plot_scale <- quantile(a_m, 0.95, na.rm = TRUE)
      
      # Prepare rasters
      r_m <- setValues(r0, tmp_m)
      r_m[r_m > plot_scale] <- plot_scale
      r_sd <- setValues(r0, tmp_sd)
      
      # Get matrix values from rasters
      rm_values <- getValues(r_m)
      rsd_values <- getValues(r_sd)
      
      # Reshape vector into a matrix with dimensions matching r_m
      mean_values <- matrix(rm_values, nrow = nlat, ncol = nlon, byrow = FALSE)
      sd_values <- matrix(rsd_values, nrow = nlat, ncol = nlon, byrow = FALSE)
      
      # Write values to NetCDF variables
      ncvar_put(nc_file, "mean_values", mean_values, start = c(1, 1, m, i), count = c(nlon, nlat, 1, 1))
      ncvar_put(nc_file, "sd_values", sd_values, start = c(1, 1, m, i), count = c(nlon, nlat, 1, 1))
    }
  }
  
  # Close the NetCDF file
  nc_close(nc_file)
}

