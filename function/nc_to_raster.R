#' =============================================================================
#' @name nc_to_raster
#' @description Extracts environmental data from a .nc file and stores them in a 
#' SpatRaster (from the terra package) per month, allowing easier movement of data 
#' in the pipeline.
#' @param MONTH month to extract the data from
#' @param NC environmental .nc file to extract the data from
#' @param MIN_DEPTH minimum depth range for extraction
#' @param MAX_DEPTH maximum depth range for extraction

nc_to_raster <- function(MONTH, NC, MIN_DEPTH, MAX_DEPTH){
  
  # --- 1. Read .nc file information
  # --- 1.1. Open the netCDF file
  nc <- nc_open(NC)
  
  # --- 1.2. Extract the variable raw data
  # Assumes that the first variable in the .nc file is the one of interest
  varname <- names(nc$var) # Get the variable name
  varname <- varname %>% .[!grepl("bnds|crs", .)] # Remove bounds or crs if the names are in the wrong order (i.e., the .nc not in the right format)
  ncvar <- ncvar_get(nc, varname[1]) # Extract the variable data
  
  # --- 1.3. Extract dimension sizes and reorder according to the variable's shape
  # This ensures that the dimensions match the shape of the extracted variable data
  nc_dimsize0 <- lapply(nc$dim, function(x) length(x$vals)) # Get dimension lengths
  id <- match(dim(ncvar), unlist(nc_dimsize0)) # Match dimensions with variable
  nc_dimsize0 <- nc_dimsize0[id] # Reorder dimensions
  
  # --- 1.4. Extract latitude, longitude, and resolution information
  lon <- as.numeric(nc$dim$lon$vals) # Longitude values
  lat <- as.numeric(nc$dim$lat$vals) # Latitude values
  res <- abs(lon[1] - lon[2]) # Resolution is the difference between two consecutive longitudes
  
  # --- 2. Reorder the variable data and dimensions based on expected dimension names
  id <- match(c("lon", "lat", "time", "depth"), names(nc_dimsize0)) # Match dimensions
  id <- id[!is.na(id)] # Remove any NA values from the matching process
  ncvar <- aperm(ncvar, id) # Transpose the variable data to reorder dimensions
  nc_dimsize <- nc_dimsize0[id] # Reorder the dimension sizes to match
  
  # --- 3. Extract time dimension if it exists
  if (!is.na(names(nc_dimsize)[3]) && names(nc_dimsize)[3] == "time") {
    # --- 3.1. By default, extract the first time layer
    t_id <- 1
    if (nc_dimsize[[3]] > 1) {
      # --- 3.2. If multiple time layers exist, extract the layer corresponding to the specified month
      t_id <- MONTH
      # --- 3.3. Subset the variable data by time
      if (length(nc_dimsize) == 3) {
        ncvar <- ncvar[,,t_id]
      } else {
        ncvar <- ncvar[,,t_id,]
      }
    }
    # --- 3.4. Remove the time dimension from the list of dimensions
    nc_dimsize <- nc_dimsize[-3]
  }
  
  # --- 4. Extract depth dimension if it exists
  if (!is.na(names(nc_dimsize)[3]) && names(nc_dimsize)[3] == "depth") {
    depth <-  nc$dim$depth$vals # Extract the depth available
    d_id <- 1 # Default to first depth layer
    if (nc_dimsize[[3]] > 1) {
      # --- 4.2. Select depth range within the specified minimum and maximum depth
      depth_bnds <- abs(depth[-length(depth)] - depth[-1]) # Calculate depth intervals
      depth <- depth[-length(depth)] # Remove last element of depth vector
      depth_id <- which(depth >= MIN_DEPTH & depth <= MAX_DEPTH) # Find depth layers within the specified range
      
      # --- 4.3. Compute weighted mean across depth bounds for selected layers
      ncvar <- apply(ncvar[,,depth_id], c(1, 2), function(x) sum(x * depth_bnds[depth_id] / sum(depth_bnds[depth_id])))
    }
  }
  
  # --- 5. Create a SpatRaster object using the terra package
  r <- rast(nrow = length(lat), ncol = length(lon), 
            xmin = min(lon) - 0.5 * res, xmax = max(lon) + 0.5 * res, 
            ymin = min(lat) - 0.5 * res, ymax = max(lat) + 0.5 * res)
  
  # Assign the extracted variable data to the SpatRaster object
  values(r) <- as.vector(ncvar)
  
  # --- 6. Reverse latitude and/or longitude if needed
  # Reverse latitude if it is in descending order
  if (lat[1] - lat[2] < 0) {
    r <- terra::flip(r, direction = "vertical")
  }
  # Reverse longitude if it is in ascending order
  if (lon[1] - lon[2] > 0) {
    r <- terra::flip(r, direction = "horizontal")
  }
  
  # --- 7. Close the netCDF file to free up resources
  nc_close(nc)
  
  # --- 8. Return the SpatRaster object and the original dimension sizes
  return(list(r, nc_dimsize0))
}
