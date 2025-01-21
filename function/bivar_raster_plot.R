#' Functions for bivariate raster plotting in R
#' 
#' @function colmat :
#' @return : two-dimensional color matrix
#' 
#' @param pal : a color palette used for the first dimension variable
#' @param saturation : the degree of saturation related to the second variable
#' @param xlab @param ylab : second (x) and first (y) variable names
#' 
#' @function colmat_plot
#' @return : a plot of the color matrix
#' 
#' @param colourmatrix : output of colmat
#' @param xlab @param ylab : second (x) and first (y) variable names
#' 
#' @function bivar_map
#' @return a raster with the cell values corresponding to the index of colmat
#' 
#' @param rasterx : raster corresponding to variable 1, same format as rastery
#' @param rastery : raster corresponding to variable 2, same format as rasterx
#' @param colourmatrix : output of colmat
#' @param cutx : user defined interval for the color classes of x
#' @param cuty : user defined interval for the color classes of y
#' 
#' @references
#' Simplified version of :
#' https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188#file-bivarrasterplot-r

# Function that produces the color matrix
colmat <- function(xmax = "#1F867B",
                   ymax = "#B64A60",
                   nbreaks = 100){
  
  # --- 1. Required Packages
  require(RColorBrewer)
  require(colorspace)
  
  # --- 2. Build the two extremes color columns
  left_pal <- colorRampPalette(c("white",ymax))(nbreaks)
  right_pal <- colorRampPalette(c(xmax, "black"))(nbreaks)
  
  # --- 3. Interpolate colors in between
  col_matrix <- NULL
  for(i in 1:nbreaks){
    col_matrix <- rbind(col_matrix, colorRampPalette(colors = c(left_pal[i], right_pal[i]))(nbreaks))
  }
  return(col_matrix)
}

# Plotting the color matrix

colmat_plot <- function(colormatrix, xlab, ylab){
  image(x = matrix(seq(1,length(colormatrix)), nrow = nrow(colormatrix), ncol = ncol(colormatrix), byrow = TRUE),
        col = colormatrix,
        main = "",
        xlab = xlab, ylab = ylab,
        axes = FALSE)
  box()
  grid(lty = "solid", col = "black")
}

# Function to assign color-code to raster data
bivar_map <- function(rasterx, rastery, colormatrix, cutx = NULL, cuty = NULL){
  require(terra)
  require(virtualspecies)
  
  # Matrix of color IDs, where rows represent splits of rastery and columns represent splits of rasterx
  colorid <- matrix(seq(1:length(colormatrix)), nrow = nrow(colormatrix), ncol = ncol(colormatrix), byrow = TRUE)
  
  # Use classify to split raster values into bins based on the number of rows/columns of the colormatrix
  splitx <- classify(rasterx, nrow(colormatrix))  # Split into bins based on the number of rows in the colormatrix
  splity <- classify(rastery, ncol(colormatrix))  # Split into bins based on the number of columns
  
  # Handle custom breaks if provided
  if(!is.null(cutx) & !is.null(cuty)){
    if(length(cutx) == (ncol(colormatrix) + 1) & length(cuty) == (nrow(colormatrix) + 1)){
      # Ensure raster values are within cutx and cuty bounds
      rasterx[rasterx < min(cutx) | rasterx > max(cutx)] <- NA
      rastery[rastery < min(cuty) | rastery > max(cuty)] <- NA
      
      splitx <- classify(rasterx, cutx)
      splity <- classify(rastery, cuty)
    } else {
      stop("cutx and cuty dimensions should be respectively nrow(colormatrix) + 1 and ncol(colormatrix) + 1")
    }
  }
  
  # Create an empty raster z based on rasterx, initialize with NA
  z <- terra::rast(rasterx)
  values(z) <- NA  # Initialize all cells as NA
  
  # Get the valid (non-NA) cells from splitx and splity
  valid_cells <- !is.na(values(splitx)) & !is.na(values(splity))
  
  # Get the indices from splitx and splity for the valid cells
  valid_splitx <- values(splitx+1)[valid_cells] # +1 because it gets the lower bound (0 value = index 1 in color_id)
  valid_splity <- values(splity+1)[valid_cells] # +1 because it gets the lower bound (0 value = index 1 in color_id)
  
  # Use the bin indices from valid_splitx and valid_splity to get the correct color index from colorid
  values(z)[valid_cells] <- colorid[cbind(valid_splitx, valid_splity)]  # Note: Rows (splity) and Columns (splitx)
  
  # Define color plot based on the range of values in z
  col_plot <- colormatrix[min(values(z), na.rm = TRUE):max(values(z), na.rm = TRUE)]
  
  return(list(z, col_plot))
}
