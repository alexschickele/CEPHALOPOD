
#' @name target_transformation
#' @param x vector of target values
#' @param REVERSE boolean, reverse the transformation or not ?
#' @param PARAM parameters used for the transformation if reverse TRUE
#' @return the transformed values of x

target_transformation <- function(x, REVERSE = FALSE, PARAM = NULL){
  library(bestNormalize)
  
  # --- 1. Forward transformation /w. REVERSE = FALSE
  if(REVERSE == FALSE){
    # --- 1.2. Apply the transformation
    quantiles <- quantile(x, seq(0,1,0.01), na.rm = TRUE)
    x_out <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
    
    # --- 1.3. Wrap up
    return(list(out = x_out, q_obj = quantiles))
  } # end reverse false
  
  # --- 2. Reverse transformation
  if(REVERSE == TRUE){
    # --- 2.1. Load transformation parameters from query object
    quantiles <- PARAM$q_obj
    
    # --- 2.2. Apply the transformation
    x_out <- quantiles[x] %>% as.numeric()
    
    # --- 2.3. Wrap up
    return(x_out)
    
  } # end reverse true
  
} # end function
