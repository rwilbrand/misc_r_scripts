#------------------------------------------------------------------------------
# Collection of various convenience functions written by Robert Wilbrand
# Reproduction permitted with attribution
# 
# Required packages: tidyverse, raster
#------------------------------------------------------------------------------
library(raster)
library(tidyverse)
#------------------------------------------------------------------------------

# Wrapper function to time execution of a function
speedtest <- function(func, ..., returnOutput = TRUE){
  # returns time elapsed alongside the actual output
  # set returnOutput=FALSE to output ONLY time elapsed
  starttime <- Sys.time()
  out <- func(...)
  elapsed <- Sys.time() - starttime
  if (returnOutput){
    output <- list("Output" = out,
                   "Elapsed" = elapsed)
  } else output <- elapsed
  return(output)
}

# Wrapper function to set a specified seed prior to function execution
seedFunc <- function(func, seed, ...){
  set.seed(seed)
  func(...)
}

# Given a data frame, generate exhaustive formula (interactions optional)
exhFormula <- function(df, i = 1, exclude = NULL, interx = FALSE){
  # Dependent variable is assumed to be in column 1 by default
  # Columns not to be used can be specified by index with exclude
  sep <- ifelse(interx, " * ", " + ")
  j <- c(i, exclude)
  paste(colnames(df)[i], " ~ ",
        paste(colnames(df)[-j],
              collapse = sep)) %>%
    as.formula
}

# Function to split geometry of sfc_POINT objects into x any y columns
xyPoints <- function(sfcp){
  # intended for use with gridded sfc_POINT objects
  # such that geom_raster() produces nicer plots than geom_sf()
  sfcp %>%
    as_tibble %>%
    mutate(across(geometry, ~as.character(.x) %>%
                    str_remove_all("[:alpha:]|[(,)]"))) %>%
    separate(geometry, into = c("x","y"), sep = "[[:blank:]]", convert = T)
}

# Determine lead digit of a number for Benford's Law testing
# using string manipulation rather than a mathematical approach
leadDigit <- function(n){
  as.character(n) %>% 
    str_remove_all("[-.0]") %>% 
    str_sub(1,1) %>% 
    as.integer
}

# Function to subset multiple elements from a list into a new list
# takes numbered or named lists as input with a corresponding vector
pickFromList <- function(ls, vc) map(vc, ~pluck(ls, .x))

# Functions to test whether queried bits are on (bits start at 1)
bitTest    <- function(n, bit) intToBits(n)[bit] == T
bitTestAll <- function(n, ...) all(as.logical(intToBits(n)[c(...)]))
bitTestAny <- function(n, ...) any(as.logical(intToBits(n)[c(...)]))

# Functions to test whether queried bits are on (bits start at 0)
bitTest0    <- function(n, bit) intToBits(n)[bit + 1] == T
bitTestAll0 <- function(n, ...) all(as.logical(intToBits(n)[c(...) + 1]))
bitTestAny0 <- function(n, ...) any(as.logical(intToBits(n)[c(...) + 1]))

#------------------------------------------------------------------------------
# This section is for functions to work with rasters
# Note that R is inefficient at dealing with large rasters
# Consider a different language (e.g. Python) or GIS program if necessary

# Random generation functions are mostly to generate artificial test data
# Generate random float value raster with dimensions n*m
genFloatRaster <- function(n, m, low = 0, high = 1, seed = NULL){
  # optional: use specified seed
  if (!is.null(seed)) set.seed(seed)
  r <- raster(xmn = 0,
              ymn = 0,
              ymx = m,
              xmx = n,
              vals = runif(n*m, min = low, max = high),
              crs = NULL,
              resolution = 1)
  return(r)
}

# Generate random binary or unary raster with dimensions n*m
genBinaryRaster <- function(n, m, seed = NULL, unary = FALSE){
  # set unary = TRUE to produce unary raster instead
  # optional: use specified seed
  if (!is.null(seed)) set.seed(seed)
  r <- raster(xmn = 0,
              ymn = 0,
              ymx = m,
              xmx = n,
              vals = sample(0:1, n*m, replace = T),
              crs = NULL,
              resolution = 1)
  if (unary) r[r == 0] <- NA
  return(r)
}

# Generate random integer raster with dimensions n*m
genIntegerRaster <- function(n, m, high, seed = NULL, includeZero = FALSE,
                             ValNA = NULL, categorical = FALSE){
  # set categorical = TRUE to produce categorical raster
  # optional: use specified seed
  if (!is.null(seed)) set.seed(seed)
  low <- ifelse(includeZero, 0, 1)
  r <- raster(xmn = 0,
              ymn = 0,
              ymx = m,
              xmx = n,
              vals = sample(low:high, n*m, replace = T),
              crs = NULL,
              resolution = 1)
  if (!is.null(ValNA)) r[r == ValNA] <- NA # set NA Value if specified
  if (categorical) r <- ratify(r)
  return(r)
}


# For a binary raster, returns majority pixel value for specified filter size 
majority <- function(rast, fsize){
  # fsize has to be an odd integer value
  m <- matrix(1,fsize,fsize)
  focal(rast, m, mean, pad = T, padValue = 0.5) %>% 
    round
}

# For a binary raster, fill in pixels surrounded by foreground
fillGaps <- function(rast){
  m <- matrix(1,3,3)
  m[2,2] <- 9 # set center pixel weight to 9 to differentiate foreground
  rcl <- cbind(0:17, c(rep(0,8), rep(1,10))) # reclassification matrix
  focal(rast, m, sum, pad = T, padValue = 1) %>% 
    reclassify(rcl)
}

# For a binary raster, remove pixels without neighbors
killIsolated <- function(rast){
  m <- matrix(1,3,3)
  m[2,2] <- 9 # set center pixel weight to 9 to differentiate foreground
  rcl <- cbind(0:17, c(rep(0,10), rep(1,8))) # reclassification matrix
  focal(rast, m, sum, pad = T, padValue = 0) %>% 
    reclassify(rcl)
}

# For a binary raster, remove isolated pixels and fill surrounded pixels
cleanupRaster <- function(rast){
  rast %>% fillGaps %>% killIsolated
}

# For a binary raster, detect edges
edgeFinder <- function(rast){
  # For default case of edgewidth 1 pixel and queen contiguity
  # Edge = foreground with at least 1 non-foreground neighbor
  # Outer limits of the raster also count as edge
  m <- matrix(1,3,3)
  m[2,2] <- 9 # set center pixel weight to 9 to differentiate foreground
  rcl <- cbind(0:17, c(rep(0,9), rep(1,8), 0)) # reclassification matrix
  focal(rast, m, sum, pad = T, padValue = 0) %>% 
    reclassify(rcl)
}

# For a binary or unary raster, return patch label and size
accounting <- function(rast, contiguity = "rook"){
  # beware of using with large rasters
  # speed test for completely random rasters:
  # 100*100 : ca. 23 sec
  # 200*200 : > 2 mins
  # 500*500 : > 14 mins
  # improvements of computing time may be possible with a more sophisticated
  # algorithm, e.g. computing minima row- and column-wise first
  
  # turn into unary raster
  rast[rast == 0] <- NA
  # create auxiliary raster with values = cell number
  aux.ras <- raster(rast)
  values(aux.ras) <- 1:ncell(aux.ras)
  lab.ras <- aux.ras*rast
  # set weights matrix according to contiguity case
  if (contiguity == "rook"){
    weights <- matrix(c(NA,1,NA,1,1,1,NA,1,NA), nrow = 3, byrow = T)
  } else weights <- matrix(1,3,3)
  # first iteration of consolidating patches
  new.lab <- focal(lab.ras, weights, min, na.rm = T, pad = T)*rast
  n_iter <- 1
  # replace raster with updated raster until they are identical
  while(!identical(lab.ras, new.lab)){
    lab.ras <- new.lab
    new.lab <- focal(lab.ras, weights, min, na.rm = T, pad = T)*rast
    n_iter <- n_iter + 1
  }
  # relabel sequentially
  new.ids <- unique(new.lab[!is.na(new.lab)])
  rcm2 <- cbind(new.ids, 1:length(new.ids))
  new.lab <- reclassify(new.lab, rcm2)
  # perform pixel count per patch
  zone_mat <- zonal(rast, new.lab, 'count')
  n_pixel <- reclassify(new.lab, zone_mat)
  # prepare output
  out_stack <- stack(new.lab, n_pixel, n_iter*rast)
  names(out_stack) <- c("PatchLabels",
                        "PatchSizes",
                        "Original") # Value = #iterations performed
  return(out_stack)
}

# Plot selected layer of a raster (stack) with location of min and max values
minmaxplot <- function(rast, lyr = 1){
  rstlyr <- rast %>% subset(lyr)
  xymin <- which.min(rstlyr) %>% xyFromCell(rast, .)
  xymax <- which.max(rstlyr) %>% xyFromCell(rast, .)
  plot(rast, lyr)
  points(xymin, pch = 3, cex = 0.5)
  points(xymax, pch = 18, col = 'red')
}
