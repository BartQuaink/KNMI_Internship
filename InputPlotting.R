#' ---------------------------
#'
#' Script name: 
#'
#' Purpose of script:
#'
#' Author: Bart Quaink
#'
#' Date Created: 2021-05-12
#'
#' Email: bartquaink@live.nl
#'
#' ---------------------------
#'
#' Notes:
#'   
#'
#' ---------------------------

#' set working directory


#' ---------------------------

#options(scipen = 6, digits = 4) # for no scientific notation

#' ---------------------------

#' load up the packages we will need:

require(pacman) # function to load all libraries

libs <- c("ggplot2", "scales", "reshape2")
pacman::p_load(char = libs)


# Quick helper function: Get reverse y log scale function
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}  

# Pathing -----------------------------------------------------------------

#'
#'
#'

# Read paths from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2){
  stop("Two arguments are required: input_directory--path output_directory--path", call. = F)
}

input_directory <- args[1]
output_directory <- args[2]
