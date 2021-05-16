#' FILE: 
#'
#' This file has the required (pre)processing functions for working with IASI
#' data in 1DVar. Each function has a header which explains the usage, outpu-
#' t and the resulting output.


# set up required libs
libs <- c("raster", "ggplot2", "ncdf4", "rgdal", "ggspectra", "BBmisc", 
          "elevatr", "lubridate", "readr", "matrixcalc", "reshape2", "Hmisc",
          "dplyr", "tidyr", "zoo", "knitr")
pacman::p_load(char = libs)




# Create Input Files ------------------------------------------------------
#' This function creates all the required input files at once, calling the 
#' functions written below. Only the input IASI file and the specified output
#' path for all the files to be put in are required.


IASI_1Dvar_Inputs <- function(nc, output_path) {
  
  
  
  
}




# Spectra scaling ---------------------------------------------------------
#' Returns the correctly scaled radiances with the correct wavelengths too 
#' with just the spectral information present (8461 channels)

IASI_scaling_factors <- function(nc) {
  
  # Get amount of scaling bands used
  scaling_bands <- ncvar_get(nc, "i_def_scale_sond_nb_scale")
  
  # Take starting observations, need to take off for correct index and values
  start_index <- ncvar_get(nc, "i_def_ns_first_1b")[1] - 1 # take first one, all values are the same, -1 for 1 indexing in R
  
  # Get indexes
  startbins <- ncvar_get(nc, "i_def_scale_sond_ns_first")[1:scaling_bands] - start_index
  endbins <- ncvar_get(nc, "i_def_scale_sond_ns_last")[1:scaling_bands] - start_index
  
  # Get scaling values
  factors <- 10^-(ncvar_get(nc, "i_def_scale_sond_scale_factor")[1:scaling_bands])
  
  # create dataframe
  df <- data.frame(factor = factors, start = startbins, end = endbins)
  
  # Return
  return(df)
  
}



# Apply scaling -----------------------------------------------------------

IASI_Spectra_Scaled <- function(spectra, scaling_factors) {
  
  return(spectra[1:8461] * rep(scaling_factors$factor, scaling_factors$end - scaling_factors$start + 1))
  
}




# Wavenumber calculation --------------------------------------------------
#' This function calculates the wavenumbers based on formula as supplied by
#' the IASI level 1C document
#' wavenumber_1b(k) = 'i_def_spect_dwn1b' × (i_def_ns_first_1b + k - 2) / 100


calculate_Wavenumbers <- function(dwn1b, first) {
  
  # preallocate vector, is factor 10> faster
  wavenumber = c(length(8461))
  
  for(i in 1:8461) {
    wavenumber[i] <- dwn1b * (first + i - 2) / 100
  }
  
  return(wavenumber)
}

# Spectra loop ------------------------------------------------------------

#'
#'
#'
#'
#'
#'


write_1DVar_input <- function(basefile, spect, lat, lon, date, satzen, solzen, filen, BT = FALSE) { 
  
  # Read in basefile information used to write the converted IASI spectra in
  basef <- readLines(basefile)
  
  # Format the lines as shown above
  lin1 <- paste0("Obs ID:              1 Obs Type:          3 Satellite ID:     4")
  lin2 <- paste0("Latitude: ", sprintf("%.3f", lat), " Longitude: ", sprintf("%.3f", lon), " Elevation: 0.0")
  lin3 <- paste0("Surface type:   1", " Sat Zen Angle:   ", sprintf("%.3f", satzen), " Solar Zen. Ang.: ", sprintf("%.3f", solzen))
  
  # Get day/month/year info
  lind <- paste0("Year: ", year(date), 
                 " Month: ", sprintf("%02d", month(date)), 
                 " Day: ", sprintf("%02d", day(date)))
  
  # Append to basefile, and If BT then write Brigthness Temperatures.
  if(BT){
    basef2 <- append(basef, c(lin1, lind, lin2, lin3, "Brightness Temperatures:")) 
  }
  else{
    basef2 <- append(basef, c(lin1, lind, lin2, lin3, "Radiances:")) 
  }
  
  # Do next line if width of each element should be the same
  spect_format <- formatC(as.numeric(spect), format = 'f', flag=' ', digits = 3, width = 13)

  # Write to file
  cat(basef2, paste0(tapply(spect_format, ceiling(seq_along(spect)/6), 
                    paste0, collapse = ''), collapse = '\n'), file = filen, sep = "\n")
  
  
}

# Writing -----------------------------------------------------------------

## For testing  ####################################################################

nc <- nc_open("E:/_InternshipData/IASI/Level1c/2011/August/20110801004800_24810.nc")
basefile <- "E:/_InternshipData/1DVar_InputFiles/IASI_basefile.txt"
filen_path <- "E:/_InternshipData/1DVar_InputFiles/Outputs/" # folder to put in the data

####################################################################################

spectra_Prep <- function(nc, basefile, path) {
  
  # Set prefix
  filen_prefix <- "IASI_Formatted_"
  
  # Read all spectra in one file
  spect_all <- ncvar_get(nc, "gs_1c_spect")
  
  # Get dimension of the netcdf file, needed for all looping
  # dimensions are shown as [full spectra, snot, px, along track pixel]
  dimensions <- dim(spect_all)
  
  # Get scaling factors
  scaling_factors <- IASI_scaling(nc)
  
  # loop over all pixels, [k, j, i], so it starts with the first along track observations
  # then loops to the 30 side scanning pixels, finally the 4 pixels observations in the IASI
  # observation matrix. [along track, snot, PN]
  for(i in 1:dimensions[4]) {
    for(j in 1:dimensions[3]) {
      for(k in 1:dimensions[2]) {
        
        # Get required values for lat/lon/zenith angles
        lat <- ncvar_get(nc, "lat")[k,j,i]
        lon <- ncvar_get(nc, "lon")[k,j,i]
        satzen <- ncvar_get(nc, "pixel_zenith_angle")[k,j,i]
        solzen <- ncvar_get(nc, "pixel_solar_zenith_angle")[k,j,i]
        
        # Get correct date
        date <- as.POSIXct(ncvar_get(nc, "measurement_date"), origin="2000-01-01")[k, j]
        
        # select the spectra, the first 8461 channels contain the information
        spect <- spect_all[1:8461, k, j, i]
        
        # make new array to store, remove later for 
        spect_scaled <- array()
        
        # convert to 1D_Var required scaling factors
        spect_scaled[1:3340] <- spect[1:3340] * (10^-2)
        spect_scaled[3341:6428] <- spect[3341:6428] * (10^-3)
        spect_scaled[6429:6960] <- spect[6429:6960] * (10^-4)
        spect_scaled[6961:8140] <- spect[6961:8140] * (10^-3)
        spect_scaled[8141:8461] <- spect[8141:8461] * (10^-4)
        
        # Set filename
        filen <- paste0(path, filen_prefix, "AT", i, "_SN", j, "_PN", k, ".dat")
        
        # Write to file in correct format
        write_1DVar_input(basefile, spect = spect_scaled, lat = lat, lon = lon, date = date, solzen = solzen, satzen = satzen, filen = filen)
        
        
      }
    
    }   
    
  }
    
} # end function



# Run all -----------------------------------------------------------------

spectra_Prep(nc, basefile)



# Function: Convert -------------------------------------------------------

## This function converts mass mixing ratios in kg/kg (MMR) to ozone ppmv (VMR) in 
## standard dry air conditions. Conversion taken from ECMF site:
## link: https://confluence.ecmwf.int/pages/viewpage.action?pageId=153391710
## standard is FALSE: from MMR to VMR, if reverse = TRUE then VMR to MMR is
## Another method is from wdc.dlr.de : MMR = VMR * 1.66 * 1e-6

convert_MMR_VMR <- function(O3, rev = FALSE) {
  if(!rev) {
    return(28.9644 / 47.9982 * 1e6 * O3)
  }
  else if(rev){
    return((O3 * 47.9982) / (28.9644 * 1e6))
  }
}





# Function: Apriori -------------------------------------------------------

## this function reads the input IASI spectra and creates a correct corresponding
## background profile to be used as input. Example of such an input is in BACKGR-
## OUND_54L.dat. Headers are used to make each format too. Ozone is formatted in 
## kg/kg, so the convert_MMR_VMR is required too. Required variables are: 
##        lat:      latitude taken from IASI pixel (ncvar = "lat")
##        t:        Seconds since jan 1st 2000 (ncvar = "measurement_date")

create_BACKGROUND <- function(lat, obs_date, filen) {
  
  # set required paths
  path_ML <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_table.dat"
  path_Background <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/test_data/IASI_54L/InputData/BACKGROUND_54L.dat"

  # Base information --------------------------------------------------------
  #' This section takes the base information on which the background file will
  #' be built on. Can be expanded in the future with personal observation data
  #' but for the moment it takes the base 54 layer file and adds own ozone data
  #' to it.
  
  # Set header lines, 10 in total, this below is taken from standard file and is subject 
  # to change
  header <- c("Background file generated from
/home/h01/frpw/1DVar/sim_bg/US_Standard.dat
54 fixed layers
With simulated errors based on B-matrix
O3 converted to kg/kg (assuming dry conversion for simplicity)\n\n\n\n")
  
  # Set top information
  info <- c("x------------------- End of Header------------------------------------x
No. Background Profiles:            1
No. of Levels/Profile:             54
Unit for q:                         1  (1=ppmv, 2=kg/kg, 3=RH)
x---------------------------------------------------------------------x
Profile #    1 Follows
x---------------------------------------------------------------------x")
  
  # Read in Pressure, Temperature and Water vapour profiles
  bg_cols <- read.table(path_Background, skip = 16, nrows = 54, header = FALSE)
  
  # GCopy surface data and parameters
  surfpams <- read_lines(path_Background, skip = 70)

  # Create ozone column -----------------------------------------------------
  #' In this section the ozone column is being made. Based on the location
  #' information taken from the total retrieval file for each observation.
  #' 
  
  # round to nearest full 10-latitude zone, and get index by /10
  latidx <- round(lat, digits = -1) / 10
  
  # Get dataframe of full climatology, accounting for negative values beginning by southern hemisphere
  clim_df <- read.table(path_ML, skip = ((latidx + 9) * 70) + 1, nrow = 66, header = TRUE)
  
  # Date of IASI measure (on board UTC): Number of seconds since 1 January 2000 00:00
  # So convert this to a date, and get the month from this info. Standard is UTC. 
  # Month corresponds to the column of the clim_df + 1
  month <- month(as.POSIXct(obs_date, origin="2000-01-01"))

  # Get vertical levels and climatology
  z_levels <- clim_df[, 1]
  clim <- clim_df[, month + 1]
  
  # Convert to kg/kg
  clim_kg <- convert_MMR_VMR(clim, rev = TRUE)
  
  # Convert z_levels to hPa, following standard barometric formula
  pressure <- zstar_press(z_levels)
  
  # Create sub dataframe for interpolation
  df1 <- data.frame(pressure = pressure, kgkg = clim_kg)
  df2 <- data.frame(pressure = bg_cols$V1, kgkg = 0)
  
  # Interpolate to 54 layers
  df3 <- df1 %>% 
    complete(pressure = union(pressure, df2$pressure)) %>%
    mutate(kgkg = na.approx(kgkg, maxgap = Inf, rule = 2))
  
  # Return dataframe with values based on pressure levels from base background file
  background_vals <- merge(df2, df3, by = "pressure")[, 3]
  
  # Create full dataframe with merged ozone
  backgr_df <- data.frame(bg_cols[,1:3], ozone = rev(background_vals))
  
  # Format columns
  backgr_df$V1 <- sprintf("%.5f",backgr_df$V1)
  backgr_df$V2 <- sprintf("%.5f",backgr_df$V2)
  backgr_df$V3 <- sprintf("%.6E", backgr_df$V3)
  backgr_df$ozone <- sprintf("%.6E", backgr_df$ozone)
  
  # Write to file
  backgr_formatted <- kable(backgr_df, digits = 12, format = "simple", align = "r")[c(-1,-2)]
  cat(header, info, backgr_formatted, surfpams, sep = '\n', file = filen)
  

}

# Function: Create Bmatrix ------------------------------------------------
# lat: integer, date: Po

create_Bmatrix <- function(lat, date, c_dampen = 1E-9, fac = 0.5) {
  
  # Pathing
  path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/test_data/IASI_54L/InputData/Bmatrix_54L"
  path_ml <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_stats.dat"
  
  # Read as table
  bm_df_lines <- scan(path, quiet = TRUE, skip = 3, nlines = 1548)
  
  # Make equal size matrix
  bm_df <- matrix(bm_df_lines, nrow = 86, ncol = 86)
  
  # Create new total matrix, fill with 0 and fill 1:86 with base BMatrix
  bm_ozone <- matrix(0, nrow = 140, ncol = 140)
  bm_ozone[1:86, 1:86] <- bm_df
  
  # Get index latitude in 10-latitudinal zone
  latidx <- round(lat, digits = -1) / 10
  
  # Get dataframe of full climatology, accounting for negative values beginning by southern hemisphere
  ppmv_stats <- read.table(path_ml, skip = ((latidx + 9) * 70) + 1, nrow = 66, header = TRUE)
  
  # Convert to kg/kg
  ppmv_stats[, -1] <- convert_MMR_VMR(ppmv_stats[, -1], rev = TRUE)

  # Get correct column of the month
  table <- data.frame(zlevel = ppmv_stats$Z.level, ppmv_stats = ppmv_stats[, month(date) + 1])
  
  # convert zstar altitudes to pressure, following supplied formula in README.txt in ML climatology
  pressure_conv <- zstar_press(table$zlevel)
  
  # Interpolate to 54 layers
  stats_54L <- data.frame(
                  pressure = approx(pressure_conv, n=54)$y, 
                  values = approx(table$ppmv_stats, n=54)$y
                  )
  
  # Make covariance
  cov_matrix <- cov_profile(stats_54L$pressure, stats_54L$values, fac)
  
  # Mirror matrix so TOA is first
  cov_mirrored <- cov_matrix[ncol(cov_matrix):1, nrow(cov_matrix):1]
  
  # merge matrices
  bm_ozone[87:140, 87:140] <- cov_mirrored
  
  # Dampen with factor c (standard is c = e-09)
  diag(bm_ozone) <- diag(bm_ozone) + c_dampen
  
  # check if singular
  if(is.non.singular.matrix(bm_ozone)) {
    stop("Matrix is not positive definite, add larger dampening factor")
  }
  
  # Return BMatrix
  return(bm_ozone)

}



# Function: Alt2Pressure --------------------------------------------------

alt2pres <- function(altitude){

"
Determine site pressure from altitude.

Parameters
----------
Altitude : scalar or Series
Altitude in meters above sea level 

Returns
-------
Pressure : scalar or Series
Atmospheric pressure (Pascals)

Notes
------
The following assumptions are made

============================   ================
Parameter                      Value
============================   ================
Base pressure                  101325 Pa
Temperature at zero altitude   288.15 K
Gravitational acceleration     9.80665 m/s^2
Lapse rate                     -6.5E-3 K/m
Gas constant for air           287.053 J/(kgK)
Relative Humidity              0%
============================   ================

References
-----------

'A Quick Derivation relating altitude to air pressure' from Portland
State Aerospace Society, Version 1.03, 12/22/2004.
"

    #press <- ((44331.514 - altitude) / 11880.516) ** (1 / 0.1902632) # original
    #press <- 101325 * ( (1 - 2.25577 * 10**-5 * (altitude)) ** 5.25588 ) # second method
    
    # method following : https://www.math24.net/barometric-formula
    press <- 1013.25 * exp(-0.00012 * altitude)

    return(press)
}




# Function: Zstar to pressure ---------------------------------------------

#' This function calculates the pressure based on the zstar vertical level
#' as supplied with the McPeters-Labow Climatology. Input is a profile or list
#' of altitudes given in zstar levels.

zstar_press <- function(zstar) {
  return( 1013 / ( 10^(zstar/16) ) )
}




# Function: covariance of profile -----------------------------------------

#' This function calaculates the covariance matrix based on the error covariance
#' from the MCPeters labow climatology. Following the inverse to st.dev to covariance
#' method

cov_profile_v1 <- function(cov) {
  
  # Set dividing factor
  fac <- 1
  
  # Create equal sized empty matrix for calculations
  cov_matrix <- matrix(0, nrow = length(cov), ncol = length(cov))
                       
  # Put covariance list on the diagonal
  diag(cov_matrix) <- cov
  
  # Calculate error covariances
  for(i in 1:nrow(cov_matrix)) {
    for(j in 1:ncol(cov_matrix)) {
      # calculate covariance in i,j
      cov_matrix[i,j] <- ( sqrt(cov_matrix[i,i]) * sqrt(cov_matrix[j,j]) )
    }
  }
  
  # Return full matrix
  return(cov_matrix)
}


# Function: covariance of profile -----------------------------------------

#' This function calaculates the covariance matrix based on the error covariance
#' from the MCPeters labow climatology. Following the inverse to st.dev to covariance
#' method. Second testing method

cov_profile_v2 <- function(cov) {
  
  # Set dividing factor
  fac <- 0.5
  
  # Create equal sized empty matrix for calculations
  cov_matrix <- matrix(0, nrow = length(cov), ncol = length(cov))
  
  # Put covariance list on the diagonal
  diag(cov_matrix) <- cov
  
  # Calculate error covariances
  for(i in 1:nrow(cov_matrix)) {
    for(j in 1:ncol(cov_matrix)) {
      # get P
      P <- exp( (-abs(log10(cov_matrix[i,i]) - log10(cov_matrix[j,j]) ) )/ fac  )
      
      # calculate covariance in i,j
      cov_matrix[i,j] <- ( sqrt(cov_matrix[i,i]) * sqrt(cov_matrix[j,j]) ) * P
                         
    }
  }
  
  # Return full matrix
  return(cov_matrix)
}


# Function: covariance of profile -----------------------------------------

#' This function calaculates the covariance matrix based on the error covariance
#' from the MCPeters labow climatology. Following the inverse to st.dev to covariance
#' method. 

cov_profile <- function(pressure, st_dev, fac = 0.5) {
  
  # Create equal sized empty matrix for calculations
  cov_matrix <- matrix(0, nrow = length(st_dev), ncol = length(st_dev))
  
  # Set diagonal
  diag(cov_matrix) <- st_dev**2
  
  # Calculate error covariances
  for(i in 1:nrow(cov_matrix)) {
    for(j in 1:ncol(cov_matrix)) {
      # get Prior
      P <- exp( (-abs(log10(pressure[i]) - log10(pressure[j]) ) ) / fac  )
      
      # calculate covariance in cell i,j
      cov_matrix[i,j] <- ( sqrt(cov_matrix[i,i]) * sqrt(cov_matrix[j,j]) ) * P
      
    }
  }
  
  # Return full matrix
  return(cov_matrix)
}



# Function: Write BMatrix -------------------------------------------------

#' This function writes the BMatrix in the correct format to a filepath specified
#' 
#' 
#' 
#' 

write_BMatrix_to_file <- function(BMatrix, BM_path) {
  
  # Check if file already exists
  # Remove file if exists
  if(file.exists(BM_path)){
    prompt <- readline("BMatrix already exists! Delete file and create new? [y/n]: ")
    if(prompt == "y"){
      file.remove(BM_path)
    }
    else if(prompt =="n"){
      stop("Function Aborted")
    }
  }
  
  # Header information
  lin1 <- c("Peter Weston's NH 54 RTcoef level B matrix from L70 Covstats for Sea. Ozone error covariance added from McPeters-Labow climatology")
  lin2 <- c("1-54=temp (K), 55-83=lnq (g/kg) (bottom 29 levels only), 84=Tsurf, 85=lnq surf, 86=Tskin, 87-140= Ozone(kg/kg)")
  lin3 <- c("140") # The total amount of nrows/ncolumns for this matrix
  
  # Append to eachother
  header <- append(lin1, c(lin2, lin3))
  
  # Set correct notation for the BMatrix, with the right significance
  BMatrix_Correct_Notation <- toupper(sprintf("%16.8e", BMatrix))
  
  # Write to file
  # Do 4 times, as with the example matrix, not exactly sure why
  for(i in 1:4){
    cat(header, apply(matrix(BMatrix_Correct_Notation, ncol=5, byrow=TRUE), 1, paste, collapse="") ,
        sep="\n", file = BM_path, append = TRUE)
  }
  
  
}





# ggplot2 Functions -------------------------------------------------------

#' This function merges the plotting of the several input files to not clutter
#' the whole file with large ggplotting tools. It should adapt to the input data
#' and create the correct graphs too.
#' 

ggplot_Input_Files <- function(input_data, type) {
  
  # Plot BMatrix
  if(type == "BMatrix"){ 
    # Melt values
    data_melted <- melt(input_data)
    
    # Find breakpoints for legend
    breaks <- quantile(data_melted$value, (0:6)/6)
    
    # Plot
    ggplot(data_melted, 
           aes(x= Var1, 
               y = Var2, 
               fill = as.numeric(value))) +
      geom_tile(color = "gray90") +
      scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                           values = scales::rescale(breaks), 
                           name="Covariance error (kg/kg)") +
      labs(title = paste0("For latitudinal region: ", 
                          latidx * 10, 
                          "-", 
                          (latidx*10)+10, 
                          "N"))  +
      xlab("Z Level") +
      ylab("Z level") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.key.height=unit(2,"cm"),
            legend.position = "right") +
      coord_fixed(expand = FALSE)
    
  }
  
  
  #' Now for IASI spectra
  else if(type == "Spectra") {
   # TODO 
  }
  
}



# Function: Inter/Extrapolate values --------------------------------------



# Spectra loop ------------------------------------------------------------
#'THIS IS A BACKUP OF A PREVIOUSLY MADE FUNCTION
#'MOSTLY TO TEST OPTIMIZING WRITING TO FILE AND BETTER FORMATTING
#'USE THE FUNCTION ABOVE FOR ACTUAL USAGE

write_1DVar_input_Backup <- function(basefile, spect, lat, lon, solzen, satzen, filen) { 
  
  # Format the lines as shown above
  lin1 <- paste0("Obs ID:              1 Obs Type:          3 Satellite ID:     4")
  lin2 <- paste0("Latitude: ", sprintf("%.3f", lat), " Longitude: ", sprintf("%.3f", lon), " Elevation: 0.0")
  lin3 <- paste0("Surface type:   1", " Sat Zen Angle:   ", sprintf("%.3f", satzen), " Solar Zen. Ang.: ", sprintf("%.3f", satzen))
  
  # now append to basefile
  basef2 <- append(basef, c(lin1, lin2, lin3, "Radiances:"))
  
  # Do next line if width of each element should be the same
  spect_format <- formatC(as.numeric(spect), format = 'f', flag=' ', digits = 3, width = 12)
  
  # Write to file
  cat(basef2, paste0(tapply(spect_format, ceiling(seq_along(spect)/5), 
                            paste0, collapse = '     '), collapse = '\n'), file = filen, sep = "\n")
  
  
  
  # Merge basefile, lines and the formatted spectra correctly
  full_text <- append(basef, cat(paste0(tapply(spect_format, ceiling(seq_along(spect)/5), paste0, collapse = "      "), collapse = '\n')))
  
  # instead of cat use writeChar, tested system.time about factor 10 faster
  fileconn <- file(filen)
  writeChar(object = full_text, con = fileconn, nchar(full_text, type = "chars"))
  close(fileconn)
  
  
  # Test with fwrite too
  fwrite(full_text, file = filen, append = TRUE)
  
  
}


# Function: covariance of profile -----------------------------------------

#' This function calaculates the covariance matrix based on the error covariance
#' from the MCPeters labow climatology. Following the inverse to st.dev to covariance
#' method. 

cov_profile_ratio <- function(pressure, st_dev, fac = 1) {
  
  # Create equal sized empty matrix for calculations
  cov_matrix <- matrix(0, nrow = length(st_dev), ncol = length(st_dev))
  
  # Set diagonal
  diag(cov_matrix) <- st_dev**2
  
  # Calculate error covariances
  for(i in 1:nrow(cov_matrix)) {
    for(j in 1:ncol(cov_matrix)) {
      # get Prior
      #P <- exp( (-abs(log10(pressure[i]) - log10(pressure[j]) ) ) / fac  )
      
      # calculate covariance in cell i,j
      cov_matrix[i,j] <- ( sqrt(cov_matrix[i,i]) * sqrt(cov_matrix[j,j]) )
      
    }
  }
  
  # Return full matrix
  return(cov_matrix)
}





# ERA5 Background File ----------------------------------------------------

#' This section creates background files based on ERA5 reanalysis data
#' 

# ERA5_filen <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2d_ERA5/20110801/AUS_Sea/atmospheric_data.txt"

ERA5_Background <- function(ERA5_filen, filen) {
  
  # Read conversion coefficients, a & b
  conv_filen <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2e_NPWSaf_Profiles/ConversionFactors.txt"
  conv_factors <- read.table(conv_filen, header = T)
  
  # Read single line variables
  surf_p <- as.numeric(read.table(ERA5_filen, skip = 12, nrows = 1)[2]) # Pressure in Pa
  s_2t <- as.numeric(read.table(ERA5_filen, skip = 14, nrows = 1)[2]) # 2 Metre temperature in K
  s_2dt <- as.numeric(read.table(ERA5_filen, skip = 16, nrows = 1)[2]) # 2 Metre dewpoint temperature in K
  s_skint <- as.numeric(read.table(ERA5_filen, skip = 18, nrows = 1)[2]) # Skin Temperature in K
  cc <- as.numeric(read.table(ERA5_filen, skip = 20, nrows = 1)[2] * 100) # Cloud Cover in %

  # read in 137 level variables
  t <- read.table(ERA5_filen, skip = 22, nrows = 137, header = F) # temperature
  q <- read.table(ERA5_filen, skip = 160, nrows = 137, header = F) # water vapor
  o <- read.table(ERA5_filen, skip = 298, nrows = 137, header = F) # ozone
  
  # Get pressure levels
  p <- conv_factors$a + conv_factors$b * surf_p # in Pa
  
  # Add small value to ensure TOA 1= o
  p[1] <- p[1] + 0.005
 
  # Interpolate to 54 layers
  t_54L <- approx(t, n = 54)$y
  q_54L <- approx(q, n = 54)$y
  o_54L <- approx(o, n = 54)$y
  p_54L <- format( (approx(p, n = 54)$y / 100), digits = 4, scientific = F) # in hPa
  
  # combine into one dataframe
  columns <- data.frame(p_54L, 
                        t_54L, 
                        q_54L, 
                        o_54L)
  
  
  # Interpolate to original pressure levels ###########################
  
  df1 <- transform(columns, p_54L = as.numeric(p_54L), 
                   t_54L = as.numeric(t_54L),
                   o_54L = as.numeric(o_54L),
                   q_54L = as.numeric(q_54L))
  
  df2 <- transform(backgr_df, p_54L = as.numeric(V1))
  
  # merge values
  # df1 = columns, df2 = backgr_df
  df3 <- df1 %>% 
    complete(p_54L = union(as.numeric(p_54L), df2$p_54L)) %>%
    mutate(o_54L = na.approx(o_54L, maxgap = Inf, rule = 2),
           q_54L = na.approx(q_54L, maxgap = Inf, rule = 2),
           t_54L = na.approx(t_54L, maxgap = Inf, rule = 2))
  
  df4 <- merge(df2, df3, by = "p_54L")[, c(1, 6, 7, 8)]
  
  #####################################################################
  
  # Invert
  columns <- data.frame(map_df(df4, rev))
  
  # Create headers and lines
  header <- "Background file generated from
ERA reanalysis data
54 fixed layers
pressure in hPa
O3 in kg/kg




x------------------- End of Header------------------------------------x
No. Background Profiles:            1
No. of Levels/Profile:             54
Unit for q:                         2  (1=ppmv, 2=kg/kg, 3=RH)
x---------------------------------------------------------------------x
Profile #    1 Follows
x---------------------------------------------------------------------x"

  # Columns
  columns_formatted <- kable(columns, digits = 12, format = "simple", align = "r")[c(-1,-2)]
  
  # Surface information
  surfT <- paste0("Surface Temperature (K):     ", s_2t)
  surfH <- paste0("Surface Humidity (kg/kg):     ", sprintf("%.6E", q[137, 2]))
  skinT <- paste0("Skin Temperature (K):     ", s_skint)
  surfP <- paste0("Surface Pressure (hPa):     ", surf_p / 100)
  uwind <- paste0("10m U-Wind (m/s):     0.00000")
  
  
  # Write to file
  cat(header, 
      columns_formatted, 
      surfT,
      surfH,
      skinT,
      surfP,
      uwind, uwind, 
      sep = '\n', 
      file = filen)
  
  
}




# Rad to BT convert -------------------------------------------------------

#' This section converts the full IASI spectrum radiance values to Brightness
#' Temperatures by means of Planck function assuming perfect black body. 
#' Input spectra in radiance,s with added wavelength too
#' 
#' Output : BT in K
#' 

RAD_to_BT <- function(spectra) {
  
  # Constants from 
  c1 <- 1.191042E-5
  c2 <- 1.4387752
  spectra$radiance <- spectra$radiance * 1E5 # to correct units (mW/m2-sr-m-1)

  # From ncc.nesdis.noaa.gov/dta/planck
  spectra$BT <- (c2 * spectra$wavenumber) /
                (log(c1 * spectra$wavenumber^3 / spectra$radiance + 1 ))
  
  # Return spectra with all information
  return(spectra)
}




# Rmatrix add error -------------------------------------------------------

#' read in RMatrix: sigma (by sqrt() of values). Add 1% of observation value
#' at that point. Then square and rewrite RMatrix.
#' 2nd: add error value of IASI observation in NETCDF4 file.
#' 

add_noise_to_Rmatrix <- function(output_path, spectra) {
  
  # Set path to base Rmatrix in radiances
  path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR\NWPSAF_1DVar_1.3/IASI_COEFFS_DIR_RAD/Rmatrix_8461_instnoise_PCRTTOV.out"
  
  # Read in values of original RMatrix
  rmatrix_values <- scan(path, skip = 849)
  
  # add 1% noise
  rmatrix_values_noise <- (sqrt(rmatrix_values) + (spectra$radiance / 100)) ^ 2
  
  # Read base information
  base <- readLines(path, n = 849)
  
  # Format values correctly
  rmatrix_values_formatted <- formatC(as.numeric(rmatrix_values_noise), format = 'f', flag=' ', digits = 6, width = 13)
  
  # Write to file
  cat(base, paste0(tapply(rmatrix_values_formatted, ceiling(seq_along(rmatrix_values_formatted)/6), 
                            paste0, collapse = ''), collapse = '\n'), file = output_path, sep = "\n")
  
}


