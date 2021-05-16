# install required libs
libs <- c("raster", "ggplot2", "rgdal", "lubridate", "readr", "sf", "dplyr", "tidyverse", "ncdf4", "BBmisc")
pacman::p_load(char = libs)


#' This script runs an exmaple profile all the way to setting up the correct background files
#' and BMatrices based on the functions created. 
#' 


# Source functions --------------------------------------------------------
source("C:/Users/bartq/Documents/WUR/Internship/F_ScriptsandCode/F3_R/F3b_SupportFunctions/1DVar_IASI_Input_Preprocessing.R")

# - ----------------------------------------------------------------------

#' Main section: reading in an example spectra above the USA [4,4,640]
## IASI spectra
retrieval_path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E4_FullRundownOutputs/v5_OzoneHole/" # put them all in this one for an example

# Required basefile for the data
basefile <- "E:/_InternshipData/1DVar_InputFiles/IASI_basefile.txt"

## BMatrix pathing and files
BM_path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/E3a_BMatrices/BMatrix_54L" # output
ML_path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_stats.dat" # MLL input error covariance



## Locations 2011-8-01:
#nc_path <- "E:/_InternshipData/IASI/Level1c/2011/August/20110801004800_24810.nc"

#idx <- c(4, 4, 640) # 38.83, -86.45 (located above Indiana) cc: 0%
#idx <- c(4, 20, 120) # 25.06, 128.10 (located above PH Sea) cc: 11%
#idx <- c(4, 20, 200) # above sea between Indonesia and Au, 12°34'12.0"S 119°36'00.0"E
#idx <- c(1, 30, 500) # -52.49, -22.32 (located above Parana, Brasil above SAA)

## Locations 2011-09-25:
nc_path <- "E:/_InternshipData/IASI/Level1c/2011/September/W_XX-EUMETSAT-Darmstadt,HYPERSPECT+SOUNDING,METOPA+IASI_C_EUMP_20110925140708_25599_eps_o_l1c.nc"

# Location above ozone hole, following TEMIS
idx <- c(4, 5, 88) # lat/lon : -56.65572/-75.58702



# Get spectra -------------------------------------------------------------

# Read in dataset
nc <- nc_open(nc_path)

## The following are the required parameters to create everything. The full spectra
## is needed, location information, time of retrieval and the zenith angles for now
## are the required input files

# Read all 8700 values for index 4,4,650. Start at first radiance observation
spect <- ncvar_get(nc, "gs_1c_spect", start = c(1, idx), count = c(8700, 1, 1, 1))

# and the lat/lon/time information too
lat <- ncvar_get(nc, "lat", start = idx, count = c(1, 1, 1))
lon <- ncvar_get(nc, "lon", start = idx, count = c(1, 1, 1))
date_s <- ncvar_get(nc, "measurement_date", start = idx[2:3], count = c(1, 1))
date <- as.POSIXct(date_s, origin="2000-01-01", tz="UTC") # convert to normal YYYY/MM/DD S format

# Plus zenith angles
satzen_angle<- ncvar_get(nc, "pixel_zenith_angle", start = idx, count = c(1, 1, 1))
solzen_angle <- ncvar_get(nc, "pixel_solar_zenith_angle", start = idx, count = c(1, 1, 1))

# Cloud fraction information
cloud_fraction <- ncvar_get(nc, "geum_avhrr_1b_cloud_fraction", start = idx, count = c(1, 1, 1))

# Get scaling factors
scaling <- ncvar_get(nc, "i_def_scale_sond_scale_factor")



# Spectra conversion and plotting -----------------------------------------

# Get wavenumbers
wavenumbers <- calculate_Wavenumbers(ncvar_get(nc, 'i_def_spect_dwn1b')[1], 
                                     ncvar_get(nc, 'i_def_ns_first_1b')[1])

# Get scaled spectra following the level 1C user guide
scaling_factors <- IASI_scaling_factors(nc)

# get spectra scaled
spectra_scaled <- IASI_Spectra_Scaled(spect, scaling_factors) # in w/m2/sr/m-1

# get spectra scaled for 1DVar (mW/m2-sr-m-1)
spect_scaled <- spectra_scaled * 1E5

# merge for plotting
spectra_scaled_wn <- data.frame(wavenumber = wavenumbers, radiance = spectra_scaled)

# Convert to Brightness temperatures
spectra_BT <- RAD_to_BT(spectra_scaled_wn) # results in wavenumbers in cm-1, rad in mw/m2/sr/cm-1, BT in K

## Again as above, the following line combines all the abovementioned steps and outputs for the
## whole spectra and netcdf file all files into a folder
#spectra_Prep(nc, basefile, retrieval_path)


# Create BMatrix ----------------------------------------------------------

#' Example: This is also put in the complete file, but required to show the
#' cov_profile step when creating an error covariance based on the MLL cli-
#' matology. 

# Get index latitude in 10-latitudinal zone
latidx <- round(lat, digits = -1) / 10

# Get dataframe of full climatology, accounting for negative values beginning by southern hemisphere
ppmv_stats <- read.table(ML_path, skip = ((latidx + 9) * 70) + 1, nrow = 66, header = TRUE)

# Get correct column of the month
table <- data.frame(zlevel = ppmv_stats$Z.level, ppmv_stats = ppmv_stats[, month(date) + 1])

# convert zstar altitudes to pressure, following supplied formula in README.txt in ML climatology
pressure_conv <- zstar_press(table$zlevel) # also a required function

# Interpolate to 54 layers
stats_54L <- data.frame(
  pressure = approx(pressure_conv, n=54)$y, 
  values = approx(table$ppmv_stats, n=54)$y
)

# Then run the function and retrieve the error covariance profile in pressure
cov_matrix <- cov_profile(stats_54L$pressure, stats_54L$values, fac = 0.5)

# Test if solvable, FALSE if there is no error so there is no issue
is.error(solve(cov_matrix))

# Plot example
ggplot_Input_Files(cov_matrix, type = "BMatrix")


#' Abovementioned steps are put into one file to recieve the total BMatrix
#' as shown here below.

# This function creates the full BMatrix
BMatrix <- create_Bmatrix(lat, date, fac = 0.5, c_dampen = 0)

# Add different damping values (c)
BMatrix_8 <- create_Bmatrix(lat, date, c = 1E-8)
BMatrix_7 <- create_Bmatrix(lat, date, c = 1E-7)
BMatrix_6 <- create_Bmatrix(lat, date, c = 1E-6)

# Create plot
ggplot_Input_Files(BMatrix, type = "BMatrix")



# Create background -------------------------------------------------------

create_BACKGROUND(lat, date, filen = paste0(retrieval_path, "BACKGROUND_54L.dat"))
ERA5_Background(ERA5_filen = "C:/Users/bartq/Documents/WUR/Internship/E_Data/E4_FullRundownOutputs/v5_OzoneHole/input/atmospheric_data.txt", 
                filen = paste0(retrieval_path, "ERA5_BACKGROUND_54L.dat"))


# Write all to file -------------------------------------------------------

# Write to file
write_BMatrix_to_file(BMatrix, paste0(retrieval_path, "BMatrix_54L"))

# Write spectra to file
write_1DVar_input(basefile, spect_scaled, lat, lon, date, satzen_angle, solzen_angle, filen = paste0(retrieval_path, "Observation_IASI_54L.dat"))

# Write BT to file
write_1DVar_input(basefile, spectra_BT$BT, lat, lon, date, satzen_angle, solzen_angle, filen = paste0(retrieval_path, "Observation_IASI_54L_BT.dat"), BT = TRUE)

# This function loops all the abovementioned steps into one:
# 

