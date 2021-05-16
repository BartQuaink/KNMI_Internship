# set up required libs
libs <- c("ggplot2", "lubridate", "readr")
pacman::p_load(char = libs)



# Pathing -----------------------------------------------------------------

# Folder
path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/E3b_TestOutputs/0204_ControldataComparison/"

# Read in all different profiles created using the different options
profile_all <- read.table(paste0(path, "IASI_54L_All/Retrieved_profiles.dat"), header = F, skip = 3, nrows = 54)
profile_neither <- read.table(paste0(path, "IASI_54L_Neither/Retrieved_profiles.dat"), header = F, skip = 3, nrows = 54)
profile_noSurfq <- read.table(paste0(path, "IASI_54L_NoSurfq/Retrieved_profiles.dat"), header = F, skip = 3, nrows = 54)
profile_noSurfT <- read.table(paste0(path, "IASI_54L_NoSurfT/Retrieved_profiles.dat"), header = F, skip = 3, nrows = 54)

# Set column headers
headers <- c("Pressure (hPa)", "T (K)", "q (ppmv)", "Ozone (ppmv)", "T (K) - Background", "q (ppmv) - Background", "Ozone (ppmv) - Background")

# Name all headers
colnames(profile_all) <- headers
colnames(profile_neither) <- headers
colnames(profile_noSurfq) <- headers
colnames(profile_noSurfT) <- headers


# Comparison --------------------------------------------------------------

# Ratio
ratio <- (profile_all$`T (K)` - profile_neither$`T (K)`) / profile_neither$`T (K)`






# Plotting ----------------------------------------------------------------
ggplot() +
  geom_path(data = profile_all, aes(x = `Ozone (ppmv)`, y = `Pressure (hPa)`, colour = "red")) +
  geom_path(data = profile_neither, aes(x = `Ozone (ppmv)`, y = `Pressure (hPa)`, colour = "blue")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  xlab("Ozone Mixing Ratio (ppmv)")+
  theme_bw()

