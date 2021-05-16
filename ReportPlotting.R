
# RMatrix -----------------------------------------------------------------

require(ggplot2)
require(scales)
require(reshape2)

#'This script reads in the RMatrix and plots the corresponding error covariances
#'


path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/IASI_COEFFS_DIR_RAD/Rmatrix_8461_instnoise_PCRTTOV.out"

# Read in
rmatrix <- data.frame(channel = scan(path, skip = 2, nmax = 8461),
                      error = scan(path, skip = 849, nmax = 8461))

# Plot
ggplot(rmatrix, aes(x = channel, y = error)) +
  geom_path() +
  theme_bw() +
  labs(title = "R-Matrix (Observation Error)",
       x ="Channel number",
       y = "Observation error, radiometric noise in W/m-2/sr-1/m-1",
       caption = "R-Matrix supplied with 1D-Var")





# Minimisation ------------------------------------------------------------

#' Plot how miminisation over time occurs
#' Starts at iteration 0, with the Cost Function value, not gradient added
#' 

# Iterations
iter <- 0:10

min_normal <- data.frame(iter = iter, cf_normal = c(133228.05, 2924.31, 1752.60, 1419.17, 1394.22, 1392.48, NA, NA, NA, NA, NA))
min_IASI <- data.frame(iter = iter, cf_IASI = c(481242.85, 19350.67, 24292.15, 16411.33, 11497.50, 11345.02, 12239.72, 14166.27, 14257.93, NA, NA))
min_BMatrix <- data.frame(iter = iter, cf_BMatrix = c(133228.05, 15114.07, 3243.19, 13230.43, 2571.36, 17590.98, 2466.54, 16661.00, 2879.13, 13653.94, 2752.98))
min_Personal <- data.frame(iter = iter, cf_Personal = c(472837.89, 61767.47, 41346.77, 89332.09, 44477.08, 93719.20, 40861.04, 91786.60, 74126.14, 78750.11, 15726.75))

# Min personal 2
min_P2 <- data.frame(iter = iter, cf_P = c(6289.30, 101674.87, 15922.63, 84726.652, 12164.219, 89947.18, 7037.59, 85479.967, 7126.012, 86486.782, 7395.858))

# hardcode in the values for now
minimisation_df <- Reduce(function(x, y) merge(x, y, all = T), list(min_normal, min_IASI, min_BMatrix, min_Persona))

# 2nd method
minimisation_df <- Reduce(function(x, y) merge(x, y, all = T), list(min_normal, min_P2))

minimisation_melted <- melt(minimisation_df, id = "iter")

# Plot
ggplot(minimisation_melted, aes(x = iter, y = value, colour = variable)) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_y_log10(labels = scales::label_comma()) +
  annotation_logticks(sides = "l") +
  labs(title = "Comparison between minimisation cost function values",
       x = "Iteration",
       y = "Cost Function value",
       color = "Legend") +
  scale_color_manual(values = c("red", "green", "blue", "orange"),
                     labels = c("Sample Data", "Personal Input", "BMatrix", "Personal Data")) +
  theme_bw()


# And without the first 0th iteration
minimisation_melted_1 <- subset(minimisation_melted, iter!=0)

# Plot
ggplot(minimisation_melted_1, aes(x = iter, y = value, colour = variable)) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_y_log10(labels = scales::label_comma()) +
  annotation_logticks(sides = "l") +
  labs(title = "Comparison between minimisation cost function values",
       subtitle = "Initial state filtered out",
       x = "Iteration",
       y = "Cost Function value",
       color = "Legend") +
  scale_color_manual(values = c("red", "green", "blue", "orange"),
                     labels = c("Sample Data", "IASI Spectra", "BMatrix", "Personal Data")) +
  theme_bw()




# IASI spectra ------------------------------------------------------------

# plot
ggplot(data = spectra_scaled_wn, aes(x = wavenumber, y = radiance)) +
  geom_path(colour ="blue") +
  labs(title = "Example of an IASI spectra",
       subtitle = paste0("Location: 38.78, -86.45. Observation time: ", date),
       x = "Wavenumber cm-1",
       y = "Radiance in W/m2/sr/m-1 ") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_bw()




# Climatology -------------------------------------------------------------

# melt values
clim_df_melted <- melt(clim_df, id = "Z.level")

ggplot(data = clim_df_melted, aes(y = Z.level, x = value, colour = variable)) +
  geom_path() +
  labs(title = "Ozone profiles from MLL climatology for 30N-40N",
       x = "Z level",
       y = "Ozone values in ppmv") +
  guides(colour=guide_legend(title="Month")) +
  theme_bw()



# Background profiles -----------------------------------------------------

# copy df
df1 <- columns
df2 <- transform(backgr_df, p_54L = as.numeric(V1))

# rename
colnames(df2) <- c("p_54L")

df1 <- transform(df1, p_54L = as.numeric(p_54L), 
                 t_54L = as.numeric(t_54L),
                 o_54L = as.numeric(o_54L),
                 q_54L = as.numeric(q_54L))

# merge values
# df1 = columns, df2 = backgr_df
df3 <- df1 %>% 
  complete(p_54L = union(as.numeric(p_54L), df2$p_54L)) %>%
  mutate(o_54L = na.approx(o_54L, maxgap = Inf, rule = 2),
         q_54L = na.approx(q_54L, maxgap = Inf, rule = 2),
         t_54L = na.approx(t_54L, maxgap = Inf, rule = 2))

# Get values of 
df4 <- merge(df2, df3, by = "p_54L")[, c(1, 6, 7, 8)]


# GGPLOT
df4_melted <- melt(df4, id = "p_54L")

## Differences after interpolation
ggplot() +
  geom_path(data = df1, aes(x = o_54L, y = p_54L), colour = "red") +
  geom_path(data = df3, aes(x = o_54L, y = p_54L), colour = "blue") +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
  annotation_logticks(sides = "l") +
  labs(title = "Comparison before/after interpolation",
       subtitle = "blue = after, red = before",
       x= "Ozone (kg/kg)",
       y = "Pressure (hPa)") +
  theme_bw()


# Get reverse y log scale function
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}   
  
  
  

# Output files ------------------------------------------------------------

#' Read all the output files as outputted in a standard folder, just link folder and the 
#' graphs will do the rest
#' 
#'   

output_dir <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E4_FullRundownOutputs/v4_ERA5_Background/Outputs/"  
  
# Read files
files <- list.files(output_dir, recursive = T, full.names = T)


########################## PROFILES.DAT ##################################
# Plot 1st one
retrieved_profiles <- read.table(files[10], skip = 3, nrows = 54, header = F)
colnames(retrieved_profiles) <- c("Pressure (hPa)", "T (K) (r)", "q (kg/kg) (r)", "Ozone (r)", "T (K) (b)", "q (kg/kg) (b)", "Ozone (b)")


# Plot Temperature
ggplot(data = retrieved_profiles) +
  geom_path(aes(x = `T (K) (r)`, y = `Pressure (hPa)`), colour = "red") +
  geom_path(aes(x = `T (K) (b)`, y = `Pressure (hPa)`), colour = "blue")  +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
  annotation_logticks(sides = "l") +
  labs(title = "Comparison retrieved/background for Temperature",
       subtitle = "red: retrieved, blue: background. Surf T (r): 300.629, Surf T (b): 298.835",
       x= "T (K)",
       y = "Pressure (hPa)") +
  theme_bw()


# Plot q
ggplot(data = retrieved_profiles) +
  geom_path(aes(x = `q (kg/kg) (r)`, y = `Pressure (hPa)`), colour = "red") +
  geom_path(aes(x = `q (kg/kg) (b)`, y = `Pressure (hPa)`), colour = "blue")  +
  scale_x_continuous(labels = scales::label_comma(), trans = "log10") +
  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
  annotation_logticks(sides = "lb") +
  labs(title = "Comparison retrieved/background for Water Vapour",
       subtitle = "red: retrieved, blue: background. Surf q (r): 0.011, Surf q (b): 0.014",
       x= "q (kg/kg)",
       y = "Pressure (hPa)") +
  theme_bw()


# Plot Ozone
ggplot(data = retrieved_profiles) +
  geom_path(aes(x = `Ozone (r)`, y = `Pressure (hPa)`), colour = "red") +
  geom_path(aes(x = `Ozone (b)`, y = `Pressure (hPa)`), colour = "blue")  +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
  annotation_logticks(sides = "lb") +
  labs(title = "Comparison retrieved/background for Ozone",
       subtitle = "red: retrieved, blue: background.",
       x= "q (kg/kg)",
       y = "Pressure (hPa)") +
  theme_bw()



########################## MINIMISATION ##################################



########################## JACOBIANS #####################################
#' BgJacobian contains the jacobian corresponding to the profile elements
#' selected in the BMatrix calclated from background profile, and retjacob-
#' ians is the output at the end of the minimisation. In format 10E12.4 in
#' size (NumChans, NumElements). NumChans = number of channels used in the
#' 1Dvar, and NumElements are elements in BMatrix used (140) 

bgJacobians_df <- read.table(files[4], skip = 3, header = F)
retJacobians_df <- read.table(files[8], skip = 3, header = F)

# Reformat to correct size
bgjacobians <- matrix(unlist(t(bgJacobians_df)), byrow = F, 314, 140) # 314 normal size, 298 excluding SWIR
retjacobians <- matrix(unlist(t(retJacobians_df)), byrow = F, 314, 140)

## Plot
# Melt values
bgjacobians_melted <- melt(bgjacobians)
retjacobians_melted <- melt(retjacobians)

# Find breakpoints for legend
breaks_bg <- quantile(bgjacobians_melted$value, (0:10)/10)
breaks_ret <- quantile(retjacobians_melted$value, (0:10)/10)

# Plot background jacobians
ggplot(bgjacobians_melted, 
       aes(x= Var1, 
           y = Var2, 
           fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                       values = scales::rescale(breaks_bg), 
                       name="") +
  labs(title = paste0("Background Jacobians for observation: ")) +
  xlab("Number of Channels") +
  ylab("Number of Elements") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height=unit(2,"cm"),
        legend.position = "right") +
  coord_fixed(expand = FALSE)

# Plot retrieved jacobians
ggplot(retjacobians_melted, 
       aes(x= Var1, 
           y = Var2, 
           fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                       values = scales::rescale(breaks_ret), 
                       name="") +
  labs(title = paste0("Retrieved Jacobians for observation: ")) +
  xlab("Number of Channels") +
  ylab("Number of Elements") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height=unit(2,"cm"),
        legend.position = "right") +
  coord_fixed(expand = FALSE)



## Some channels as a lineplot





########################### BRIGHTNESS T #################################
BT <- read.table(files[9], skip = 2, header = T)

# Create new diff dataframe
BT_diff <- data.frame(Channel = BT$Channel, `B-0` = BT$Background - BT$Observed, `B-T` = BT$Background - BT$Retrieved)

# Melt for plotting
BT_Melt <- melt(BT, id = "Channel")
BT_diff_melt <- melt(BT_diff, id = "Channel")

# Plot Values
ggplot(data = BT_Melt, aes(x = value, y = Channel, colour = variable)) +
  geom_path()+
  theme_bw()


ggplot(data = BT_diff_melt, aes(x = value, y = Channel, colour = variable)) +
  geom_path()+
  theme_bw()






########################## AVERAGING KERNEL ##############################
# Written out in format 10E12.4, size (NumElements, NumElements)
# NumElements is determined by the BMatrix, for our 54L case this is 140*140

AVK <- read.table(files[3], skip = 3, nrows = 1960, header = F)

# Reformat to correct size
AVK <- matrix(unlist(t(AVK)), byrow = T, 140, 140)

# Test plot
matplot(t(AVK), 
        type = "l",
        main="Averaging Kernel",
        ylab = "dy,dp",
        xlab = "Element number in BMatrix,
        ")

# Reshape to long
AVK_melt <- melt(AVK) # Var2 is row number

# plot complete
ggplot(AVK_melt, aes(x = value, y = Var1)) +
  geom_path(aes(color = Var2, group = Var2)) +
  scale_colour_gradient(low = "blue", high = "green") +
  labs(title = "Averaging Kernel",
       x = "Value",
       y = "Element number",
       color = "element nr") +
  theme_bw()


# Plot ozone only
AVK_Ozone <- AVK_melt[AVK_melt[,1] >= 87 & AVK_melt[,2] >= 87, ]

ggplot(AVK_Ozone, aes(x = value, y = Var1)) +
  geom_path(aes(color = Var2, group = Var2)) +
  scale_colour_gradient(low = "blue", high = "green") +
  theme_bw()





# RMatrix + Observation Error ---------------------------------------------
# Diag Rmatrix is in ^2, so square for true covariance value

path_RMatrix_RAD <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/IASI_COEFFS_DIR_RAD/Rmatrix_8461_instnoise_PCRTTOV.out"
path_RMatrix_BT <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/IASI_COEFFS_DIR_BT/Rmatrix_8461_instnoise_PCRTTOV.out"

# Read all elements one by one
rmatrix_RAD <- sqrt(scan(path_RMatrix_RAD, skip = 849))
rmatrix_BT <- sqrt(scan(path_RMatrix_BT, skip = 849))

# Plot Rad
ggplot(spectra_BT, aes(x = wavenumber, y = radiance)) +
  geom_ribbon(aes(ymin = radiance - rmatrix_RAD, 
                  ymax = radiance + rmatrix_RAD), 
              linetype = 2, alpha = 0.1, colour = "blue") +
  geom_line() +
  xlim(c(1000,1070)) +
  theme_bw()

# Plot BT
ggplot(spectra_BT, aes(x = wavenumber, y = BT)) +
  geom_ribbon(aes(ymin = BT - rmatrix_BT, 
                  ymax = BT + rmatrix_BT), 
              linetype = 2, alpha = 0.1, colour = "blue") +
  geom_line() +
  xlim(c(1000,1070)) +
  theme_bw()


# Relative difference
rmatrix_ret_BT <- rmatrix_BT / spectra_BT$BT
rmatrix_ret_RAD <- rmatrix_RAD / spectra_BT$radiance

# plot
ggplot()+
  geom_path(rmatrix_ret_RAD)






# Jacobians as linegraph --------------------------------------------------

ggplot(data = bgjacobians_melted, aes(x = Var1, y = as.numeric(value), color = Var2)) +
  geom_path() +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_colour_viridis() +
  labs(title = "Background Jacobians per element",
       x = "channel selection number",
       y = "Value",
       legend = "BMatrix element number") +
  labs(color="BMatrix element number") +
  xlim(c(146,162))+
  theme_bw()


