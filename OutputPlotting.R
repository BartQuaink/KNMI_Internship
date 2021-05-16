
# Script Information ------------------------------------------------------

#'
#'
#'
#'
#'

# set required libraries
require(ggplot2)
require(scales)
require(reshape2)
require(viridis)

#' Quick helper function: Get reverse y log scale function
#' Makes it a single function so we can reverse the hPa going from the surface
#' to the top of the atmosphere while also making a log scale of the values

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
args             <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2){
  stop("Two arguments are required: input_directory--path output_directory--path", call. = F)
}

input_directory  <- args[1]
output_directory <- args[2]


# Read in files -----------------------------------------------------------

files <- list.files(input_directory, recursive = T, full.names = T)

# Read in all matrices/profiles
# A-Matrix.out
A_Matrix_values <- scan(files[1], what = "numeric") # Read all values
A_Matrix        <- suppressWarnings(na.omit(as.numeric(A_Matrix_values))) # Remove characters and text, suppress warnings
A_Matrix        <- A_Matrix[!A_Matrix == 1 & !A_Matrix == 140] # Remove iteration number and matrix size
A_Matrix_m      <- matrix(A_Matrix[c(19601:length(A_Matrix))], nrow = 140, ncol = 140) # Take last iteration

data_melted_A   <- melt(A_Matrix_m) # melt values
breaks_A        <- quantile(data_melted_A$value, (0:10)/10) # find breakpoints for plotting


# AM-Matrix.out
Am_Matrix_values <- scan(files[2], what = "numeric") # Read all values
Am_Matrix        <- suppressWarnings(na.omit(as.numeric(Am_Matrix_values))) # Remove characters and text, suppress warnings
Am_Matrix        <- Am_Matrix[!Am_Matrix == 1 & !Am_Matrix == 140] # Remove iteration number and matrix size
Am_Matrix_m      <- matrix(Am_Matrix[c(1:length(Am_Matrix))], nrow = 140, ncol = 140) # make matrix

data_melted_Am   <- melt(Am_Matrix_m) # melt values
breaks_Am        <- quantile(data_melted_Am$value, (0:10)/10) # find breakpoints for plotting


# AveragingKernel.out
AVK       <- read.table(files[3], skip = 3, nrows = 1960, header = F)
AVK       <- matrix(unlist(t(AVK)), byrow = T, 140, 140) # Reformat to correct size

AVK_melt  <- melt(AVK) # melt values, Var2 is row number
AVK_Ozone <- AVK_melt[AVK_melt[,1] >= 87 & AVK_melt[,2] >= 87, ] # take ozone only


# BgJacobian.out
bgJacobians_df     <- read.table(files[4], skip = 3, header = F)
bgjacobians        <- matrix(unlist(t(bgJacobians_df)), byrow = F, 298, 140) # 314 normal size, 298 excluding SWIR

bgjacobians_melted <- melt(bgjacobians) # melt values
breaks_bg          <- quantile(bgjacobians_melted$value, (0:10)/10) # find breakpoints for plotting


# Minimisation.log
minimisation        <- readLines(files[5])
# some data extraction: first extract only rows with cost function lines and values
extract_values      <- lapply(minimisation, function(x) {if (grepl(" Cost Function =   ", x)) unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))})
minimisation_empty  <- which(sapply(extract_values, is.null)) # get indexes of empty return rows
minimisation        <- minimisation[-c(minimisation_empty)] # remove from initial dataframe
minimisation_values <- data.frame(iteration = 0:(length(minimisation)-1), # add iteration numbers from previous line
                                  cost_function = read.table(text = minimisation, fill = TRUE)[[4]]) # get just the cost function values



# Minimisation_BT.log
# not sure what to plot from this yet

# ProfileQC.dat
# There is not much in this file, so leave it for now.


# RetJacobian.out
retJacobians_df     <- read.table(files[8], skip = 3, header = F)
retjacobians        <- matrix(unlist(t(retJacobians_df)), byrow = F, 298, 140)

retjacobians_melted <- melt(retjacobians)
breaks_ret          <- quantile(retjacobians_melted$value, (0:10)/10)


# Retrieved_BTs.dat
BT           <- read.table(files[9], skip = 2, header = T)
BT_diff      <- data.frame(Channel = BT$Channel, `B-O` = BT$Background - BT$Observed, `B-R` = BT$Background - BT$Retrieved) # calculate differences

BT_Melt      <- melt(BT, id = "Channel")
BT_diff_melt <- melt(BT_diff, id = "Channel")


# Retrieved_Profiles.dat
retrieved_profiles_values           <- read.table(files[10], skip = 3, nrows = 54, header = F)
retrieved_surface_values            <- read.table(files[10], skip = 57, nrows = 4, header = F)

colnames(retrieved_profiles_values) <- c("Pressure (hPa)", "T (K) (r)", "q (kg/kg) (r)", "Ozone (r)", "T (K) (b)", "q (kg/kg) (b)", "Ozone (b)")


# Plotting ----------------------------------------------------------------


#' This section plots the A- and Am-Matrix.out files. These are the analysis
#' error covaraince matrix and the propagated measurement noise matrix resp-
#' ectively. These are the expected retrieval error and the expected contri-
#' bution to the retrieval error from the observation noise. 
#' 
#' Output is in the same order as the elements from the input BMatrix and a-
#' re outputted for each iteration. 


# A-Matrix.out
A_Matrix_ggplot         <- ggplot(data_melted_A, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
                                  geom_tile(color = "gray90") +
                                  scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                                                       values = scales::rescale(breaks_A), 
                                                       name="Analysis error") +
                                  labs(title = "Analysis Covariance error, latest iteration results ")  +
                                  xlab("BMatrix element number") +
                                  ylab("BMatrix element number") +
                                  theme(plot.title = element_text(hjust = 0.5), 
                                        legend.key.height=unit(2,"cm"),
                                        legend.position = "right") +
                                  coord_fixed(expand = FALSE)


#Am-Matrix.out
Am_Matrix_ggplot        <- ggplot(data_melted_Am, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
                                  geom_tile(color = "gray90") +
                                  scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                                                       values = scales::rescale(breaks_Am), 
                                                       name="Noise") +
                                  labs(title = "Propagated measurement noise")  +
                                  xlab("BMatrix element number") +
                                  ylab("BMatrix element number") +
                                  theme(plot.title = element_text(hjust = 0.5), 
                                        legend.key.height=unit(2,"cm"),
                                        legend.position = "right") +
                                  coord_fixed(expand = FALSE)



#' This section plots the averaging kernel from the files in a vertical
#' line graph. All levels will be plotted. 


AVK_ggplot              <- ggplot(AVK_melt, aes(x = value, y = Var1)) +
                                  geom_path(aes(color = Var2, group = Var2)) +
                                  scale_colour_gradient(low = "blue", high = "green") +
                                  labs(title = "Averaging Kernel",
                                       x = "Value",
                                       y = "Element number",
                                       color = "Element nr") +
                                  theme_bw()

AVK_ozone_ggplot        <- ggplot(AVK_Ozone, aes(x = value, y = Var1)) +
                                  geom_path(aes(color = Var2, group = Var2)) +
                                  scale_colour_gradient(low = "blue", high = "green") +
                                  theme_bw()



#' This section plots both the background and the retrieval jacobian
#' values.
#' 
#' 

bgJacobians_ggplot      <- ggplot(bgjacobians_melted, 
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

bgJac_Lines_ggplot      <- ggplot(data = bgjacobians_melted, aes(x = Var1, y = as.numeric(value), color = Var2)) +
                                  geom_path() +
                                  scale_y_continuous(labels = scales::label_comma()) +
                                  scale_colour_viridis() +
                                  labs(title = "Background Jacobians per element",
                                       x = "channel selection number",
                                       y = "Value",
                                       legend = "BMatrix element number") +
                                  labs(color="BMatrix element number") +
                                  theme_bw() 


retJacobians_ggplot     <- ggplot(retjacobians_melted, 
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


retJac_Lines_ggplot     <- ggplot(data = retjacobians_melted, aes(x = Var1, y = as.numeric(value), color = Var2)) +
                                  geom_path() +
                                  scale_y_continuous(labels = scales::label_comma()) +
                                  scale_colour_viridis() +
                                  labs(title = "Background Jacobians per element",
                                       x = "channel selection number",
                                       y = "Value",
                                       legend = "BMatrix element number") +
                                  labs(color="BMatrix element number") +
                                  theme_bw()





#' This section plots the values of the minimisation cost function
#' per iteration. This is done for both the standard minimisation
#' and the brightness temperatures. 

minimisation_ggplot     <- ggplot(minimisation_values, aes(x = iteration, y = cost_function)) +
                                  geom_smooth(se = FALSE) +
                                  scale_x_continuous(breaks = c(0:length(minimisation))) +
                                  scale_y_log10(labels = scales::label_comma()) +
                                  annotation_logticks(sides = "l") +
                                  labs(title = "Cost function values per iteration",
                                       x = "Iteration",
                                       y = "Cost Function value") +
                                  theme_bw()


minimisation_BT_ggplot  <- ggplot()

#' This section plots the output Brightness Temperatures, for the direct
#' retrieval and the difference on each vertical level.

profiles_BT_ggplot      <- ggplot(data = BT_Melt, aes(x = value, y = Channel, colour = variable)) +
                                  geom_path() +
                                  labs(title = "Brightness Temperatures absolute values",
                                       x = "Temperature (K)",
                                       y = "Channel nr.",
                                       color = "Legend") +
                                  scale_y_log10() +
                                  annotation_logticks(sides = "l") +
                                  theme_bw()


profiles_BT_diff_ggplot <- ggplot(data = BT_diff_melt, aes(x = value, y = Channel, colour = variable)) +
                                  geom_vline(xintercept = 0, linetype = "dashed") +
                                  geom_path() +
                                  labs(title = "Brightness Temperatures differences",
                                       subtitle = "B = Background, O = Observation, R = Retrieved",
                                       x = "Temperature (K)",
                                       y = "Channel nr.",
                                       color = "Legend") +
                                  scale_y_log10() +
                                  annotation_logticks(sides = "l") +
                                  scale_colour_manual(values = c("red", "blue"), labels = c("B-O", "B-T")) +
                                  theme_bw()

#' This section plots the output profiles from Retrieved_profiles.dat 
#' and assigns them to variables. Three profiles are being retrieved:
#' Temperature, Water Vapour, Ozone.

colors                  <- c("Retrieved" = "red", "Background" = "blue")

# Water vapour q in kg/kg
profiles_q_ggplot       <- ggplot(data = retrieved_profiles_values) +
                                  geom_path(aes(x = `q (kg/kg) (r)`, y = `Pressure (hPa)`, colour = "Retrieved")) +
                                  geom_path(aes(x = `q (kg/kg) (b)`, y = `Pressure (hPa)`, colour = "Background"))  +
                                  scale_x_continuous(labels = scales::label_comma(), trans = "log10") +
                                  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
                                  annotation_logticks(sides = "lb") +
                                  labs(title = "Comparison retrieved/background for Water Vapour",
                                       subtitle = sprintf("Surface q (retrieved): %sK\nSurface q (background): %sK", 
                                                          retrieved_surface_values[2, 4], 
                                                          retrieved_surface_values[2, 5]),
                                       x= "q (kg/kg)",
                                       y = "Pressure (hPa)",
                                       colour = "Legend") +
                                  scale_color_manual(values = colors) +
                                  theme_bw()


# Temperature difference
profiles_T_ggplot       <- ggplot(data = retrieved_profiles_values) +
                                  geom_path(aes(x = `T (K) (r)`, y = `Pressure (hPa)`, colour = "Retrieved")) +
                                  geom_path(aes(x = `T (K) (b)`, y = `Pressure (hPa)`, colour = "Background"))  +
                                  scale_x_continuous(labels = scales::label_comma()) +
                                  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
                                  annotation_logticks(sides = "l") +
                                  labs(title = "Comparison retrieved/background for Temperature",
                                     subtitle = sprintf("Surface T (retrieved): %sK\nSurface T (background): %sK", 
                                                        retrieved_surface_values[1, 4], 
                                                        retrieved_surface_values[1, 5]),
                                     x= "T (K)",
                                     y = "Pressure (hPa)",
                                     colour = "Legend") +
                                  scale_color_manual(values = colors) +
                                  theme_bw()
  
# Ozone
profiles_O3_ggplot      <- ggplot(data = retrieved_profiles_values) +
                                  geom_path(aes(x = `Ozone (r)`, y = `Pressure (hPa)`, colour = "Retrieved")) +
                                  geom_path(aes(x = `Ozone (b)`, y = `Pressure (hPa)`, colour = "Background"))  +
                                  scale_x_continuous(labels = scales::label_comma()) +
                                  scale_y_continuous(trans=reverselog_trans(10), labels = scales::label_comma()) +
                                  annotation_logticks(sides = "lb") +
                                  labs(title = "Comparison retrieved/background for Ozone",
                                       x= "q (kg/kg)",
                                       y = "Pressure (hPa)",
                                       colour = "Legend") +
                                  scale_color_manual(values = colors)+
                                  theme_bw()



# Save to file ------------------------------------------------------------

# Matrices
ggsave(plot     = A_Matrix_ggplot, 
       filename = "AMatrix.pdf", 
       device   = "pdf",
       path     = output_directory)

ggsave(plot     = Am_Matrix_ggplot, 
       filename = "AM_Matrix.pdf", 
       device   = "pdf",
       path     = output_directory)

# Averaging kernel
ggsave(plot     = AVK_ggplot, 
       filename = "AVK.pdf", 
       device   = "pdf",
       path     = output_directory)

# Jacobians
ggsave(plot     = bgJacobians_ggplot, 
       filename = "bgJacobians.pdf", 
       device   = "pdf",
       height   = 21,
       width    = 29.7,
       units    = "cm",
       path     = output_directory)

ggsave(plot     = bgJac_Lines_ggplot, 
       filename = "bgJacobians_lines.pdf", 
       device   = "pdf",
       height   = 21,
       width    = 29.7,
       units    = "cm",
       path     = output_directory)

ggsave(plot     = retJacobians_ggplot, 
       filename = "retJacobians.pdf", 
       device   = "pdf",
       height   = 21,
       width    = 29.7,
       units    = "cm",
       path     = output_directory)

ggsave(plot     = retJac_Lines_ggplot, 
       filename = "retJacobians_lines.pdf", 
       device   = "pdf",
       height   = 21,
       width    = 29.7,
       units    = "cm",
       path     = output_directory)

# Retrieved BTs
ggsave(plot     = profiles_BT_ggplot,
       filename = "retBT.pdf",
       device   = "pdf",
       path     = output_directory
       ) 

ggsave(plot     = profiles_BT_diff_ggplot,
       filename = "BT_Difference.pdf",
       device   = "pdf",
       path     = output_directory)

# Retrieved Profiles
ggsave(plot     = profiles_T_ggplot, 
       filename = "Temperature.pdf", 
       device   = "pdf",
       path     = output_directory)

ggsave(plot     = profiles_q_ggplot, 
       filename = "WaterVapour.pdf", 
       device   = "pdf",
       path     = output_directory)

ggsave(plot     = profiles_O3_ggplot, 
       filename = "Ozone.pdf", 
       device   = "pdf",
       path     = output_directory)



# End of file -------------------------------------------------------------