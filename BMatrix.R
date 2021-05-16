## Header




# BMatrix  ----------------------------------------------------------------

## Libraries
libs <- c("rhdf5", "reshape2", "ggplot2", "ncdf4", "soundgen", "corrplot")
pacman::p_load(char = libs)


# File paths
path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/NWPSAF_1DVar_1.3/test_data/IASI_54L/InputData/Bmatrix_54L"

# Read in example BMatrix
# Short workaround, probably not optimized
rows_BM <- readLines(path)
rows_BM <- tail(rows_BM, -3)
rows_BM <- rows_BM[1:1548] # first 1548 rows are the matrix, there's duplicates for some reason, need to check

# Merge together in a 86*86 matrix
BMatrix <- do.call(rbind, 
                   tapply(rows_BM, ceiling(seq_along(rows_BM)/18), function(x) 
                     strsplit(trimws(paste0(x, collapse = ' ')), '\\s+')[[1]]))

# Make numeric matrix
BMatrix_num <- data.matrix(BMatrix)
mode(BMatrix_num) <- "numeric"


# Cut out individual components
BMatrix_T <- BMatrix_num[1:54, 1:54]
BMatrix_lnq <- BMatrix_num[55:83, 55:83]
BMatrix_T_Surf <- BMatrix_num[84, ]
BMatrix_lnq_Surf <- BMatrix_num[85, ]
BMatrix_T_Skin <- BMatrix_num[86, ]

# Plot different values
corrplot(BMatrix_T, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 150, palette = "spectral", rev = TRUE))
corrplot(BMatrix_lnq, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 150, palette = "spectral", rev = TRUE))

# Individual values for the final three covariances
T_Surf <- BMatrix_T_Surf[length(BMatrix_T_Surf) - 3]
lnq_Surf <- BMatrix_lnq_Surf[length(BMatrix_lnq_Surf) -2]
T_Skin <- tail(BMatrix_T_Skin, 1)



# Ozone BMatrix -----------------------------------------------------------

## This section makes our own BMatrix/covariance matrix for observed ozone values
## used from the example .ncdf4 file, an apriori file

apriori_path <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/E3a_BMatrices/S-O3M_GOME_NHP_02_M01_20210311013857Z_20210311014157Z_N_O_20210311023812Z.hdf5"
apriori_nc <- nc_open(apriori_path)

# Get an error covariance profile
apriori_error <- ncvar_get(apriori_nc, "DATA/AprioriErrorCovariance") # dimension 720x42x42
apriori_error_covariance <- apriori_error[1, , ]


# Interpolate to 54 layers
BMatrix_ozone_interp <- soundgen:::interpolMatrix(apriori_error_covariance, nr = 54, nc = 54, interpol = "approx")

# Mirror matrix to get from TOA to surface instead of currently surface to TOA
BMatrix_ozone_interp_mirrored <- BMatrix_ozone_interp[ncol(BMatrix_ozone_interp):1, nrow(BMatrix_ozone_interp):1]

# Plot for example
corrplot(BMatrix_ozone_interp_mirrored, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 150, palette = "spectral", rev = TRUE))

# e macht om correlatie lagen aan te geven
# deel door interpolaties



# Merge all together ------------------------------------------------------

n <- 140 # dimension of new BMatrix

# Make empty matrix
BMatrix_w_Ozone <- matrix(0, ncol = n, nrow = n)

# Fill BMatrix with values
BMatrix_w_Ozone[1:54,1:54] <- BMatrix_T
BMatrix_w_Ozone[55:83, 55:83] <- BMatrix_lnq
BMatrix_w_Ozone[84, 84] <- T_Surf
BMatrix_w_Ozone[85, 85] <- lnq_Surf
BMatrix_w_Ozone[86, 86] <- T_Skin
BMatrix_w_Ozone[87:140, 87:140] <- BMatrix_ozone_interp_mirrored


# Add Complete BMatrix
BMatrix_w_Ozone[1:86, 1:86] <- BMatrix_num

# Just diagonal of ozone matrix
# Add to matrix
ozone_diag_m <- matrix(0, ncol = 54, nrow = 54)
diag(ozone_diag_m) <- diag(BMatrix_ozone_interp)

# Change previous matrix
BMatrix_w_Ozone[84:137, 84:137] <- ozone_diag_m
  
# Plot resulting BMatrix
corrplot(BMatrix_w_Ozone, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 150, palette = "spectral", rev = TRUE))



# Damping -----------------------------------------------------------------

## This section is to dampen the values of the BMatrix as currently it's ill-fitted
## and non invertible. Copy matrix an do some test with adding a c value
BM_tst <- BMatrix_ozone_interp

# Set Marquardt-Levenberg coefficient
# after trying with all kinds of values, 1e-9 is the smallest value that makes it work
# 1e-10 still gives  an error, makes sense as the width of all decimal places is 10, so 
# 1e-10 does not alter the values anymore.
c <- 1e-9

# Add some values
diag(BM_tst) <- diag(BM_tst) + c

# Test if error shows, FALSE if solvable and invertible
is.error(solve(BM_tst))

# then solve and see actual values
solve(BM_tst)

# Then add to matrix
BMatrix_w_Ozone[84:137, 84:137] <- BM_tst # still non invertible



## Second method: damping after merging the total matrix
BMatrix_w_Ozone[84:137, 84:137] <- BMatrix_ozone_interp # Set back original ozone matrix

# Dampen again
diag(BMatrix_w_Ozone) <- diag(BMatrix_w_Ozone) + c

# Plus test
is.error(solve(BMatrix_w_Ozone)) # FALSE -> Good

# Write to file, as 1DVar wants -------------------------------------------

# File location
BM_File <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/E3a_BMatrices/BMatrix_54L"

# First header files
lin1 <- c("Peter Weston's NH 54 RTcoef level B matrix from L70 Covstats for Sea. Ozone error covariance added from GOME O3 observation")
lin2 <- c("1-54=temp (K), 55-83=lnq (g/kg) (bottom 29 levels only), 84=Tsurf, 85=lnq surf, 86=Tskin, 87-140= Ozone(ppmv)")
lin3 <- c("140") # The total amount of nrows/ncolumns for this matrix

# Append to eachother
header <- append(lin1, c(lin2, lin3))

# Write to file
BMatrix_Correct_Notation <- toupper(sprintf("%16.8e", BMatrix_w_Ozone))

# Do 4 times, as with the example matrix, not exactly sure why TODO
for(i in 1:4){
cat(header, apply(matrix(BMatrix_Correct_Notation, ncol=5, byrow=TRUE), 1, paste, collapse="") ,
    sep="\n", file = BM_File, append = TRUE)
}




# Write to file with function ---------------------------------------------

# Write to file
BMatrix_Correct_Notation <- toupper(sprintf("%16.8e", BMatrix))

# Do 4 times, as with the example matrix, not exactly sure why TODO
for(i in 1:4){
  cat(header, apply(matrix(BMatrix_Correct_Notation, ncol=5, byrow=TRUE), 1, paste, collapse="") ,
      sep="\n", file = BM_File, append = TRUE)
}





# Abovementioned method causes a singular matrix to aoccur: Matrix is not positive definite
# so we now  try to just use the diagonal and not the complete 54*54 interpolated ozone matrix
# First look at which values are not colinear
colnames(model.matrix(BMatrix_w_Ozone))[-which(is.na(coef(BMatrix_w_Ozone))==TRUE)][-1]
which(is.na(coef(BMatrix_w_Ozone))==TRUE)

# short if singular function, if false then not invertible
is.singular <- function(m) "matrix" %in% class(try(solve(m),silent=TRUE))
is.singular(BMatrix_w_Ozone)



# Write all to file -------------------------------------------------------

# this section is for the email to James
write.table(BMatrix_num, "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/BMatrix_54L.csv", 
            sep = ",", row.names = FALSE, col.names = FALSE)


corrplot(BMatrix_num, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 50, palette = "RdYlBu", rev = TRUE))
corrplot(BMatrix_num, is.corr = FALSE, method = "color", addgrid.col = "grey", col = hcl.colors(n = 50, palette = "Zissou 1", rev = TRUE))

corrplot(BMatrix_w_Ozone, is.corr = FALSE, 
         method = "color", 
         sig.level = 0.01, insig = "blank", 
         addgrid.col = "grey", 
         tl.col="black", tl.srt=45,
         col = hcl.colors(n = 50, palette = "Zissou 1", rev = TRUE))


# reshape
bm_test <- melt(BMatrix)
bm_oz <- melt(BMatrix_w_Ozone)


# Normal
ggplot(bm_test, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red" , "white", "yellow", "green", "blue", "purple"),
                       values = scales::rescale(c(-0.4609482, 0, 0.05, 0.5, 2, 4.109625)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)

# With ozone
ggplot(bm_oz, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red" , "white", "yellow", "green", "blue", "purple"),
                       values = scales::rescale(c(-0.4609482, 0, 0.1, 1, 3, 5.41319)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)


# With ozone, highlight values between 
ggplot(bm_oz, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red" , "white", "yellow"),
                       values = scales::rescale(c(-0.4609482, 0, 0.5)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)




# Make Bmatrix from MLL Climatology ---------------------------------------

# # -----------------------------------------------------------------------


# Testing
path_stats <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_stats.dat"
ppmv_stats <- read.table(path_stats, skip = ((latidx + 9) * 70) + 1, nrow = 66, header = TRUE)

# run it
newcov <- cov_profile(ppmv_stats[,2])
newcov_v2 <- cov_profile_v2(ppmv_stats[,2])

# See where the differences are
cov_differences <- newcov - newcov_v2

# Plot
cov_melted <- melt(newcov)
cov_melted_v2 <- melt(newcov_v2)
cov_differences_melted <- melt(cov_differences)

# Corrplot
corrplot(newcov, is.corr = F, method = "color")
corrplot(newcov_v2, is.corr = F, method = "color")

# Normal
ggplot(cov_melted, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red", "yellow", "green", "blue", "purple"),
                       values = scales::rescale(c(0.007, 0.2, 0.4, 0.9, 1.437)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)

# V2
ggplot(cov_melted_v2, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red" , "yellow", "green", "blue", "purple"),
                       values = scales::rescale(c(0.000980,0.2, 0.4, 0.9, 1.437)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)


# differences
ggplot(cov_differences_melted, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("red" , "yellow", "green", "blue", "purple"),
                       values = scales::rescale(c(0, 0.05, 0.2, 0.35, 0.511)), 
                       name="Covariance error") +
  theme_bw() +
  coord_fixed(expand = FALSE)


# Plot Climatology Stats --------------------------------------------------

path_stats <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_stats.dat"
path_ppmv <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E2_Auxiliary/E2c_Climatologies/ML_Climatology/ML_ppmv_table.dat"

# read
ppmv_stats <- read.table(path_stats, skip = ((latidx + 9) * 70) + 1, nrow = 66, header = TRUE)
ppmv_clim <- read.table(path_ppmv, skip = ((latidx +9 ) * 70) + 1, nrow = 66, header = TRUE)

# Ratio
ppmv_percentage <- (ppmv_stats / ppmv_clim) * 100
ppmv_percentage[, 1] <- ppmv_stats[, 1]

# melt
ppmv_melt_percentage <- melt(ppmv_percentage, id ="Z.level")
ppmv_melt <- melt(ppmv_stats, id = "Z.level")

# plot percentages for one latitudinal region
ggplot(ppmv_melt_percentage, aes(x = value, y = `Z.level`, colour = variable)) +
  geom_path() +
  xlab("Standard Deviation (%)") +
  labs(colour = "Month", title = paste0("For latitudinal region: ", latidx * 10, "-", (latidx*10)+10, "N"))  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# Plot standard deviations too
ggplot(ppmv_melt, aes(x = value, y = `Z.level`, colour = variable)) +
  geom_path() +
  xlab("Standard Deviation") +
  labs(colour = "Month", title = paste0("For latitudinal region: ", latidx * 10, "-", (latidx*10)+10, "N"))  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  
  


# # -----------------------------------------------------------------------

# Get parameters
pressure_lvls <- zstar_press(ppmv_stats[,1])
kgkg_values <- convert_MMR_VMR(ppmv_stats[,2], rev = TRUE)

# Intrpolate to 54 layers
pressure_interp <- approx(pressure_lvls, n=54)$y
kgkg_interp <- approx(kgkg_values, n=54)$y

# Final version
newcov_v3 <- cov_profile(pressure_interp, kgkg_interp)

# melt
newcov_v3_melted <- melt(newcov_v3)

# Find breakpoints for legend
breaks <- quantile(newcov_v3_melted$value, (0:6)/6)

# Plot
ggplot(newcov_v3_melted, aes(x= Var1, y = Var2, fill = as.numeric(value))) +
  geom_tile(color = "gray90") +
  scale_fill_gradientn(colours = c("purple", "blue", "green", "yellow", "orange", "red"),
                       values = scales::rescale(breaks), 
                       name="Covariance error (kg/kg)") +
  labs(title = paste0("For latitudinal region: ", latidx * 10, "-", (latidx*10)+10, "N"))  +
  xlab("Z Level") +
  ylab("Z level") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height=unit(2,"cm"),
        legend.position = "right") +
  coord_fixed(expand = FALSE)


# Write to file -----------------------------------------------------------


# Mirror matrix to get from TOA to surface instead of currently surface to TOA 
# for when used as input in the original BMatrix
newcov_v3 <- newcov_v3[ncol(newcov_v3):1, nrow(newcov_v3):1]
    
# Fill in 
BMatrix_w_Ozone[84:137, 84:137] <- newcov_v3

# Dampen
diag(BMatrix_w_Ozone) <- diag(BMatrix_w_Ozone) + c

# Test if solvable
is.error(solve(BMatrix_w_Ozone)) # FALSE -> Good

# File location
BM_File <- "C:/Users/bartq/Documents/WUR/Internship/E_Data/E3_1DVAR/E3a_BMatrices/BMatrix_54L"

# First header files
lin1 <- c("Peter Weston's NH 54 RTcoef level B matrix from L70 Covstats for Sea. Ozone error covariance added from GOME O3 observation")
lin2 <- c("1-54=temp (K), 55-83=lnq (g/kg) (bottom 29 levels only), 84=Tsurf, 85=lnq surf, 86=Tskin, 87-140= Ozone(ppmv)")
lin3 <- c("140") # The total amount of nrows/ncolumns for this matrix

# Append to eachother
header <- append(lin1, c(lin2, lin3))

# Write to file
BMatrix_Correct_Notation <- toupper(sprintf("%16.8e", BMatrix_w_Ozone))

# Remove file if exists
if(file.exists(BM_File)){
  file.remove(BM_File)
}

# Do 4 times, as with the example matrix, not exactly sure why TODO
for(i in 1:4){
  cat(header, apply(matrix(BMatrix_Correct_Notation, ncol=5, byrow=TRUE), 1, paste, collapse="") ,
      sep="\n", file = BM_File, append = TRUE)
}



# Test ratios -------------------------------------------------------------

# diag: kgkg / stdev voor ratio
# niet diag: 0.25**2

# in 66 levels:
kgkg <- convert_MMR_VMR(clim_df$AUG, rev = T)
stdev <- convert_MMR_VMR(ppmv_stats$AUG, rev = T)
ratio <- stdev / kgkg

# Interpolate to 54 layers
ratio_matrix <- data.frame(
  pressure = approx(pressure_conv, n=54)$y, 
  ratio = approx(ratio, n=54)$y,
  values = approx(kgkg_ratio, n=54)$y
)

# calculate values
ratio_matrix_covariance <- cov_profile_ratio(ratio_matrix$pressure, ratio_matrix$ratio, fac = 1)

# Get covariance
rmc <- cov(ratio_matrix_covariance)

# plot
ggplot_Input_Files(ratio_matrix_covariance, type = "BMatrix")

# error plots
ggplot(ratio_matrix, aes(x = pressure, y = values)) +
         geom_line() +
         geom_point() +
         geom_errorbar(aes(ymin = (values - (ratio * values)), 
                           ymax = (values + (ratio * values)),
                           width = 0.2)) +
         coord_flip() +
         theme_bw() +
         scale_y_reverse(limits = c(1015, 0), 
                         expand = c(0,0)) +
         scale_x_log10()


# cov
rmc_melt <- melt(rmc)

# plot matrix  
data_melted <- melt(ratio_matrix_covariance)

# get breaks
breaks <- quantile(data_melted$value, (0:6)/6)

ggplot(rmc_melt, 
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

