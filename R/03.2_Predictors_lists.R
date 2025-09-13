#_____________________________________________________________________________
#                        Predictors: Monthly Stacks 
#_____________________________________________________________________________




library(tidyverse)
library(raster)
library(terra)
library(ncdf4)
library(rerddap)
library(rerddapXtracto)
library(stringr)
library(lubridate)
source("R/00_Helper_Functions.R") 


# Import DATA -------------------------------------------------------------



mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1.tif")

wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1.tif")

sst_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1.tif")

depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")

slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1.tif")

rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1.tif")

dist_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1.tif")

id_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/id_stack_0.1.tif")

month_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/month_stack_0.1.tif")

year_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/year_stack_0.1.tif")

lat_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lat_stack_0.1.tif")

lon_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lon_stack_0.1.tif")




# Convert to factor
month_stack_0.1 <- as.factor(month_stack_0.1)
year_stack_0.1 <- as.factor(year_stack_0.1)
# Check again (should be TRUE)
is.factor(month_stack_0.1)

# Plot the categorical raster
plot(month_stack_0.1)




pred_stack_lst <- list(
  depth_stack_0.1, 
  slope_stack_0.1, 
  rough_stack_0.1, 
  dist_stack_0.1,
  sst_stack_0.1,
  mld_stack_0.1,
  uv_stack_0.1,
  wz_stack_0.1,
  chl_stack_0.1,
  lon_stack_0.1,
  lat_stack_0.1,
  id_stack_0.1,
  month_stack_0.1,
  year_stack_0.1
)


pred_stack_lst


# Name the list elements
names(pred_stack_lst) <- c("depth", 
                           "slope", 
                           "roughness", 
                           "dist2000",
                           "thetao", 
                           "mltost",  
                           "uv",
                           "wz",
                           "chl",
                           "lon", 
                           "lat", 
                           "id",
                           "month",
                           "year")

names(pred_stack_lst)

names(pred_stack_lst$thetao[[1]])
names(pred_stack_lst$thetao[[186]])



# Create an empty list to store monthly stacks
monthly_stacks_lst <- list()

# Loop through each month (assuming 168 layers for each variable)
for (month in 1:186) {
  # Extract the month name from the original stack
  month_name <- names(pred_stack_lst$thetao)[month]
  
  # Create a stack for the current month
  monthly_stack <- rast(list(
    pred_stack_lst$depth[[month]], 
    pred_stack_lst$slope[[month]], 
    pred_stack_lst$roughness[[month]], 
    pred_stack_lst$dist2000[[month]],
    pred_stack_lst$thetao[[month]], 
    pred_stack_lst$mltost[[month]], 
    pred_stack_lst$uv[[month]], 
    pred_stack_lst$wz[[month]],
    pred_stack_lst$chl[[month]],
    pred_stack_lst$lon[[month]], 
    pred_stack_lst$lat[[month]], 
    pred_stack_lst$id[[month]],
    pred_stack_lst$month[[month]],
    pred_stack_lst$year[[month]]
  ))
  
  # Set names for the layers
  names(monthly_stack) <- c("depth", 
                            "slope", 
                            "roughness", 
                            "dist2000",
                            "thetao", 
                            "mltost",  
                            "uv",
                            "wz",
                            "chl",
                            "lon", 
                            "lat", 
                            "id",
                            "month",
                            "year")
  
  # Assign the monthly stack to the corresponding month name
  monthly_stacks_lst[[month_name]] <- monthly_stack
}







monthly_stacks_lst_0.1 <- monthly_stacks_lst

monthly_stacks_lst_0.1$`2025-06-01`
monthly_stacks_lst_0.1[[186]]
names(monthly_stacks_lst_0.1[[186]])





# Extebded and transformed  -----------------------------------------------



mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1_ext.tif")

wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1_ext.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1_ext.tif")

sst_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1_ext.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1_ext.tif")

depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1_ext.tif")

slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1_ext.tif")

rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1_ext.tif")

dist_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1_ext.tif")
dist_stack_0.1


lat_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lat_stack_0.1_ext.tif")

lon_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lon_stack_0.1_ext.tif")

month_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/month_stack_0.1_ext.tif")

year_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/year_stack_0.1_ext.tif")


id_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/id_stack_0.1_ext.tif")

dist_seamount_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_seamounts_month_0.1.tif")
dist_seamount_stack_0.1
dist_knolls_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_knolls_month_0.1.tif")


# Convert to factor
month_stack_0.1 <- as.factor(month_stack_0.1)
year_stack_0.1 <- as.factor(year_stack_0.1)
# Check again (should be TRUE)
is.factor(month_stack_0.1)

## Transformations:

# depth_stack_0.1 <- -depth_stack_0.1
# plot(depth_stack_0.1[[1]])
# depth_stack_0.1 <- log1p(depth_stack_0.1)
print(depth_stack_0.1)
plot(depth_stack_0.1[[1]])
# terra::units(depth_stack_0.1) <- "log1p(m)"
terra::units(depth_stack_0.1) <- "m"


# slope_stack_0.1 <- log1p(slope_stack_0.1)
# terra::units(slope_stack_0.1) <- "log1p(degrees)"
terra::units(slope_stack_0.1) <- "degrees"
slope_stack_0.1
# 
# 
# rough_stack_0.1 <- log1p(rough_stack_0.1)
# terra::units(rough_stack_0.1) <- "log1p(roughness)"
terra::units(rough_stack_0.1) <- "roughness"
rough_stack_0.1



# dist_log_stack_0.1 <- log1p(dist_stack_0.1)
# terra::units(dist_log_stack_0.1) <- "log1p(m)"
# dist_log_stack_0.1
terra::units(dist_stack_0.1) <- "km"
dist_stack_0.1


# uv_stack_0.1 <- log1p(uv_stack_0.1)
# terra::units(uv_stack_0.1) <- "log1p(m s-1)"
terra::units(uv_stack_0.1) <- "m s-1"
uv_stack_0.1


chl_stack_0.1 <- log(chl_stack_0.1)
terra::units(chl_stack_0.1) <- "log(milligram m-3)"
chl_stack_0.1

# mld_stack_0.1 <-log(mld_stack_0.1)
# terra::units(mld_stack_0.1) <- "log(m)"
terra::units(mld_stack_0.1) <- "m"
mld_stack_0.1

terra::units(dist_knolls_stack_0.1) <- "km"
terra::units(dist_seamounts_stack_0.1) <- "km"

pred_stack_lst <- list(
  depth_stack_0.1, 
  slope_stack_0.1, 
  rough_stack_0.1, 
  dist_stack_0.1,
  sst_stack_0.1,
  mld_stack_0.1,
  uv_stack_0.1,
  wz_stack_0.1,
  chl_stack_0.1,
  lon_stack_0.1,
  lat_stack_0.1,
  id_stack_0.1,
  month_stack_0.1,
  year_stack_0.1,
  dist_seamount_stack_0.1, 
  dist_knolls_stack_0.1
)


pred_stack_lst


# Name the list elements
names(pred_stack_lst) <- c("depth", 
                           "slope", 
                           "roughness", 
                           "dist2000",
                           "thetao", 
                           "mltost",  
                           "uv",
                           "wz",
                           "chl",
                           "lon", 
                           "lat", 
                           "id",
                           "month",
                           "year",
                           "dist_seamount",
                           "dist_knoll")

names(pred_stack_lst)

names(pred_stack_lst$thetao[[1]])
names(pred_stack_lst$thetao[[186]])



# Create an empty list to store monthly stacks
monthly_stacks_lst <- list()

# Loop through each month (assuming 168 layers for each variable)
for (month in 1:186) {
  # Extract the month name from the original stack
  month_name <- names(pred_stack_lst$thetao)[month]
  
  # Create a stack for the current month
  monthly_stack <- rast(list(
    pred_stack_lst$depth[[month]], 
    pred_stack_lst$slope[[month]], 
    pred_stack_lst$roughness[[month]], 
    pred_stack_lst$dist2000[[month]],
    pred_stack_lst$thetao[[month]], 
    pred_stack_lst$mltost[[month]], 
    pred_stack_lst$uv[[month]], 
    pred_stack_lst$wz[[month]],
    pred_stack_lst$chl[[month]],
    pred_stack_lst$lon[[month]], 
    pred_stack_lst$lat[[month]], 
    pred_stack_lst$id[[month]],
    pred_stack_lst$month[[month]],
    pred_stack_lst$year[[month]],
    pred_stack_lst$dist_seamount[[month]],
    pred_stack_lst$dist_knoll[[month]]
  ))
  
  # Set names for the layers
  names(monthly_stack) <- c("depth", 
                            "slope", 
                            "roughness", 
                            "dist2000",
                            "thetao", 
                            "mltost",  
                            "uv",
                            "wz",
                            "chl",
                            "lon", 
                            "lat", 
                            "id",
                            "month",
                            "year",
                            "dist_seamount",
                            "dist_knoll")
  
  # Assign the monthly stack to the corresponding month name
  monthly_stacks_lst[[month_name]] <- monthly_stack
}







monthly_stacks_lst_0.1_trans <- monthly_stacks_lst

monthly_stacks_lst_0.1_trans$`2025-06-01`
monthly_stacks_lst_0.1_trans[[186]]
names(monthly_stacks_lst_0.1_trans[[186]])





      