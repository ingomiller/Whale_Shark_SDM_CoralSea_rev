#_____________________________________________________________________________
#                        Prediction Rasters
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




# Vertical Current Velocity -----------------------------------------------

#____________ Bluelink BRAN2020 - wo
## Note, we start with Bluelink as these data come in 0.1 degrees resolution, which is our target res, so we will use this as reference layer for the other predictors


nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Environmental_Varibales_Downloaded/Bluelink/Monthly"
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

nc_bran <- terra::rast(nc_files)
crs(nc_bran) <- "EPSG:4326"
print(nc_bran)
names(nc_bran)

plot(nc_bran[[1]])


# Extract the indices of the first layer of each month
first_layer_indices <- seq(1, nlyr(nc_bran), by = 51)

# Subset the SpatRaster to keep only the first layer of each month
Wo_monthly_0.1 <- subset(nc_bran, first_layer_indices)

Wo_monthly_0.1

# monthly dates from 2005-01-16 to 2023-12-16
dates_monthly <- base::seq.Date(
  from = base::as.Date("2005-01-01"),
  to   = base::as.Date("2023-12-01"),
  by   = "month"
)

dates_monthly

# change layer names to dates 
names(Wo_monthly_0.1) <- dates_monthly

names(Wo_monthly_0.1)
print(Wo_monthly_0.1)

plot(Wo_monthly_0.1[[165]], range  = c(-0.00005, 0.00005))

## crop to study extent 

# aoi <- terra::ext(134, 180, -45, 10)
aoi <- terra::ext(110, 180, -50, 50)
aoi

Wo_monthly_0.1_crop <- crop(Wo_monthly_0.1, aoi)

## crop to temporal extent of data  
Wo_monthly_0.1_crop <- Wo_monthly_0.1_crop[[61:228]]

range <- c(-0.00003, 0.00003)
plot(Wo_monthly_0.1_crop[[168]], range = range)


## save a reference raster for later conversions 
# writeRaster(Wo_monthly_0.1_crop[[1]], filename = "data/processed/reference_raster_0.1.tif", overwrite = TRUE)
writeRaster(Wo_monthly_0.1_crop[[1]], filename = "data/processed/reference_raster_0.1_ext.tif", overwrite = TRUE)

#____________ Adding CMEMS Wo for >=2024-01

output_path <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files"
nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/wo_after2024/"
depth_layer <- 5  # User-defined depth layer (e.g., 5m)

reference_raster_0.1 <- terra::rast("data/processed/reference_raster_0.1_ext.tif")
reference_raster_0.1

create_monthly_rasters_from_files_CM(nc_folder, depth_layer, reference_raster = reference_raster_0.1, output_path = output_path, output_filename = "Copernicus_W_5m_Monthly_Stack_2024to2025_10km")

Wz_CM_monthly <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Copernicus_W_5m_Monthly_Stack_2024to2025_10km.tif")
  

print(Wz_CM_monthly)
print(Wo_monthly_0.1_crop)

Wz_Stack_BRAN_CM <- c(Wo_monthly_0.1_crop, Wz_CM_monthly)
Wz_Stack_BRAN_CM[[168]]
print(Wz_Stack_BRAN_CM)

# writeRaster(Wz_Stack_BRAN_CM, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1.tif", overwrite = TRUE)
writeRaster(Wz_Stack_BRAN_CM, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1_ext.tif", overwrite = TRUE)



# Horizontal Current Velocity: uv -----------------------------------------

nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/mld.thetao.uo.vo.so/"

nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

# inspect files tructure
nc <- ncdf4::nc_open(nc_files[18])
nc
nc_close(nc)

# Read the NetCDF file separately for "vo" and "uo"
vo_rast <- terra::rast(nc_files, subds = "vo")
uo_rast <- terra::rast(nc_files, subds = "uo")

names(vo_rast) <- time(vo_rast)
names(uo_rast) <- time(uo_rast)

# only use dates needed to save time 
vo_rast <- vo_rast[[61:246]]
uo_rast <- uo_rast[[61:246]]

print(vo_rast)
print(uo_rast)
names(uo_rast)

plot(vo_rast[[1]])

currents_stack <- c(vo_rast, uo_rast)
currents_stack

## crop to target area because this caluclations takes a lot of time 
pad_deg <- 2
e <- terra::ext(reference_raster_0.1)
e_big <- terra::ext(
  terra::xmin(e) - pad_deg,
  terra::xmax(e) + pad_deg,
  terra::ymin(e) - pad_deg,
  terra::ymax(e) + pad_deg
)

reference_raster_0.1_plus <- terra::extend(reference_raster_0.1, e_big)
reference_raster_0.1_plus
currents_stack_crop <- terra::crop(currents_stack, reference_raster_0.1_plus)
currents_stack_crop

## Calculate mean velocity uv

currents_uv <- calculate_uv(currents_stack)

print(currents_uv)
names(currents_uv) <- time(currents_uv)

plot(currents_uv[[1]])




# Resample to match extent and resolution
currents_uv_0.1 <- terra::resample(currents_uv, reference_raster_0.1, method = "bilinear", threads = TRUE)

print(currents_uv_0.1[[186]])

varnames(currents_uv_0.1) <- "uv (mean velocity)"
plot(currents_uv_0.1[[168]])



# writeRaster(currents_uv_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1.tif", overwrite = TRUE)

writeRaster(currents_uv_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1_ext.tif", overwrite = TRUE)





# Chlorophyll -------------------------------------------------------------

nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/chl/"

nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

# inspect files tructure
nc <- ncdf4::nc_open(nc_files[1])
nc
nc_close(nc)

chl_stack <- terra::rast(nc_files, subds = "CHL")
print(chl_stack[[259]])
names(chl_stack) <- time(chl_stack)

#fix duplicates (2015 was duplicated)
dates <- names(chl_stack)
dups <- dates[base::duplicated(dates)] |> base::unique()
dups

# keep first occurrence of each name
keep <- !base::duplicated(dates)
keep
chl_stack_nodup <- chl_stack[[keep]]
chl_stack_nodup

chl_stack <- chl_stack_nodup[[61:247]]
print(chl_stack)

# Resample to match extent and resolution
chl_stack_0.1 <- terra::resample(chl_stack, reference_raster_0.1, method = "bilinear", threads = TRUE)
chl_stack_0.1 <- terra::mask(chl_stack_0.1, reference_raster_0.1)

print(chl_stack_0.1)
plot(chl_stack_0.1[[13]])

# writeRaster(chl_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1.tif", overwrite = TRUE)

writeRaster(chl_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1_ext.tif", overwrite = TRUE)


# SST Data ----------------------------------------------------------------

# SST_info <- rerddap::info('jplMURSST41mday', url = 'https://coastwatch.pfeg.noaa.gov/erddap/')
# SST_info
# str(SST_info)
# 
# 
# # Lon_range <- c(-179.99, 180)
# # Lat_range <- c(-50, 50)
# # time_range <- c("2010-01-01", 
# #                 "2010-01-31")
# # sst_nc1 <- rerddap::griddap(SST_info,
# #                             longitude = Lon_range,
# #                             latitude = Lat_range,
# #                             fields = "sst",
# #                             time =  time_range,
# #                             read = FALSE
# # )
# 
# 
# 
# 
# # Download global data 
# res <- download_mursst41_monthly(
#   from_year = 2005L,
#   out_dir   = "/Volumes/Ingo_PhD_2/MUR_SST_Global/",
#   progress  = "cli",   # or "text" if you don’t use {cli}
#   write_log = TRUE
# )

nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/mld.thetao.uo.vo.so/"

nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

# inspect files tructure
nc <- ncdf4::nc_open(nc_files[18])
nc
nc_close(nc)

sst_stack <- terra::rast(nc_files, subds = "thetao")
print(sst_stack)
names(sst_stack) <- time(sst_stack)

sst_stack


sst_stack <- sst_stack[[61:246]]
print(sst_stack)
plot(sst_stack[[1]])



# from nc attrs
fill   <- -32767
scale  <- 0.000732444226741791
offset <- 21


terra::crs(sst_stack) <- "EPSG:4326"


#set fill to NA (before any interpolation)
sst_clean <- sst_stack
sst_clean[sst_clean == fill] <- NA
sst_clean



# Resample to match extent and resolution
sst_stack_0.1 <- terra::resample(sst_clean, reference_raster_0.1, method = "bilinear", threads = TRUE)

print(sst_stack_0.1)
plot(sst_stack_0.1[[1]])

# writeRaster(sst_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1.tif", overwrite = TRUE)

writeRaster(sst_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1_ext.tif", overwrite = TRUE)



# Mixed layer depth -------------------------------------------------------



nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/mld.thetao.uo.vo.so/"

nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

# inspect files tructure
nc <- ncdf4::nc_open(nc_files[18])
nc
nc_close(nc)

mld_stack <- terra::rast(nc_files, subds = "mlotst")
print(mld_stack[[1]])
names(mld_stack) <- time(mld_stack)

mld_stack

mld_stack <- mld_stack[[61:246]]
print(mld_stack)
plot(mld_stack[[1]])



# from nc attrs
fill   <- -32767

terra::crs(mld_stack) <- "EPSG:4326"

#set fill to NA (before any interpolation)
mld_clean <- mld_stack
mld_clean[mld_clean == fill] <- NA
mld_clean

plot(mld_clean[[1]])


# Resample to match extent and resolution
mld_stack_0.1 <- terra::resample(mld_clean, reference_raster_0.1, method = "bilinear", threads = TRUE)

print(mld_stack_0.1)
plot(mld_stack_0.1[[3]])

# writeRaster(mld_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1.tif", overwrite = TRUE)

writeRaster(mld_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1_ext.tif", overwrite = TRUE)



#  STATIC RASTERS ---------------------------------------------------------

gebco_nc <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GEBCO_2023//GEBCO_2023.nc"

bathy_raster <- terra::rast(gebco_nc)
bathy_raster

terra::crs(bathy_raster) <- "EPSG:4326"
plot(bathy_raster, range = c(-11000, 0))


# set values >0 to NA in deoth raster 
bathy_raster <- terra::ifel(bathy_raster > 0, NA_real_, bathy_raster)
plot(bathy_raster)
print(bathy_raster)


Slope <- terra::terrain(bathy_raster, v = "slope", unit = "degrees", neighbors = 8)
Rough <- terra::terrain(bathy_raster, v = "roughness", neighbors = 8)

names(bathy_raster) <- "Depth"
names(Slope) <- "Slope"
names(Rough) <- "Roughness"


# Resample to match extent and resolution
reference_raster_0.1 <- terra::rast("data/processed/reference_raster_0.1_ext.tif")
reference_raster_0.1
bathy_raster_0.1 <- terra::resample(bathy_raster, reference_raster_0.1, method = "bilinear", threads = TRUE)
slope_raster_0.1 <- terra::resample(Slope, reference_raster_0.1, method = "bilinear", threads = TRUE)
rough_raster_0.1 <- terra::resample(Rough, reference_raster_0.1, method = "bilinear", threads = TRUE)

plot(bathy_raster_0.1)
print(bathy_raster_0.1)

land_mass <- terra::vect("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Basemap_files//ne_10m_land.shp")
sf::st_crs(land_mass)
# Transform CRS  to match the CRS of the raster stack
land_mass_transf <- terra::project(land_mass, crs(bathy_raster_0.1))


# Mask out land masses
bathy_raster_0.1_mask <- terra::mask(bathy_raster_0.1, land_mass_transf, inverse = TRUE)
slope_raster_0.1_mask <- terra::mask(slope_raster_0.1, land_mass_transf, inverse = TRUE)
rough_raster_0.1_mask <- terra::mask(rough_raster_0.1, land_mass_transf, inverse = TRUE)
plot(bathy_raster_0.1_mask)
print(bathy_raster_0.1_mask)
plot(slope_raster_0.1_mask)
plot(rough_raster_0.1_mask)



## Replicate 
# monthly dates from 2010-01 to 2025-08
dates_monthly <- base::seq.Date(
  from = base::as.Date("2010-01-01"),
  to   = base::as.Date("2025-08-01"),
  by   = "month"
)

depth_stack_0.1 <- c(replicate(188, bathy_raster_0.1_mask)) |>  terra::rast() 
depth_stack_0.1
slope_stack_0.1 <- c(replicate(188, slope_raster_0.1_mask)) |>  terra::rast()
rough_stack_0.1 <- c(replicate(188, rough_raster_0.1_mask)) |>  terra::rast()



names(depth_stack_0.1) <- dates_monthly
names(slope_stack_0.1) <- dates_monthly
names(rough_stack_0.1) <- dates_monthly

varnames(depth_stack_0.1) <- "Depth"
varnames(slope_stack_0.1) <- "Slope"
varnames(rough_stack_0.1) <- "Roughness"

print(depth_stack_0.1)
print(rough_stack_0.1)
print(slope_stack_0.1)


# writeRaster(depth_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif", overwrite = TRUE)
# writeRaster(slope_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1.tif", overwrite = TRUE)
# writeRaster(rough_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1.tif", overwrite = TRUE)

writeRaster(depth_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1_ext.tif", overwrite = TRUE)
writeRaster(slope_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1_ext.tif", overwrite = TRUE)
writeRaster(rough_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1_ext.tif", overwrite = TRUE)


depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")
depth_stack_0.1

# Distance to Isobath  ----------------------------------------------------

# create contours
cont_200 <- terra::as.contour(bathy_raster_0.1_mask, levels = -200)
cont_1000 <- terra::as.contour(bathy_raster_0.1_mask, levels = -1000)
cont_2000 <- terra::as.contour(bathy_raster_0.1_mask, levels = -2000)


terra::plot(bathy_raster_0.1_mask, main = "Isobaths")
terra::plot(cont_200, col = "black", add = TRUE)
terra::plot(cont_1000, col = "red", add = TRUE)
terra::plot(cont_2000, col = "blue", add = TRUE)


## calculate distance between the two layers
dist_200_0.1 <- terra::distance(bathy_raster_0.1_mask, cont_200, unit = "m", method = "haversine") # haversine = great circle distance 
dist_1000_0.1 <- terra::distance(bathy_raster_0.1_mask, cont_1000, unit = "m", method = "haversine")
dist_2000_0.1 <- terra::distance(bathy_raster_0.1_mask, cont_2000, unit = "m", method = "haversine")


plot(dist_2000_0.1)

## mask to just within the marine area
dist_200_mask_0.1 <- terra::mask(dist_200_0.1, land_mass_transf, inverse = TRUE)
dist_1000_mask_0.1 <- terra::mask(dist_1000_0.1, land_mass_transf, inverse = TRUE)
dist_2000_mask_0.1 <- terra::mask(dist_2000_0.1, land_mass_transf, inverse = TRUE)


plot(dist_200_mask_0.1)
plot(dist_2000_mask_0.1)





# Replicate
dates_monthly <- base::seq.Date(
  from = base::as.Date("2010-01-01"),
  to   = base::as.Date("2025-08-01"),
  by   = "month"
)

dist200_stack_0.1 <- c(replicate(188, dist_200_mask_0.1)) |>  terra::rast() 
dist1000_stack_0.1 <- c(replicate(188, dist_1000_mask_0.1)) |>  terra::rast() 
dist2000_stack_0.1 <- c(replicate(188, dist_2000_mask_0.1)) |>  terra::rast() 

names(dist200_stack_0.1) <- dates_monthly
names(dist1000_stack_0.1) <- dates_monthly
names(dist2000_stack_0.1) <- dates_monthly

dist200_stack_0.1 <- dist200_stack_0.1 /1000
dist1000_stack_0.1 <- dist1000_stack_0.1 /1000
dist2000_stack_0.1 <- dist2000_stack_0.1 /1000

varnames(dist200_stack_0.1) <- "dist200"
varnames(dist1000_stack_0.1) <- "dist1000"
varnames(dist2000_stack_0.1) <- "dist2000"

units(dist200_stack_0.1) <- "km"
units(dist1000_stack_0.1) <- "km"
units(dist2000_stack_0.1) <- "km"



print(dist200_stack_0.1[[1]])
print(dist1000_stack_0.1)
print(dist2000_stack_0.1)

terra::plot(dist200_stack_0.1[[1]])
terra::plot(dist2000_stack_0.1[[1]])



writeRaster(dist200_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist200_month_0.1.tif", overwrite = TRUE)
writeRaster(dist1000_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist1000_month_0.1.tif", overwrite = TRUE)
writeRaster(dist2000_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1.tif", overwrite = TRUE)


writeRaster(dist200_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist200_month_0.1_ext.tif", overwrite = TRUE)
writeRaster(dist1000_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist1000_month_0.1_ext.tif", overwrite = TRUE)
writeRaster(dist2000_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1_ext.tif", overwrite = TRUE)



# Distance to Seamouints --------------------------------------------------

seamounts   <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Seamounts_data/Global_2011_ModelledSeamounts_ZSL/DownloadPack-14_001_ZSL002_ModelledSeamounts2011_v1/01_Data/seamounts/Seamounts.shp")            # seamount polygons or points
knolls      <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Seamounts_data/Global_2011_ModelledSeamounts_ZSL/DownloadPack-14_001_ZSL002_ModelledSeamounts2011_v1/01_Data/Knolls/Knolls.shp")               # knoll polygons or points

# Convert to same CRS (important!)
seamounts   <- sf::st_transform(seamounts, 4326)
knolls      <- sf::st_transform(knolls, 4326)


reference_raster_0.1 <- terra::rast("data/processed/reference_raster_0.1_ext.tif")
reference_raster_0.1

sea_ras_low <- terra::rasterize(vect(seamounts), reference_raster_0.1, field = 1)
sea_ras_low

dist_seamount_m <- terra::distance(sea_ras_low)

# Convert to km
dist_seamount_km <- dist_seamount_m / 1000
dist_seamount_km
plot(dist_seamount_km)

varnames(dist_seamount_km) <- "dist_seamount"
terra::units(dist_seamount_km) <- "km"


land_mass <- terra::vect("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Basemap_files//ne_10m_land.shp")
sf::st_crs(land_mass)
# Transform CRS  to match the CRS of the raster stack
land_mass_transf <- terra::project(land_mass, crs(reference_raster_0.1))


# Mask out land masses
dist_seamount_km_mask <- terra::mask(dist_seamount_km, land_mass_transf, inverse = TRUE)
plot(dist_seamount_km_mask)






sea_knolls_low <- terra::rasterize(vect(knolls), reference_raster_0.1, field = 1)
sea_knolls_low

dist_knolls_m <- terra::distance(sea_knolls_low)




# Convert to km
dist_knolls_km <- dist_knolls_m / 1000
dist_knolls_km
plot(dist_knolls_km)

varnames(dist_knolls_km) <- "dist_knoll"
terra::units(dist_knolls_km) <- "km"

# Mask out land masses
dist_knolls_km_mask <- terra::mask(dist_knolls_km, land_mass_transf, inverse = TRUE)
plot(dist_knolls_km_mask)



# Replicate
dates_monthly <- base::seq.Date(
  from = base::as.Date("2010-01-01"),
  to   = base::as.Date("2025-08-01"),
  by   = "month"
)

dist_seamounts_stack_0.1 <- c(replicate(188, dist_seamount_km_mask)) |>  terra::rast() 
dist_knolls_stack_0.1 <- c(replicate(188, dist_knolls_km_mask)) |>  terra::rast() 

names(dist_seamounts_stack_0.1) <- dates_monthly
names(dist_knolls_stack_0.1) <- dates_monthly



writeRaster(dist_seamounts_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_seamounts_month_0.1.tif", overwrite = TRUE)
writeRaster(dist_knolls_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_knolls_month_0.1.tif", overwrite = TRUE)







# Temporal rasters --------------------------------------------------------


month_stack <- depth_stack_0.1

values(month_stack) <- 0
varnames(month_stack) <- "Month"
month_stack


assign_month <- function(layer, layer_name) {
  month <- as.integer(substr(layer_name, 6, 7))  # Extract and convert month to integer
  layer[] <- month  # Assign month value to all cells
  return(layer)
}

# Apply the function to each layer in the stack
month_stack_0.1 <- app(month_stack, function(layer, layer_name) {
  assign_month(layer, layer_name)
}, layer_name = names(month_stack))


print(month_stack_0.1)
print(month_stack_0.1[[50]])
plot(month_stack_0.1[[2]])



year_stack <- depth_stack_0.1
values(year_stack) <- 0
varnames(year_stack) <- "Year"
year_stack


# Define the function to extract and assign the year value to each layer
assign_year <- function(layer, layer_name) {
  year <- as.integer(substr(layer_name, 1, 4))  
  layer[] <- year
  return(layer)
}

# Apply the function to each layer in the stack
year_stack_0.1 <- app(year_stack, function(layer, layer_name) {
  assign_year(layer, layer_name)
}, layer_name = names(year_stack))


year_stack_0.1



### LAT & LON raster

template_raster <- depth_stack_0.1[[1]]

# Create a raster for longitude and latitude
lon_rast <- terra::init(template_raster, fun  = "x")
lat_rast <- terra::init(template_raster, fun  = "y")
lon_rast
lat_rast
plot(lon_rast)
plot(lat_rast)

lon_stack_0.1 <- terra::rast(replicate(188, lon_rast))
names(lon_stack_0.1) <- dates_monthly
varnames(lon_stack_0.1) <- "lon"
lon_stack_0.1

lat_stack_0.1 <- terra::rast(replicate(188, lat_rast))
names(lat_stack_0.1) <- dates_monthly
varnames(lat_stack_0.1) <- "lat"
lat_stack_0.1


# writeRaster(month_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/month_stack_0.1.tif", overwrite = TRUE)
# 
# writeRaster(year_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/year_stack_0.1.tif", overwrite = TRUE)
# 
# writeRaster(lat_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lat_stack_0.1.tif", overwrite = TRUE)
# 
# writeRaster(lon_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lon_stack_0.1.tif", overwrite = TRUE)


writeRaster(month_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/month_stack_0.1_ext.tif", overwrite = TRUE)

writeRaster(year_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/year_stack_0.1_ext.tif", overwrite = TRUE)

writeRaster(lat_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lat_stack_0.1_ext.tif", overwrite = TRUE)

writeRaster(lon_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lon_stack_0.1_ext.tif", overwrite = TRUE)







# dummy raster for ID -----------------------------------------------------


id_rast <- depth_stack_0.1[[1]]
names(id_rast) <- "id"
id_rast
values(id_rast)[!is.na(values(id_rast))] = 0
plot(id_rast)
id_stack_0.1 <- terra::rast(replicate(188, id_rast))
names(id_stack_0.1) <- dates_monthly
varnames(id_stack_0.1) <- "id"
id_stack_0.1

writeRaster(id_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/id_stack_0.1.tif", overwrite = TRUE)

writeRaster(id_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/id_stack_0.1_ext.tif", overwrite = TRUE)





# Transformations ---------------------------------------------------------


depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")
slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1.tif")
rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1.tif")
dist2000_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1.tif")

mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1.tif")

wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1.tif")



# preserve negative values 
depth_cbrt_stack_0.1 <- terra::app(depth_stack_0.1, \(v) base::sign(v) * base::abs(v)^(1/3))
terra::units(depth_cbrt_stack_0.1) <- "cbrt(m)"

depth_stack_0.1 <- -depth_stack_0.1
plot(depth_stack_0.1)

depth_cbrt_stack_0.1


slope_cbrt_stack_0.1 <- (slope_stack_0.1)^(1/3)
terra::units(slope_cbrt_stack_0.1) <- "cbrt(degrees)"
slope_cbrt_stack_0.1


rough_cbrt_stack_0.1 <- rough_stack_0.1^(1/3)
terra::units(rough_cbrt_stack_0.1) <- "cbrt(rough)"
rough_cbrt_stack_0.1

dist2000_stack_0.1
dist2000_cbrt_stack_0.1 <- (dist2000_stack_0.1)^(1/3)
terra::units(dist2000_cbrt_stack_0.1) <- "cbrt(m)"
dist2000_cbrt_stack_0.1
dist2000_stack_0.1


uv_cbrt_stack_0.1 <- (uv_stack_0.1)^(1/3)
terra::units(uv_cbrt_stack_0.1) <- "cbrt(m s-1)"
uv_cbrt_stack_0.1
uv_stack_0.1


chl_log_stack_0.1 <- log(chl_stack_0.1)
terra::units(chl_log_stack_0.1) <- "log(milligram m-3)"
chl_stack_0.1
chl_log_stack_0.1

wz_cbrt_stack_0.1 <- terra::app(wz_stack_0.1, \(v) base::sign(v) * base::abs(v)^(1/3))
terra::units(wz_cbrt_stack_0.1) <- "cbrt(m s-1)"
wz_stack_0.1
wz_cbrt_stack_0.1



mld_log_stack_0.1 <-log(mld_stack_0.1)
terra::units(mld_log_stack_0.1) <- "log(m)"
mld_stack_0.1
mld_log_stack_0.1


writeRaster(depth_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_transf_month_0.1.tif", overwrite = TRUE)
writeRaster(slope_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_transf_month_0.1.tif", overwrite = TRUE)
writeRaster(rough_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(dist2000_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(uv_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(chl_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(wz_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(mld_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_transf_month_0.1.tif", overwrite = TRUE)








# Creating a mean raster stack for spatial blocking method (blockC --------


mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1.tif")

wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1.tif")

sst_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1.tif")

depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")

slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1.tif")

rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1.tif")

dist_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1.tif")


sst_mean <- terra::app(sst_stack_0.1, fun = mean, na.rm = TRUE)
mld_mean <- terra::app(mld_stack_0.1, fun = mean, na.rm = TRUE)
wz_mean <- terra::app(wz_stack_0.1, fun = mean, na.rm = TRUE)
uv_mean <- terra::app(uv_stack_0.1, fun = mean, na.rm = TRUE)
chl_mean <- terra::app(chl_stack_0.1, fun = mean, na.rm = TRUE)
chl_mean <- terra::mask(chl_mean, depth_rast)
plot(chl_mean)
depth_rast <- depth_stack_0.1[[1]]
slope_rast <- slope_stack_0.1[[1]]
rough_rast <- rough_stack_0.1[[1]]
dist_rast <- dist_stack_0.1[[1]]


names(sst_mean) <- "thetao"
names(mld_mean) <- "mltost"
names(wz_mean) <- "wz"
names(uv_mean) <- "uv"
names(chl_mean) <- "chl"
names(depth_rast) <- "depth"
names(slope_rast) <- "slope"
names(rough_rast) <- "roughness"
names(dist_rast) <- "dist2000"

# temporal dummy 

month <- depth_rast
names(month) <- "month"
month

base::set.seed(42L)

# Create a random month (1..12) raster from your template
month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() |>
  terra::as.factor() # form amchine learning 

month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() 

month


# Quick sanity checks
terra::freq(month)
month

mean_month_predistor_stack <- c(sst_mean,
                                chl_mean,
                                uv_mean,
                                mld_mean,
                                wz_mean,
                                depth_rast,
                                slope_rast,
                                rough_rast,
                                dist_rast,
                                month)

mean_month_predistor_stack


writeRaster(mean_month_predistor_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1.tif", overwrite = TRUE)


writeRaster(mean_month_predistor_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_gam.tif", overwrite = TRUE)






# creating a mean raster for GAMMs ----------------------------------------

# Creating a mean raster stack for spatial blocking method (blockC --------


mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_transf_month_0.1.tif")

wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_transf_month_0.1.tif")

sst_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_transf_month_0.1.tif")

depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_transf_month_0.1.tif")

slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_transf_month_0.1.tif")

rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_transf_month_0.1.tif")

dist_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_transf_month_0.1.tif")



sst_mean <- terra::app(sst_stack_0.1, fun = mean, na.rm = TRUE)
mld_mean <- terra::app(mld_stack_0.1, fun = mean, na.rm = TRUE)
wz_mean <- terra::app(wz_stack_0.1, fun = mean, na.rm = TRUE)
uv_mean <- terra::app(uv_stack_0.1, fun = mean, na.rm = TRUE)
depth_rast <- depth_stack_0.1[[1]]
chl_mean <- terra::app(chl_stack_0.1, fun = mean, na.rm = TRUE)
chl_mean <- terra::mask(chl_mean, depth_rast)
plot(chl_mean)

slope_rast <- slope_stack_0.1[[1]]
rough_rast <- rough_stack_0.1[[1]]
dist_rast <- dist_stack_0.1[[1]]


names(sst_mean) <- "thetao"
names(mld_mean) <- "mltost"
names(wz_mean) <- "wz"
names(uv_mean) <- "uv"
names(chl_mean) <- "chl"
names(depth_rast) <- "depth"
names(slope_rast) <- "slope"
names(rough_rast) <- "roughness"
names(dist_rast) <- "dist2000"

# temporal dummy 

month <- depth_rast
names(month) <- "month"
month

base::set.seed(42L)

# Create a random month (1..12) raster from your template

month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() 

month


# Quick sanity checks
terra::freq(month)
month

mean_month_predistor_stack <- c(sst_mean,
                                chl_mean,
                                uv_mean,
                                mld_mean,
                                wz_mean,
                                depth_rast,
                                slope_rast,
                                rough_rast,
                                dist_rast,
                                month)

mean_month_predistor_stack


writeRaster(mean_month_predistor_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_transf_month_predictor_stack_0.1_gam.tif", overwrite = TRUE)




# extended raster stacks to cover lrager area for simulated tracks  -------


mld_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_month_0.1_ext.tif")


wz_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_month_0.1_ext.tif")

chl_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_month_0.1_ext.tif")

sst_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/thetao_sst_month_0.1_ext.tif")

uv_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_month_0.1_ext.tif")

depth_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1_ext.tif")

slope_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_month_0.1_ext.tif")

rough_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_month_0.1_ext.tif")

dist200_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist200_month_0.1_ext.tif")

dist1000_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist1000_month_0.1_ext.tif")

dist2000_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_month_0.1_ext.tif")
dist2000_stack_0.1


lat_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lat_stack_0.1_ext.tif")

lon_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/lon_stack_0.1_ext.tif")

dist_seamounts_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_seamounts_month_0.1.tif")

dist_knolls_stack_0.1 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist_knolls_month_0.1.tif")

sst_mean <- terra::app(sst_stack_0.1, fun = mean, na.rm = TRUE)
sst_mean
mld_mean <- terra::app(mld_stack_0.1, fun = mean, na.rm = TRUE)
wz_mean <- terra::app(wz_stack_0.1, fun = mean, na.rm = TRUE)
uv_mean <- terra::app(uv_stack_0.1, fun = mean, na.rm = TRUE)
depth_rast <- depth_stack_0.1[[1]]
slope_rast <- slope_stack_0.1[[1]]
rough_rast <- rough_stack_0.1[[1]]
dist200_rast <- dist200_stack_0.1[[1]]
dist1000_rast <- dist1000_stack_0.1[[1]]
dist2000_rast <- dist2000_stack_0.1[[1]]
lat_rast <- lat_stack_0.1[[1]]
lon_rast <- lon_stack_0.1[[1]]
dist_seamounts <- dist_seamounts_stack_0.1[[1]]
dist_knolls <- dist_knolls_stack_0.1[[1]]

chl_mean <- terra::app(chl_stack_0.1, fun = mean, na.rm = TRUE)
chl_mean <- terra::mask(chl_mean, depth_rast)
plot(chl_mean)

id_rast <- depth_stack_0.1[[1]]
names(id_rast) <- "id"
id_rast
values(id_rast)[!is.na(values(id_rast))] = 0
plot(id_rast)




names(sst_mean) <- "thetao"
names(mld_mean) <- "mltost"
names(wz_mean) <- "wz"
names(uv_mean) <- "uv"
names(chl_mean) <- "chl"
names(depth_rast) <- "depth"
names(slope_rast) <- "slope"
names(rough_rast) <- "roughness"
names(dist200_rast) <- "dist200"
names(dist1000_rast) <- "dist1000"
names(dist2000_rast) <- "dist2000"
names(id_rast) <- "id"
names(lat_rast) <- "lat"
names(lon_rast) <- "lon"
names(dist_seamounts) <- "dist_seamount"
names(dist_knolls) <- "dist_knoll"

# temporal dummy 

month <- depth_rast
names(month) <- "month"
month

base::set.seed(42)

# Create a random month (1..12) raster from your template
month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() |>
  terra::as.factor() # form amchine learning 

month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() 

month






mean_month_predistor_stack <- c(sst_mean,
                                chl_mean,
                                uv_mean,
                                mld_mean,
                                wz_mean,
                                depth_rast,
                                slope_rast,
                                rough_rast,
                                dist200_rast,
                                dist1000_rast,
                                dist2000_rast,
                                month,
                                id_rast,
                                lat_rast,
                                lon_rast,
                                dist_seamounts,
                                dist_knolls)

mean_month_predistor_stack


writeRaster(mean_month_predistor_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext.tif", overwrite = TRUE)

writeRaster(mean_month_predistor_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_ML.tif", overwrite = TRUE)





## Transformations (use log1p)


# depth_stack_0.1 <- -depth_stack_0.1
# plot(depth_stack_0.1[[1]])
# depth_log_stack_0.1 <- log1p(depth_stack_0.1)
# print(depth_log_stack_0.1)
# plot(depth_log_stack_0.1[[1]])
# terra::units(depth_log_stack_0.1) <- "log1p(m)"


# slope_log_stack_0.1 <- log1p(slope_stack_0.1)
# terra::units(slope_log_stack_0.1) <- "log1p(degrees)"
# slope_log_stack_0.1

# 
# rough_log_stack_0.1 <- log1p(rough_stack_0.1)
# terra::units(rough_log_stack_0.1) <- "log1p(roughness)"
# rough_log_stack_0.1



# dist_log_stack_0.1 <- log1p(dist_stack_0.1)
# terra::units(dist_log_stack_0.1) <- "log1p(m)"
# dist_log_stack_0.1

dist200_stack_0.1



# uv_log_stack_0.1 <- log1p(uv_stack_0.1)
# terra::units(uv_log_stack_0.1) <- "log1p(m s-1)"
# uv_log_stack_0.1
# uv_stack_0.1


chl_log_stack_0.1 <- log(chl_stack_0.1)
terra::units(chl_log_stack_0.1) <- "log(milligram m-3)"
chl_stack_0.1
chl_log_stack_0.1

# wz_cbrt_stack_0.1 <- terra::app(wz_stack_0.1, \(v) base::sign(v) * base::abs(v)^(1/3))
# terra::units(wz_cbrt_stack_0.1) <- "cbrt(m s-1)"
# wz_stack_0.1
# wz_cbrt_stack_0.1



# mld_log_stack_0.1 <-log(mld_stack_0.1)
# terra::units(mld_log_stack_0.1) <- "log(m)"
# mld_stack_0.1
# mld_log_stack_0.1


writeRaster(depth_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_transf_month_0.1.tif", overwrite = TRUE)
writeRaster(slope_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/slope_transf_month_0.1.tif", overwrite = TRUE)
writeRaster(rough_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/rough_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(dist_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/dist2000_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(uv_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/uv_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(chl_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/chl_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(wz_cbrt_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/Wz_transf_month_0.1.tif", overwrite = TRUE)

writeRaster(mld_log_stack_0.1, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mltost_mld_transf_month_0.1.tif", overwrite = TRUE)



## Mean raster transformed and extended
sst_mean <- terra::app(sst_stack_0.1, fun = mean, na.rm = TRUE)
# mld_mean <- terra::app(mld_log_stack_0.1, fun = mean, na.rm = TRUE)
mld_mean <- terra::app(mld_stack_0.1, fun = mean, na.rm = TRUE)
wz_mean <- terra::app(wz_stack_0.1, fun = mean, na.rm = TRUE)
# uv_mean <- terra::app(uv_log_stack_0.1, fun = mean, na.rm = TRUE)
uv_mean <- terra::app(uv_stack_0.1, fun = mean, na.rm = TRUE)
# depth_rast <- depth_log_stack_0.1[[1]]
depth_rast <- depth_stack_0.1[[1]]
# slope_rast <- slope_log_stack_0.1[[1]]
# rough_rast <- rough_log_stack_0.1[[1]]
slope_rast <- slope_stack_0.1[[1]]
rough_rast <- rough_stack_0.1[[1]]
dist200_rast <- dist200_stack_0.1[[1]]
dist1000_rast <- dist1000_stack_0.1[[1]]
dist2000_rast <- dist2000_stack_0.1[[1]]
dist_seamounts_rast <- dist_seamounts_stack_0.1[[1]]
dist_knolls_rast <- dist_knolls_stack_0.1[[1]]


chl_mean <- terra::app(chl_log_stack_0.1, fun = mean, na.rm = TRUE)
chl_mean <- terra::mask(chl_mean, depth_rast)
plot(chl_mean)



names(sst_mean) <- "thetao"
names(mld_mean) <- "mltost"
names(wz_mean) <- "wz"
names(uv_mean) <- "uv"
names(chl_mean) <- "chl"
names(depth_rast) <- "depth"
names(slope_rast) <- "slope"
names(rough_rast) <- "roughness"
names(dist200_rast) <- "dist200"
names(dist1000_rast) <- "dist1000"
names(dist2000_rast) <- "dist2000"
names(id_rast) <- "id"
names(dist_seamounts_rast) <- "dist_seamount"
names(dist_knolls_rast) <- "dist_knoll"

# temporal dummy 

month <- depth_rast
names(month) <- "month"
month

base::set.seed(42)

# Create a random month (1..12) raster from your template
month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() |>
  terra::as.factor() # form amchine learning 

plot(month)

month <- month |>
  terra::rast() |>
  (\(r) {
    vals <- base::sample.int(12L, size = terra::ncell(r), replace = TRUE)
    terra::setValues(r, vals)
  })() 

month


sst_mean
chl_mean
uv_mean
mld_mean
wz_mean
depth_rast
slope_rast
rough_rast
dist200_rast
dist1000_rast
dist2000_rast
month
id_rast
lat_rast
lon_rast

mean_month_predistor_stack_log <- c(sst_mean,
                                chl_mean,
                                uv_mean,
                                mld_mean,
                                wz_mean,
                                depth_rast,
                                slope_rast,
                                rough_rast,
                                dist200_rast,
                                dist1000_rast,
                                dist2000_rast,
                                month,
                                id_rast,
                                lat_rast,
                                lon_rast,
                                dist_seamounts_rast,
                                dist_knolls_rast)

mean_month_predistor_stack_log


writeRaster(mean_month_predistor_stack_log, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log.tif", overwrite = TRUE)

writeRaster(mean_month_predistor_stack_log, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log_ML.tif", overwrite = TRUE)



