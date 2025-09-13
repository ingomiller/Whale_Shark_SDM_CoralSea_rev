#_____________________________________________________________________________
#                        Extractions: Distance to Drop-off
#_____________________________________________________________________________

## to remove old rast temp files: terra::tmpFiles(current = FALSE, old = TRUE, remove = TRUE)


library(tidyverse)
library(raster)
library(terra)
library(ncdf4)
library(sf)
library(fasterRaster)
grassDir <- "/Applications/GRASS-8.4.app/Contents/Resources"
fasterRaster::faster(grassDir = grassDir)

stopifnot(dir.exists(grassDir))
fasterRaster::faster(grassDir = grassDir, verbose = TRUE, debug = TRUE)

#Some Grass system checks
Sys.setenv(GISBASE = grassDir)
Sys.setenv(GDAL_DATA = file.path(grassDir, "share", "gdal"))
Sys.setenv(GRASS_PROJSHARE = file.path(grassDir, "share", "proj"))
Sys.setenv(
  PATH = paste(
    file.path(grassDir, "bin"),
    file.path(grassDir, "scripts"),
    Sys.getenv("PATH"),
    sep = ":"
  )
)
Sys.which(c("g.gisenv","g.version"))
system2("g.version", "-g", stdout = TRUE)






# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy.rds")

sight <- readRDS("data/work_files/Sightings_Validation_data_bathy.rds")
# str(sight)

tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy.rds")
tracks <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy.rds")
tracks <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds") 

tracks |>  dplyr::summarise(min = min(date),
                            max = max(date))

str(tracks)

# Load Bathy/Geography Files ---------------------------------------------------

gebco_nc <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GEBCO_2023//GEBCO_2023.nc"

bathy <- terra::rast(gebco_nc)

land_mass <- terra::vect("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Basemap_files//ne_10m_land.shp")
sf::st_crs(land_mass)

# Transform CRS  to match the CRS of the raster stack
land_mass_transf <- project(land_mass, crs(bathy))

plot(bathy, range = c(-11000, 0))
plot(land_mass_transf)

tracks |> 
  sf::st_drop_geometry() |>
  rstatix::get_summary_stats(lat, lon)


# File transformations/Preperations ---------------------------------------

# Transform CRS  to match the CRS of the raster stack
land_mass_transf <- project(land_mass, crs(bathy))

plot(bathy, range = c(-11000, 0))
plot(land_mass_transf)


## crop raseter to save moemory:
# Define the extent of the region of interest
aoi_extent <- ext(110, 180, -50, 50)  
aoi_extent


land <- crop(land_mass_transf, aoi_extent)

bathy_crop <- crop(bathy, aoi_extent)
print(bathy_crop)

plot(bathy_crop, axes = TRUE, box = FALSE, range = c(-11000, 0))
plot(land, add = TRUE, border = "black")



## transform to fasterRaster object 
terra::crs(bathy_crop) <- "EPSG:4326"
bathy_crop <- fasterRaster::fast(bathy_crop)
land <- fasterRaster::fast(land)


## Mask out land
land_bathy <- fasterRaster::mask(bathy_crop, land, inverse = FALSE)
ocean_bathy <- fasterRaster::mask(bathy_crop, land_bathy, inverse = TRUE)

plot(ocean_bathy, axes = TRUE, box = FALSE)
plot(land, add = TRUE, border = "black")



# Get Contours & Distances------------------------------------------------------------
cont_200 <- fasterRaster::as.contour(ocean_bathy, levels = c(-200))
plot(ocean_bathy, main  = "200m isobath")
plot(cont_200, col = "blue", add = TRUE)


dist_to_cont200 <- fasterRaster::distance(ocean_bathy, cont_200, unit = "m")


## mask to just within the marine area
dist_to_cont_marine200 <- mask(dist_to_cont200, ocean_bathy)

dist_to_cont_marine200 <- terra::rast(dist_to_cont_marine200)
terra::crs(dist_to_cont_marine200) <- "EPSG:4326"
dist_to_cont_marine200

plot(dist_to_cont_marine200)




cont_1000 <- fasterRaster::as.contour(ocean_bathy, levels = c(-1000))
plot(ocean_bathy, main  = "1000m isobath")
plot(cont_1000, col = "hotpink", add = TRUE)



dist_to_cont1000 <- fasterRaster::distance(ocean_bathy, cont_1000, unit = "m")


## mask to just within the marine area
dist_to_cont_marine1000 <- mask(dist_to_cont1000, ocean_bathy)

dist_to_cont_marine1000 <- terra::rast(dist_to_cont_marine1000)
terra::crs(dist_to_cont_marine1000) <- "EPSG:4326"
dist_to_cont_marine1000

plot(dist_to_cont_marine1000)







cont_2000 <- fasterRaster::as.contour(ocean_bathy, levels = c(-2000))
plot(ocean_bathy, main  = "2000m isobath")
plot(cont_2000, col = "red", add = TRUE)



dist_to_cont2000 <- fasterRaster::distance(ocean_bathy, cont_2000, unit = "m")


## mask to just within the marine area
dist_to_cont_marine2000 <- mask(dist_to_cont2000, ocean_bathy)

dist_to_cont_marine2000 <- terra::rast(dist_to_cont_marine2000)
terra::crs(dist_to_cont_marine2000) <- "EPSG:4326"
dist_to_cont_marine2000

plot(dist_to_cont_marine2000)




plot(ocean_bathy, main  = "Isobath")
plot(cont_2000, col = "red", add = TRUE)
plot(cont_1000, col = "black", add = TRUE)
plot(cont_200, col = "blue", add = TRUE)





terra::writeRaster(dist_to_cont_marine200, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance200_for_extractions.tif", overwrite = TRUE)


terra::writeRaster(dist_to_cont_marine1000, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance1000_for_extractions.tif", overwrite = TRUE)


terra::writeRaster(dist_to_cont_marine2000, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance2000_for_extractions.tif", overwrite = TRUE)

# Extraction --------------------------------------------------------------

dist_to_cont_marine2000 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance2000_for_extractions.tif")
dist_to_cont_marine1000 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance1000_for_extractions.tif")
dist_to_cont_marine200 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance200_for_extractions.tif")


plot(dist_to_cont_marine2000)

dt <- sight
dt <- tracks |> as.data.frame()

dt_v <- vect(dt, geom = c("lon", "lat"), crs = crs(dist_to_cont_marine2000))

dt$dist200 <- terra::extract(dist_to_cont_marine200, dt_v)
dt$dist1000 <- terra::extract(dist_to_cont_marine1000, dt_v)
dt$dist2000 <- terra::extract(dist_to_cont_marine2000, dt_v)



dt_dist <- dt |>
  dplyr::mutate(dist2000 = dist2000$distance,
                dist1000 = dist1000$distance,
                dist200 = dist200$distance)

glimpse(dt_dist)

sight_dist <- dt_dist
tracks_dist <- dt_dist

# Save  -------------------------------------------------------------------

tracks_dist |> dplyr::summarise(min = min(date),
                                max = max(date))

saveRDS(sight, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist.rds")

saveRDS(sight_dist, "data/work_files/Sightings_Validation_data_bathy_dist.rds")

saveRDS(tracks_dist, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist.rds")

saveRDS(tracks_dist, "data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist.rds")
saveRDS(tracks_dist, "data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist.rds")

saveRDS(tracks_dist, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds") 



# Distance to Seamounts ---------------------------------------------------

seamounts   <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Seamounts_data/Global_2011_ModelledSeamounts_ZSL/DownloadPack-14_001_ZSL002_ModelledSeamounts2011_v1/01_Data/seamounts/Seamounts.shp")            # seamount polygons or points
knolls      <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Seamounts_data/Global_2011_ModelledSeamounts_ZSL/DownloadPack-14_001_ZSL002_ModelledSeamounts2011_v1/01_Data/Knolls/Knolls.shp")               # knoll polygons or points

tracks_ready <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
str(tracks_ready)

# Convert to same CRS (important!)
seamounts   <- sf::st_transform(seamounts, 4326)
knolls      <- sf::st_transform(knolls, 4326)


base <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Bathymetry_Distance2000_for_extractions.tif")   # template raster, CRS 4326


# define bounding box around your occurrences
ext_occ <- terra::ext(vect(tracks_ready))
ext_occ
# ext_buffered <- terra::extend(ext_occ, c(5,5,5,5))  # add ~5° buffer
# ext_buffered

# crop base raster
# base_crop <- terra::crop(base, ext_buffered)
# base_crop

base_lowres <- terra::aggregate(base, fact = 10)  # 4× coarser (~0.04°)
base_lowres
sea_ras_low <- terra::rasterize(vect(seamounts), base_lowres, field = 1)

crs_proj <- "EPSG:6933"  # World Cylindrical Equal Area
# or something regional like "ESRI:102033" (Asia South Equidistant Conic)

# Reproject your rasterized seamount layer
sea_ras_proj <- terra::project(sea_ras_low, crs_proj)

# Compute distances in metres (since projection units = metres)
dist_seamount_m <- terra::distance(sea_ras_proj)

# Convert to km
dist_seamount_km <- dist_seamount_m / 1000
plot(dist_seamount_km)

dist_seamount_km_4326 <- terra::project(dist_seamount_km, "EPSG:4326", method = "bilinear")
plot(dist_seamount_km_4326)
dist_seamount_km_4326

terra::writeRaster(dist_seamount_km_4326, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Seamounts_Distance.tif", overwrite = TRUE)

tracks_ready$dist_seamount <- terra::extract(dist_seamount_km_4326, vect(tracks_ready))[,2]
str(tracks_ready)



knolls_ras_low <- terra::rasterize(vect(knolls), base_lowres, field = 1)

crs_proj <- "EPSG:6933"  # World Cylindrical Equal Area
# or something regional like "ESRI:102033" (Asia South Equidistant Conic)

# Reproject your rasterized seamount layer
knoll_ras_proj <- terra::project(knolls_ras_low, crs_proj)

# Compute distances in metres (since projection units = metres)
dist_knoll_m <- terra::distance(knoll_ras_proj)

# Convert to km
dist_knoll_km <- dist_knoll_m / 1000
plot(dist_knoll_km)

dist_knoll_km_4326 <- terra::project(dist_knoll_km, "EPSG:4326", method = "bilinear")
plot(dist_knoll_km_4326)
dist_knoll_km_4326

terra::writeRaster(dist_knoll_km_4326, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/Knolls_Distance.tif", overwrite = TRUE)

tracks_ready$dist_knoll <- terra::extract(dist_knoll_km_4326, vect(tracks_ready))[,2]
str(tracks_ready)



saveRDS(tracks_ready, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_seamounts.rds")





