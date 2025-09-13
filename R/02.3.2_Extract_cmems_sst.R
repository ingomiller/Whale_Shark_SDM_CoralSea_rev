#_____________________________________________________________________________
#                        Extractions: Copernicus Marine (CMEMS)
#                               - sst -
#_____________________________________________________________________________
### *Data downlaaded usoing Copernicus Marine Toolbox via Python

library(raster)
library(sf)
library(sp)
library(doParallel)
library(foreach)
library(terra)
source("R/00_Helper_Functions.R") 





# Get Raster Data from CMEMS ----------------------------------------------

# Set the directory containing the NetCDF files
nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/daily/mld.thetao.uo.vo.so/"


# List all NetCDF files in the directory
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

## Create a raster stack containing all the daily Chlorophyll data

sst_stack <- terra::rast(nc_files[6:22], subds = "thetao")
sst_stack

# add dates as names from time stamps
names(sst_stack) <- terra::time(sst_stack)

# add missing CRS information 
terra::crs(sst_stack) <- "EPSG:4326"

plot(sst_stack[[168]])



# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist.rds")

sight <- readRDS("data/work_files/Sightings_Validation_data_bathy_dist.rds")

tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist.rds")

tracks <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist.rds")
tracks <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist.rds")


dt <- sight
dt <- tracks
dt.sf <- dt |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  dplyr::mutate(Date = base::as.Date(date))

str(dt.sf)


# Extractions -------------------------------------------------------------

# quick plot to make sure everyhting is aligned!
dev.off()
v <- terra::vect(dt.sf)
bb <- terra::ext(v)
buff <- 2  # ~2° buffer added

terra::plot(
  sst_stack[[168]],
  xlim = c(terra::xmin(bb) - buff, terra::xmax(bb) + buff),
  ylim = c(terra::ymin(bb) - buff, terra::ymax(bb) + buff)
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")


dt_thetao <- CMEMS_extract_by_date_surf(sst_stack, dt.sf, date_col = "Date", varname = "thetao", buffer_m = 25000)

glimpse(dt_thetao)


dt_thetao |> 
  as.data.frame() |> 
  rstatix::get_summary_stats(thetao, type = "common")

# check NAs
dt_thetao |>
  dplyr::filter(is.na(thetao)) |> 
  dplyr::select(date, lat, lon, thetao) |> 
  print(n=50)

# Save  -------------------------------------------------------------------
sight_thetao <- dt_thetao
tracks_thetao <- dt_thetao

saveRDS(sight_thetao, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst.rds")

saveRDS(sight_thetao, "data/work_files/Sightings_Validation_data_bathy_dist_sst.rds")


saveRDS(tracks_thetao, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst.rds")


saveRDS(tracks_thetao, "data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst.rds")
saveRDS(tracks_thetao, "data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst.rds")



