#_____________________________________________________________________________
#                        Extractions: Copernicus Marine (CMEMS)
#                               - chlorophyll -
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
nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/daily/chl/"


# List all NetCDF files in the directory
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files





## Create a raster stack containing all the daily Chlorophyll data

chl_stack <- terra::rast(nc_files, subds = "CHL")
chl_stack

# add dates as names from time stamps
names(chl_stack) <- terra::time(chl_stack)

# add missing CRS information 
terra::crs(chl_stack) <- "EPSG:4326"
chl_stack

plot(mld_stack[[168]])



# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst_uv.curr_mld.rds")

sight <- readRDS("data/work_files/Sightings_Validation_data_bathy_dist_sst_uv_mld.rds")


tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst_uv.curr_mld.rds")

tracks <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld.rds")
tracks <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld.rds")

tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist_sst_uv_mld.rds")


tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist_sst_uv_mld.rds")


tracks <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist_sst_uv_mld.rds")

dt <- sight
dt <- tracks

str(dt)

dt.sf <- dt |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  dplyr::mutate(Date = base::as.Date(date))



str(dt.sf)



# Extractions -------------------------------------------------------------

# quick plot to make sure everyhting is aligned!
dev.off()
v <- terra::vect(dt.sf)
v
bb <- terra::ext(v)
buff <- 2  # ~2° buffer added

terra::plot(
  chl_stack[[168]],
  xlim = c(terra::xmin(bb) - buff, terra::xmax(bb) + buff),
  ylim = c(terra::ymin(bb) - buff, terra::ymax(bb) + buff),
  range = c(0, 2)
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")


str(dt.sf)


dt_chl <- CMEMS_extract_by_date_surf(chl_stack, dt.sf, date_col = "Date", varname = "chl", buffer_m = 25000)

glimpse(dt_chl)


dt_chl |> 
  as.data.frame() |> 
  rstatix::get_summary_stats(chl, type = "common")


# check NAs
dt_chl |>
  dplyr::filter(is.na(chl)) |> 
  dplyr::select(date, lat, lon, chl) |> 
  print(n=50)


# Save  -------------------------------------------------------------------
sight_chl <- dt_chl
tracks_chl <- dt_chl

mapview::mapview(tracks_chl |> dplyr::filter(PA ==1) |>  dplyr::select(-Date))


saveRDS(sight_chl, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

saveRDS(sight_chl, "data/work_files/Sightings_Validation_data_bathy_dist_sst_uv_mld_chl.rds")

saveRDS(tracks_chl, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

saveRDS(tracks_chl, "data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
saveRDS(tracks_chl, "data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

saveRDS(tracks_chl, "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist_sst_uv_mld_chl.rds")




saveRDS(tracks_chl, "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist_sst_uv_mld_chl.rds")


saveRDS(tracks_chl, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist_sst_uv_mld_chl.rds")


