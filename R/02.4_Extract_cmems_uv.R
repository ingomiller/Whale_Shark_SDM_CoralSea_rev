#_____________________________________________________________________________
#                        Extractions: Copernicus Marine (CMEMS)
#                               - current velocity -
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

## Create a raster stack containing all the daily currents data
curr_stack_u <- terra::rast(nc_files, subds = "uo")
curr_stack_v <- terra::rast(nc_files, subds = "vo")

# add dates as names from time stamps
names(curr_stack_u) <- terra::time(curr_stack_u)
names(curr_stack_v) <- terra::time(curr_stack_v)

# add missing CRS information 
terra::crs(curr_stack_u) <- "EPSG:4326"
terra::crs(curr_stack_v) <- "EPSG:4326"
# terra::crs(curr_stack_v) <- "OGC:CRS84"

print(curr_stack_u)
print(curr_stack_v)


# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst.rds")

sight <- readRDS("data/work_files/Sightings_Validation_data_bathy_dist_sst.rds")


tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst.rds")

tracks <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst.rds")
tracks <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst.rds")

tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist_sst.rds")

tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist_sst.rds")


tracks <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist_sst.rds")


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
  curr_stack_u[[168]],
  xlim = c(terra::xmin(bb) - buff, terra::xmax(bb) + buff),
  ylim = c(terra::ymin(bb) - buff, terra::ymax(bb) + buff)
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")




str(dt.sf)


dt_uo <- CMEMS_extract_by_date_surf(curr_stack_u, dt.sf, date_col = "Date", varname = "uo", buffer_m = 25000)
dt_uo

dt_uo_vo <- CMEMS_extract_by_date_surf(curr_stack_v, dt_uo, date_col = "Date", varname = "vo", buffer_m = 25000)

glimpse(dt_uo_vo)




# Calculations  -----------------------------------------------------------


dt_uv <- dt_uo_vo |> 
  dplyr::mutate(
    uv = sqrt(uo^2 + vo^2), # overall velocity
    curr_dir = atan2(uo, vo) * 180 / pi, # adding current direction 
    curr_dir = ifelse(curr_dir < 0, curr_dir + 360, curr_dir) #ensuring angles are in 360
  )


str(dt_uv)

dt_uv |> 
  as.data.frame() |> 
  rstatix::get_summary_stats(uv, type = "common")

# check NAs
dt_uv |>
  dplyr::filter(is.na(uv)) |> 
  dplyr::select(date, lat, lon, uv) |> 
  print(n=50)


# Save  -------------------------------------------------------------------

sight_uv <- dt_uv
track_uv <- dt_uv

saveRDS(sight_uv, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst_uv.curr_.rds")


saveRDS(sight_uv, "data/work_files/Sightings_Validation_data_bathy_dist_sst_uv.rds")

saveRDS(track_uv, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst_uv.curr_.rds")


saveRDS(track_uv, "data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr.rds")
saveRDS(track_uv, "data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr.rds")



saveRDS(track_uv, "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist_sst_uv.rds")



saveRDS(track_uv, "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist_sst_uv.rds")


saveRDS(track_uv, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist_sst_uv.rds")
