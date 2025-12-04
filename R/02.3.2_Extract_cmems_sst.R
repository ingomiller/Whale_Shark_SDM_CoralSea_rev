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

sst_stack <- terra::rast(nc_files, subds = "thetao")
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

tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist.rds")

tracks <- readRDS( "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_final_extA.rds")
tracks <- readRDS( "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_final_extA_5d.rds")
tracks <- readRDS( "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist.rds")


tracks <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist.rds")

dt <- sight
dt <- tracks
dt.sf <- dt |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  dplyr::mutate(Date = base::as.Date(date))

# fix date issue 
dt.sf <- dt.sf |> dplyr::mutate(Date = dplyr::case_when(
  Date == as.Date("2024-02-29") ~ as.Date("2024-03-01"),
  TRUE                          ~ Date)
)

str(dt.sf)

min(dt.sf$Date)
max(dt.sf$Date)


stack_dates <- terra::time(sst_stack) |>
  base::as.Date(tz = "UTC")

# 2) Dates requested by your sf
dt_dates <- base::as.Date(dt.sf$Date, tz = "UTC")

# 3) Which requested dates are missing from the stack?
missing_mask <- !(dt_dates %in% stack_dates)

# 4) Quick summary of missing dates (which and how many rows each)
missing_dates_tbl <- dt.sf |>
  dplyr::mutate(Date_UTC = base::as.Date(.data$Date)) |>
  dplyr::filter(!( .data$Date_UTC %in% stack_dates )) |>
  dplyr::count(.data$Date_UTC, sort = TRUE)

missing_dates_tbl


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
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(thetao, type = "common")

dt_thetao |> 
  as.data.frame() |> 
  dplyr::filter(thetao < 25)


## extending 2 degrees not much impact 


# check NAs
dt_thetao |>
  dplyr::filter(is.na(thetao)) |> 
  dplyr::select(date, lat, lon, thetao, Depth) |> 
  print(n=50)

dt_thetao |>
  ggplot2::ggplot(ggplot2::aes(x = thetao, fill = factor(PA))) +
  ggplot2::geom_density(alpha = 0.3, position = "identity", trim = FALSE)

# Save  -------------------------------------------------------------------
sight_thetao <- dt_thetao
tracks_thetao <- dt_thetao

saveRDS(sight_thetao, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst.rds")

saveRDS(sight_thetao, "data/work_files/Sightings_Validation_data_bathy_dist_sst.rds")


saveRDS(tracks_thetao, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst.rds")


saveRDS(tracks_thetao, "data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst.rds")
saveRDS(tracks_thetao, "data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst.rds")

saveRDS(tracks_thetao, "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy_dist_sst.rds")

saveRDS(tracks_thetao, "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy_dist_sst.rds")


saveRDS(tracks_thetao, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy_dist_sst.rds")













# SST gradient - Fronts  --------------------------------------------------



ext_reg <- terra::ext(135, 175, -45, 5)  # add a few degrees buffer around your true region

sst_reg <- sst_stack |>
  terra::crop(ext_reg)



sst_reg


t_sst <- sst_reg |> 
  terra::time()

idx <- (as.Date(t_sst) >= as.Date("2019-01-01")) &
  (as.Date(t_sst) <= as.Date("2025-06-30"))

# 3) Subset the raster by those layers
sst_reg <- sst_reg[[idx]]

w <- matrix(1, nrow = 3, ncol = 3)

# 1) 3x3 moving max and min of SST
sst_max <- sst_reg |>
  terra::focal(
    w         = w,
    fun       = max,
    na.policy = "omit"
    # filename  = "sst_max_3x3.tif",
    # wopt      = list(datatype = "FLT4S", overwrite = TRUE)
  )

sst_min <- sst_reg |>
  terra::focal(
    w         = w,
    fun       = min,
    na.policy = "omit"
    # filename  = "sst_min_3x3.tif",
    # wopt      = list(datatype = "FLT4S", overwrite = TRUE)
  )

# 2) Difference centre -> neighbours
diff_max <- sst_max - sst_reg
diff_min <- sst_reg - sst_min

# 3) Max absolute difference (this is Δ°C per cell, at 0.083°)
sst_slope_max <- max(diff_max, diff_min)

varnames(sst_slope_max) <- "sst_slope"

sst_slope_max[[360]]
plot(sst_slope_max[[150]])



tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")


dt <- tracks
dt.sf <- dt |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  dplyr::mutate(Date = base::as.Date(date))

# fix date issue 
dt.sf <- dt.sf |> dplyr::mutate(Date = dplyr::case_when(
  Date == as.Date("2024-02-29") ~ as.Date("2024-03-01"),
  TRUE                          ~ Date)
)

str(dt.sf)

min(dt.sf$Date)
max(dt.sf$Date)


stack_dates <- terra::time(sst_slope_max) |>
  base::as.Date(tz = "UTC")

# 2) Dates requested by your sf
dt_dates <- base::as.Date(dt.sf$Date, tz = "UTC")

# 3) Which requested dates are missing from the stack?
missing_mask <- !(dt_dates %in% stack_dates)

# 4) Quick summary of missing dates (which and how many rows each)
missing_dates_tbl <- dt.sf |>
  dplyr::mutate(Date_UTC = base::as.Date(.data$Date)) |>
  dplyr::filter(!( .data$Date_UTC %in% stack_dates )) |>
  dplyr::count(.data$Date_UTC, sort = TRUE)

missing_dates_tbl


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


dt_sst_slope <- CMEMS_extract_by_date_surf(sst_slope_max, dt.sf, date_col = "Date", varname = "sst_slope", buffer_m = 25000)

glimpse(dt_sst_slope)


dt_sst_slope |> 
  as.data.frame() |> 
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(sst_slope, type = "common")


## extending 2 degrees not much impact 


# check NAs
dt_sst_slope |>
  dplyr::filter(is.na(thetao)) |> 
  dplyr::select(date, lat, lon, sst_slope, depth) |> 
  print(n=50)

dt_sst_slope |>
  ggplot2::ggplot(ggplot2::aes(x = sst_slope, fill = factor(PA))) +
  ggplot2::geom_density(alpha = 0.3, position = "identity", trim = FALSE)

saveRDS(dt_sst_slope, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")



