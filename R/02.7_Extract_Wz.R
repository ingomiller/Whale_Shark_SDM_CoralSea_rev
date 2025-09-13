#_____________________________________________________________________________
#                   Extractions: Bluelink BRAN2020 & CMEMS 
#                       - vertical current velocity Wz -
#_____________________________________________________________________________library(raster)



library(sf)
library(sp)
library(doParallel)
library(foreach)
library(terra)
library(remora)
source("R/00_Helper_Functions.R") 




# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

sight <- readRDS("data/work_files/Sightings_Validation_data_bathy_dist_sst_uv_mld_chl.rds")

tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

# tracks <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
# tracks <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")

# after thinning and already all other variables extracted
tracks <- readRDS("data/work_files/Tracks_sims_50_thinned_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
tracks <- readRDS("data/work_files/Tracks_mp_sims_30_thinned_2018_2025_bathy_dist_sst_uv.curr_mld_chl.rds")





dt <- sight
dt <- tracks

dt.sf <- dt |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  dplyr::mutate(Date = base::as.Date(date))


mapview::mapview(dt.sf |> dplyr::filter(PA ==1) |>  dplyr::select(-Date))









# Bluelink BRAN2020 for data until Dec 2023 -------------------------------

max(dt$date)

str(dt.sf)

input_df <- dt.sf |> as.data.frame() |> 
  dplyr::filter(Date <= as.Date("2023-12-31"))
  # dplyr::slice(1:100)
str(input_df)

max(input_df$Date)



# rm(input.lst)

tictoc::tic("Bluelink extraction took: ")
res <- extractWz(
  df = input_df,
  X = "lon",
  Y = "lat",
  datetime = "Date",
  folder_name = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Environmental_Varibales_Downloaded/Bluelink/Daily/",
  export_path = "data/temp/wz_results",
  max_depth = -200,
  fill_gaps = TRUE,
  buffer = 0.25, # buffer in degree
  verbose = TRUE,
  export_step = TRUE
)
tictoc::toc()

res

res |> 
  as.data.frame() |> 
  rstatix::get_summary_stats(Wz, type = "common")

max(res$date)

# check NAs
res |>
  dplyr::filter(is.na(Wz)) |> 
  dplyr::select(date, lat, lon, uv) |> 
  print(n=50)


mapview::mapview(res |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> dplyr::filter(PA ==1) |>  dplyr::select(-Date))



# CMEMS for data >2024 -----------------------------------------------
# nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Data_Global/daily/wo"
# 
# 
# # List all NetCDF files in the directory
# nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
# nc_files



input_df_2 <- dt.sf |> as.data.frame() |> 
  dplyr::filter(Date >= as.Date("2024-01-01")) 

# manipulating the 29/02/24 issue by making 1st of March, but then reversing it back to old date by just adding new colunmn which will then be removed 
input_df_3 <- input_df_2 |> dplyr::mutate(Date2 = dplyr::case_when(
  Date == as.Date("2024-02-29") ~ as.Date("2024-03-01"),
  TRUE                          ~ Date)
)

input_df_3 |>  dplyr::filter(Date2 == as.Date("2024-02-29"))

min(input_df_2$Date)

tictoc::tic("CMEMS wz extraction took: ")
res_2 <- extractWz_CM(df = input_df_3,
                                   X = "lon", 
                                   Y = "lat", 
                                   datetime = "Date2", 
                                   folder_name = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/daily/wo/",
                                   max_depth = -200,
                                   fill_gaps = TRUE, buffer = 0.25,
                                   export_path =  "data/temp/wz_results_CM",
                      keep_nc_files = TRUE)
tictoc::toc()

glimpse(res_2)






# check NAs
res_2 |>
  dplyr::filter(is.na(Wz)) |> 
  dplyr::select(date, lat, lon, uv)


mapview::mapview(res_2 |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> dplyr::filter(PA ==1) |>  dplyr::select(-Date))

# Merging and savres_2# Merging and saving ------------------------------------------------------


str(res)
str(res_2)

dt_wz <- dplyr::bind_rows(res, res_2 |> dplyr::select(-Date2)) |> 
  dplyr::select(-aux.date, -geometry) |>
  # dplyr::arrange(id, rep, date ) |> 
  dplyr::arrange(id, date ) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
str(dt_wz)

sight_wz <- dt_wz
track_wz <- dt_wz

mapview::mapview(track_wz  |> dplyr::filter(PA ==1) |>  dplyr::select(-Date))


saveRDS(sight_wz, "data/work_files/Sightings_PA_w_dynSDM_100_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")

saveRDS(sight_wz, "data/work_files/Sightings_Validation_data_bathy_dist_sst_uv_mld_chl_wz.rds")
saveRDS(sight_wz, "data/work_files/Sightings_Validation_data_bathy_dist_sst_uv_mld_chl_wz_full.rds")

saveRDS(track_wz, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")
saveRDS(track_wz, "data/work_files/Tracks_sims_30_thinned_2010_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")
saveRDS(track_wz, "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")

# aslo save as processed file:
saveRDS(sight_wz, "data/processed/Sightings_PA_w_dynSDM_100_2010_2025_extract.rds")

saveRDS(sight_wz, "data/work_files/Sightings_Validation_data_extract.rds")
saveRDS(sight_wz, "data/work_files/Sightings_Validation_data_extract_full.rds")

saveRDS(track_wz, "data/processed/Tracks_PA_w_dynSDM_10_2010_2025_extract.rds")
saveRDS(track_wz, "data/processed/Tracks_PA_w_3to7days_dynSDM_10_2010_2025_extract.rds")
saveRDS(track_wz, "data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract.rds")
saveRDS(track_wz, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract.rds")



