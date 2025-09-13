#_____________________________________________________________________________
#                        Extractions: SST (Erddap)
#_____________________________________________________________________________



library(rerddap)
library(rerddapXtracto)
library(ncdf4)
library(parsedate)
library(plotdap)
library(gganimate)
library(tictoc)
library(doParallel)
library(foreach)
library(tidyverse)
library(future.apply)
library(future)
library(progressr)
source("R/00_Helper_Functions.R") 




# Import location data  ---------------------------------------------------

sight <- readRDS("data/work_files/Sightings_PA_w_dynSDM_10_raw_2010_2025_bathy_dist.rds")
str(sight)

tracks <- readRDS("data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy_dist.rds")




# Extractions -------------------------------------------------------------

dt <- sight
dt <- tracks

# define input file 
input.df <- dt |> as.data.frame() |> 
  dplyr::mutate(Date = base::as.Date(date))
str(input.df)


## Set Erddap dataset 
SST_info <- rerddap::info('jplMURSST41', url = 'https://coastwatch.pfeg.noaa.gov/erddap/')
SST_info
str(SST_info)
my_interp <- c('Mean',  '16')



SST.df <- extract_ERDDAP_data(
  input.df  = input.df,
  info      = SST_info,
  id_col    = "id",
  lat_col   = "lat",
  lon_col   = "lon",
  date_col  = "date",
  rep_col   = "rep",
  variable  = "analysed_sst",
  xlen = 0.02, ylen = 0.02,
  parallel  = TRUE,
  pb_style  = "cli"
)



# Extract releavnt columns, and merge names with original df, so you can merge them
str(SST.df)
glimpse(SST.df)

SST_df <- SST.df |> 
  dplyr::as_tibble() |>
  dplyr::select(id, requested.date, lat, lon, rep, mean.analysed_sst) |>
  dplyr::rename(Date = requested.date, 
                MUR_SST = mean.analysed_sst) |>
  dplyr::mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d", tz="UTC"),
                Date = base::as.Date(Date))

str(SST_df)
str(input.df)


# Join MER_SST col to original df 
sight_SST <- dplyr::left_join(input.df, SST_df, by = c("id", "Date", "lat", "lon", "rep")) 

str(sight_SST)
print(sight_SST)


sight_SST |> 
  as.data.frame() |> 
  rstatix::get_summary_stats(MUR_SST, type = "common")


# Save  -------------------------------------------------------------------

saveRDS(sight_SST, "data/work_files/Sightings_PA_w_dynSDM_10_raw_2010_2025_bathy_dist_sst.rds")


