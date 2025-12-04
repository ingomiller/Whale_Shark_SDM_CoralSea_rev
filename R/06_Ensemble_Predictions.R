#_____________________________________________________________________________
#                        Models: Ensemble
#_____________________________________________________________________________


# remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
# citation("SDMtune")
library(tidyverse)
library(raster)
library(terra)
library(ncdf4)
library(rerddap)
library(rerddapXtracto)
library(stringr)
library(lubridate)
library(SDMtune)
library(patchwork)
source("R/00_Helper_Functions.R") 
library(rasterVis)
library(tidyterra)
library(ggspatial)
library(sf)
library(basemaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(grid)
library(patchwork)




## Themes/Palletes

#Generate a color palette using cmocean
cm_ocean_palette <- cmocean::cmocean(name = "balance", alpha = 1, start = 0.2, end = 0.8, direction = 1)

#for levelplot background
myTheme <- rasterVis::BTCTheme()
myTheme$panel.background$col = 'gray' 





# Tracking Data -----------------------------------------------------------




### Maxent Models 

# mean_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_Maxent_mean_rev_mp.tif")
# 
# months_all_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_2018_2025_rev_mp.tif")
# 
# monthly_max <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_means_rev.tif")
# 
# season_4_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons4_means_rev_mp.tif")
# 
# season_2_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp.tif")

mean_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_Maxent_mean_rev_mp_crwPA.tif")

months_all_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_2018_2025_rev_mp_crwPA.tif")

monthly_max <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_means_rev_crwPA.tif")

season_4_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons4_means_rev_mp_crwPA.tif")

season_2_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp_crwPA.tif")

names(months_all_brt)

### GAMM Models


# mean_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_mean_rev_mp.tif")
# 
# months_all_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_2018_2025_rev_mp.tif")
# 
# monthly_gam <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_means_rev_mp.tif")
# 
# season_4_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons4_means_rev_mp.tif")
# 
# season_2_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp.tif")

mean_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_mean_rev_mp_crw_PA.tif")

months_all_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_2018_2025_rev_mp_crwPA.tif")

monthly_gam <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_rev_mp_crwPA.tif")

season_4_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons4_means_rev_mp_crwPA.tif")

season_2_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp_crwPA.tif")



### BRT

# mean_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_BRT_mean_rev_mp.tif")
# 
# months_all_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_2019_2025_rev_mp.tif")
# 
# monthly_brt <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_means_rev_mp.tif")
# 
# season_4_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp.tif")
# 
# season_2_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp.tif")

mean_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_BRT_mean_rev_mp_crwPA.tif")

months_all_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_2019_2025_rev_mp_crwPA.tif")

monthly_brt <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_means_rev_mp_crwPA.tif")

season_4_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp_crwPA.tif")

season_2_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA.tif")



plot(season_2_brt, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))

plot(season_2_max, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))

plot(season_2_gam, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))



# gams need to be cropped
template <- mean_max

mean_gam <- mean_gam |>
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  


months_all_gam <- months_all_gam |>
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  

monthly_gam <- monthly_gam |>
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  

season_4_gam <- season_4_gam |>
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear") 

season_2_gam <- season_2_gam |>
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear") 



### MESS maps

# mess_mean <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_mean_stack.tif")
# 
# mess_monthly_stack <-terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_monthly_stack.tif")
# 
# mess_mean_calendar <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_monthly_means.tif")
# 
# seasonal_mess_stack_2 <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_seasons_2.tif")
# 
# seasonal_mess_stack_4 <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_seasons_4.tif")

mess_mean <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_mean_stack_crwPA.tif")

mess_monthly_stack <-terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_monthly_stack_crwPA.tif")

mess_mean_calendar <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_monthly_means_crwPA.tif")

seasonal_mess_stack_2 <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_seasons_2_crwPA.tif")

seasonal_mess_stack_4 <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_seasons_4_crwPA.tif")

mess_mean <- mess_mean |> 
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  

mess_monthly_stack <- mess_monthly_stack |> 
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  

mess_mean_calendar <- mess_mean_calendar |> 
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  

seasonal_mess_stack_2 <- seasonal_mess_stack_2 |> 
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  


seasonal_mess_stack_4 <- seasonal_mess_stack_4 |> 
  terra::crop(terra::ext(template)) |>
  terra::resample(template, method = "bilinear")  



mess_mean

plot(seasonal_mess_stack_2, range = c(-200, 0))
seasonal_mess_stack_2

mess_extrap <- terra::ifel(mess_mean < 0, 1, NA_real_)
mess_extrap

mess_extrap_df <- mess_extrap |>
  terra::as.data.frame(xy = TRUE, na.rm = TRUE) |> 
  tidyr::pivot_longer(
    cols = -c(x, y),
    # names_to  = "Season",
    values_to = "flag"
  ) |>
  dplyr::filter(!is.na(flag))

# cols: x, y, mean   (mean = the layer name)
# names(mess_extrap_df) <- c("x", "y", "flag")
mess_extrap_df



names(seasonal_mess_stack_2) <- c("Monsoon Season (Nov - Apr)", 
                                  "Dry Season (May - Oct)")
seasonal_mess_stack_2_extrap <- terra::ifel(seasonal_mess_stack_2 < 0, 1, NA_real_)

seasonal_mess_stack_2_extrap_df <- seasonal_mess_stack_2_extrap|>
  terra::as.data.frame(xy = TRUE) |>
  tidyr::pivot_longer(
    cols = -c(x, y),
    names_to  = "lyr",
    values_to = "flag"
  ) |>
  dplyr::filter(!is.na(flag))



names(seasonal_mess_stack_4) <- c("Jan - Mar", 
                                  "Apr - Jun",
                                  "Jul - Sep",
                                  "Oct - Dec")
seasonal_mess_stack_4_extrap <- terra::ifel(seasonal_mess_stack_4 < 0, 1, NA_real_)

seasonal_mess_stack_4_extrap_df <- seasonal_mess_stack_4_extrap|>
  terra::as.data.frame(xy = TRUE) |>
  tidyr::pivot_longer(
    cols = -c(x, y),
    names_to  = "lyr",
    values_to = "flag"
  ) |>
  dplyr::filter(!is.na(flag))




# Ensemble  ---------------------------------------------------------------

mean_max
mean_gam
mean_brt

models_mean <- c(mean_gam, mean_max, mean_brt)
names(models_mean) <- c("GAM", "MaxEnt", "BRT")

plot(models_mean, range = c(0,1))

## using TSS test weights 
weights <- c(0.33, 0.33, 0.28)

## Now with a weighted mean model (weighted using model mean test AUCs)
mean_climate_ensemble <- terra::weighted.mean(models_mean, w = weights)
mean_climate_ensemble_thres <- mean_climate_ensemble > 0.5

plot(mean_climate_ensemble)
plot(mean_climate_ensemble_thres)



season2 <- c(season_2_gam, season_2_max, season_2_brt)
season2
names(season2) <- c("GAM: Monsoon Season", "GAM: Dry Season ", "MaxEnt: Monsoon Season", "MaxEnt: Dry Season", "BRT: Monsoon Season", "BRT: Dry Season ")


plot(season2, range = c(0, 1))

monsoon <- season2[[c(1, 3, 5)]]
dry <- season2[[c(2, 4, 6)]]
monsoon

plot(monsoon, range = c(0, 1))
plot(dry, range = c(0, 1))


monsoon_ensemble <- terra::weighted.mean(monsoon, w = weights)
dry_ensemble <- terra::weighted.mean(dry, w = weights)



season_2_ensemble <- c(monsoon_ensemble, dry_ensemble)
names(season_2_ensemble) <- c("Monsoon_Season_Ensemble_Tracking", "Dry_Season_Ensemble_Tracking")
season_2_ensemble


season_2_ensemble_thres <- season_2_ensemble > 0.5

plot(season_2_ensemble, range = c(0,1))
plot(season_2_ensemble_thres)

plot(season_2_ensemble, range = c(0, 1), xlim = c(140, 170), ylim = c(-40, 0), axes =FALSE, legend =FALSE)



high = 1

season_2_ensemble_map <-  rasterVis::levelplot(season_2_ensemble,
                                                 layout=c(2, 1),
                                                 par.settings = myTheme,
                                                 main = "Seasonal whale shark habitat suitability (Tracking Ensemble)",
                                                 names.attr = c("Monsoon (Nov - Apr)", 
                                                                "Dry (May - Oct)"),
                                                 zlim = c(0, high),
                                                 at = seq(0, high, length.out = 100),
                                                 # col.regions =  topo.colors(100),
                                                 col.regions = cm_ocean_palette,
                                                 # col.regions =  viridisLite::inferno(100),
                                                 colorkey = list(space = "right", 
                                                                 length = 0.5, 
                                                                 height = 0.75,
                                                                 labels = list(
                                                                   at = c(0, high),  # Positions for "Low" and "High"
                                                                   labels = c("Low", "High"))))  # Labels to use

season_2_ensemble_map





## Quartals 

quartals <- c(season_4_gam, season_4_max, season_4_brt)
quartals
names(quartals) <- c("GAM: Q1", "GAM: Q2", "GAM: Q3", "GAM: Q4",
                       "MaxEnt: Q1", "MaxEnt: Q2", "MaxEnt: Q3", "MaxEnt: Q4", 
                       "BRT: Q1", "BRT: Q2", "BRT: Q3", "BRT: Q4")

plot(quartals)

q1 <- quartals[[c(1, 5, 9)]]
q2 <- quartals[[c(2, 6, 10)]]
q3 <- quartals[[c(3, 7, 11)]]
q4 <- quartals[[c(4, 8, 12)]]
q1

plot(q1)
plot(q4)


q1_ensemble <- terra::weighted.mean(q1, w = weights)
q2_ensemble <- terra::weighted.mean(q2, w = weights)
q3_ensemble <- terra::weighted.mean(q3, w = weights)
q4_ensemble <- terra::weighted.mean(q4, w = weights)



quartals_ensemble <- c(q1_ensemble, q2_ensemble, q3_ensemble, q4_ensemble)
quartals_ensemble

names(quartals_ensemble) <- c("Jan - Mar", 
                             "Apr - Jun",
                             "Jul - Sep",
                             "Oct - Dec")
quartals_ensemble
quartals_ensemble_thres <- quartals_ensemble > 0.5

plot(quartals_ensemble, range = c(0, 1))
plot(quartals_ensemble_thres)



## Monthly 

monthly_brt
weights

monthly <- c(monthly_brt, monthly_max, monthly_gam)
monthly

w_brt  <- weights[1]   # your actual weights
w_max  <- weights[2]
w_gam  <- weights[3]

w_sum <- w_brt + w_max + w_gam

monthly_ensemble <- (
  monthly_brt * w_brt +
    monthly_max * w_max +
    monthly_gam * w_gam
) / w_sum

names(monthly_ensemble) <- base::month.abb  

monthly_ensemble

plot(monthly_ensemble, range = c(0,1))




# fill near coastal NAs for mapping  --------------------------------------
w3 <- matrix(1, 3, 3)
w5 <- matrix(1, 5, 5)  # fallback window


mean_climate_ensemble <-
  mean_climate_ensemble |>
  # Pass 1: fill NAs with 3x3 local mean
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 2: repeat 3x3 to grow into slightly larger gaps
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 3 (optional): one 5x5 sweep for stubborn holes
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w5, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Clamp to [0, 1]
  (\(r) terra::ifel(r < 0, 0, terra::ifel(r > 1, 1, r)))()

plot(mean_climate_ensemble)


season_2_ensemble <-
  season_2_ensemble |>
  # Pass 1: fill NAs with 3x3 local mean
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 2: repeat 3x3 to grow into slightly larger gaps
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 3 (optional): one 5x5 sweep for stubborn holes
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w5, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Clamp to [0, 1]
  (\(r) terra::ifel(r < 0, 0, terra::ifel(r > 1, 1, r)))()

plot(season_2_ensemble)



quartals_ensemble <-
  quartals_ensemble |>
  # Pass 1: fill NAs with 3x3 local mean
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 2: repeat 3x3 to grow into slightly larger gaps
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 3 (optional): one 5x5 sweep for stubborn holes
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w5, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Clamp to [0, 1]
  (\(r) terra::ifel(r < 0, 0, terra::ifel(r > 1, 1, r)))()

plot(quartals_ensemble)

monthly_ensemble <-
  monthly_ensemble |>
  # Pass 1: fill NAs with 3x3 local mean
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 2: repeat 3x3 to grow into slightly larger gaps
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 3 (optional): one 5x5 sweep for stubborn holes
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w5, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Clamp to [0, 1]
  (\(r) terra::ifel(r < 0, 0, terra::ifel(r > 1, 1, r)))()

plot(monthly_ensemble)


# unchanged_ok <- (mean_climate_ensemble[!is.na(x)] == x[!is.na(x)])
# all(unchanged_ok[])  # should be TRUE









# Pretty Maps -------------------------------------------------------------

high = 1

cities <- data.frame(Loc = c("Cooktown", "Cairns",  "Mackay", "Brisbane", "Sydney"),
                     Group = c("Town", "Town",  "Town", "Town", "Town"),
                     # Season = c("Monsoon Season (Nov - Apr)"),
                     # lyr = 1,
                     lat = c(-15.4758, -16.918246, -21.1434, -27.4705, -33.911609),
                     lon = c(145.2471, 145.771359, 149.1868, 153.026, 151.186715)) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities2 <- data.frame(Loc = c( "Port\nMoresby"),
                      lat = c(  -9.4790),
                      lon = c( 147.1494)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB <- data.frame(Loc = c( "Wreck\nBay"),
                 lat = c(  -12.132504),
                 lon = c( 143.893818)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

str(cities)





australia_map_data <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name == "Australia" & name == "Papua")
print(australia_map_data)

australia_map_data

world <- ne_countries(scale = 10, returnclass = "sf")

AUS_PNG <- world |> filter(name == "Australia" | name == "Papua New Guinea" | name == "Indonesia")
plot(AUS_PNG)

# Define the extent of the study area 
ext <- c(xmin = 140, xmax = 170, ymin = -40, ymax = 0)

# Create the inset map of Australia + Papua New Guinea and add the study area rectangle
AUS_map <- ggplot() +
  geom_sf(data = world, fill = "lightgray", colour = NA, size = 0.1) +
  coord_sf(xlim = c(112, 170), ylim = c(-44, 1), expand = FALSE) +
  geom_rect(aes(xmin = ext["xmin"], xmax = ext["xmax"], ymin = ext["ymin"], ymax = ext["ymax"]), 
            fill = NA, color = "black", 
            linetype = "solid", size = 0.5) +  # Study area rectangle
  theme_void() +  
  theme(panel.border = element_blank(),
        panel.background = element_blank())  

AUS_map




cm_ocean_palette <- cmocean::cmocean(name = "balance", alpha = 1, start = 0, end = 1, direction = 1)(256)

sz = 2.5


P1 <- ggplot2::ggplot() +
  # basemaps::basemap_gglayer(ext_bbox, map_service = "esri", map_type = "world_ocean_reference") +
  # scale_fill_identity() + 
  tidyterra::geom_spatraster(data = mean_climate_ensemble) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  geom_sf(data = cities,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = cities2,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities2,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 1, 
                                nudge_y = 0, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = WB,
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 2,
                                nudge_y = -0.5,
                                size = 2,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  scale_fill_gradientn(
    # colours = cm_ocean_palette,
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(100, "pt"),
      barwidth  = grid::unit(10,  "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability", x = "Longitude", y = "Latitude", title = "") +
  
  scale_y_continuous(limits = c(-40, 0), breaks = seq(-35, -5, by = 5), expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170), breaks = seq(145, 165, by = 5),expand = c(0,0)) +
  
  # ggspatial::annotation_scale(location = "bl", 
  #                  width_hint = 0.25,
  #                  pad_x = unit(.5, "cm"),
  #                  pad_y = unit(.5, "cm")) +
  
  ggspatial::annotation_north_arrow(location = "bl",
                                    which_north = "true", 
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"),
                                    pad_x = unit(1, "cm"),
                                    pad_y = unit(1.5, "cm"),
                                    style =  north_arrow_fancy_orienteering) +
  
  guides(alpha = "none") +
  
  # Add the expanded inset map as an annotation
  annotation_custom(
    grob = ggplotGrob(AUS_map),  # Convert the inset map to a grob
    xmin = 141, xmax = 150, ymin = -33, ymax = -25  # Adjust position & size of inset map
  ) +
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 151, y = -13.5, label = "Coral Sea", color = "grey5", size = sz, fontface = "italic") +
  annotate("text", x = 153, y = -7.5, label = "Solomon\nSea", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "grey5", size = sz, fontface = "italic")+
  # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "grey5", size = 2, fontface = "italic")+ 
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    # legend.key.width = rel(0.25),
    # legend.key.height = rel(0.75),
    legend.title = element_text(size=8),
    legend.text = element_text(size =8),
    axis.text = element_text(size = 8),
    #axis.title = element_text(size = 10),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) 




P1

mess_extrap_df

P1.m <- P1 +
  # adding MESS map
  ggplot2::geom_tile(
    data = mess_extrap_df,
    mapping = ggplot2::aes(x = x, y = y),
    fill = "black",
    alpha = 0.3,
    width = 0.1,  # match raster resolution
    height = 0.1,
    inherit.aes = FALSE
  ) +
  tidyterra::geom_spatraster_contour(
    data = mess_mean, breaks = 0,
    colour = "white", linewidth = 0.5, linetype = "dotted",
    show.legend = FALSE
  )
  
  

P1.m



ggsave("HSM_Map_Mean_Climate_Ensemble.png", plot = P1.m, path = "outputs/figures", scale =1, width = 14, height = 17, units = "cm", dpi = 600)



path = "/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/"

ggsave("Figure_3.png", plot = P1.m, path = path, scale =1, width = 14, height = 17, units = "cm", dpi = 600)
ggsave("Figure_3.pdf", plot = P1.m, path = path, scale =1, width = 14, height = 17, units = "cm", dpi = 600, device = "pdf")
ggsave("Figure_3.eps", plot = P1.m, path = path, scale =1, width = 14, height = 17, units = "cm", dpi = 600)



## 2 seasons 

print(season_2_ensemble)
names(season_2_ensemble) <- c("Monsoon Season (Nov - Apr)", 
                                  "Dry Season (May - Oct)")

season_2_ensemble[[1]]




cities_A <- data.frame(Loc = c(
                             "Cairns",  
                             "Brisbane", 
                             "Sydney"), 
                     lat = c(
                             -16.918246, 
                             -27.4705, 
                             -33.911609),
                     lon = c(
                             145.771359, 
                             153.026, 
                             151.186715),
                     lyr = "Monsoon Season (Nov - Apr)") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities_B <- data.frame(Loc = c(
                               "Cairns",  
                               "Brisbane", 
                               "Sydney"), 
                       lat = c(
                               -16.918246, 
                               -27.4705, 
                               -33.911609),
                       lon = c( 
                               145.771359, 
                               153.026, 
                               151.186715),
                       lyr = "Dry Season (May - Oct)") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


cities_2 <- rbind(cities_A, cities_B)
cities_2$lyr <- factor(cities_2$lyr, levels = c("Monsoon Season (Nov - Apr)", "Dry Season (May - Oct)"))


PM_A <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Monsoon Season (Nov - Apr)") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


PM_B <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Dry Season (May - Oct)") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


PM <- rbind(PM_A, PM_B)
PM$lyr <- factor(PM$lyr, levels = c("Monsoon Season (Nov - Apr)", "Dry Season (May - Oct)"))

WB_A <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Monsoon Season (Nov - Apr)") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB_B <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Dry Season (May - Oct)") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


WB <- rbind(WB_A, WB_B)
WB$lyr <- factor(WB$lyr, levels = c("Monsoon Season (Nov - Apr)", "Dry Season (May - Oct)"))




P2 <- ggplot2::ggplot() +
  # basemaps::basemap_gglayer(ext_bbox, map_service = "esri", map_type = "world_ocean_reference") +
  # scale_fill_identity() + 
  tidyterra::geom_spatraster(data = season_2_ensemble) +
  
  facet_wrap(~lyr) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  geom_sf(data = subset(cities_2, lyr == "Monsoon Season (Nov - Apr)"),
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = subset(cities_2, lyr == "Monsoon Season (Nov - Apr)"),
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = subset(PM, lyr == "Dry Season (May - Oct)"),
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = subset(PM, lyr == "Dry Season (May - Oct)"),
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 0.5, 
                                nudge_y = 0, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = subset(WB, lyr == "Monsoon Season (Nov - Apr)"),
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 3,
                                nudge_y = -0.5,
                                size = 2,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  scale_fill_gradientn(
    # colours = cm_ocean_palette,
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(10, "pt"),
      barwidth  = grid::unit(100,  "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability", x = "Longitude", y = "Latitude", title = "") +
  
  scale_y_continuous(limits = c(-40, 0), breaks = seq(-35, -5, by = 5), expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170), breaks = seq(145, 165, by = 5),expand = c(0,0)) +
  
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = 3, fontface = "bold") +
  # annotate("text", x = 151, y = -13.5, label = "Coral Sea", color = "grey95", size = sz, fontface = "italic") +
  # annotate("text", x = 153, y = -7.5, label = "Solomon\nSea", color = "grey95", size = sz, fontface = "italic")+
  # annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "grey95", size = sz, fontface = "italic")+
  # # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  # annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "grey95", size = sz, fontface = "italic")+
  # annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "grey95", size = 2, fontface = "italic")+ 
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    # legend.key.width = rel(0.25),
    # legend.key.height = rel(0.75),
    legend.title = element_text(size=8),
    legend.text = element_text(size =8),
    axis.text = element_text(size = 8),
    #axis.title = element_text(size = 10),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    strip.text.x = element_text(size = 10, color = "black", face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) 


P2

seasonal_mess_stack_2_extrap_df

seasonal_mess_stack_2

# Contours at level 0 for each layer
c_monsoon <- terra::as.contour(seasonal_mess_stack_2[["Monsoon Season (Nov - Apr)"]], levels = 0)
c_dry     <- terra::as.contour(seasonal_mess_stack_2[["Dry Season (May - Oct)"]],     levels = 0)

# Tag with Season and combine
c_monsoon$lyr <- "Monsoon Season (Nov - Apr)"
c_dry$lyr     <- "Dry Season (May - Oct)"
c_all <- rbind(c_monsoon, c_dry)

lvl <- c("Monsoon Season (Nov - Apr)", "Dry Season (May - Oct)")

seasonal_mess_stack_2_extrap_df$lyr <- factor(seasonal_mess_stack_2_extrap_df$lyr, levels = lvl)

c_all
c_all$lyr  <- factor(c_all$lyr,  levels = lvl)

P2.m <- P2 +
  # adding MESS map
  ggplot2::geom_tile(
    data = seasonal_mess_stack_2_extrap_df,
    mapping = ggplot2::aes(x = x, y = y),
    fill = "black",
    alpha = 0.3,
    width = 0.1,  # match raster resolution
    height = 0.1,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_sf(
    data = sf::st_as_sf(c_all),
    colour = "white",
    linewidth = 0.5,
    linetype = "dotted",
    inherit.aes = FALSE
  )



P2.m


ggsave("HSM_Map_2seasons_Climate_Ensemble.png", plot = P2.m, path = "outputs/figures", scale =1, width = 18, height = 15, units = "cm", dpi = 600)



path = "/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/"

ggsave("Figure_4.png", plot = P2.m, path = path, scale =1, width = 18, height = 15, units = "cm", dpi = 600)
ggsave("Figure_4.pdf", plot = P2.m, path = path, scale =1, width = 18, height = 15, units = "cm", dpi = 600, device = "pdf")
ggsave("Figure_4.eps", plot = P2.m, path = path, scale =1, width = 18, height = 15, units = "cm", dpi = 600)









## 4 seasons 
quartals_ensemble
print(quartals_ensemble)
names(quartals_ensemble) <- c("Jan - Mar", 
                              "Apr - Jun",
                              "Jul - Sep",
                              "Oct - Dec")

quartals_ensemble[[1]]




cities_A <- data.frame(Loc = c(
  "Cairns",  
  "Townville",
  "Brisbane", 
  "Sydney"), 
  lat = c(
    -16.918246, 
    -19.2599,
    -27.4705, 
    -33.911609),
  lon = c(
    145.771359, 
    146.8137,
    153.026, 
    151.186715),
  lyr = "Jan - Mar") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities_B <- data.frame(Loc = c(
  "Cairns",  
  "Townville",
  "Brisbane", 
  "Sydney"), 
  lat = c(
    -16.918246, 
    -19.2599,
    -27.4705, 
    -33.911609),
  lon = c(
    145.771359, 
    146.8137,
    153.026, 
    151.186715),
  lyr = "Apr - Jun") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities_C <- data.frame(Loc = c(
  "Cairns",  
  "Townville",
  "Brisbane", 
  "Sydney"), 
  lat = c(
    -16.918246, 
    -19.2599,
    -27.4705, 
    -33.911609),
  lon = c(
    145.771359, 
    146.8137,
    153.026, 
    151.186715),
  lyr = "Jul - Sep") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities_D <- data.frame(Loc = c(
  "Cairns",  
  "Townville",
  "Brisbane", 
  "Sydney"), 
  lat = c(
    -16.918246, 
    -19.2599,
    -27.4705, 
    -33.911609),
  lon = c(
    145.771359, 
    146.8137,
    153.026, 
    151.186715),
  lyr = "Oct - Dec") |>   # This must match the facet name) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


cities_4 <- rbind(cities_A, cities_B, cities_C, cities_D)
cities_4$lyr <- factor(cities_4$lyr, levels = c("Jan - Mar", "Apr - Jun", "Jul - Sep", "Oct - Dec"))
cities_4

PM_A <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Jan - Mar") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


PM_B <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Apr - Jun") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

PM_C <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Jul - Sep") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

PM_D <- data.frame(Loc = c( "Port\nMoresby"),
                   lat = c(  -9.4790),
                   lon = c( 147.1494),
                   lyr = "Oct - Dec") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

PM <- rbind(PM_A, PM_B, PM_C, PM_D)
PM$lyr <- factor(PM$lyr, levels = c("Jan - Mar", "Apr - Jun", "Jul - Sep", "Oct - Dec"))

WB_A <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Jan - Mar") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB_B <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Apr - Jun") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB_C <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Jul - Sep") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB_D <- data.frame(Loc = c( "Wreck\nBay"),
                   lat = c(  -12.132504),
                   lon = c( 143.893818),
                   lyr = "Oct - Dec") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB <- rbind(WB_A, WB_B, WB_C, WB_D)
WB$lyr <- factor(WB$lyr, levels = c("Jan - Mar", "Apr - Jun", "Jul - Sep", "Oct - Dec"))









sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
str(sight)

val_sights <- sight |>
  dplyr::filter(PA == 1) |> 
  dplyr::mutate(
    lyr = dplyr::case_when(
      month %in% 1:3   ~ "Jan - Mar",
      month %in% 4:6   ~ "Apr - Jun",
      month %in% 7:9   ~ "Jul - Sep",
      month %in% 10:12 ~ "Oct - Dec"
    ),
    lyr = factor(
      lyr,
      levels = c("Jan - Mar", "Apr - Jun", "Jul - Sep", "Oct - Dec")
    )
  )

str(val_sights)

set.seed(1)
sight_sample <- val_sights |>
  dplyr::slice_sample(prop = 0.7)   # or n = 200

sight_sample

P3 <- ggplot2::ggplot() +
  # basemaps::basemap_gglayer(ext_bbox, map_service = "esri", map_type = "world_ocean_reference") +
  # scale_fill_identity() + 
  tidyterra::geom_spatraster(data = quartals_ensemble) +
  
  facet_wrap(~lyr, ncol = 2) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  # seasonal sightings – will appear in the correct panel via `lyr`
  ggplot2::geom_sf(
    # data = val_sights,
    data = sight_sample,
    shape = 21,
    fill  = "hotpink1",
    colour = "black",
    alpha = 0.5,
    size = 0.5,
    show.legend = FALSE
  ) +
  
  # geom_sf(
  #   data = sight_sample,
  #   shape = 21,
  #   fill  = NA,
  #   colour = "white",
  #   alpha = 0.5,
  #   size  = 0.5,
  #   stroke = 0.25,
  #   show.legend = FALSE
  # ) +
  
  geom_sf(data = subset(cities_4, lyr == "Jan - Mar"),
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = subset(cities_4, lyr == "Jan - Mar"),
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 1.5, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  geom_sf(data = subset(cities_4, lyr == "Oct - Dec"),
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = subset(cities_4, lyr == "Oct - Dec"),
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 1.5, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = subset(PM, lyr == "Apr - Jun"),
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = subset(PM, lyr == "Apr - Jun"),
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 0.5, 
                                nudge_y = 0, 
                                size = 1.5, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = subset(WB, lyr == "Oct - Dec"),
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 4,
                                nudge_y = -0.5,
                                size = 1.5,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  scale_fill_gradientn(
    # colours = cm_ocean_palette,
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(10, "pt"),
      barwidth  = grid::unit(100,  "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability", x = "Longitude", y = "Latitude", title = "") +
  
  scale_y_continuous(limits = c(-40, 0), breaks = seq(-35, -5, by = 10), expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170), breaks = seq(145, 165, by = 10),expand = c(0,0)) +
  
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = 3, fontface = "bold") +
  # annotate("text", x = 151, y = -13.5, label = "Coral Sea", color = "grey95", size = sz, fontface = "italic") +
  # annotate("text", x = 153, y = -7.5, label = "Solomon\nSea", color = "grey95", size = sz, fontface = "italic")+
  # annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "grey95", size = sz, fontface = "italic")+
  # # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  # annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "grey95", size = sz, fontface = "italic")+
  # annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "grey95", size = 2, fontface = "italic")+ 
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    # legend.key.width = rel(0.25),
    # legend.key.height = rel(0.75),
    legend.title = element_text(size=8),
    legend.text = element_text(size =8),
    axis.text = element_text(size = 8),
    #axis.title = element_text(size = 10),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.background = element_blank(),
    strip.text.x = element_text(size = 8, color = "black", face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) 


P3


seasonal_mess_stack_4

# Contours at level 0 for each layer
c_q1 <- terra::as.contour(seasonal_mess_stack_4[["Jan - Mar"]], levels = 0)
c_q2 <- terra::as.contour(seasonal_mess_stack_4[["Apr - Jun"]], levels = 0)
c_q3 <- terra::as.contour(seasonal_mess_stack_4[["Jul - Sep"]], levels = 0)
c_q4 <- terra::as.contour(seasonal_mess_stack_4[["Oct - Dec"]], levels = 0)

c_q1
# Tag with Season and combine
c_q1$lyr <- "Jan - Mar"
c_q2$lyr <- "Apr - Jun"
c_q3$lyr <- "Jul - Sep"
c_q4$lyr <- "Oct - Dec"
c_q <- rbind(c_q1, c_q2, c_q3, c_q4)

lvl <- c("Jan - Mar", 
         "Apr - Jun",
         "Jul - Sep",
         "Oct - Dec")

seasonal_mess_stack_4_extrap_df$lyr <- factor(seasonal_mess_stack_4_extrap_df$lyr, levels = lvl)

c_q
c_q$lyr  <- factor(c_q$lyr,  levels = lvl)

P3.m <- P3 +
  # adding MESS map
  ggplot2::geom_tile(
    data = seasonal_mess_stack_4_extrap_df,
    mapping = ggplot2::aes(x = x, y = y),
    fill = "black",
    alpha = 0.3,
    width = 0.1,  # match raster resolution
    height = 0.1,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_sf(
    data = sf::st_as_sf(c_q),
    colour = "white",
    linewidth = 0.2,
    linetype = "dotted",
    inherit.aes = FALSE
  )



P3.m


ggsave("HSM_Map_4seasons_Climate_Ensemble.png", plot = P3.m, path = "outputs/figures", scale =1, width = 18, height = 10, units = "cm", dpi = 600)



path = "/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/"

ggsave("Figure_5.png", plot = P3.m, path = path, scale =1, width = 18, height = 21, units = "cm", dpi = 600)
ggsave("Figure_5.pdf", plot = P3.m, path = path, scale =1, width = 18, height = 15, units = "cm", dpi = 600, device = "pdf")
ggsave("Figure_5.eps", plot = P3.m, path = path, scale =1, width = 18, height = 15, units = "cm", dpi = 600)








# Monthly maps ------------------------------------------------------------


val_sights_month <- sight |>
  dplyr::filter(PA == 1) |>
  dplyr::mutate(
    lyr = factor(
      base::month.abb[month],
      levels = base::month.abb          # Jan, Feb, ..., Dec
    )
  )

# set.seed(1)
# sight_sample_month <- val_sights_month |>
#   dplyr::slice_sample(prop = 0.7)



P_monthly <- ggplot2::ggplot() +
  # Ensemble SDM for each month
  tidyterra::geom_spatraster(data = monthly_ensemble) +
  
  # 12 calendar months
  ggplot2::facet_wrap(~lyr, ncol = 4) +
  
  # Land
  ggplot2::geom_sf(
    data   = world,
    fill   = "grey50",
    colour = "grey20",
    linewidth = 0.1
  ) +
  
  # Monthly sightings (in matching month panels via `lyr`)
  ggplot2::geom_sf(
    data   = val_sights_month,
    shape  = 21,
    fill   = "hotpink1",
    colour = "black",
    alpha  = 0.5,
    size   = 0.5,
    show.legend = FALSE
  ) +
  
  # # Cities – appear in all panels
  # ggplot2::geom_sf(
  #   data   = cities_all,
  #   shape  = 21,
  #   colour = "black",
  #   fill   = "yellow",
  #   alpha  = 0.5,
  #   size   = 1,
  #   show.legend = FALSE
  # ) +
  # ggsflabel::geom_sf_text_repel(
  #   data   = cities_all,
  #   aes(label = Loc),
  #   colour = "black",
  #   nudge_x = -2.5,
  #   nudge_y = 0.25,
  #   size    = 1.5,
  #   force   = 1,
  #   force_pull = 10,
  #   seed    = 10
  # ) +
  
  # Spatial extent
  ggplot2::coord_sf(
    xlim   = c(140, 170),
    ylim   = c(-40, 0),
    expand = FALSE
  ) +
  
  # Colourbar / palette
  ggplot2::scale_fill_gradientn(
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = ggplot2::guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight    = grid::unit(10, "pt"),
      barwidth     = grid::unit(100, "pt")
    )
  ) +
  
  ggplot2::labs(
    fill = "Rel. Habitat\nSuitability",
    x    = "Longitude",
    y    = "Latitude",
    title = ""
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(-40, 0),
    breaks = seq(-35, -5, by = 10),
    expand = c(0, 0)
  ) +
  ggplot2::scale_x_continuous(
    limits = c(140, 170),
    breaks = seq(145, 165, by = 10),
    expand = c(0, 0)
  ) +
  
  # Country abbreviations
  # ggplot2::annotate("text", x = 144,   y = -29.5, label = "AUS", colour = "black", size = 3, fontface = "bold") +
  # ggplot2::annotate("text", x = 144,   y = -6,    label = "PNG", colour = "black", size = 3, fontface = "bold") +
  # ggplot2::annotate("text", x = 165,   y = -21,   label = "NC",  colour = "black", size = 3, fontface = "bold") +
  # ggplot2::annotate("text", x = 166.75,y = -15.5, label = "VU",  colour = "black", size = 3, fontface = "bold") +
  # ggplot2::annotate("text", x = 158.5, y = -8,    label = "SI",  colour = "black", size = 3, fontface = "bold") +
  
  ggplot2::theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = ggplot2::element_text(size = 8),
    legend.text       = ggplot2::element_text(size = 8),
    axis.text         = ggplot2::element_text(size = 8),
    axis.title        = ggplot2::element_blank(),
    panel.border      = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.background  = ggplot2::element_blank(),
    strip.text.x      = ggplot2::element_text(size = 8, colour = "black", face = "bold"),
    strip.background  = ggplot2::element_blank(),
    plot.margin       = grid::unit(c(0, 0, 0, 0), "cm")
  )

P_monthly



mess_mean_calendar_backup <- mess_mean_calendar
mess_mean_calendar <- mess_mean_calendar_backup
names(mess_mean_calendar) <- base::month.abb  # "Jan", ..., "Dec" if not already

mess_mean_calendar <- terra::ifel(mess_mean_calendar < 0, 1, NA_real_)
mess_mean_calendar

monthly_mess_extrap_df <- mess_mean_calendar|>
  terra::as.data.frame(xy = TRUE) |>
                         tidyr::pivot_longer(
                           cols = -c(x, y),
                           names_to  = "lyr",
                           values_to = "flag"
                         ) |>
  dplyr::filter(!is.na(flag)) |> 
  dplyr::filter(flag < 0) |>
  dplyr::mutate(
    lyr = factor(lyr, levels = base::month.abb)
  )
monthly_mess_extrap_df 

# Create contours at level 0 for each calendar month
cont_list <- lapply(base::month.abb, function(m) {
  r <- mess_mean_calendar[[m]]
  c <- terra::as.contour(r, levels = 0)
  c$lyr <- m
  c
})

c_m <- do.call(rbind, cont_list)

c_m$lyr <- factor(c_m$lyr, levels = base::month.abb)


P_monthly_m <- P_monthly +
  # dark shaded extrapolated MESS cells
  ggplot2::geom_tile(
    data = monthly_mess_extrap_df,
    mapping = ggplot2::aes(x = x, y = y),
    fill   = "black",
    alpha  = 0.3,
    width  = 0.1,
    height = 0.1,
    inherit.aes = FALSE
  ) +
  # white dotted MESS = 0 contour
  ggplot2::geom_sf(
    data        = sf::st_as_sf(c_m),
    colour      = "white",
    linewidth   = 0.2,
    linetype    = "dotted",
    inherit.aes = FALSE
  )

P_monthly_m



ggsave("HSM_Map_monthly_Climate_Ensemble.png", plot = P_monthly_m, path = "outputs/figures", scale =1, width = 18, height = 21, units = "cm", dpi = 600)



path = "/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/"

ggsave("Figure_S9.png", plot = P_monthly_m, path = path, scale =1, width = 18, height = 21, units = "cm", dpi = 600)
ggsave("Figure_6.pdf", plot = P_monthly_m, path = path, scale =1, width = 18, height = 21, units = "cm", dpi = 600, device = "pdf")
ggsave("Figure_6.eps", plot = P_monthly_m, path = path, scale =1, width = 18, height = 21, units = "cm", dpi = 600)






# Animations --------------------------------------------------------------

## animate the predictions:

library(animation)
library(gganimate)






GBR_zone <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/GBR_Zones/Great_Barrier_Reef_Marine_Park_Boundary.shp")


GBR_boundary_dt <- st_boundary(GBR_zone) %>%
  st_cast("LINESTRING") %>%
  st_coordinates() %>%
  as.data.frame()

GBR_boundary_dt



GBR_features <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/GBR_Zones/Great_Barrier_Reef_Features.shp")

GBR_features





# Set up animation parameters
ani.options(interval = 2)  # Set duration of each frame in seconds
#

month_labels <- names(monthly_track_en)



p <- ggplot() +
  geom_spatraster(data = monthly_ensemble) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  geom_sf(data = cities,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = cities2,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities2,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 1, 
                                nudge_y = 0, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = WB,
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 2,
                                nudge_y = -0.5,
                                size = 2,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  scale_fill_gradientn(
    # colours = cm_ocean_palette,
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(100, "pt"),
      barwidth  = grid::unit(10,  "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability", x = "Longitude", y = "Latitude", title = "{current_frame}") +
  
  scale_y_continuous(limits = c(-40, 0), breaks = seq(-35, -5, by = 5), expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170), breaks = seq(145, 165, by = 5),expand = c(0,0)) +
  
  # ggspatial::annotation_scale(location = "bl", 
  #                  width_hint = 0.25,
  #                  pad_x = unit(.5, "cm"),
  #                  pad_y = unit(.5, "cm")) +
  
  ggspatial::annotation_north_arrow(location = "bl",
                                    which_north = "true", 
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"),
                                    pad_x = unit(1, "cm"),
                                    pad_y = unit(1.5, "cm"),
                                    style =  north_arrow_fancy_orienteering) +
  
  guides(alpha = "none") +
  
  # Add the expanded inset map as an annotation
  annotation_custom(
    grob = ggplotGrob(AUS_map),  # Convert the inset map to a grob
    xmin = 141, xmax = 150, ymin = -33, ymax = -25  # Adjust position & size of inset map
  ) +
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 151, y = -13.5, label = "Coral Sea", color = "grey5", size = sz, fontface = "italic") +
  annotate("text", x = 153, y = -7.5, label = "Solomon\nSea", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "grey85", size = sz, fontface = "italic")+
  # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "grey5", size = 2, fontface = "italic")+ 
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = rel(1.5),
    legend.key.width = rel(0.75),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_blank(),
    plot.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    strip.text.x = element_text(size = 10, color = "black", face = "bold", hjust = 0),
    strip.background = element_blank()
  )


p


# Animate the raster layers using gganimate
animation <- p +
  gganimate::transition_manual(lyr) +
  gganimate::ease_aes('linear')



gganimate::animate(
  animation, 
  width = 800,           # Increase width (adjust as needed)
  height = 800,          # Increase height (adjust as needed)
  res = 120,             # Increase resolution (dots per inch)
  nframes = 100, 
  fps = 5, 
  end_pause = 5,            # Adds a 20-frame pause at the end
  renderer = gganimate::gifski_renderer(loop = TRUE)  
)


# Save the animation
gganimate::anim_save("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/HSM_Map_mean_Monthly_Ensemble_crw_mp.gif", animation = last_animation())


# Animate and save as an MP4 video
gganimate::animate(
  animation, 
  width = 800,           # Width in pixels
  height = 800,          # Height in pixels
  res = 120,             # Resolution (DPI)
  nframes = 100,         # Total frames for smooth animation
  fps = 5,              # Frames per second
  end_pause = 5,        # Pause at the end of the animation
  renderer = gganimate::av_renderer("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/HSM_Map_mean_Monthly_Ensemble_crw_mp.mp4")  # MP4 renderer
)










## all months 


compute_weighted_means_all_dates <- function(monthly_gam, monthly_max, monthly_brt, weights) {
  # Ensure that all rasters have the same number of layers (dates)
  if (nlyr(monthly_gam) != nlyr(monthly_max) || nlyr(monthly_gam) != nlyr(monthly_brt)) {
    stop("All model layers must have the same number of layers (dates).")
  }
  
  # Get the dates from one of the layers (assumed to match across rasters)
  date_labels <- names(monthly_gam)
  
  # Initialize a list to store the weighted mean rasters for each date
  weighted_means_all_dates <- list()
  
  # Loop through each date layer
  for (i in seq_along(date_labels)) {
    # Extract the corresponding layers for the current date from each model
    date_layers <- c(monthly_gam[[i]], monthly_max[[i]], monthly_brt[[i]])
    
    # Compute the weighted mean for the current date
    weighted_mean <- terra::weighted.mean(date_layers, w = weights)
    
    # Store the result
    weighted_means_all_dates[[i]] <- weighted_mean
  }
  
  # Combine the weighted mean rasters into a SpatRaster stack
  weighted_mean_stack <- rast(weighted_means_all_dates)
  
  # Set names for each layer (date labels)
  names(weighted_mean_stack) <- date_labels
  
  return(weighted_mean_stack)
}


# Use the function to calculate the weighted means for all dates
months_track_en <- compute_weighted_means_all_dates(monthly_gam = months_all_gam, monthly_max = months_all_max, months_all_brt, weights)

months_track_en

months_track_en <-
  months_track_en |>
  # Pass 1: fill NAs with 3x3 local mean
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 2: repeat 3x3 to grow into slightly larger gaps
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w3, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Pass 3 (optional): one 5x5 sweep for stubborn holes
  (\(r) terra::cover(
    r,
    terra::focal(r, w = w5, fun = base::mean, na.rm = TRUE)
  ))() |>
  # Clamp to [0, 1]
  (\(r) terra::ifel(r < 0, 0, terra::ifel(r > 1, 1, r)))()


date_labels <- format(as.Date(names(months_track_en)), "%Y-%b")
date_labels



p <- ggplot() +
  geom_spatraster(data = months_track_en) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  geom_sf(data = cities,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = cities2,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 0.5, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities2,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 1, 
                                nudge_y = 0, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = WB,
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 2,
                                nudge_y = -0.5,
                                size = 2,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  scale_fill_gradientn(
    # colours = cm_ocean_palette,
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(100, "pt"),
      barwidth  = grid::unit(10,  "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability", x = "Longitude", y = "Latitude", title = "{date_labels[as.integer(frame)]}") +
  
  scale_y_continuous(limits = c(-40, 0), breaks = seq(-35, -5, by = 5), expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170), breaks = seq(145, 165, by = 5),expand = c(0,0)) +
  
  # ggspatial::annotation_scale(location = "bl", 
  #                  width_hint = 0.25,
  #                  pad_x = unit(.5, "cm"),
  #                  pad_y = unit(.5, "cm")) +
  
  ggspatial::annotation_north_arrow(location = "bl",
                                    which_north = "true", 
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"),
                                    pad_x = unit(1, "cm"),
                                    pad_y = unit(1.5, "cm"),
                                    style =  north_arrow_fancy_orienteering) +
  
  guides(alpha = "none") +
  
  # Add the expanded inset map as an annotation
  annotation_custom(
    grob = ggplotGrob(AUS_map),  # Convert the inset map to a grob
    xmin = 141, xmax = 150, ymin = -33, ymax = -25  # Adjust position & size of inset map
  ) +
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = 3, fontface = "bold") +
  annotate("text", x = 151, y = -13.5, label = "Coral Sea", color = "grey5", size = sz, fontface = "italic") +
  annotate("text", x = 153, y = -7.5, label = "Solomon\nSea", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "grey85", size = sz, fontface = "italic")+
  # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "grey5", size = sz, fontface = "italic")+
  annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "grey5", size = 2, fontface = "italic")+ 
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = rel(1.5),
    legend.key.width = rel(0.75),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_blank(),
    plot.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    strip.text.x = element_text(size = 10, color = "black", face = "bold", hjust = 0),
    strip.background = element_blank()
  )


p


# Animate the raster layers using gganimate
animation <- p +
  gganimate::transition_manual(lyr) +
  gganimate::ease_aes('linear')



gganimate::animate(
  animation, 
  width = 800,           # Increase width (adjust as needed)
  height = 800,          # Increase height (adjust as needed)
  res = 120,             # Increase resolution (dots per inch)
  nframes = 200, 
  fps = 10, 
  end_pause = 5,            # Adds a 20-frame pause at the end
  renderer = gganimate::gifski_renderer(loop = TRUE)  
)


# Save the animation
gganimate::anim_save("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/HSM_Map_Monthly_all_Ensemble_crw_mp.gif", animation = last_animation())


# Animate and save as an MP4 video
gganimate::animate(
  animation, 
  width = 800,           # Width in pixels
  height = 800,          # Height in pixels
  res = 120,             # Resolution (DPI)
  nframes = 200,         # Total frames for smooth animation
  fps = 10,              # Frames per second
  end_pause = 5,        # Pause at the end of the animation
  renderer = gganimate::av_renderer("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/HSM_Map_Monthly_all_Ensemble_crw_mp.mp4")  # MP4 renderer
)




# SPATIAL BLOCK CV VS RANDOM FOLDS CV -------------------------------------

mean_brt.rf <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_BRT_mean_rev_mp_crwPA_RANDOM_FOLDS_CV.tif")

months_all_brt.rf <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_2019_2025_rev_mp_crwPA_RANDOM_FOLDS_CV.tif")

monthly_brt.rf <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_means_rev_mp_crwPA_RANDOM_FOLDS_CV.tif")

season_4_brt.rf <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp_crwPA_RANDOM_FOLDS_CV.tif")

season_2_brt.rf <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA_RANDOM_FOLDS_CV.tif")



sb_rf_comp <- c(season_2_brt, season_2_brt.rf)
sb_rf_comp 

names(sb_rf_comp) <- c(
  "Monsoon – Block CV",
  "Dry – Block CV",
  "Monsoon – Random folds",
  "Dry – Random folds"
)

plot(sb_rf_comp, range = c(0,1))


sb_df <- sb_rf_comp |>
  tidyterra::as_tibble(xy = TRUE) |>
  tidyr::pivot_longer(
    cols      = -c(x, y),
    names_to  = "lyr",
    values_to = "value"
  )


sb_df <- sb_df |>
  dplyr::mutate(
    Season = dplyr::case_when(
      stringr::str_detect(lyr, "Monsoon") ~ "Monsoon Season (Nov - Apr)",
      stringr::str_detect(lyr, "Dry")     ~ "Dry Season (May - Oct)",
      TRUE ~ NA_character_
    ),
    CV = dplyr::case_when(
      stringr::str_detect(lyr, "Block")   ~ "Block CV",
      stringr::str_detect(lyr, "Random")  ~ "Random CV",
      TRUE ~ NA_character_
    )
  )

sb_df <- sb_df |>
  dplyr::mutate(
    Season = factor(
      Season,
      levels = c(
        "Monsoon Season (Nov - Apr)",
        "Dry Season (May - Oct)"
      )
    ),
    CV = factor(
      CV,
      levels = c(
        "Block CV",
        "Random CV"
      )
    )
  )

str(sb_df)






P_folds <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = sb_df,
    mapping = ggplot2::aes(x = x, y = y, fill = value)
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(CV),
    cols = ggplot2::vars(Season)
  ) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  

  scale_fill_gradientn(
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(10, "pt"),
      barwidth  = grid::unit(100, "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability",
       x = "Longitude",
       y = "Latitude",
       title = "") +
  
  scale_y_continuous(limits = c(-40, 0),
                     breaks = seq(-35, -5, by = 5),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170),
                     breaks = seq(145, 165, by = 5),
                     expand = c(0,0)) +
  
  # your annotate("text", ...) calls here
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size=8),
    legend.text  = element_text(size =8),
    axis.text    = element_text(size = 8),
    axis.title   = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    strip.text.x = element_text(size = 10, color = "black", face = "bold"),
    strip.text.y = element_text(size = 10, color = "black", face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )


P_folds

ggsave("HSM_BRT_Seasonal_Block_vs_RandomCV.png", plot = P_folds, path = "outputs/final_figures", scale =1, width = 18, height = 21, units = "cm", dpi = 600)





# Comparing the Model types -----------------------------------------------




season2_models <- c(season_2_brt, season_2_max, season_2_gam)
season2_models



names(season2_models) <- c(
  "Monsoon – BRT",
  "Dry – BRT",
  "Monsoon – MaxEnt",
  "Dry – MaxEnt",
  "Monsoon – GAMM",
  "Dry – GAMM"
)

plot(season2_models, range = c(0,1))


season2_models.df <- season2_models |>
  tidyterra::as_tibble(xy = TRUE) |>
  tidyr::pivot_longer(
    cols      = -c(x, y),
    names_to  = "lyr",
    values_to = "value"
  )


season2_models.df <- season2_models.df |>
  dplyr::mutate(
    Season = dplyr::case_when(
      stringr::str_detect(lyr, "Monsoon") ~ "Monsoon Season (Nov - Apr)",
      stringr::str_detect(lyr, "Dry")     ~ "Dry Season (May - Oct)",
      TRUE ~ NA_character_
    ),
    Model = dplyr::case_when(
      stringr::str_detect(lyr, "BRT")   ~ "BRT",
      stringr::str_detect(lyr, "MaxEnt")  ~ "MaxEnt",
      stringr::str_detect(lyr, "GAMM")  ~ "GAMM",
      TRUE ~ NA_character_
    )
  )

season2_models.df <- season2_models.df |>
  dplyr::mutate(
    Season = factor(
      Season,
      levels = c(
        "Monsoon Season (Nov - Apr)",
        "Dry Season (May - Oct)"
      )
    ),
    Model = factor(
      Model,
      levels = c(
        "BRT", "MaxEnt", "GAMM"
      )
    )
  )

str(season2_models.df)






P_models <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = season2_models.df,
    mapping = ggplot2::aes(x = x, y = y, fill = value)
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(Model),
    rows = ggplot2::vars(Season)
  ) +
  
  ggplot2::geom_sf(data = world, fill = "grey50", colour = "grey20", linewidth = 0.1) +
  
  
  scale_fill_gradientn(
    colours = viridisLite::turbo(n = 100, direction = 1, begin = 0, end = 1),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::number_format(accuracy = 0.1),
    oob     = scales::squish,
    name    = "Rel. Habitat\nSuitability",
    guide   = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = grid::unit(10, "pt"),
      barwidth  = grid::unit(100, "pt")
    )
  ) +
  
  labs(fill = "Rel. Habitat\nSuitability",
       x = "Longitude",
       y = "Latitude",
       title = "") +
  
  scale_y_continuous(limits = c(-40, 0),
                     breaks = seq(-35, -5, by = 5),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(140, 170),
                     breaks = seq(145, 165, by = 5),
                     expand = c(0,0)) +
  
  # your annotate("text", ...) calls here
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size=8),
    legend.text  = element_text(size =8),
    axis.text    = element_text(size = 8),
    axis.title   = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.background = element_blank(),
    strip.text.x = element_text(size = 10, color = "black", face = "bold"),
    strip.text.y = element_text(size = 10, color = "black", face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )


P_models

ggsave("HSM_Seasonal_Model_Comparison.png", plot = P_models, path = "outputs/final_figures", scale =1, width = 18, height = 17, units = "cm", dpi = 600)




