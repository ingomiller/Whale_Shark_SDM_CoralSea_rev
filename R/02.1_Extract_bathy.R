#_____________________________________________________________________________
#                        Extractions: Bathymetry
#_____________________________________________________________________________

## to remove old rast temp files: terra::tmpFiles(current = FALSE, old = TRUE, remove = TRUE)


library(raster)
library(terra)
library(tidyverse)
library(sf)
library(marmap)
library(stars)



# Import location data  ---------------------------------------------------

sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_raw_2010_2025.rds")
str(sight)

# track <- readRDS

locs <- readRDS("data/processed/Tracking_PA_w_dynSDM_10_raw_2010_2025.rds")
locs <- readRDS("data/processed/Tracking_PA_3days_w_dynSDM_10_raw_2010_2025.rds")


# final dataset
locs <- readRDS("data/work_files/Tracks_mp_sims_30_thinned_2018_2025_final.rds")

locs <- readRDS( "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_final_extA_10d.rds")

locs <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final.rds")

str(locs)

glimpse(locs)

unique(locs$Tag_type)

tracks <- locs

track_p <- locs |>
  dplyr::filter(PA == 1) |>
  dplyr::filter(Tag_type %in% c("SPLASH10", "PSAT", "PSAT_SPOT", "SPOT")) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

v <- terra::vect(track_p)
bb <- terra::ext(v)
buff <- 1  # ~2° buffer added




# Load Bathymetry Files ---------------------------------------------------

## GBR 30m bathymetry 
bathy30_a <-  terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GBR30//Great_Barrier_Reef_A_2020_30m_MSL_cog.tif")
bathy30_b <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GBR30//Great_Barrier_Reef_B_2020_30m_MSL_cog.tif")
bathy30_c <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GBR30//Great_Barrier_Reef_C_2020_30m_MSL_cog.tif")
bathy30_d <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GBR30//Great_Barrier_Reef_D_2020_30m_MSL_cog.tif")


## PNG 100m bathymetry

bathy_Papua_100 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/PNG//GulfOfPapua_100m_bathymetry.tif")


## GEBCO 2023 (global to fill gaps)
gebco_nc <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Geography/Bathymetry/GEBCO_2023//GEBCO_2023.nc"


bathy_Gebco <- terra::rast(gebco_nc)


terra::plot(bathy_Gebco, range = c(-11000, 0))
cols <- terra::map.pal("bgyr", 100) 
cols_full <- c(cols, "grey10")
brks <- c(seq(-9000, 0, length.out = 101), Inf)


terra::plot(
  bathy_Gebco,
  xlim   = c(terra::xmin(bb) - buff, terra::xmax(bb) + buff),
  ylim   = c(terra::ymin(bb) - buff, terra::ymax(bb) + buff),
  range  = c(-9000, 0),   # fix scale; keeps legend as a gradient
  col    = cols,
  #breaks = 10,
  legend = TRUE
)
# land_mask <- bathy_Gebco >= 0
terra::plot(land_mask, add = TRUE, col = c("#00000000", "grey10"), legend = FALSE)
terra::points(v, pch = 21, cex = 0.25, col = "white", bg = NA)



# Calculate Slope and Roughness -------------------------------------------

slope30_a <- terra::terrain(bathy30_a, v = "slope", unit = "degrees", neighbors = 8)
slope30_b <- terra::terrain(bathy30_b, v = "slope", unit = "degrees", neighbors = 8)
slope30_c <- terra::terrain(bathy30_c, v = "slope", unit = "degrees", neighbors = 8)
slope30_d <- terra::terrain(bathy30_d, v = "slope", unit = "degrees", neighbors = 8)

Slope_Papua100 <- terra::terrain(bathy_Papua_100, v = "slope", unit = "degrees", neighbors = 8)

Slope_Gebco <- terra::terrain(bathy_Gebco, v = "slope", unit = "degrees", neighbors = 8)


rough30_a <- terra::terrain(bathy30_a, v = "roughness", neighbors = 8)
rough30_b <- terra::terrain(bathy30_b, v = "roughness", neighbors = 8)
rough30_c <- terra::terrain(bathy30_c, v = "roughness", neighbors = 8)
rough30_d <- terra::terrain(bathy30_d, v = "roughness", neighbors = 8)

Rough_Papua100 <- terra::terrain(bathy_Papua_100, v = "roughness", neighbors = 8)

Rough_Gebco <- terra::terrain(bathy_Gebco, v = "roughness", neighbors = 8)







# Extractions -------------------------------------------------------------

# input data
# dt <- sight
dt <- tracks
dt
dt.v <- terra::vect(dt)

# stepwise extractions 
dt$bathy30_a <- terra::extract(bathy30_a, dt.v, method = "simple", ID = FALSE)[,1]
dt$bathy30_b <- terra::extract(bathy30_b, dt.v, method = "simple", ID = FALSE)[,1]
dt$bathy30_c <- terra::extract(bathy30_c, dt.v, method = "simple", ID = FALSE)[,1]
dt$bathy30_d <- terra::extract(bathy30_d, dt.v, method = "simple", ID = FALSE)[,1]

dt$slope30_a <- terra::extract(slope30_a, dt.v, method = "simple", ID = FALSE)[,1]
dt$slope30_b <- terra::extract(slope30_b, dt.v, method = "simple", ID = FALSE)[,1]
dt$slope30_c <- terra::extract(slope30_c, dt.v, method = "simple", ID = FALSE)[,1]
dt$slope30_d <- terra::extract(slope30_d, dt.v, method = "simple", ID = FALSE)[,1]


dt$rough30_a <- terra::extract(rough30_a, dt.v, method = "simple", ID = FALSE)[,1]
dt$rough30_b <- terra::extract(rough30_b, dt.v, method = "simple", ID = FALSE)[,1]
dt$rough30_c <- terra::extract(rough30_c, dt.v, method = "simple", ID = FALSE)[,1]
dt$rough30_d <- terra::extract(rough30_d, dt.v, method = "simple", ID = FALSE)[,1]


dt$bathy_Papua_100 <- terra::extract(bathy_Papua_100, dt.v, method = "simple", ID = FALSE)[,1]
dt$slope_Papua_100 <- terra::extract(Slope_Papua100, dt.v, method = "simple", ID = FALSE)[,1]
dt$rough_Papua_100 <- terra::extract(Rough_Papua100, dt.v, method = "simple", ID = FALSE)[,1]

dt$bathy_Gebco <- terra::extract(bathy_Gebco, dt.v, method = "simple", ID = FALSE)[,1]
dt$slope_Gebco <- terra::extract(Slope_Gebco, dt.v, method = "simple", ID = FALSE)[,1]
dt$rough_Gebco <- terra::extract(Rough_Gebco, dt.v, method = "simple", ID = FALSE)[,1]




# Consolidate  ------------------------------------------------------------

dt_merged <- dt |>
  dplyr::mutate(Depth = dplyr::coalesce(bathy30_a, bathy30_b, bathy30_c, bathy30_d, bathy_Papua_100, bathy_Gebco)) |> 
  dplyr::mutate(Slope = dplyr::coalesce(slope30_a, slope30_b, slope30_c, slope30_d, slope_Papua_100, slope_Gebco)) |> 
  dplyr::mutate(Roughness = dplyr::coalesce(rough30_a, rough30_b, rough30_c, rough30_d, rough_Papua_100, rough_Gebco)) |> 
  dplyr::mutate(Bathy_Data_Source = dplyr::case_when(
    !is.na(bathy30_a) ~ "GBR_30",
    !is.na(bathy30_b) ~ "GBR_30",
    !is.na(bathy30_c) ~ "GBR_30",
    !is.na(bathy30_d) ~ "GBR_30",
    !is.na(bathy_Papua_100) ~ "PAPUA_100",
    !is.na(bathy_Gebco) ~ "GEBCO_2023",
    TRUE ~ NA_character_)) |>
  dplyr::select(-bathy30_a, -bathy30_b, -bathy30_c, -bathy30_d, -slope30_a, -slope30_b, -slope30_c, -slope30_d, -rough30_a, -rough30_b, -rough30_c, -rough30_d, -bathy_Papua_100, -slope_Papua_100, -rough_Papua_100, -bathy_Gebco, -slope_Gebco, -rough_Gebco)


dplyr::glimpse(dt_merged)


# sight_bathy <- dt_merged
track_bathy <- dt_merged




# Save --------------------------------------------------------------------

## save temporary working file for next extractions

# saveRDS(sight_bathy, "data/work_files/Sightings_PA_w_dynSDM_10_raw_2010_2025_bathy.rds")
# saveRDS(track_bathy, "data/work_files/Tracks_PA_w_dynSDM_10_raw_2010_2025_bathy.rds")
saveRDS(track_bathy, "data/work_files/Tracks_PA_w_dynSDM_30_raw_2018_2025_final_bathy.rds")

saveRDS(track_bathy, "data/work_files/Tracks_PA_w_dynSDM_30_th_2018_2025_final_bathy.rds")



saveRDS(track_bathy, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final_bathy.rds")








