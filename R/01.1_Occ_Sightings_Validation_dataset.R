#_____________________________________________________________________________
#                        Sightings: Validation Dataset
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
library(dynamicSDM)
source("R/00_Helper_Functions.R")


# here we collect sightings for which we know that sharks were likley in a "foirgaging" or "low movement persistence" state. This dataset is then used as external independent dataset to validate the models performance. 
# the input data set is based on Miller et al. 2025 and further screened for information on sharks either being observed feeding, or groups of sharks being reported. This process has been doen for the ISRA proposal. We also explore GBIF and OBIS records for any further potential occurences to supplement our data set -> only one additional occurrence that undpoubdatly ppoint to a gropup of sharks wwas found in the OBIS database. Any other occurrecnes were more than one sharks were recorded on the same day/same location, could not be clearly identified as actual multiple occurrences (they could have been recorded multiple time but the same animal)

# import Data  -----------------------------------------------------------



df <- readxl::read_xlsx("data/input/ISRA_whaleshark_groups_feeding.xlsx", sheet = "Sheet1")

str(df)



library(rgbif)




bbox <- sf::st_bbox(c(xmin = 140, xmax = 180, ymin = -50, ymax = 10), crs = 4326)

# Convert the bounding box to an sf polygon object
bbox_polygon <- sf::st_as_sfc(bbox)

# Convert the polygon geometry to WKT format
bbox_wkt <- bbox_polygon |> 
  st_geometry() |> 
  st_as_text()

# occ_search(scientificName = "Rhincodon typus")
# occ_data(scientificName = "Rhincodon typus")

# taxonKey <- name_backbone("Rhincodon typus")$usageKey
# taxonKey
# occ_search(taxonKey = taxonKey)




occ_download(pred("taxonKey", 2417522),
             pred("hasGeospatialIssue", FALSE),
             pred("hasCoordinate", TRUE),
             pred("occurrenceStatus","PRESENT"), 
             pred_gte("year", 1950),
             pred_not(pred_in("basisOfRecord",c("LIVING_SPECIMEN"))),
             pred_within(bbox_wkt),
             user = 'ingo_miller', 
             pwd = 'asmingo2022', 
             email = 'ingo.miller@my.jcu.edu.au', 
             format = "SIMPLE_CSV") 


occ_download_cancel(key="0055703-251009101135966",
                    user = 'ingo_miller',
                    pwd = 'asmingo2022')

occ_download_wait('0055703-251009101135966') # checks if download is finished

gbif <- occ_download_get('0055703-251009101135966') |> 
  occ_download_import()

gbif_citation('0055703-251009101135966')

glimpse(gbif)

gbif_2 <- gbif |>  
  dplyr::mutate(date = lubridate::make_date(year, month, day),
                lat = decimalLatitude,
                lon = decimalLongitude,
                n_sharks = individualCount,
                flag = issue) |> 
  dplyr::filter(!is.na(date))


dups <- gbif_2 |> 
  dplyr::filter(duplicated(date) | duplicated(date, fromLast = TRUE)) |> 
  print(n=100)


dups |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  mapview::mapview()

gbif_3 <- gbif_2 |>  
  dplyr::filter(n_sharks > 1) |> 
  dplyr::transmute(date = date,
                   lat = decimalLatitude,
                   lon = decimalLongitude,
                   n_sharks = individualCount,
                   flag = issue) |> 
  dplyr::filter(!is.na(date))

str(gbif_3)




gbif_3|>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  mapview::mapview()


library(robis)
citation('robis')

obis_df <- occurrence(scientificname = "Rhincodon typus", 
                      startdate = as.Date("1950-01-01"),
                      geometry = bbox_wkt)


obis_df
str(obis_df)



obis_df |>
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE) |> 
  mapview::mapview()


obis_df_2 <- obis_df |> 
  #dplyr::filter(!coordinateUncertaintyInMeters >20000) |> 
  dplyr::mutate(
    ID = paste0("OBIS_", occurrenceID),
    date = as.Date(eventDate, format = "%Y-%m-%d")) |> 
  dplyr::transmute(date = date,
                   id = ID,
                   lat = decimalLatitude,
                   lon = decimalLongitude,
                   Occ_Source = "Sighting")


obis_df_2


dups <- obis_df_2 |> 
  dplyr::filter(duplicated(date) | duplicated(date, fromLast = TRUE)) |> 
  dplyr::arrange(date) |> 
  print(n=100)





# combine -----------------------------------------------------------------

str(df)
str(gbif_3)


gbif_3_renamed <- gbif_3 |>
  dplyr::rename(
    Date     = date,
    Lat      = lat,
    Lon      = lon,
    n_sharks = n_sharks
  ) |>
  # add missing columns to match df
  dplyr::mutate(
    Time = NA_character_,
    `night/day` = NA_character_,
    Location = NA_character_,
    feeding = FALSE,
    Observations = NA_character_,
    Observer = NA_character_,
    Aerial = NA_character_,
    Comments = NA_character_,
    Footage = NA,
    Link = NA_character_
  ) |>
  # reorder columns to match df
  dplyr::select(names(df))

# now bind
df_all <- dplyr::bind_rows(df, gbif_3_renamed)

df_all

df_all <- df_all |> 
  dplyr::rename(date = Date,
                time = Time,
                lat = Lat,
                lon = Lon) |> 
  dplyr::select(-time)

df_final <- df_all |>
  dplyr::filter(is.na(flag)) |> 
  dplyr::mutate(x = lon, y = lat,
                year = year(date),
                month = month(date),
                day = as.numeric(day(date)))


# df_expanded <- df_all |> 
#   tidyr::uncount(weights = n_sharks, .remove = FALSE)
# 
# 
# 
# df_final <- df_expanded |> 
#   dplyr::filter(is.na(flag))




# # Biases ------------------------------------------------------------------
# 
# 
# occ_sight <- df_final |> 
#   dplyr::mutate(x = lon, y = lat,
#                 year = year(date),
#                 month = month(date),
#                 day = as.numeric(day(date)))
# 
# 
# terra::plot(terra::vect(occ_sight[, c("x", "y")],
#                         geom = c("x", "y"),
#                         crs = "+proj=longlat +datum=WGS84 +no_defs"))
# 
# 
# # Biases ------------------------------------------------------------------
# ## we use dynamic SDM to deal spatiotemporal biases 
# 
# occ_sight_chck <- dynamicSDM::spatiotemp_check(
#   occ_sight,
#   na.handle = "exclude", # NAs in coordinates or dates 
#   #duplicate.handle = "exclude", #will deal with those during thinning preocedures
#   coord.handle =  "exclude",
#   date.handle =  "exclude",
#   coordclean = FALSE, # becasue it will flag sea locations -> you bloody terrestrial turds
#   date.res = "day",
#   coordclean.species = "Rhincodon typus",
#   coordclean.handle =  "report"
# )
# 
# 
# 
# 
# 
# occ_sight_cropped <- spatiotemp_resolution(occ.data = occ_sight_chck,
#                                            spatial.res = 4,
#                                            temporal.res = "day")
# 
# terra::plot(terra::vect(occ_sight_cropped[, c("x", "y")],
#                         geom = c("x", "y"),
#                         crs = "+proj=longlat +datum=WGS84 +no_defs"))
# 
# 
# 
# # filter out locations on land 
# occ_sight_cropped_sf <- occ_sight_cropped |>
#   sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# 
# # get land polygons (Natural Earth, 1:50m scale)
# sf::sf_use_s2(FALSE)
# land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |> 
#   sf::st_transform(4326) |>
#   sf::st_make_valid() |> 
#   sf::st_union()
# land
# 
# # keep only points NOT on land
# on_land <- sf::st_intersects(occ_sight_cropped_sf, land, sparse = FALSE)[,1]
# # keep only those NOT on land
# occ_filtered_1 <- occ_sight_cropped_sf[!on_land, ]
# 
# mapview::mapview(occ_filtered_1)
# 
# 
# occ_filtered_1 <- as.data.frame(occ_filtered_1)
# 
# 
# # a first check of biases
# bias_results <- spatiotemp_bias(occ.data = occ_filtered_1,
#                                 temporal.level = c("day", "month", "year"),
#                                 plot = TRUE,
#                                 spatial.method = "simple",
#                                 prj = "+proj=longlat +datum=WGS84")
# 
# bias_results
# 
# 
# 
# 
# 
# ## some temporal filtering
# 
# occ_filtered_2 <- spatiotemp_extent(
#   occ.data = occ_filtered_1,
#   temporal.ext = c("2010-01-01", "2024-12-30"),
#   spatial.ext = c(136, 179, -40, 0),
#   prj = "+proj=longlat +datum=WGS84"
# )
# 
# 
# 
# # create WB core area for spatial bias check 
# cent_sf <- sf::st_as_sf(data.frame(x = 143.87, y = -12.18),
#                         coords = c("x","y"), crs = 4326)
# 
# 
# bias_results <- spatiotemp_bias(occ.data = occ_filtered_2,
#                                 temporal.level = c("day", "month", "year"),
#                                 plot = TRUE,
#                                 spatial.method = "simple",
#                                 centroid = cent_sf,
#                                 radius = 1,
#                                 prj = "EPSG:3577"
#                                 #prj = "+proj=longlat +datum=WGS84"
# )
# 
# bias_results
# 
# 
# 
# 
# 
# 
# 
# ### data thinning to avoid spatioremporal correlation 
# occ_sight_thin <- spatiotemp_thin(occ.data = occ_filtered_2,   
#                                   temporal.method = "day",
#                                   temporal.dist = 0,        # ignore temporal thinning altogether 
#                                   spatial.split.degrees = 1,
#                                   spatial.dist = 500,
#                                   iterations = 5
# )
# 
# 
# 
# 
# 
# 
# 
# 
# # Thinned
# par(mfrow = c(1, 2))
# 
# terra::plot(terra::vect(occ_sight_thin[, c("x", "y")],
#                         geom = c("x", "y"),
#                         crs = "+proj=longlat +datum=WGS84 +no_defs"), col = "red", main = "thinned")
# terra::plot(terra::vect(occ_filtered_2[, c("x", "y")],
#                         geom = c("x", "y"),
#                         crs = "+proj=longlat +datum=WGS84 +no_defs"), add = FALSE, main = "sightings")
# 
# dev.off()
# 
# set.seed(999)
# 
# 
# 
# 
# # focus on spatial bias centered around Wreck Bay with 100 km buffer -> because this is the area we will use the randomly sample absence locations from, thus, random sampling on whole region (i.e., simple) for spatial bias doe snot make sense here, as we will adress this bias by using a buffer 
# bias_results <- spatiotemp_bias(occ.data = occ_sight_thin,
#                                 temporal.level = c("day", "month", "year"),
#                                 plot = TRUE,
#                                 spatial.method = "core",
#                                 centroid = cent_sf,
#                                 radius = 0.5,
#                                 prj = "EPSG:3577"
#                                 #prj = "+proj=longlat +datum=WGS84"
# )
# 
# 
# bias_results
# 
# # now focus on temporal biases (most importan there is month as we will use monht in the models)
# bias_results <- spatiotemp_bias(occ.data = occ_sight_thin,
#                                 temporal.level = c("day", "month", "year"),
#                                 plot = TRUE,
#                                 spatial.method = "simple",
#                                 # centroid = cent_sf,
#                                 # radius = 100000,
#                                 prj = "EPSG:3577"
#                                 #prj = "+proj=longlat +datum=WGS84"
# )
# 
# bias_results
# ## temporal bias based on day is not biased! good; go ahead
# 



## remove any sightings in same grid cell and month of trainign data 

training_df <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")
sight_df <- df_final
str(sight_df)

predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log_ML.tif")

r <- predictors_mean[[1]]

train_cells <- cellFromXY(r, cbind(training_df$lon, training_df$lat))
sight_cells <- cellFromXY(r, cbind(sight_df$lon, sight_df$lat))

train_key  <- paste0(train_cells, "_", training_df$month)
sight_key  <- paste0(sight_cells, "_", sight_df$month)
keep_ct    <- !(sight_key %in% train_key)
sight_clean_ct <- sight_df[keep_ct, ]


# Pseudoabsence sampling  -------------------------------------------------



## creating ocean buffer and land barrier for pseudo-absences 

input <- sight_df

occ_sf <- sf::st_as_sf(input, coords = c("x", "y"), crs = 4326)

aoi <- occ_sf |>
  sf::st_as_sf(coords = c("x","y"), crs = 4326) |>
  sf::st_transform(3577) |>
  sf::st_buffer(500000) |>          # 500 km
  sf::st_union() |>                 # works on all sf versions
  sf::st_make_valid()

land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |>
  sf::st_transform(3577) |>
  sf::st_make_valid() |>
  sf::st_crop(sf::st_bbox(aoi))

plot(land$geometry)

land_buf <- land |>
  sf::st_buffer(1000) |>
  sf::st_union() |>
  sf::st_make_valid()

plot(land_buf)

ocean <- sf::st_difference(aoi, land_buf) |>
  sf::st_make_valid()
plot(ocean)









#Pseudo-absences generated within spatial and temporal buffer

#######################################
## ATTENTION: dynamicSDM has some issue swith the buffer as it does niot really create the distances specified! Thus, a modified version of the code is used so that acriual buffer ditances provided are realised.
#######################################

# pseudo_abs_buff <- dynamicSDM::spatiotemp_pseudoabs(occ.data = occ_sight_thin,
#                                         spatial.method = "buffer",
#                                         temporal.method = "buffer",
#                                         spatial.ext = ocean,
#                                         spatial.buffer = c(25000, 100000),
#                                         temporal.buffer = 10, # sample 10 days before and after 
#                                         n.pseudoabs = nrow(occ_sight_thin) * 50)



set.seed(42)                # for reproducibility
k <- 5                     # pseudo-absences per occurrence

pseudo_abs_buff <- purrr::map_dfr(seq_len(nrow(occ_sf)), function(i) {
  out_i <- spatiotemp_pseudoabs_fix(
    spatial.method  = "buffer",
    temporal.method = "buffer",
    occ.data        = occ_sf[i, , drop = FALSE],     # <- one occurrence at a time
    spatial.ext     = ocean,                         # ocean mask to avoid absences sampled on land 
    spatial.buffer  = c(10000, 100000),              # 10-100 km (10 to avoid absence in same grod than presnece and 100km to consider their potential movement in a day)
    temporal.buffer = 5,                             
    n.pseudoabs     = k,                             # exactly k per occurrence
    prj             = "+proj=longlat +datum=WGS84"   # same as package default
  )
  
  out_i |>
    dplyr::mutate(occ_index = i, rep = dplyr::row_number()) |>
    dplyr::select(x, y, year, month, day, occ_index, rep)
})


str(pseudo_abs_buff)



# plot presences and absences 
aus_sf <- rnaturalearth::ne_countries( scale = 10, returnclass = "sf") |>
  sf::st_transform(4326) 
aus <- terra::vect(aus_sf)

dev.off()
terra::plot(terra::vect(pseudo_abs_buff[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), pch = 21, col = "red", add=F)
plot(ocean, border = "green4", add=T)
terra::plot(terra::vect(input[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"),pch = 21, col = "blue",add=T)
terra::plot(aus, border = "black", col = "grey20", add = TRUE)









# adding to occurences 

pseudo_abs_buff$lat <- pseudo_abs_buff$y
pseudo_abs_buff$lon <- pseudo_abs_buff$x
pseudo_abs_buff$PA <- rep(0, nrow(pseudo_abs_buff))


input <- input |> 
  dplyr::mutate(PA = 1,
                rep = 0,
                source = "Sightings_record")

str(input)
str(pseudo_abs_buff)


pseudo_abs_buff_1 <- pseudo_abs_buff |> 
  dplyr::mutate(lon = x, lat = y,
                Occ_Source = "Sighting",
                date = as.POSIXct(
                  paste(year, month, day, sep = "-"),
                  format = "%Y-%m-%d"),
                source = "dynamicSDM",
                id = paste0("bg_", occ_index)) |> 
  dplyr::select(-occ_index)





str(input)
str(pseudo_abs_buff_1)

sight_prep <- input |>
  dplyr::transmute(
    lon        = lon,
    lat        = lat,
    x          = x,
    y          = y,
    date       = date,
    year       = year,
    month      = month,
    day        = day,
    PA         = 1,                        
    source     = source,
    Occ_Source = "Sighting",
    id         = NA_character_            
  )


pseudo_prep <- pseudo_abs_buff_1 |>
  dplyr::transmute(
    lon        = lon,
    lat        = lat,
    x          = x,
    y          = y,
    date       = date,
    year       = year,
    month      = month,
    day        = day,
    PA         = 0,                         
    source     = source,
    Occ_Source = Occ_Source,
    id         = id
  )




sight_val_dt <- as.data.frame(rbind(sight_prep, pseudo_prep))
str(sight_val_dt)


sight_val_dt |> dplyr::summarise(start = min(date),
                           end = max(date))


saveRDS(sight_val_dt, "data/processed/Sightings_Validation_data.rds")



