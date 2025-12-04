
#_____________________________________________________________________________
#                         SIGHTINGS DATA PREP
#_____________________________________________________________________________


# remotes::install_github("r-a-dobson/dynamicSDM", force = TRUE, build_vignettes = TRUE)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------

sightings_2025 <- readRDS("data/input/WhaleSharks_Sightings_Final_CONSERVATIVE_SDM_2010_2025.rds")
str(sightings_2025)


occ_sight <- sightings_2025 |> 
  dplyr::mutate(x = lon, y = lat,
                year = year(date),
                month = month(date),
                day = as.numeric(day(date)))


terra::plot(terra::vect(occ_sight[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))



# Biases ------------------------------------------------------------------
## we use dynamic SDM to deal spatiotemporal biases 

occ_sight_chck <- dynamicSDM::spatiotemp_check(
  occ_sight,
  na.handle = "exclude", # NAs in coordinates or dates 
  #duplicate.handle = "exclude", #will deal with those during thinning preocedures
  coord.handle =  "exclude",
  date.handle =  "exclude",
  coordclean = FALSE, # becasue it will flag sea locations -> you bloody terrestrial turds
  date.res = "day",
  coordclean.species = "Rhincodon typus",
  coordclean.handle =  "report"
)





occ_sight_cropped <- spatiotemp_resolution(occ.data = occ_sight_chck,
                                           spatial.res = 4,
                                           temporal.res = "day")

terra::plot(terra::vect(occ_sight_cropped[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))



# filter out locations on land 
occ_sight_cropped_sf <- occ_sight_cropped |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# get land polygons (Natural Earth, 1:50m scale)
sf::sf_use_s2(FALSE)
land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |> 
  sf::st_transform(4326) |>
  sf::st_make_valid() |> 
  sf::st_union()
land

# keep only points NOT on land
on_land <- sf::st_intersects(occ_sight_cropped_sf, land, sparse = FALSE)[,1]
# keep only those NOT on land
occ_filtered_1 <- occ_sight_cropped_sf[!on_land, ]

mapview::mapview(occ_filtered_1)


occ_filtered_1 <- as.data.frame(occ_filtered_1)


# a first check of biases
bias_results <- spatiotemp_bias(occ.data = occ_filtered_1,
                                temporal.level = c("day", "month", "year"),
                                plot = TRUE,
                                spatial.method = "simple",
                                prj = "+proj=longlat +datum=WGS84")

bias_results

# --> based on this we should drop 2010 t0 2013, but lets try without doing so 



## some temporal filtering

occ_filtered_2 <- spatiotemp_extent(
  occ.data = occ_filtered_1,
  temporal.ext = c("2010-01-01", "2024-12-30"),
  spatial.ext = c(136, 179, -40, 0),
  prj = "+proj=longlat +datum=WGS84"
)

# map to manually filter duplicates, land locations etc
occ_filtered_2 |>
  dplyr::mutate(SourceType = case_when(
    grepl("^GBIF", id) ~ "GBIF",
    grepl("^OBIS", id) ~ "OBIS",
    TRUE               ~ "BOF"
  )) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  #dplyr::filter(grepl("^(GBIF|OBIS)", id)) |>
  mapview::mapview(zcol = "SourceType")

occ_filtered_3 <- occ_filtered_2 |> 
  dplyr::filter(!id %in% c("OBIS_QLD-Wildnet-7028891", 
                           "OBIS_QLD-Wildnet-6316639", 
                           "OBIS_QLD-Wildnet-7080661", 
                           "GBIF_4465758841")) |> 
  dplyr::mutate(lat =  dplyr::case_when(id == "GBIF_4102881317" ~ -9.9323,
                                        TRUE ~ lat))

# double check it worked
occ_filtered_3 |>
  dplyr::mutate(SourceType = case_when(
    grepl("^GBIF", id) ~ "GBIF",
    grepl("^OBIS", id) ~ "OBIS",
    TRUE               ~ "BOF"
  )) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  #dplyr::filter(year == 2011) |> 
  mapview::mapview(zcol = "SourceType")


# create WB core area for spatial bias check 
cent_sf <- sf::st_as_sf(data.frame(x = 143.87, y = -12.18),
                        coords = c("x","y"), crs = 4326)


bias_results <- spatiotemp_bias(occ.data = occ_filtered_3,
                                temporal.level = c("day", "month", "year"),
                                plot = TRUE,
                                spatial.method = "simple",
                                centroid = cent_sf,
                                radius = 500000,
                                prj = "EPSG:3577"
                                #prj = "+proj=longlat +datum=WGS84"
)

bias_results





# removing undersampled months 
occ_filtered_4 <- occ_filtered_3 |> 
  dplyr::filter(!month %in% c(5, 6, 7, 8))




bias_results <- spatiotemp_bias(occ.data = occ_filtered_4,
                                temporal.level = c("day", "month", "year"),
                                plot = TRUE,
                                spatial.method = "core",
                                centroid = cent_sf,
                                radius = 1,
                                prj = "EPSG:3577"
                                #prj = "+proj=longlat +datum=WGS84"
)

bias_results



str(occ_filtered_4)


### data thinning to avoid spatioremporal correlation 
occ_sight_thin <- spatiotemp_thin(occ.data = occ_filtered_3,   # we use the all years except 2025 and all months -> can still be culled later if models turn out to be struggling with capturing seasonality due to this 
                                  temporal.method = "day",
                                  temporal.dist = 15,
                                  spatial.split.degrees = 1,
                                  spatial.dist = 25000,
                                  iterations = 100
)








# Thinned
par(mfrow = c(1, 2))
terra::plot(terra::vect(occ_sight_thin[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), col = "red", main = "thinned")
terra::plot(terra::vect(occ_filtered_3[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), add = FALSE, main = "sightings")

dev.off()

set.seed(999)




# focus on spatial bias centered around Wreck Bay with 100 km buffer -> because this is the area we will use the randomly sample absence locations from, thus, random sampling on whole region (i.e., simple) for spatial bias doe snot make sense here, as we will adress this bias by using a buffer 
bias_results <- spatiotemp_bias(occ.data = occ_sight_thin,
                                temporal.level = c("day", "month", "year"),
                                plot = TRUE,
                                spatial.method = "core",
                                centroid = cent_sf,
                                radius = 1,
                                prj = "EPSG:3577"
                                #prj = "+proj=longlat +datum=WGS84"
)


bias_results

# now focus on temporal biases (most importan there is month as we will use monht in the models)
bias_results <- spatiotemp_bias(occ.data = occ_sight_thin,
                                temporal.level = c("day", "month", "year"),
                                plot = TRUE,
                                spatial.method = "simple",
                                # centroid = cent_sf,
                                # radius = 100000,
                                prj = "EPSG:3577"
                                #prj = "+proj=longlat +datum=WGS84"
)

bias_results
## temporal bias based on day is not biased! good; go ahead







# Pseudoabsence sampling  -------------------------------------------------



## creating ocean buffer and land barrier for pseudo-absences 
occ_sf <- sf::st_as_sf(occ_sight_thin, coords = c("x", "y"), crs = 4326)

aoi <- occ_sight_thin |>
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
k <- 100                     # pseudo-absences per occurrence

pseudo_abs_buff <- purrr::map_dfr(seq_len(nrow(occ_sf)), function(i) {
  out_i <- spatiotemp_pseudoabs_fix(
    spatial.method  = "buffer",
    temporal.method = "buffer",
    occ.data        = occ_sight_thin[i, , drop = FALSE],     # <- one occurrence at a time
    spatial.ext     = ocean,                         # ocean mask to avoid absences sampled on land 
    spatial.buffer  = c(25000, 500000),              # 10-100 km (10 to avoid absence in same grod than presnece and 100km to consider their potential movement in a day)
    temporal.buffer = 5,                             # sample only the same dates as presence 
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
terra::plot(terra::vect(occ_sight_thin[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"),pch = 21, col = "blue",add=T)
terra::plot(aus, border = "black", col = "grey20", add = TRUE)



## check the distances 
# occurrences with index used by your function
occ_sf <- occ_sight_thin |>
  dplyr::mutate(occ_index = dplyr::row_number()) |>
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326)

# pseudo-absences back in lon/lat
pseudo_sf <- pseudo_abs_buff |>
  sf::st_as_sf(coords = c("x","y"), crs = 4326)

# reorder occurrences to align with each pseudo's occ_index
occ_match <- occ_sf[pseudo_abs_buff$occ_index, ]

# geodesic distance (vectorised, by element)
d_m <- sf::st_distance(pseudo_sf, occ_match, by_element = TRUE)
d_km <- units::set_units(d_m, "km") |> units::drop_units()

# summarise
dplyr::tibble(dist_km = d_km) |>
  dplyr::summarise(min = min(dist_km), q25 = quantile(dist_km, 0.25),
                   med = median(dist_km), q75 = quantile(dist_km, 0.75),
                   max = max(dist_km))



#--> makes sense; distances are betweeen 25 and 500 km 





# adding to occurences 

pseudo_abs_buff$lat <- pseudo_abs_buff$y
pseudo_abs_buff$lon <- pseudo_abs_buff$x
pseudo_abs_buff$PA <- rep(0, nrow(pseudo_abs_buff))


occ_sight_thin <- occ_sight_thin |> 
  dplyr::mutate(PA = 1,
                rep = 0,
                source = "Sightings_record")

str(occ_sight_thin)
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



occ_sight_thin_1 <- occ_sight_thin |> 
  dplyr::select(-geometry)

str(occ_sight_thin_1)
str(pseudo_abs_buff_1)

complete.dataset <- as.data.frame(rbind(occ_sight_thin_1, pseudo_abs_buff_1))
str(complete.dataset)

str(occ_sight)

occ_sight_clean <- occ_sight |>
  dplyr::mutate(
    date = as.POSIXct(date),
    x    = as.numeric(lon),
    y    = as.numeric(lat)
  ) |>
  dplyr::filter(!is.na(date), is.finite(x), is.finite(y)) |>
  dplyr::distinct() |>
  as.data.frame()

str(occ_sight_clean)

anyNA(complete.dataset[, c("date","x","y", "year", "month", "day" )])
anyNA(occ_sight_clean[, c("date","x","y", "year", "month", "day" )])



input <- complete.dataset |> 
  dplyr::select(x, y, year, month, day)
str(input)






# Weights for occurrences based on sampling -------------------------------
## create a (pseudo-)sampling effort file 
### Note, we don't have sampling effort records we can use here, but instead we create a pseudo-effort file baswd on how much was sampled in each monhts, resulting of occurrences in undersampled months getting a higher weight to counter-balance this bias 


crs_m <- 3577     # EPSG:3577 Australian Albers (metres)
cell_km <- 10     # ~25 km cells (tie to your spatial.dist = 100 km; 25–50 km works)
cell_m  <- cell_km * 1000


sight_sf <- occ_sight |>
  dplyr::transmute(
    date  = as.Date(date),
    year  = lubridate::year(date),
    month = lubridate::month(date),
    day   = lubridate::day(date),
    x     = as.numeric(lon),
    y     = as.numeric(lat)
  ) |>
  dplyr::filter(is.finite(x), is.finite(y), !is.na(date)) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 4326) |>
  sf::st_transform(crs_m)

# ---- Build a hex grid over the data extent ----
bbox <- sf::st_as_sfc(sf::st_bbox(sight_sf))
# Expand bbox a bit to avoid edge clipping
bbox_buf <- sf::st_buffer(sf::st_transform(bbox, crs_m), dist = cell_m)
grid_hex <- sf::st_make_grid(bbox_buf, cellsize = cell_m, square = FALSE) |>
  sf::st_as_sf() |>
  dplyr::mutate(cell_id = dplyr::row_number())

names(grid_hex)

# ---- Snap each sighting to a hex cell ----
sight_join <- sf::st_join(sight_sf, grid_hex, left = FALSE)  # drop points outside grid, if any

# ---- Define time unit for effort collapsing (month shown; swap to day/week if needed) ----
sight_join <- sight_join |>
  dplyr::mutate(ym = sprintf("%04d-%02d", year, month))

# ---- One effort event per (cell_id, ym) using cell centroid ----
grid_centroids <- grid_hex |>
  dplyr::mutate(centroid = sf::st_centroid(sf::st_geometry(grid_hex))) |>
  sf::st_set_geometry("centroid")

effort_events_sf <- sight_join |>
  dplyr::distinct(cell_id, ym, year, month) |>
  dplyr::left_join(grid_centroids |> dplyr::select(cell_id, centroid),
                   by = "cell_id") |>
  sf::st_as_sf()


effort_events_sf

# Transform to WGS84
ee_ll <- sf::st_transform(effort_events_sf, 4326)

# Extract coordinates from the *active* geometry (whatever its column is named)
coords <- sf::st_coordinates(ee_ll)

# Add lon/lat, drop geometry, and keep the minimal columns
effort_events_lonlat <- ee_ll |>
  dplyr::mutate(
    lon = coords[, 1],
    lat = coords[, 2]
  ) |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    x = lon,
    y = lat,
    year  = as.numeric(year),
    month = as.numeric(month),
    day   = 15L,        # or the modal/median day per cell–month if you prefer
    day = as.numeric(day)
  ) |>
  as.data.frame()


effort_events_lonlat


effort_events_lonlat |> dplyr::group_by(month) |> 
  dplyr::summarise(N = n())



complete_dataset_weights <- dynamicSDM::spatiotemp_weights(occ.data = input,
                                                           samp.events = effort_events_lonlat,
                                                           spatial.dist = 10000,
                                                           temporal.dist = 30,
                                                           prj = "EPSG:3577")


# merge
complete_dataset_w <- dplyr::left_join(complete.dataset, complete_dataset_weights)
str(complete_dataset_w)

complete_dataset_w <-  complete_dataset_w |> dplyr::filter(lon > -100)


saveRDS(complete_dataset_w, "data/processed/Sightings_PA_w_dynSDM_100_raw_2010_2025.rds")



# # check if still duplicates detected
# 
# chck <- dynamicSDM::spatiotemp_check(
#   complete_dataset_w,
#   duplicate.handle = "exclude", #will deal with those during thinning preocedures
#   date.res = "day"
# )
# 
# # none, so all good!


