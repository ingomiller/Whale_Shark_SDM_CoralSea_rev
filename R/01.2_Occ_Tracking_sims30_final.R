#_____________________________________________________________________________
#                         TRACKING DATA PREP
#_____________________________________________________________________________


# remotes::install_github("r-a-dobson/dynamicSDM", force = TRUE, build_vignettes = TRUE)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------
# thsi si for the final dataset of 30 simulation; everyhting befofee was with 50 simulations to findout how many simualtions are needed to proceed
# locs <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")

locs <- readRDS("data/input/locs_ani_mp_customintrp_daily_PA_sim30_final.rds")


str(locs)
sf::st_crs(locs)

locs |> 
  dplyr::group_by(PA) |> 
  dplyr::summarise(n())

max(locs$date)

locs |> 
  dplyr::group_by(id) |> 
  dplyr::summarise(
    n_pres = sum(PA == 1, na.rm = TRUE),
    n_abs  = sum(PA == 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_pres)) |> 
  print(n=100)



locs <- locs |> 
  # limit to June 2025 due to availability of remote sensing data 
  dplyr::mutate(Date = as.Date(date)) |> 
  dplyr::filter(Date < as.Date("2025-07-01")) |>
  dplyr::mutate(year = year(date),
                month = month(date),
                day = as.numeric(day(date))) |> 
  dplyr::filter(x >= 0 & y <= 180)  |>
  dplyr::select(-Date)





# crop absences to extend of presences +/- 1 degree buffer 

pres <- locs |> dplyr::filter(PA == 1)
if (nrow(pres) == 0) stop("No presences found (PA == 1).")

# 1) Bounding box from presences (ignore NAs)
pad <- 10  # degrees
lon_rng <- range(pres$lon, na.rm = TRUE)
lat_rng <- range(pres$lat, na.rm = TRUE)

lon_min <- lon_rng[1] - 1
lon_max <- lon_rng[2] + pad
lat_min <- lat_rng[1] - pad
lat_max <- lat_rng[2] + pad

# 2) Filter absences by lon/lat ranges
abs_crop <- locs |>
  dplyr::filter(PA == 0) |>
  dplyr::filter(dplyr::between(lon, lon_min, lon_max),
                dplyr::between(lat, lat_min, lat_max))

# 3) Combine back
locs_crop <- dplyr::bind_rows(pres, abs_crop)

# Optional quick check
message("Original: ", nrow(locs), 
        " | Pres: ", nrow(pres),
        " | Abs kept: ", nrow(abs_crop),
        " | New total: ", nrow(locs_crop))

ggplot() +
  geom_point(data = locs, aes(lon, lat), color = "grey80", size = 1) +
  geom_point(data = abs_crop, aes(lon, lat), color = "blue", size = 1) +
  geom_point(data = pres, aes(lon, lat), color = "red", size = 1.5) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude",
       title = "Presences (red) and Cropped Absences (blue)") +
  theme_minimal()



locs <- locs_crop


locs_crop |>
  dplyr::filter(PA == 0) |>
  ggplot2::ggplot(ggplot2::aes(lon, lat)) +
  ggplot2::geom_point(alpha = 0.2, size = 0.5)

# fix ccords
coords <- sf::st_coordinates(locs)  # matrix with X,Y
locs$x <- coords[,1]
locs$y <- coords[,2]




glimpse(locs)

unique(locs$id)
unique(locs$id_split)

max(locs$date)
max(locs$lon)
min(locs$lon)


locs |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(
    P = sum(PA == 1, na.rm = TRUE),
    A = sum(PA == 0, na.rm = TRUE)
  ) |> 
  print(n=100)


str(locs)

track_p <- locs |> 
  # sf::st_drop_geometry() |>
  dplyr::filter(PA == 1)


sims_a <- locs |> 
  # sf::st_drop_geometry() |> 
  dplyr::filter(PA == 0)



terra::plot(terra::vect(track_p))
terra::plot(terra::vect(sims_a))


# 
# ## clip absences to extent of presences
# 
# pres_sf <- st_as_sf(track_p, coords = c("lon","lat"), crs = 4326, remove = FALSE)
# abs_sf  <- st_as_sf(sims_a,  coords = c("lon","lat"), crs = 4326, remove = FALSE)
# 
# pres_sf
# 
# hull <- track_p |>
#   sf::st_union() |>
#   sf::st_convex_hull() |>
#   sf::st_make_valid()
# 
# 
# # 1) Validity check (often fixes "Loop … crosses edge …")
# hull <- sf::st_make_valid(hull)
# 
# # 2) Pick a suitable projected CRS for your area (Australia-wide equal area)
# #    EPSG:3577 (GDA94 / Australian Albers) works well for QLD/Coral Sea
# hull_proj <- hull |>
#   sf::st_transform(3577)
# 
# # 3) Buffer by 100 km (dist is in metres in projected CRS)
# hull_buf_proj <- hull_proj |>
#   sf::st_buffer(dist = 100000)
# 
# # 4) Back to WGS84 for plotting/joins if needed
# hull_buf <- hull_buf_proj |>
#   sf::st_transform(4326)
# 
# 
# # 5) Spatial subset
# inside <- sf::st_within(sims_a, hull_buf, sparse = FALSE)[, 1]
# sims_a_clip <- sims_a[inside, , drop = FALSE]
# 
# 
# 
# 
# terra::plot(terra::vect(track_p))
# terra::plot(terra::vect(sims_a))
# terra::plot(terra::vect(sims_a_clip))
# 
# 
# str(track_p)
# 





# Biases ------------------------------------------------------------------
## we use dynamic SDM to deal spatiotemporal biases 


track_p_df <- track_p |> 
  sf::st_drop_geometry() |> 
  as.data.frame() 
  
sims_a_df <- sims_a |> 
    sf::st_drop_geometry() |> 
  as.data.frame() 
  

str(track_p_df)

track_p_chck <- dynamicSDM::spatiotemp_check(
  track_p_df,
  na.handle = "exclude", # NAs in coordinates or dates 
  duplicate.handle = "exclude", #will deal with those during thinning preocedures
  coord.handle =  "exclude",
  date.handle =  "exclude",
  coordclean = FALSE, # becasue it will flag sea locations -> you bloody terrestrial turds
  date.res = "day",
  coordclean.species = "Rhincodon typus",
  coordclean.handle =  "report"
)

track_p_chck

sims_a_chck <- dynamicSDM::spatiotemp_check(
  sims_a_df,
  na.handle = "exclude", # NAs in coordinates or dates 
  duplicate.handle = "exclude", #will deal with those during thinning preocedures
  coord.handle =  "exclude",
  date.handle =  "exclude",
  coordclean = FALSE, # becasue it will flag sea locations -> you bloody terrestrial turds
  date.res = "day",
  coordclean.species = "Rhincodon typus",
  coordclean.handle =  "report"
)



track_p_cropped <- dynamicSDM::spatiotemp_resolution(occ.data = track_p_chck,
                                                     spatial.res = 4,
                                                     temporal.res = "day")

sims_a_cropped <- dynamicSDM::spatiotemp_resolution(occ.data = sims_a_chck,
                                                    spatial.res = 4,
                                                    temporal.res = "day")




# filter out locations on land 
track_p_cropped_sf <- track_p_cropped |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

sims_a_cropped_sf <- sims_a_cropped |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)



# get land polygons (Natural Earth, 1:50m scale)
sf::sf_use_s2(FALSE)
land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |> 
  sf::st_transform(4326) |>
  sf::st_make_valid() |> 
  sf::st_union()
land

# keep only points NOT on land
on_land <- sf::st_intersects(track_p_cropped_sf, land, sparse = FALSE)[,1]
track_p_filtered_1 <- track_p_cropped_sf[!on_land, ]

str(track_p_filtered_1)

mapview::mapview(track_p_filtered_1)
track_p_filtered_1 <- as.data.frame(track_p_filtered_1)

on_land <- sf::st_intersects(sims_a_cropped_sf, land, sparse = FALSE)[,1]
sims_a_filtered_1 <- sims_a_cropped_sf[!on_land, ]

mapview::mapview(sims_a_filtered_1)
sims_a_filtered_1 <- as.data.frame(sims_a_filtered_1)









# a first check of biases
bias_results <- dynamicSDM::spatiotemp_bias(occ.data = track_p_filtered_1,
                                            temporal.level = c("day", "month", "year"),
                                            plot = TRUE,
                                            spatial.method = "core",
                                            prj = "+proj=longlat +datum=WGS84")

bias_results

# some annual bias but we add this as random effect; spatial bias to be dealt with by ID and later on by cpatial blocking






# Thinning ----------------------------------------------------------------

## data thinning to avoid spatioremporal correlation
track_p_thin <- dynamicSDM::spatiotemp_thin(occ.data = track_p_filtered_1,
                                            temporal.method = "day",
                                            temporal.dist = 1, # or 1 if daily data is used
                                            spatial.split.degrees = 1,
                                            spatial.dist = 10000,
                                            iterations = 5
)



# Thinned
par(mfrow = c(1, 2))
terra::plot(terra::vect(track_p_filtered_1[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), add = FALSE, main = "before spatiotemp thinning", pch =21, cex =0.5)
terra::plot(terra::vect(track_p_thin[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), col = "red", main = "thinned", add = FALSE, pch =21, cex = 0.5)

dev.off()

set.seed(999)


# create WB core area for spatial bias check 
cent_sf <- sf::st_as_sf(data.frame(x = 143.87, y = -12.18),
                        coords = c("x","y"), crs = 4326)

bias_results <- dynamicSDM::spatiotemp_bias(occ.data = track_p_thin,
                                            temporal.level = c("day", "month", "year"),
                                            plot = TRUE,
                                            spatial.method = "core",
                                            centroid = cent_sf,
                                            radius = 1,
                                            prj = "EPSG:3577"
                                            #prj = "+proj=longlat +datum=WGS84"
)

bias_results
bias_results$Plots[4]





str(track_p_filtered_1 )


# thinning per id
track_pres_thin_id <- track_p_filtered_1 |>
  dplyr::group_split(id_split) |>
  purrr::map_dfr(\(.df) {
    occ_in <- .df |>
      dplyr::transmute(
        x     = x, 
        y     = y,
        year  = year,
        month = month,
        day   = day
      )
    
    keep <- dynamicSDM::spatiotemp_thin(
      occ.data               = occ_in,
      temporal.method        = "day",
      temporal.dist          = 1,        
      spatial.split.degrees  = 1,     
      spatial.dist           = 10000,    
      iterations             = 5
    )
    
    # robust join back to original rows 
    keep_key <- keep |>
      dplyr::mutate(.key = paste(x, y, year, month, day, sep = "_")) |>
      dplyr::select(.key)
    
    .df |>
      dplyr::mutate(.key = paste(x, y, year, month, day, sep = "_")) |>
      dplyr::semi_join(keep_key, by = ".key") |>
      dplyr::select(-.key)
  })




bias_results2 <- dynamicSDM::spatiotemp_bias(occ.data = track_pres_thin_id,
                                             temporal.level = c("day", "month", "year"),
                                             plot = TRUE,
                                             spatial.method = "core",
                                             centroid = cent_sf,
                                             radius = 0.5,
                                             prj = "EPSG:3577"
                                             #prj = "+proj=longlat +datum=WGS84"
)


bias_results2
bias_results2$Plots[4]


saveRDS(bias_results2, "outputs/tests/Bias_Results_Presences_sim30_mp_dynamicSDM_final.rds")

# Option a: Thinning absences by removing same days removed from presences --------

# keys of ID×day rows we REMOVED from presences --> to remove the asscoaiated absences 
# pres_drop_keys <- track_p_filtered_1 |>
#   dplyr::filter(PA == 1, rep == 0) |>
#   dplyr::anti_join(track_pres_thin_id, by = c("id","year","month","day")) |>
#   dplyr::distinct(id, year, month, day)




# Bias reduction for Absences  --------------------------------------------


# .get_k_target <- function(df) {
#   k <- df |>
#     dplyr::pull(rep) |>
#     max(na.rm = TRUE)
#   if (is.finite(k)) as.integer(k) else NA_integer_
# }



# filter out locations on land 
sims_a_sf <- sims_a_filtered_1 |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)



# keep only points NOT on land
on_land <- sf::st_intersects(sims_a_sf, land, sparse = FALSE)[,1]
sims_a_filtered <- sims_a_sf[!on_land, ]

mapview::mapview(sims_a_filtered)
sims_a_filtered <- as.data.frame(sims_a_filtered)

# # remove dates that were deleted in presences 
# sims_input <- sims_a_filtered |> dplyr::anti_join(pres_drop_keys, by = c("id","year","month","day"))


# bias_results <- spatiotemp_bias(occ.data = sims_input |> as.data.frame(),
#                                 temporal.level = c("day", "month", "year"),
#                                 plot = TRUE,
#                                 spatial.method = "core",
#                                 centroid = cent_sf,
#                                 radius = 1,
#                                 prj = "EPSG:3577"
#                                 #prj = "+proj=longlat +datum=WGS84"
# )
# 
# 
# bias_results
# 

# str(sims_input)
# str(track_pres_thin_id)
# 
# 
# locs_thinned_PA_1 <- dplyr::bind_rows(track_pres_thin_id, sims_input) |> 
#   dplyr::arrange(id, rep, date)
# 
# glimpse(locs_thinned_PA_1)
# 
# locs_thinned_PA_1 |> dplyr::summarise(min = min(date),
#                                       max = max(date))
# 
# 
# locs_thinned_PA_1 |> 
#   dplyr::group_by(id) |> 
#   dplyr::summarise(
#     n_pres = sum(PA == 1, na.rm = TRUE),
#     n_abs  = sum(PA == 0, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   arrange(desc(n_pres)) |> 
#   print(n=100)

# remove 176407 as only one presence left after thinning procesdures
# locs_thinned_PA_2 <- locs_thinned_PA_1 |> 
#   dplyr::filter(!id %in% c("176407"))
# 
# 
# # saveRDS(locs_thinned_PA_1, "data/processed/Tracking_PA_w_dynSDM_10_raw_2010_2025.rds")
# saveRDS(locs_thinned_PA_2, "data/processed/Tracking_PA_3to10days_w_dynSDM_10_raw_2010_2025.rds")


# Option B: Thinning absences by ID the same way than presences -  --------


# Fast helper: intended k per ID×day from your 'rep' column
get_k <- function(d) {
  k <- suppressWarnings(max(d$replicate, na.rm = TRUE))
  if (is.finite(k)) as.integer(k) else NA_integer_
}

sims_abs_thin_id <- sims_a_filtered |>
  dplyr::group_split(id_split) |>
  # dplyr::group_split(id) |>
  purrr::map_dfr(\(d) {
    
    # 1) collapse exact duplicates (same x,y,year,month,day)
    d_nodup <- d |>
      dplyr::distinct(x, y, year, month, day, .keep_all = TRUE)
    
    # 2) dynamicSDM thinning
    occ <- d_nodup |>
      dplyr::transmute(x = x, y = y, year = year, month = month, day = day)
    
    keep <- dynamicSDM::spatiotemp_thin(
      occ.data              = occ,
      temporal.method        = "day",
      temporal.dist          = 1,
      spatial.split.degrees  = 1,
      spatial.dist           = 10000,
      iterations             = 5
    )
    
    # 3) map back to retained rows
    keep_key <- dplyr::tibble(
      x = round(keep$x, 6),
      y = round(keep$y, 6),
      year = keep$year,
      month = keep$month,
      day = keep$day
    )
    d_key <- d_nodup |>
      dplyr::mutate(x = round(x, 6), y = round(y, 6))
    
    idx <- vctrs::vec_in(
      dplyr::select(d_key, x, y, year, month, day),
      keep_key
    )
    kept <- d_nodup[idx, , drop = FALSE]
    
    # 4) enforce per-day limit based on replicate count
    kept |>
      dplyr::group_by(id_split, year, month, day) |>
      # dplyr::group_by(id, year, month, day) |>
      dplyr::group_modify(\(g, key) {
        k <- get_k(g)
        if (is.finite(k) && nrow(g) > k) dplyr::slice_sample(g, n = k) else g
      }) |>
      dplyr::ungroup()
  })





bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = sims_abs_thin_id,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "core",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)





bias_results_a1


saveRDS(bias_results_a1, "outputs/tests/Bias_Results_Absences_sim30_mp_dynamicSDM_final.rds")
bias_results_a1 <- readRDS("outputs/tests/Bias_Results_Absences_sim30_mp_dynamicSDM_final.rds")


locs_thinned_PA <- dplyr::bind_rows(track_pres_thin_id, sims_abs_thin_id) |> 
  dplyr::arrange(id, rep, date)

glimpse(locs_thinned_PA)


locs_thinned_PA |> 
  dplyr::group_by(id) |> 
  dplyr::summarise(
    n_pres = sum(PA == 1, na.rm = TRUE),
    n_abs  = sum(PA == 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_pres)) |> 
  print(n=100)

locs_thinned_PA |> 
  dplyr::group_by(PA) |> 
  dplyr::summarise(N = n())




## --> we go with absences thinned per ID 




### remove absences within 25 km of a presence 

sp_thinning_PA_overlap <- function(input_df, grid_res = 0.1) {
  library(future.apply)
  library(progressr)
  
  # Drop geometry if input is sf
  if ("sf" %in% class(input_df)) {
    input_df <- sf::st_drop_geometry(input_df)
  }
  
  # Separate presence and absence points
  presence_points <- input_df %>% filter(PA == "1")
  absence_points <- input_df %>% filter(PA == "0")
  
  # Initialize progressor
  p <- progressor(along = seq_len(nrow(presence_points)))
  
  # Function to process each presence point and find overlapping absences
  process_presence_point <- function(i) {
    # Define the grid cell for the current presence point
    lat_range <- c(presence_points$lat[i] - grid_res / 2, presence_points$lat[i] + grid_res / 2)
    lon_range <- c(presence_points$lon[i] - grid_res / 2, presence_points$lon[i] + grid_res / 2)
    
    # Identify absences in the same grid cell and on the same date
    absences_to_remove <- absence_points %>%
      dplyr::filter(date == presence_points$date[i] & 
                      lat >= lat_range[1] & lat <= lat_range[2] & 
                      lon >= lon_range[1] & lon <= lon_range[2])
    
    # Update the progress bar
    p()
    
    # Return the rows to be removed
    return(absences_to_remove)
  }
  
  # Run the thinning process in parallel to identify all absences to remove
  absences_to_remove_list <- future_lapply(seq_len(nrow(presence_points)), process_presence_point)
  
  # Combine all absences to remove into a single dataframe
  absences_to_remove <- do.call(rbind, absences_to_remove_list)
  
  # Remove these absences from the absence data frame
  absence_points <- dplyr::anti_join(absence_points, absences_to_remove, by = c("id", "rep", "date", "lon", "lat"))
  
  # Combine the remaining absences with all presence points
  thinned_data <- dplyr::bind_rows(presence_points, absence_points)
  
  return(thinned_data)
}

future::plan("sequential") # in sequence
progressr::handlers(global = TRUE)
# progressr::handlers("progress")
progressr::handlers("cli")
str(locs_crop)


locs_ol <- sp_thinning_PA_overlap(locs_thinned_PA, grid_res = 0.5) 

locs_ol <- locs_ol |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

locs_ol |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())

locs <- locs_ol



# remove undersampled months
locs |> dplyr::filter(PA==1) |> dplyr::group_by(month) |> 
  dplyr::summarise(n = n()) |>  print(n=50)

# remove undersampled months

mapview::mapview(
  locs  |> dplyr::filter(PA == 1) |> dplyr::filter(month %in% c(7, 8, 9)),
  col.regions = "steelblue1", layer.name = "underrepresented months"
) +
  mapview::mapview(
    locs  |> dplyr::filter(PA == 1) |> dplyr::filter(!month %in% c(7, 8, 9)),
    col.regions = "firebrick", layer.name = "good months"
  ) 

locs_reduced <- locs |> 
  dplyr::filter(!month %in% c(7, 8, 9)) 



locs <- locs_reduced



saveRDS(locs, "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_final_extA_10d.rds")








# Adding buffered random pseudo absences ----------------------------------

## creating ocean buffer and land barrier for pseudo-absences 
occ_sf <- sf::st_as_sf(track_pres_thin_id, coords = c("x", "y"), crs = 4326)

aoi <- track_pres_thin_id |>
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

# ocean <- sf::st_difference(aoi, land_buf) |>
#   sf::st_make_valid()
# plot(ocean)
# 

access_ll <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = 142.5,
      xmax = 180,        # or a smaller xmax if you want
      ymin = -30,        # tweak lat range to your study area
      ymax =   -5
    ),
    crs = sf::st_crs(4326)
  )
)

# transform 
access_aus <- access_ll |>
  sf::st_transform(3577)

# 3) Build ocean as before, then clip it to the accessible region
ocean <- sf::st_difference(aoi, land_buf) |>
  sf::st_make_valid() |>
  sf::st_intersection(access_aus) |>   # <- apply longitude cut-off here
  sf::st_make_valid()

plot(ocean)




set.seed(42)                # for reproducibility
k <- 30                   # pseudo-absences per occurrence

pseudo_abs_buff <- purrr::map_dfr(seq_len(nrow(track_pres_thin_id)), function(i) {
  
  # 1) Get the presence row + its ID
  occ_i <- track_pres_thin_id[i, , drop = FALSE]
  pres_id <- occ_i |>
    dplyr::pull(id_split)         # <-- change "id" to your actual ID column name
  
  # 2) Generate pseudo-absences for this occurrence
  out_i <- spatiotemp_pseudoabs_fix(
    spatial.method  = "buffer",
    temporal.method = "buffer",
    occ.data        = occ_i,                    # one occurrence at a time
    spatial.ext     = ocean,                    # ocean mask to avoid land
    spatial.buffer  = c(25000, 500000),         # 25–500 km
    temporal.buffer = 5,                        # same +/- 5 days as presence
    n.pseudoabs     = k,                        # exactly k per occurrence
    prj             = "+proj=longlat +datum=WGS84"
  )
  
  # 3) Attach index, replicate number, and presence ID
  out_i |>
    dplyr::mutate(
      occ_index = i,                            # numeric index of presence
      rep       = dplyr::row_number(),         # 1..k within presence
      id        = pres_id                      # presence ID propagated
    ) |>
    dplyr::select(x, y, year, month, day, occ_index, id, rep)
})

str(track_pres_thin_id)
str(pseudo_abs_buff)

# plot presences and absences 
aus_sf <- rnaturalearth::ne_countries( scale = 10, returnclass = "sf") |>
  sf::st_transform(4326) 
aus <- terra::vect(aus_sf)

dev.off()
terra::plot(terra::vect(pseudo_abs_buff[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), pch = 21, col = "red", add=F)
plot(ocean, border = "green4", add=TRUE)
terra::plot(terra::vect(track_pres_thin_id[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"),pch = 21, col = "blue",add=T)
terra::plot(aus, border = "black", col = "grey20", add = TRUE)




## Thinning

bias_results <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff,
                                            temporal.level = c("day", "month", "year"),
                                            plot = TRUE,
                                            spatial.method = "core",
                                            prj = "+proj=longlat +datum=WGS84")

bias_results




pseudo_abs_buff_thin <- dynamicSDM::spatiotemp_thin(occ.data = pseudo_abs_buff,
                                            temporal.method = "day",
                                            temporal.dist = 1, # or 1 if daily data is used
                                            spatial.split.degrees = 1,
                                            spatial.dist = 10000,
                                            iterations = 5
)


bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "core",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)

bias_results_a1


get_k <- function(d) {
  k <- suppressWarnings(max(d$replicate, na.rm = TRUE))
  if (is.finite(k)) as.integer(k) else NA_integer_
}

pseudo_abs_buff_thin_id <- pseudo_abs_buff |>
  dplyr::group_split(id) |>
  # dplyr::group_split(id) |>
  purrr::map_dfr(\(d) {
    
    # 1) collapse exact duplicates (same x,y,year,month,day)
    d_nodup <- d |>
      dplyr::distinct(x, y, year, month, day, .keep_all = TRUE)
    
    # 2) dynamicSDM thinning
    occ <- d_nodup |>
      dplyr::transmute(x = x, y = y, year = year, month = month, day = day)
    
    keep <- dynamicSDM::spatiotemp_thin(
      occ.data              = occ,
      temporal.method        = "day",
      temporal.dist          = 1,
      spatial.split.degrees  = 1,
      spatial.dist           = 10000,
      iterations             = 5
    )
    
    # 3) map back to retained rows
    keep_key <- dplyr::tibble(
      x = round(keep$x, 6),
      y = round(keep$y, 6),
      year = keep$year,
      month = keep$month,
      day = keep$day
    )
    d_key <- d_nodup |>
      dplyr::mutate(x = round(x, 6), y = round(y, 6))
    
    idx <- vctrs::vec_in(
      dplyr::select(d_key, x, y, year, month, day),
      keep_key
    )
    kept <- d_nodup[idx, , drop = FALSE]
    
    # 4) enforce per-day limit based on replicate count
    kept |>
      dplyr::group_by(id, year, month, day) |>
      # dplyr::group_by(id, year, month, day) |>
      dplyr::group_modify(\(g, key) {
        k <- get_k(g)
        if (is.finite(k) && nrow(g) > k) dplyr::slice_sample(g, n = k) else g
      }) |>
      dplyr::ungroup()
  })





bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin_id,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "core",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)

bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin_id,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "simple",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)





bias_results_a1






pseudo_abs_buff <- pseudo_abs_buff_thin_id


# adding to occurences 

pseudo_abs_buff$lat <- pseudo_abs_buff$y
pseudo_abs_buff$lon <- pseudo_abs_buff$x
pseudo_abs_buff$PA <- rep(0, nrow(pseudo_abs_buff))

str(track_pres_thin_id)

track_pres_thin_id <- track_pres_thin_id |> 
  dplyr::mutate(PA = 1,
                rep = 0,
                source = "Sightings_record")

str(occ_sight_thin)
str(pseudo_abs_buff)


pseudo_abs_buff_1 <- pseudo_abs_buff |> 
  dplyr::mutate(lon = x, lat = y,
                date = as.POSIXct(
                  paste(year, month, day, sep = "-"),
                  format = "%Y-%m-%d"),
                model = "random",
                id_split =  id, 
                id = if_else(str_detect(id_split, "_\\d{1,2}$"),
                             str_remove(id_split, "_\\d{1,2}$"),
                             id_split),
                PA = 0
                ) |> 
  dplyr::select(-occ_index)


str(pseudo_abs_buff_1)










track_pres_thin_id_1 <- track_pres_thin_id |> 
  dplyr::select(-geometry)

str(track_pres_thin_id_1)
str(pseudo_abs_buff_1)

track_pa_all <- dplyr::bind_rows(
  track_pres_thin_id_1,
  pseudo_abs_buff_1
)
str(track_pa_all)


track_pa_all |> 
  dplyr::filter(PA==0)



anyNA(track_pa_all[, c("date","x","y", "year", "month", "day" )])




### remove absences within 25 km of a presence 

sp_thinning_PA_overlap <- function(input_df, grid_res = 0.1) {
  library(future.apply)
  library(progressr)
  
  # Drop geometry if input is sf
  if ("sf" %in% class(input_df)) {
    input_df <- sf::st_drop_geometry(input_df)
  }
  
  # Separate presence and absence points
  presence_points <- input_df %>% filter(PA == "1")
  absence_points <- input_df %>% filter(PA == "0")
  
  # Initialize progressor
  p <- progressor(along = seq_len(nrow(presence_points)))
  
  # Function to process each presence point and find overlapping absences
  process_presence_point <- function(i) {
    # Define the grid cell for the current presence point
    lat_range <- c(presence_points$lat[i] - grid_res / 2, presence_points$lat[i] + grid_res / 2)
    lon_range <- c(presence_points$lon[i] - grid_res / 2, presence_points$lon[i] + grid_res / 2)
    
    # Identify absences in the same grid cell and on the same date
    absences_to_remove <- absence_points %>%
      dplyr::filter(date == presence_points$date[i] & 
                      lat >= lat_range[1] & lat <= lat_range[2] & 
                      lon >= lon_range[1] & lon <= lon_range[2])
    
    # Update the progress bar
    p()
    
    # Return the rows to be removed
    return(absences_to_remove)
  }
  
  # Run the thinning process in parallel to identify all absences to remove
  absences_to_remove_list <- future_lapply(seq_len(nrow(presence_points)), process_presence_point)
  
  # Combine all absences to remove into a single dataframe
  absences_to_remove <- do.call(rbind, absences_to_remove_list)
  
  # Remove these absences from the absence data frame
  absence_points <- dplyr::anti_join(absence_points, absences_to_remove, by = c("id", "rep", "date", "lon", "lat"))
  
  # Combine the remaining absences with all presence points
  thinned_data <- dplyr::bind_rows(presence_points, absence_points)
  
  return(thinned_data)
}

future::plan("sequential") # in sequence
progressr::handlers(global = TRUE)
# progressr::handlers("progress")
progressr::handlers("cli")
str(locs_crop)


locs_ol <- sp_thinning_PA_overlap(track_pa_all, grid_res = 0.25) 

locs_ol <- locs_ol |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

locs_ol |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())

locs <- locs_ol



# remove undersampled months
locs |> dplyr::filter(PA==1) |> dplyr::group_by(month) |> 
  dplyr::summarise(n = n()) |>  print(n=50)

# remove undersampled months

mapview::mapview(
  locs  |> dplyr::filter(PA == 1) |> dplyr::filter(month %in% c(7, 8, 9)),
  col.regions = "steelblue1", layer.name = "underrepresented months"
) +
  mapview::mapview(
    locs  |> dplyr::filter(PA == 1) |> dplyr::filter(!month %in% c(7, 8, 9)),
    col.regions = "firebrick", layer.name = "good months"
  ) 

locs_reduced <- locs |> 
  dplyr::filter(!month %in% c(7, 8, 9)) 



locs <- locs_reduced

locs |> dplyr::filter(id_split == "172900_1") |> print(n=200)



saveRDS(locs, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final.rds")

locs <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final.rds")

str(locs)

locs |> dplyr::group_by(PA) |> dplyr::summarise(N = n())














# Ading random background  ------------------------------------------------


## creating ocean buffer and land barrier for pseudo-absences 
occ_sf <- sf::st_as_sf(track_pres_thin_id, coords = c("x", "y"), crs = 4326)

aoi <- track_pres_thin_id |>
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

# ocean <- sf::st_difference(aoi, land_buf) |>
#   sf::st_make_valid()
# plot(ocean)
# 

access_ll <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = 142.5,
      xmax = 180,        # or a smaller xmax if you want
      ymin = -30,        # tweak lat range to your study area
      ymax =   -5
    ),
    crs = sf::st_crs(4326)
  )
)

# transform 
access_aus <- access_ll |>
  sf::st_transform(3577)

# 3) Build ocean as before, then clip it to the accessible region
ocean <- sf::st_difference(aoi, land_buf) |>
  sf::st_make_valid() |>
  sf::st_intersection(access_aus) |>   # <- apply longitude cut-off here
  sf::st_make_valid()

plot(ocean)




set.seed(42)                # for reproducibility
k <- 30                   # pseudo-absences per occurrence

pseudo_abs_buff <- purrr::map_dfr(seq_len(nrow(track_pres_thin_id)), function(i) {
  
  # 1) Get the presence row + its ID
  occ_i <- track_pres_thin_id[i, , drop = FALSE]
  pres_id <- occ_i |>
    dplyr::pull(id_split)         # <-- change "id" to your actual ID column name
  
  # 2) Generate pseudo-absences for this occurrence
  out_i <- spatiotemp_pseudoabs_fix(
    spatial.method  = "buffer",
    temporal.method = "buffer",
    occ.data        = occ_i,                    # one occurrence at a time
    spatial.ext     = ocean,                    # ocean mask to avoid land
    spatial.buffer  = c(25000, 500000),         # 25–500 km
    temporal.buffer = 5,                        # same +/- 5 days as presence
    n.pseudoabs     = k,                        # exactly k per occurrence
    prj             = "+proj=longlat +datum=WGS84"
  )
  
  # 3) Attach index, replicate number, and presence ID
  out_i |>
    dplyr::mutate(
      occ_index = i,                            # numeric index of presence
      rep       = dplyr::row_number(),         # 1..k within presence
      id        = pres_id                      # presence ID propagated
    ) |>
    dplyr::select(x, y, year, month, day, occ_index, id, rep)
})

str(track_pres_thin_id)
str(pseudo_abs_buff)

# plot presences and absences 
aus_sf <- rnaturalearth::ne_countries( scale = 10, returnclass = "sf") |>
  sf::st_transform(4326) 
aus <- terra::vect(aus_sf)

dev.off()
terra::plot(terra::vect(pseudo_abs_buff[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"), pch = 21, col = "red", add=F)
plot(ocean, border = "green4", add=TRUE)
terra::plot(terra::vect(track_pres_thin_id[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"),pch = 21, col = "blue",add=T)
terra::plot(aus, border = "black", col = "grey20", add = TRUE)




## Thinning

bias_results <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff,
                                            temporal.level = c("day", "month", "year"),
                                            plot = TRUE,
                                            spatial.method = "core",
                                            prj = "+proj=longlat +datum=WGS84")

bias_results




pseudo_abs_buff_thin <- dynamicSDM::spatiotemp_thin(occ.data = pseudo_abs_buff,
                                                    temporal.method = "day",
                                                    temporal.dist = 1, # or 1 if daily data is used
                                                    spatial.split.degrees = 1,
                                                    spatial.dist = 10000,
                                                    iterations = 5
)


bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "core",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)

bias_results_a1


get_k <- function(d) {
  k <- suppressWarnings(max(d$replicate, na.rm = TRUE))
  if (is.finite(k)) as.integer(k) else NA_integer_
}

pseudo_abs_buff_thin_id <- pseudo_abs_buff |>
  dplyr::group_split(id) |>
  # dplyr::group_split(id) |>
  purrr::map_dfr(\(d) {
    
    # 1) collapse exact duplicates (same x,y,year,month,day)
    d_nodup <- d |>
      dplyr::distinct(x, y, year, month, day, .keep_all = TRUE)
    
    # 2) dynamicSDM thinning
    occ <- d_nodup |>
      dplyr::transmute(x = x, y = y, year = year, month = month, day = day)
    
    keep <- dynamicSDM::spatiotemp_thin(
      occ.data              = occ,
      temporal.method        = "day",
      temporal.dist          = 1,
      spatial.split.degrees  = 1,
      spatial.dist           = 10000,
      iterations             = 5
    )
    
    # 3) map back to retained rows
    keep_key <- dplyr::tibble(
      x = round(keep$x, 6),
      y = round(keep$y, 6),
      year = keep$year,
      month = keep$month,
      day = keep$day
    )
    d_key <- d_nodup |>
      dplyr::mutate(x = round(x, 6), y = round(y, 6))
    
    idx <- vctrs::vec_in(
      dplyr::select(d_key, x, y, year, month, day),
      keep_key
    )
    kept <- d_nodup[idx, , drop = FALSE]
    
    # 4) enforce per-day limit based on replicate count
    kept |>
      dplyr::group_by(id, year, month, day) |>
      # dplyr::group_by(id, year, month, day) |>
      dplyr::group_modify(\(g, key) {
        k <- get_k(g)
        if (is.finite(k) && nrow(g) > k) dplyr::slice_sample(g, n = k) else g
      }) |>
      dplyr::ungroup()
  })





bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin_id,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "core",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)

bias_results_a1 <- dynamicSDM::spatiotemp_bias(occ.data = pseudo_abs_buff_thin_id,
                                               temporal.level = c("day"),
                                               plot = TRUE,
                                               spatial.method = "simple",
                                               centroid = cent_sf,
                                               radius = 0.25,
                                               prj = "EPSG:3577"
                                               #prj = "+proj=longlat +datum=WGS84"
)





bias_results_a1






pseudo_abs_buff <- pseudo_abs_buff_thin_id


# adding to occurences 

pseudo_abs_buff$lat <- pseudo_abs_buff$y
pseudo_abs_buff$lon <- pseudo_abs_buff$x
pseudo_abs_buff$PA <- rep(0, nrow(pseudo_abs_buff))

str(track_pres_thin_id)

track_pres_thin_id <- track_pres_thin_id |> 
  dplyr::mutate(PA = 1,
                rep = 0,
                source = "Sightings_record")

str(occ_sight_thin)
str(pseudo_abs_buff)


pseudo_abs_buff_1 <- pseudo_abs_buff |> 
  dplyr::mutate(lon = x, lat = y,
                date = as.POSIXct(
                  paste(year, month, day, sep = "-"),
                  format = "%Y-%m-%d"),
                model = "random",
                id_split =  id, 
                id = if_else(str_detect(id_split, "_\\d{1,2}$"),
                             str_remove(id_split, "_\\d{1,2}$"),
                             id_split),
                PA = 0
  ) |> 
  dplyr::select(-occ_index)


str(pseudo_abs_buff_1)










track_pres_thin_id_1 <- track_pres_thin_id |> 
  dplyr::select(-geometry)

str(track_pres_thin_id_1)
str(pseudo_abs_buff_1)

track_pa_all <- dplyr::bind_rows(
  track_pres_thin_id_1,
  pseudo_abs_buff_1
)
str(track_pa_all)


track_pa_all |> 
  dplyr::filter(PA==0)



anyNA(track_pa_all[, c("date","x","y", "year", "month", "day" )])




### remove absences within 25 km of a presence 

sp_thinning_PA_overlap <- function(input_df, grid_res = 0.1) {
  library(future.apply)
  library(progressr)
  
  # Drop geometry if input is sf
  if ("sf" %in% class(input_df)) {
    input_df <- sf::st_drop_geometry(input_df)
  }
  
  # Separate presence and absence points
  presence_points <- input_df %>% filter(PA == "1")
  absence_points <- input_df %>% filter(PA == "0")
  
  # Initialize progressor
  p <- progressor(along = seq_len(nrow(presence_points)))
  
  # Function to process each presence point and find overlapping absences
  process_presence_point <- function(i) {
    # Define the grid cell for the current presence point
    lat_range <- c(presence_points$lat[i] - grid_res / 2, presence_points$lat[i] + grid_res / 2)
    lon_range <- c(presence_points$lon[i] - grid_res / 2, presence_points$lon[i] + grid_res / 2)
    
    # Identify absences in the same grid cell and on the same date
    absences_to_remove <- absence_points %>%
      dplyr::filter(date == presence_points$date[i] & 
                      lat >= lat_range[1] & lat <= lat_range[2] & 
                      lon >= lon_range[1] & lon <= lon_range[2])
    
    # Update the progress bar
    p()
    
    # Return the rows to be removed
    return(absences_to_remove)
  }
  
  # Run the thinning process in parallel to identify all absences to remove
  absences_to_remove_list <- future_lapply(seq_len(nrow(presence_points)), process_presence_point)
  
  # Combine all absences to remove into a single dataframe
  absences_to_remove <- do.call(rbind, absences_to_remove_list)
  
  # Remove these absences from the absence data frame
  absence_points <- dplyr::anti_join(absence_points, absences_to_remove, by = c("id", "rep", "date", "lon", "lat"))
  
  # Combine the remaining absences with all presence points
  thinned_data <- dplyr::bind_rows(presence_points, absence_points)
  
  return(thinned_data)
}

future::plan("sequential") # in sequence
progressr::handlers(global = TRUE)
# progressr::handlers("progress")
progressr::handlers("cli")
str(locs_crop)


locs_ol <- sp_thinning_PA_overlap(track_pa_all, grid_res = 0.25) 

locs_ol <- locs_ol |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

locs_ol |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())

locs <- locs_ol



# remove undersampled months
locs |> dplyr::filter(PA==1) |> dplyr::group_by(month) |> 
  dplyr::summarise(n = n()) |>  print(n=50)

# remove undersampled months

mapview::mapview(
  locs  |> dplyr::filter(PA == 1) |> dplyr::filter(month %in% c(7, 8, 9)),
  col.regions = "steelblue1", layer.name = "underrepresented months"
) +
  mapview::mapview(
    locs  |> dplyr::filter(PA == 1) |> dplyr::filter(!month %in% c(7, 8, 9)),
    col.regions = "firebrick", layer.name = "good months"
  ) 

locs_reduced <- locs |> 
  dplyr::filter(!month %in% c(7, 8, 9)) 



locs <- locs_reduced

locs |> dplyr::filter(id_split == "172900_1") |> print(n=200)



saveRDS(locs, "data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final.rds")

locs <- readRDS("data/work_files/Tracks_mp_RandomBuf_30_thinned_2018_2025_final.rds")

str(locs)

locs |> dplyr::group_by(PA) |> dplyr::summarise(N = n())







