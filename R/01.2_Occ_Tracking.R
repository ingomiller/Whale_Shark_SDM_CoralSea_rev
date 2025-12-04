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

# locs <- readRDS("data/input/locs_ani_crw_customintrp_PA_sim10.rds")
# 3 to 7 day interpolated tracks (3 for Spota and 7 for double tagged individuals with PAT tags)

locs <- readRDS("data/input/locs_ani_mp_customintrp_daily_PA_sim50.rds")

locs <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl_wz.rds")
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



# fix ccords
coords <- sf::st_coordinates(locs)  # matrix with X,Y
locs$x <- coords[,1]
locs$y <- coords[,2]



locs <- locs |> 
  # limit to June 2025 due to availability of remote sensing data 
  dplyr::filter(Date < as.Date("2025-07-01")) |> 
  dplyr::mutate(lat = y,
                lon = x,
                year = year(date),
                month = month(date),
                day = as.numeric(day(date))) |> 
  dplyr::filter(x >= 0 & y <= 180)  |>
  dplyr::select(-Date)

glimpse(locs)

unique(locs$id)
unique(locs$id_split)

max(locs$date)
max(locs$x)
min(locs$y)


locs |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(
    P = sum(PA == 1, na.rm = TRUE),
    A = sum(PA == 0, na.rm = TRUE)
  ) |> 
  print(n=100)

## looks like we somehow added some track segmentrs that shouldn't be there (i.e, were not foraging and thus have no simulated data; needs to be removed now)
valid_ids <- locs |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(A = sum(PA == 0, na.rm = TRUE)) |>
  dplyr::filter(A > 0) |>
  dplyr::pull(id_split)

valid_ids

locs_filtered <- locs |>
  dplyr::filter(id_split %in% valid_ids)

locs_filtered |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(
    P = sum(PA == 1, na.rm = TRUE),
    A = sum(PA == 0, na.rm = TRUE)
  ) |> 
  print(n=100)

locs <- locs_filtered

mapview::mapview(locs |> dplyr::filter(PA ==1))

names(locs)
table(locs$PA, useNA = "ifany")
summary(locs$rep)


# using first 30 simulations 

locs1 <- locs |>
  dplyr::group_by(id_split) |>
  dplyr::mutate(
    replicate = dplyr::case_when(
      # For pseudo-absences, re-index replicate numbers (1–200)
      PA == 0 ~ as.integer(factor(rep, levels = sort(unique(rep[PA == 0])))),
      # For presences, just keep 0 so they stay included
      PA == 1 ~ 0L
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(replicate <= 30 | replicate == 0) |>  # keep first 30 reps + presences
  as.data.frame()

locs1 |> dplyr::group_by(PA) |> 
  dplyr::summarise(N= n())

max(locs1$replicate)

locs <- locs1
str(locs)


track_p <- locs |> 
  dplyr::filter(PA == 1)

track_p |> rstatix::get_summary_stats(lon, lat)


sims_a <- locs |> 
  dplyr::filter(PA == 0)

sims_a |> rstatix::get_summary_stats(lat, lon)

terra::plot(terra::vect(track_p[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))

terra::plot(terra::vect(sims_a[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))



## clip absences to extent of presences

pres_sf <- sf::st_as_sf(track_p, coords = c("lon","lat"), crs = 4326, remove = FALSE)
abs_sf  <- sf::st_as_sf(sims_a,  coords = c("lon","lat"), crs = 4326, remove = FALSE)


hull <- pres_sf |>
  sf::st_union() |>
  sf::st_convex_hull() |>
  sf::st_make_valid()


# 1) Validity check (often fixes "Loop … crosses edge …")
hull <- sf::st_make_valid(hull)

# 2) Pick a suitable projected CRS for your area (Australia-wide equal area)
#    EPSG:3577 (GDA94 / Australian Albers) works well for QLD/Coral Sea
hull_proj <- hull |>
  sf::st_transform(3577)

# 3) Buffer by 100 km (dist is in metres in projected CRS)
hull_buf_proj <- hull_proj |>
  sf::st_buffer(dist = 100000)

# 4) Back to WGS84 for plotting/joins if needed
hull_buf <- hull_buf_proj |>
  sf::st_transform(4326)

# 5) Make your absences match CRS before spatial predicates
abs_sf <- abs_sf |>
  sf::st_transform(4326)

# 6) Spatial subset
inside <- sf::st_within(abs_sf, hull_buf, sparse = FALSE)[, 1]
sims_a_clip <- sims_a[inside, , drop = FALSE]




terra::plot(terra::vect(track_p[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))

terra::plot(terra::vect(sims_a[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))
terra::plot(terra::vect(sims_a_clip[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84 +no_defs"))


# Biases ------------------------------------------------------------------
## we use dynamic SDM to deal spatiotemporal biases 

track_p_chck <- dynamicSDM::spatiotemp_check(
  track_p,
  na.handle = "exclude", # NAs in coordinates or dates 
  duplicate.handle = "exclude", #will deal with those during thinning preocedures
  coord.handle =  "exclude",
  date.handle =  "exclude",
  coordclean = FALSE, # becasue it will flag sea locations -> you bloody terrestrial turds
  date.res = "day",
  coordclean.species = "Rhincodon typus",
  coordclean.handle =  "report"
)

sims_a_chck <- dynamicSDM::spatiotemp_check(
  sims_a_clip,
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
    
    # robust join back to original rows (avoid float equality pitfalls)
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


saveRDS(bias_results2, "outputs/tests/Bias_Results_Presences_sim30_mp_dynamicSDM.rds")

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


saveRDS(bias_results_a1, "outputs/tests/Bias_Results_Absences_sim30_mp_dynamicSDM.rds")


locs_thinned_PA <- dplyr::bind_rows(track_pres_thin_id, sims_abs_thin_id) |> 
  dplyr::arrange(id, rep, date)

glimpse(locs_thinned_PA)
str(locs_thinned_PA)

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

mapview::mapview(locs_thinned_PA |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE))


### ---> decide what approach to go with:
# a) absences thinned the same way which results in less than 10 reps per id/date 
# b) only cull absecnces based on teh oens deleted in presences and use weights instead; preserved 1:10 ratio and will be it ba;lanced 

## --> we go with absences thinned per ID 



saveRDS(locs_thinned_PA, "data/work_files/Tracks_sims_50_thinned_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
saveRDS(locs_thinned_PA, "data/work_files/Tracks_mp_sims_30_thinned_2018_2025_bathy_dist_sst_uv.curr_mld_chl.rds")





