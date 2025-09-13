#_____________________________________________________________________________
#                               Correlations / Spatial Blocks: Sightings Data 
#_____________________________________________________________________________


# remotes::install_github("r-a-dobson/dynamicSDM", force = TRUE, build_vignettes = TRUE)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


tracks <- readRDS("data/processed/Tracks_PA_w_dynSDM_10_2010_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_3days_dynSDM_10_2010_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2010_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_seamounts.rds")


str(tracks)


tracks <- tracks |>
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |> 
  sf::st_drop_geometry()

str(tracks)

tracks |> distinct(id) |>  nrow()
tracks

vars <- c("depth", "slope", "roughness", "dist2000", "thetao", "uv", "mltost", "chl", "wz", "dist_seamount", "dist_knoll")


# check for NAs
keep <- stats::complete.cases(tracks[, vars])

tracks <- tracks |>
  dplyr::filter(keep) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 


mapview::mapview(
  tracks |> dplyr::select(-Date) |> dplyr::filter(PA == 0),
  col.regions = "steelblue1", layer.name = "Absences"
) +
  mapview::mapview(
    tracks |> dplyr::select(-Date) |> dplyr::filter(PA == 1),
    col.regions = "firebrick1", layer.name = "Presences"
  )

str(tracks)
crs(tracks)



mapview::mapview(
  tracks |> dplyr::select(-Date) |> dplyr::filter(PA == 1) |> dplyr::filter(month %in% c(7, 8, 9)),
  col.regions = "steelblue1", layer.name = "underrepresented months"
) +
  mapview::mapview(
    tracks |> dplyr::select(-Date) |> dplyr::filter(PA == 1) |> dplyr::filter(!month %in% c(7, 8, 9)),
    col.regions = "firebrick", layer.name = "good months"
  ) 

mapview::mapview(
  tracks |> dplyr::select(-Date) |> dplyr::filter(PA == 1) |> dplyr::filter(dist2000 > 300000),
  col.regions = "steelblue1", layer.name = "underrepresented months"
) 
  


## removing extreme outliers

plots <- tracks |> dplyr::mutate(PA = as.factor(PA))


(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


tracks_out <- tracks |> 
  dplyr::filter(
    # between(chl, mean(chl, na.rm = TRUE) - 4 * sd(chl, na.rm = TRUE), mean(chl, na.rm = TRUE) + 4 * sd(chl, na.rm = TRUE)),
    between(uv, mean(uv, na.rm = TRUE) - 4 * sd(uv, na.rm = TRUE), mean(uv, na.rm = TRUE) + 4 * sd(uv, na.rm = TRUE)),
    between(mltost, mean(mltost, na.rm = TRUE) - 4 * sd(mltost, na.rm = TRUE), mean(mltost, na.rm = TRUE) + 4 * sd(mltost, na.rm = TRUE)),
    between(roughness, mean(roughness, na.rm = TRUE) - 4 * sd(roughness, na.rm = TRUE), mean(roughness, na.rm = TRUE) + 4 * sd(roughness, na.rm = TRUE)),
    between(slope, mean(slope, na.rm = TRUE) - 4 * sd(slope, na.rm = TRUE), mean(slope, na.rm = TRUE) + 4 * sd(slope, na.rm = TRUE)),
    # between(dist2000, mean(dist2000, na.rm = TRUE) - 4 * sd(dist2000, na.rm = TRUE), mean(dist2000, na.rm = TRUE) + 4 * sd(dist2000, na.rm = TRUE)),
    between(wz, mean(wz, na.rm = TRUE) - 4 * sd(wz, na.rm = TRUE), mean(wz, na.rm = TRUE) + 4 * sd(wz, na.rm = TRUE))
  )

plots <- tracks_out |> dplyr::mutate(PA = as.factor(PA))


(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 200) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


# tracks_filtered <- tracks |> 
#   dplyr::filter(lon < 160)
# 
# tracks <- tracks_filtered

# crop tracks to extend of raster 
raster <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")

# ext  <- terra::ext(raster)
# ext
# bbox <- sf::st_bbox(
#   c(
#     xmin = terra::xmin(ext),
#     ymin = terra::ymin(ext),
#     xmax = terra::xmax(ext),
#     ymax = terra::ymax(ext)
#   ),
#   crs = sf::st_crs(tracks)
# )
# 
# roi <- sf::st_as_sfc(bbox)
# roi
# 
# tracks <- tracks |> sf::st_crop(roi)


# crop absences to extend of presences +/- 1 degree buffer 

pres <- tracks_out |> dplyr::filter(PA == 1)
if (nrow(pres) == 0) stop("No presences found (PA == 1).")

# 1) Bounding box from presences (ignore NAs)
pad <- 1  # degrees
lon_rng <- range(pres$lon, na.rm = TRUE)
lat_rng <- range(pres$lat, na.rm = TRUE)

lon_min <- lon_rng[1] - pad
lon_max <- lon_rng[2] + pad
lat_min <- lat_rng[1] - pad
lat_max <- lat_rng[2] + pad

# 2) Filter absences by lon/lat ranges
abs_crop <- tracks_out |>
  dplyr::filter(PA == 0) |>
  dplyr::filter(dplyr::between(lon, lon_min, lon_max),
                dplyr::between(lat, lat_min, lat_max))

# 3) Combine back
tracks_crop <- dplyr::bind_rows(pres, abs_crop)

# Optional quick check
message("Original: ", nrow(tracks), 
        " | Pres: ", nrow(pres),
        " | Abs kept: ", nrow(abs_crop),
        " | New total: ", nrow(tracks_crop))

ggplot() +
  geom_point(data = tracks_out, aes(lon, lat), color = "grey80", size = 1) +
  geom_point(data = abs_crop, aes(lon, lat), color = "blue", size = 1) +
  geom_point(data = pres, aes(lon, lat), color = "red", size = 1.5) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude",
       title = "Presences (red) and Cropped Absences (blue)") +
  theme_minimal()



tracks_reduced <- tracks_crop |> 
  dplyr::filter(!month %in% c(7, 8, 9)) 




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
str(tracks_crop)

# tracks_ol <- sp_thinning_PA_overlap(tracks_crop, grid_res = 0.5)
tracks_ol <- sp_thinning_PA_overlap(tracks_reduced, grid_res = 0.5)
tracks_ol <- tracks_ol |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 
  






# autocrrelation ----------------------------------------------------------
tracks_occ <- tracks_ol |> 
  dplyr::filter(PA == 1) |> 
  as.data.frame()


ac <- dynamicSDM::spatiotemp_autocorr(tracks_occ,
                                      varname = vars,
                                      plot = TRUE,
                                      temporal.level = c("day", "month", "year")) 


# ac$Plots

## extract spatial tests results 
spatial_tests <- ac$Statistical_tests |>
  purrr::imap(~ .x$Spatial_autocorrelation |> dplyr::mutate(variable = .y)) |>
  dplyr::bind_rows() |>
  dplyr::relocate(variable, .before = dplyr::everything()) |> 
  dplyr::mutate(
    direction = dplyr::case_when(
      observed > expected ~ "positive",
      observed < expected ~ "negative",
      TRUE ~ "≈ expected"
    ),
    sig_raw = p.value < 0.05,
    p.adj   = stats::p.adjust(p.value, method = "BH"),
    sig_adj = p.adj < 0.05
  ) |>
  dplyr::arrange(p.adj)

## extract temporal test results
temporal_tests <- ac$Statistical_tests |>
  purrr::imap(function(var_list, var_name) {
    ta <- var_list$Temporal_autocorrelation
    if (base::is.null(ta)) return(NULL)
    
    c("day", "month", "year") |>
      purrr::map(function(ts) {
        h <- ta[[ts]]
        if (base::is.null(h)) return(NULL)
        
        tibble::tibble(
          variable   = var_name,
          timescale  = ts,
          statistic  = base::as.numeric(h$statistic),
          df         = base::as.integer(h$parameter),
          p.value    = h$p.value,
          estimate   = base::as.numeric(h$estimate),
          conf.low   = base::as.numeric(h$conf.int[1]),
          conf.high  = base::as.numeric(h$conf.int[2]),
          conf.level = base::attr(h$conf.int, "conf.level"),
          method     = h$method,
          data.name  = h$data.name
        )
      }) |>
      purrr::compact() |>
      dplyr::bind_rows()
  }) |>
  purrr::compact() |>
  dplyr::bind_rows() |>
  dplyr::mutate(timescale = base::factor(timescale, levels = c("day", "month", "year"))) |>
  dplyr::arrange(variable, timescale)

print(spatial_tests)
print(temporal_tests, n=50)



## negative spatial autocorrelation detcted 
## some temporal ac on daily scale to be dealt with 
## to account for this, we can use spatial blocks for the modelling rather than random CV



# Spatial Blocks ----------------------------------------------------------

#______________ use dynamic SDM method


# provide raster for spatial blocks 

raster <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")
raster <- raster[[1]]
print(raster)
plot(raster)

# crop raster to extent of occurrences 

tracks_m <- tracks_ol |>
  sf::st_make_valid() |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_transform("EPSG:3577")

str(tracks_m)


v <- terra::vect(tracks_m)
v
bb <- terra::ext(v)
buff <- 100000  # ~2° buffer added

bb <- sf::st_bbox(tracks_ol)

bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 1),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 1),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 1),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 1)
)

bbox_poly_wgs <- sf::st_as_sfc(sf::st_bbox(bb_exp, crs = sf::st_crs(4326)))
bbox_vect_wgs <- terra::vect(bbox_poly_wgs)
bbox_vect_r   <- bbox_vect_wgs |> terra::project(terra::crs(raster))

raster_crop <- raster |> terra::crop(bbox_vect_r, snap = "out")
raster_crop[raster_crop >= 0] <- NA
terra::crs(raster_crop) <- "EPSG:4326"

# reporject to metric 
raster_m <- raster_crop |>
  terra::project("EPSG:3577", method = "bilinear")


print(raster_m)
str(tracks)
plot(raster_m)


terra::plot(
  raster_m,
  range = c(-11000, 0)
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")



tracks_occ <- tracks_ol |> 
  as.data.frame()

tracks_occ_blocks <- dynamicSDM::spatiotemp_block(tracks_occ,
                                                 vars.to.block.by = vars,
                                                 spatial.layer = raster_crop,
                                                 spatial.split.degrees = 1,
                                                 temporal.block = "month",
                                                 n.blocks = 5,
                                                 iterations = 10000)



str(tracks_occ_blocks)
unique(tracks_occ_blocks$BLOCK.CATS)


tracks_occ_blocks.sf <- tracks_occ_blocks |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

tracks_occ_blocks.sf |> 
  dplyr::select(-Date) |> 
  mapview::mapview(zcol = "BLOCK.CATS", burst = TRUE)


str(tracks_occ_blocks.sf)


saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_dynSDM_10_2018_2025_extract_processed.rds")
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_3days_dynSDM_10_2018_2025_extract_processed.rds")
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract_processed.rds")
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")



#_____________ use blockCV package 
library(blockCV)
library(tmap)


tracks_m <- tracks_m |> dplyr::rename(id_obs = id)
base::row.names(tracks_m) <- as.character(seq_len(nrow(tracks_m)))

str(tracks_m)

input_occ <- tracks_m |> 
  dplyr::transmute(id_obs, 
                   date, 
                   lon, 
                   lat, 
                   PA, 
                   depth,
                   slope,
                   roughness, 
                   dist2000, 
                   thetao, 
                   uv, 
                   mltost, 
                   chl, 
                   wz,
                   dist_seamount,
                   dist_knoll)




#load monthly mean raster stack of continous variables 
input_raster_stack <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext.tif")
input_raster_stack <- input_raster_stack[[!names(input_raster_stack) %in% c("lat", "lon", "month", "id")]]

plot(input_raster_stack)

v <- terra::vect(tracks_m)
v
bb <- terra::ext(v)
buff <- 100000  # ~2° buffer added

bb <- sf::st_bbox(tracks_ol)

bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 1),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 1),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 1),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 1)
)

bbox_poly_wgs <- sf::st_as_sfc(sf::st_bbox(bb_exp, crs = sf::st_crs(4326)))
bbox_vect_wgs <- terra::vect(bbox_poly_wgs)
bbox_vect_r   <- bbox_vect_wgs |> terra::project(terra::crs(input_raster_stack))

raster_crop <- input_raster_stack |> terra::crop(bbox_vect_r, snap = "out")
# raster_crop[raster_crop >= 0] <- NA
terra::crs(raster_crop) <- "EPSG:4326"

# reporject to metric 
raster_m <- raster_crop |>
  terra::project("EPSG:3577", method = "bilinear")


print(raster_m)
str(tracks)
plot(raster_m)





terra::plot(
  raster_m[[1]]
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")


input_raster <- raster_m
input_raster
str(input_occ)


tmap::tm_shape(input_raster[[names(input_raster) != "id"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = gray.colors(10)),
    col.legend = tm_legend_hide()
  ) +
  tm_shape(input_occ) +
  tm_dots(
    fill = "PA",
    fill.scale = tm_scale_categorical(),
    size = 0.5,
    fill_alpha = 0.5
  )



# test spatial autocreelation distance
set.seed(44)
sac1 <- blockCV::cv_spatial_autocor(r = input_raster, 
                                    num_sample = 25000,
                                    plot = TRUE)
sac1$range




sac2 <- blockCV::cv_spatial_autocor(r = input_raster,
                                    x = input_occ, 
                                    column =  "PA", 
                                    num_sample = 25000,
                                    plot = TRUE)
sac2$range

str(sac2)


sac3 <- blockCV::cv_spatial_autocor(x = input_occ,
                                    column = "PA",
                                    plot = TRUE)
sac3$range

library(automap)
sac2$variograms
plot(sac1$variograms[[1]])
plot(sac2$variograms[[1]])
plot(sac3$variograms[[1]])




blockCV::cv_block_size(r = input_raster,
                       x = input_occ,
                       column = "PA",
                       min_size = 10000,
                       max_size = 1000000)








## Spatial Blocks
# random 
range = 50000
# range = 250000 # size of the blocks in metres
# range = 500000
pal <- cmocean::cmocean("balance")(201)

sb1 <- blockCV::cv_spatial(x = input_occ,
                           column = "PA", # the response column (binary or multi-class)
                           r = input_raster,
                           k = 5, # number of folds
                           size = range, 
                           hexagon = FALSE,
                           selection = "random", # random blocks-to-fold
                           iteration = 500, # find evenly dispersed folds
                           progress = TRUE,
                           seed = 666,
                           biomod2 = FALSE) # also create folds for biomod2
sb1$records

blockCV::cv_plot(cv = sb1, 
                 x = input_occ,
                 r = input_raster,
                 nrow = 2, 
                 points_alpha = 0.5)


sb1_sim <- blockCV::cv_similarity(cv = sb1,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:10,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

sb1_sim

sb1$folds_list





# # checkerboard block to CV fold assignment
# sb2 <- blockCV::cv_spatial(x = input_occ,
#                            column = "PA", # the response column (binary or multi-class)
#                            r = input_raster,
#                            hexagon = FALSE,
#                            size = range, # size of the blocks in metres
#                            selection = "checkerboard", 
#                            iteration = 500, # find evenly dispersed folds
#                            progress = TRUE,
#                            seed = 777,
#                            biomod2 = FALSE) # also create folds for biomod2
# 
# sb2$records
# tm_shape(sb2$blocks) +
#   tm_fill(
#     fill = "folds",
#     fill.scale = tm_scale_categorical()
#   )
# 
# blockCV::cv_plot(cv = sb2, 
#                  x = input_occ,
#                  r = input_raster,
#                  nrow = 2, 
#                  points_alpha = 0.5)
# 
# 
# sb2_sim <- blockCV::cv_similarity(cv = sb2,
#                                   x = input_occ,
#                                   r = input_raster,
#                                   num_plot = 1:10,
#                                   method = "MESS",
#                                   num_sample = 10000,
#                                   jitter_width = 0.2,
#                                   points_size = 1,
#                                   points_alpha = 0.5,
#                                   points_colors = pal,
#                                   progress = TRUE)
# 
# sb2_sim
# 


sb3 <- blockCV::cv_spatial(x = input_occ,
                           column = "PA", # the response column (binary or multi-class)
                           r = input_raster,
                           k = 5, # number of folds
                           hexagon = FALSE,
                           size = range, # size of the blocks in metres
                           selection = "systematic", 
                           iteration = 1000, # find evenly dispersed folds
                           seed = 888,
                           progress = TRUE,
                           biomod2 = FALSE) # also create folds for biomod2
sb3$records

tm_shape(sb3$blocks) +
  tm_fill(
    fill = "folds",
    fill.scale = tm_scale_categorical()
  )

blockCV::cv_plot(cv = sb3, 
                 x = input_occ,
                 r = input_raster,
                 nrow = 2, 
                 points_alpha = 0.5)

sb3_sim <- blockCV::cv_similarity(cv = sb3,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:10,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

sb3_sim





# spatial clustering
set.seed(6)
scv <- blockCV::cv_cluster(x = input_occ,
                           column = "PA", # optional: counting number of train/test records
                           k = 5,
                           biomod2 = FALSE)

scv


blockCV::cv_plot(cv = scv, 
                 x = input_occ,
                 r = input_raster,
                 nrow = 2, 
                 points_alpha = 0.5)


scv_sim <- blockCV::cv_similarity(cv = scv,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:10,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

scv_sim



# Nearest Neighbour Distance Matching (NNDM) LOO
nncv <- blockCV::cv_nndm(x = input_occ,
                column = "PA",
                r = input_raster,
                size = range,
                num_sample = 5000,
                sampling = "random",
                min_train = 0.25,
                plot = TRUE)

plot(nncv)



blockCV::cv_plot(cv = nncv,
                 x = input_occ,
                 r = input_raster,
                 nrow = 2,
                 points_alpha = 0.5)




nncv_sim <- blockCV::cv_similarity(cv = nncv,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:10,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

nncv_sim


# Save  -------------------------------------------------------------------

saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking.rds")
saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_3days.rds")
saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_3to7days.rds")
saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_1to4days.rds")
saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")
