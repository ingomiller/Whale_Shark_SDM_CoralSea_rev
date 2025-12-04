#_____________________________________________________________________________
#                               Correlations: Sightings Data 
#_____________________________________________________________________________


# remotes::install_github("r-a-dobson/dynamicSDM", force = TRUE, build_vignettes = TRUE)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full.rds")


str(sight)

# lets make an experiemtn and just sample east coast australia nd. bult a quick sdm from that, maybe for supplements
sight <- sight |>
  dplyr::rename(depth = Depth,
                slope = Slope,
                roughness = Roughness,
                wz = Wz) |> 
  dplyr::filter(lon > 140) |> 
  sf::st_drop_geometry()

str(sight)

vars <- c("depth", "slope", "roughness", "dist2000", "thetao", "uv", "mltost", "chl", "wz")


# check for NAs
keep <- stats::complete.cases(sight[, vars])

sight <- sight |>
  dplyr::filter(keep)


sight <- sight |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

str(sight)
crs(sight)





# autocrrelation ----------------------------------------------------------
sight_occ <- sight |> 
  #dplyr::filter(PA == 1) |> 
  as.data.frame()


ac <- dynamicSDM::spatiotemp_autocorr(sight_occ,
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



## negative spatial autocorrelation detcted for dpeth, sst, roughness, dist2000, uv, chl 
## no temporal ac exept for sst at a monthly scale 
## to account for this, we can use spatial blocks for the modelling rather than random CV



# Spatial Blocks ----------------------------------------------------------

#______________ use dynamic SDM method


# provide raster for spatial blocks 

raster <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/depth_month_0.1.tif")
raster <- raster[[1]]
print(raster)
plot(raster)

# crop raster top extent of occurrences 

sight_m <- sight |>
  sf::st_make_valid() |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_transform("EPSG:3577")

str(sight_m)


v <- terra::vect(sight_m)
v
bb <- terra::ext(v)
buff <- 200000  # ~2° buffer added

bb <- sf::st_bbox(sight)

bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 2),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 2),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 2),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 2)
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
str(sight)
plot(raster_m)





terra::plot(
  raster_m,
  range = c(-9000, 0)
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")



sight_occ_blocks <- dynamicSDM::spatiotemp_block(sight_occ,
                                       vars.to.block.by = vars,
                                       spatial.layer = raster_crop,
                                       spatial.split.degrees = 1,
                                       temporal.block = "month",
                                       n.blocks = 5,
                                       iterations = 10000)



str(sight_occ_blocks)
unique(sight_occ_blocks$REL_SAMP_EFFORT)
unique(sight_occ_blocks$BLOCK.CATS)


sight_occ_blocks.sf <- sight_occ_blocks |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
sight_occ_blocks.sf |> 
  dplyr::select(-Date) |> 
  mapview::mapview(zcol = "BLOCK.CATS", burst = TRUE)


str(sight_occ_blocks.sf)


saveRDS(sight_occ_blocks.sf, "data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
saveRDS(sight_occ_blocks.sf, "data/work_files/Sightings_Validation_data_extract_full_processed.rds")


saveRDS(sight_occ_blocks.sf, "data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed_SUPPS_MODEL.rds")
#_____________ use blockCV package 
library(blockCV)
library(tmap)


# sight_m <- sight_m |> dplyr::rename(id_obs = id)
base::row.names(sight_m) <- as.character(seq_len(nrow(sight_m)))

str(sight_m)

input_occ <- sight_m |> 
  dplyr::transmute(id, 
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
                   wz)
  

#load monthly mean raster stack of continous variables 
input_raster_stack <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1.tif")
input_raster_stack

plot(input_raster_stack)

v <- terra::vect(sight_m)
v
bb <- terra::ext(v)
buff <- 200000  # ~2° buffer added

bb <- sf::st_bbox(sight)

bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 2),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 2),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 2),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 2)
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
str(sight)
plot(raster_m)





terra::plot(
  raster_m[[1]]
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")


input_raster <- raster_m

str(input_occ)


tmap::tm_shape(input_raster) +
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
set.seed(42)
sac1 <- blockCV::cv_spatial_autocor(r = input_raster, 
                           num_sample = 50000,
                           plot = TRUE)
sac1$range




sac2 <- blockCV::cv_spatial_autocor(r = input_raster,
                           x = input_occ, 
                           column =  "PA", 
                           num_sample = 50000,
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




blockCV::cv_block_size(r = input_raster,
              x = input_occ,
              column = "PA",
              min_size = 50000,
              max_size = 1000000)




## extending the raster lightly as the points on edge are cut off otherwise
pts_v   <- terra::vect(input_occ)                    # sf -> SpatVector (same CRS as raster)
bb      <- terra::ext(pts_v)
resxy   <- terra::res(input_raster)
pad_x   <- resxy[1] * 1                              # 1 cell pad (change to 2 if you like)
pad_y   <- resxy[2] * 1

bb_pad  <- terra::ext(bb$xmin - pad_x,
                      bb$xmax + pad_x,
                      bb$ymin - pad_y,
                      bb$ymax + pad_y)

# 2) Union with the current raster extent so we only grow where needed
cur_ext <- terra::ext(input_raster)
new_ext <- terra::ext(
  min(cur_ext$xmin, bb_pad$xmin),
  max(cur_ext$xmax, bb_pad$xmax),
  min(cur_ext$ymin, bb_pad$ymin),
  max(cur_ext$ymax, bb_pad$ymax)
)

# 3) Extend the raster (new cells are NA). Works for multi-layer SpatRaster too.
input_raster_ext <- terra::extend(input_raster, new_ext)


# 4) If MESS needs actual values, drop points that fall on NA across ANY predictor layer
complete_mask <- terra::app(input_raster_ext, fun = function(...) all(!is.na(c(...))))
ok <- terra::extract(complete_mask, pts_v, ID = FALSE)[,1]
pts_ok <- pts_v[!is.na(ok) & ok, ]


str(input_occ)




## Spatial Blocks
# random 

range = 500000 # size of the blocks in metres
# input_raster <- input_raster_ext


sb1 <- blockCV::cv_spatial(x = input_occ,
                  column = "PA", # the response column (binary or multi-class)
                  r = input_raster,
                  k = 5, # number of folds
                  size = range, 
                  hexagon = FALSE,
                  selection = "random", # random blocks-to-fold
                  iteration = 100, # find evenly dispersed folds
                  progress = TRUE,
                  seed = 666,
                  biomod2 = FALSE) # also create folds for biomod2


blockCV::cv_plot(cv = sb1, 
        x = input_occ,
        r = input_raster,
        nrow = 2, 
        points_alpha = 0.5)

pal <- cmocean::cmocean("balance")(201)
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





# checkerboard block to CV fold assignment
sb2 <- blockCV::cv_spatial(x = input_occ,
                  column = "PA", # the response column (binary or multi-class)
                  r = input_raster,
                  hexagon = FALSE,
                  size = range, # size of the blocks in metres
                  selection = "checkerboard", 
                  iteration = 500, # find evenly dispersed folds
                  progress = TRUE,
                  biomod2 = FALSE) # also create folds for biomod2


tm_shape(sb2$blocks) +
  tm_fill(
    fill = "folds",
    fill.scale = tm_scale_categorical()
  )

blockCV::cv_plot(cv = sb2, 
        x = input_occ,
        r = input_raster,
        nrow = 2, 
        points_alpha = 0.5)


sb2_sim <- blockCV::cv_similarity(cv = sb2,
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

sb2_sim



sb3 <- blockCV::cv_spatial(x = input_occ,
                  column = "PA", # the response column (binary or multi-class)
                  r = input_raster,
                  k = 5, # number of folds
                  hexagon = FALSE,
                  size = range, # size of the blocks in metres
                  selection = "systematic", 
                  iteration = 500, # find evenly dispersed folds
                  seed = 666,
                  progress = TRUE,
                  biomod2 = FALSE) # also create folds for biomod2


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


sb1$records
sb3$records

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


# # Buffering LOO (also known as Spatial LOO)
# bloo <- cv_buffer(x = input_occ,
#                   column = "PA",
#                   size = 350000)
# 
# bloo
# 
# cv_plot(cv = bloo, 
#         x = input_occ,
#         r = input_raster,
#         nrow = 2, 
#         points_alpha = 0.5,
#         num_plots = c(1, 50, 100)) # only show folds 1, 50 and 100
# 
# 
# 
# 
# # Nearest Neighbour Distance Matching (NNDM) LOO
# nncv <- cv_nndm(x = input_occ,
#                 column = "PA",
#                 r = input_raster,
#                 size = 350000,
#                 num_sample = 5000, 
#                 sampling = "regular",
#                 min_train = 0.1,
#                 plot = TRUE)
# 
# nncv








# Save  -------------------------------------------------------------------

saveRDS(sb2, "data/processed/BlockCV_spatial_folds_sightings_groups.rds")


saveRDS(sb1, "data/processed/BlockCV_spatial_folds_sightings_SUPPS_MODEL.rds")
