#_____________________________________________________________________________
#                               Correlations / Spatial Blocks: Sightings Data 
#_____________________________________________________________________________


# remotes::install_github("r-a-dobson/dynamicSDM", force = TRUE, build_vignettes = TRUE)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
library(patchwork)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


# tracks <- readRDS("data/processed/Tracks_PA_w_dynSDM_10_2010_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_3days_dynSDM_10_2010_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2010_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_seamounts.rds")

# remove PAT tags from tehse and run spatial blocks again 
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")


# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final.rds")

tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_th_final.rds")


# using random buffer absences dynamicSDM
# tracks <- readRDS("data/processed/Tracks_PA_mp_RandomBuf_30_thinned_2018_2025_extract_final.rds")





str(tracks_crw)
str(tracks)


# tracks <- tracks_crw

tracks <- tracks |>
  dplyr::rename(depth = Depth,
                slope = Slope,
                roughness = Roughness,
                wz = Wz) |>
  # remove positive deoths 
  dplyr::filter(!depth >= 0) |>
  dplyr::filter(month != 10) |> 
  sf::st_drop_geometry()


tracks |> dplyr::filter(month == 10)
  


# cut off some of teh souther extent 
# tracks1 <- tracks |> 
#   dplyr::filter(lat > -27)
# 
# # tracks_old |> 
# #   as.data.frame() |>
# #   dplyr::group_by(PA) |> 
# #   rstatix::get_summary_stats(thetao, type = "common")
# 
# tracks1 |> 
#   as.data.frame() |>
#   dplyr::group_by(PA) |> 
#   rstatix::get_summary_stats(thetao, type = "common")
# 
# tracks <- tracks1

both <- dplyr::bind_rows(
  tracks_crw |> dplyr::mutate(dataset = "crw"),
  tracks     |> dplyr::mutate(dataset = "random buffered bg")
)

# Density / histogram per PA and dataset
both |>
  ggplot2::ggplot(
    ggplot2::aes(x = thetao, fill = factor(PA), colour = factor(PA))
  ) +
  ggplot2::geom_density(alpha = 0.3, trim = TRUE) +
  ggplot2::facet_wrap(~ dataset, ncol = 1)

tracks |>
  ggplot2::ggplot(ggplot2::aes(x = thetao, y = lat, colour = factor(PA))) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::scale_y_continuous()

tracks_old |>
  ggplot2::ggplot(ggplot2::aes(x = thetao, y = lat, colour = factor(PA))) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::scale_y_continuous()

tracks |>
  ggplot2::ggplot(ggplot2::aes(x = thetao, y = lon, colour = factor(PA))) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::scale_y_continuous()

tracks_old |>
  ggplot2::ggplot(ggplot2::aes(x = thetao, y = lon, colour = factor(PA))) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::scale_y_continuous()



tracks_old |> dplyr::filter(thetao < 25)
tracks |> dplyr::filter(thetao < 26)
unique(tracks$id_split)
unique(tracks_old$id_split)

tracks |>  dplyr::filter(id_split == "252215_9") |> dplyr::select(Date, lat, lon, PA, rep, thetao) 
tracks_old |>  dplyr::filter(id_split == "252215_9") |> dplyr::select(Date, lat, lon, PA, rep, thetao) 

# id252215_9 <- tracks_old |> 
#   sf::st_drop_geometry() |> 
#   as.data.frame() |> 
#   dplyr::select(-replicate) |> 
#   dplyr::filter(PA == 0 & id_split == "252215_9" & thetao < 25)
# 
# id252215_9
# 
# 
# tracks_fixed <- tracks |> 
#   dplyr::bind_rows(id252215_9) |> 
#   dplyr::arrange(id, rep, date)
# 
# tracks_fixed |>
#   ggplot2::ggplot(ggplot2::aes(x = thetao, y = lat, colour = factor(PA))) +
#   ggplot2::geom_point(alpha = 0.3) +
#   ggplot2::scale_y_continuous()
# 
# 
# 
# tracks <- tracks_fixed
# 
# str(tracks)
# 
# tracks


unique(tracks$Tag_type)
unique(tracks$month)
unique(tracks$rep)

# experiment removing PAT tags 
# tracks <- tracks |> 
#   dplyr::filter(!Tag_type %in% c("PSAT", "PSAT_SPOT"))
# 


# # crop absences to extend of presences +/- 1 degree buffer 
# 
# pres <- tracks |> dplyr::filter(PA == 1)
# if (nrow(pres) == 0) stop("No presences found (PA == 1).")
# 
# # 1) Bounding box from presences (ignore NAs)
# pad <- 1  # degrees
# lon_rng <- range(pres$lon, na.rm = TRUE)
# lat_rng <- range(pres$lat, na.rm = TRUE)
# 
# lon_min <- lon_rng[1] - pad
# lon_max <- lon_rng[2] + pad
# lat_min <- lat_rng[1] - pad
# lat_max <- lat_rng[2] + pad
# 
# # 2) Filter absences by lon/lat ranges
# abs_crop <- tracks |>
#   dplyr::filter(PA == 0) |>
#   dplyr::filter(dplyr::between(lon, lon_min, lon_max),
#                 dplyr::between(lat, lat_min, lat_max))
# 
# # 3) Combine back
# locs_crop <- dplyr::bind_rows(pres, abs_crop)
# 
# # Optional quick check
# message("Original: ", nrow(tracks), 
#         " | Pres: ", nrow(pres),
#         " | Abs kept: ", nrow(abs_crop),
#         " | New total: ", nrow(locs_crop))
# 
# ggplot() +
#   geom_point(data = tracks, aes(lon, lat), color = "grey80", size = 1) +
#   geom_point(data = abs_crop, aes(lon, lat), color = "blue", size = 1) +
#   geom_point(data = pres, aes(lon, lat), color = "red", size = 1.5) +
#   coord_equal() +
#   labs(x = "Longitude", y = "Latitude",
#        title = "Presences (red) and Cropped Absences (blue)") +
#   theme_minimal()
# 
# 
# 
# tracks <- locs_crop
# 





# tracks <- tracks |>
#   dplyr::rename(depth = Depth,
#                 slope = Slope,
#                 roughness = Roughness,
#                 wz = Wz) |>
#   sf::st_drop_geometry()

tracks <- tracks |>
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |>
  sf::st_drop_geometry()

str(tracks)

tracks |> distinct(id) |>  nrow()
tracks



# vars <- c("depth", "slope", "roughness","dist200", "dist1000", "dist2000", "thetao", "uv", "mltost", "chl", "wz", "dist_seamount", "dist_knoll")

vars <- c("depth", "slope", "roughness", "dist2000", "thetao", "uv", "mltost", "chl", "wz", "sst_slope")

rows_with_na <- tracks |>
  dplyr::filter(
    dplyr::if_any(dplyr::all_of(vars), ~ is.na(.))
  )

rows_with_na

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



# removing absence locations to achieve rougly 1:10 ratio - as it was found too many absences will make model bad 

# reduce number of replicayes 
N_reps <- 15  # try 10, 15, 20 ...

# 1) pick WHICH absence reps to keep per id_split
set.seed(99)

abs_rep_keep <- tracks |>
  dplyr::filter(PA == 0, !is.na(rep)) |>
  dplyr::distinct(id_split, rep) |>
  dplyr::group_by(id_split) |>
  dplyr::arrange(rep, .by_group = TRUE) |>
  dplyr::slice_sample(n = N_reps) |>
  dplyr::ungroup()

print(abs_rep_keep, n=200)

tracks_capped <- dplyr::bind_rows(
  tracks |> dplyr::filter(PA == 1),
  tracks |> dplyr::filter(PA == 0) |>
    dplyr::semi_join(abs_rep_keep, by = c("id_split","rep"))
)

# ?slice_sample
# 
# ## Method 2 to make sure number is the same 
# 
# set.seed(42)
# 
# pres_counts <- tracks |>
#   dplyr::filter(PA == 1) |>
#   sf::st_drop_geometry() |>
#   dplyr::count(id_split, name = "n_pres")
# 
# # 2) All candidate absences, with n_pres joined on
# abs_candidates <- tracks |>
#   dplyr::filter(PA == 0) |>
#   dplyr::inner_join(pres_counts, by = "id_split")
# 
# str(abs_candidates)
# 
# # 3) For each id_split, decide how many absences to keep:
# #    n_keep = min(#presences, #available absences)
# abs_n_keep <- abs_candidates |>
#   sf::st_drop_geometry() |>
#   dplyr::group_by(id_split) |>
#   dplyr::summarise(
#     n_keep = min(dplyr::first(n_pres), dplyr::n()),
#     .groups = "drop"
#   )
# 
# # 4) join n_keep back and randomly sample that many absences per id_split
# abs_kept <- abs_candidates |>
#   dplyr::inner_join(abs_n_keep, by = "id_split") |>
#   dplyr::group_by(id_split) |>
#   dplyr::slice_sample(n = dplyr::first(n_keep)) |>
#   dplyr::ungroup()
# 
# # 5) combine presences + balanced absences
# tracks_balanced <- dplyr::bind_rows(
#   tracks |> dplyr::filter(PA == 1),
#   abs_kept
# )
# 
# # quick check: should be ~1:1 per id_split
# dplyr::count(tracks_balanced, id_split, PA)
# sanity check

tracks |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(Pre = sum(PA==1), Abs = sum(PA==0), ratio = Abs/Pre, .groups="drop") |>
  print(n = 200)




tracks_capped |>
  dplyr::group_by(id_split) |>
  dplyr::summarise(Pre = sum(PA==1), Abs = sum(PA==0), ratio = Abs/Pre, .groups="drop") |>
  print(n = 200)

tracks_capped |> dplyr::group_by(PA) |>  dplyr::summarise(n = n())
tracks |> dplyr::group_by(PA) |>  dplyr::summarise(n = n())




tracks |> 
  as.data.frame() |>
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(thetao, type = "common")

tracks_capped |> 
  as.data.frame() |>
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(thetao, type = "common")


tracks_capped |> dplyr::filter(thetao < 25)


# Explore Data ------------------------------------------------------------

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
    # plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)

plots <- tracks_capped |> dplyr::mutate(PA = as.factor(PA))


(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


# Transformations
max(tracks$depth)
tracks |> dplyr::filter(depth >= 0)



tracks <- tracks_capped




tracks_trans <- tracks |>
    dplyr::mutate(
                  year  = factor(year),
                  # depth = case_when(depth >= 0 ~ -5, TRUE ~ depth),
                  # id = as.factor(id),
                  # id_split = as.factor(id_split),
                  dist200 = dist200/1000,
                  dist1000 = dist1000/1000,
                  dist2000 = dist2000/1000,
                  chl = log(chl)
    )

str(tracks_trans)
mapview::mapView(tracks_trans |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>  dplyr::select(-Date))

plots <- tracks_trans |> dplyr::mutate(PA = as.factor(PA))


(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 200) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


mapview::mapview(tracks_trans |> dplyr::select(-Date))

tracks_out <- tracks_trans |> 
  dplyr::filter(
    # between(uv, mean(uv, na.rm = TRUE) - 4 * sd(uv, na.rm = TRUE), mean(uv, na.rm = TRUE) + 4 * sd(uv, na.rm = TRUE)),
    # between(depth, mean(depth, na.rm = TRUE) - 4 * sd(depth, na.rm = TRUE), mean(depth, na.rm = TRUE) + 4 * sd(depth, na.rm = TRUE)),
    # between(dist2000, mean(dist2000, na.rm = TRUE) - 4 * sd(dist2000, na.rm = TRUE), mean(dist2000, na.rm = TRUE) + 4 * sd(dist2000, na.rm = TRUE)),
    # between(chl, mean(chl, na.rm = TRUE) - 4 * sd(chl, na.rm = TRUE), mean(chl, na.rm = TRUE) + 4 * sd(chl, na.rm = TRUE)),
    # between(thetao, mean(thetao, na.rm = TRUE) - 4 * sd(thetao, na.rm = TRUE), mean(thetao, na.rm = TRUE) + 4 * sd(thetao, na.rm = TRUE)),
    # between(mltost, mean(mltost, na.rm = TRUE) - 4 * sd(mltost, na.rm = TRUE), mean(mltost, na.rm = TRUE) + 4 * sd(mltost, na.rm = TRUE)),
    between(roughness, mean(roughness, na.rm = TRUE) - 4 * sd(roughness, na.rm = TRUE), mean(roughness, na.rm = TRUE) + 4 * sd(roughness, na.rm = TRUE)),
    between(slope, mean(slope, na.rm = TRUE) - 4 * sd(slope, na.rm = TRUE), mean(slope, na.rm = TRUE) + 4 * sd(slope, na.rm = TRUE)),
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

tracks |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())
tracks_trans |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())
tracks_out |> dplyr::group_by(PA) |> 
  dplyr::summarise(n = n())


tracks_out |> 
  as.data.frame() |>
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(thetao, type = "common")

mapview::mapview(tracks_out |> dplyr::filter(PA == 1) |> dplyr::select(-Date))
mapview::mapview(tracks_out |> dplyr::select(-Date))

str(tracks)

deleted_presences <- dplyr::anti_join(as.data.frame(tracks_trans), as.data.frame(tracks_out), by = names(tracks_trans)) |>
  dplyr::filter(PA == 1)


mapview::mapView(deleted_presences |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>  dplyr::select(-Date))

str(deleted_presences)

# autocrrelation ----------------------------------------------------------
tracks_occ <- tracks_out |> 
  dplyr::filter(PA == 1) |> 
  as.data.frame()


tracks_occ |> dplyr::filter(month ==10)

# tracks_occ <- tracks_trans |>
#   dplyr::filter(PA == 1) |>
#   as.data.frame()


variables <- c("depth", "slope", "roughness", "dist2000", "thetao", "uv", "mltost", "chl", "wz")

ac <- dynamicSDM::spatiotemp_autocorr(tracks_occ,
                                      varname = variables,
                                      plot = TRUE,
                                      temporal.level = c("day", "month", "year")) 


# ac$Plots

ac$Statistical_tests

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

temporal_tests |> dplyr::filter(p.value <= 0.05)

spatial_tests |> dplyr::filter(p.value <= 0.05)

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


tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")

tracks <- tracks |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# crop raster to extent of occurrences 

# tracks_m <- tracks_out |>
#   sf::st_make_valid() |>
#   sf::st_zm(drop = TRUE, what = "ZM") |>
#   sf::st_transform("EPSG:3577")


tracks_m <- tracks |>
  sf::st_make_valid() |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_transform("EPSG:3577")

# tracks_m <- tracks_trans |>
#   sf::st_make_valid() |>
#   sf::st_zm(drop = TRUE, what = "ZM") |>
#   sf::st_transform("EPSG:3577")

str(tracks_m)


v <- terra::vect(tracks_m)
v
bb <- terra::ext(v)
buff <- 100000  # ~2° buffer added

bb <- sf::st_bbox(tracks)

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



# tracks_occ <- tracks_out |>
#   as.data.frame()

tracks_occ <- tracks |>
  as.data.frame()

# tracks_occ <- tracks_trans |>
#   as.data.frame()



tracks_occ_blocks <- dynamicSDM::spatiotemp_block(tracks_occ,
                                                 vars.to.block.by = vars,
                                                 spatial.layer = raster_crop,
                                                 spatial.split.degrees = 5,
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


# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_dynSDM_10_2018_2025_extract_processed.rds")
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_3days_dynSDM_10_2018_2025_extract_processed.rds")
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract_processed.rds")
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")

# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced_woPATS.rds")
# 
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced.rds")
# 
# 
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_2.rds")
# 
# 
# # rep rediced to 20, maintainign 1:15 PA ratio approximately 
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_3.rds")
# 
# rep 15; spatial ac distance incrreased to 500 km; no lat cutoff
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")



# raddomly sampled bufeer absence locations (25-500km); rep 15; ; no lat cutoff
saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_5.rds")


# # raddomly sampled bufeer absence locations (25-500km); keep all absences 
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_6.rds")




# 
# 
# saveRDS(tracks_occ_blocks.sf, "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_removedPATS_test.rds")
# 


#_____________ use blockCV package 
library(blockCV)
library(tmap)

tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")


tracks_m <- tracks |>
  sf::st_make_valid() |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_transform("EPSG:3577")



tracks_m <- tracks_m |> dplyr::rename(id_obs = id)
row.names(tracks_m) <- as.character(seq_len(nrow(tracks_m)))

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
                   # dist200,
                   # dist1000,
                   dist2000, 
                   thetao, 
                   uv, 
                   mltost, 
                   chl, 
                   wz,
                   sst_slope
                   # dist_seamount,
                   # dist_knoll
                   )




#load monthly mean raster stack of continous variables 
input_raster_stack <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log.tif")
input_raster_stack <- input_raster_stack[[!names(input_raster_stack) %in% c("lat", 
                                                                            "lon", 
                                                                            "month", 
                                                                            "id", 
                                                                            "dist200",
                                                                            "dist1000", 
                                                                            "dist_seamount", 
                                                                            "dist_knoll"
                                                                            # "sst_slope",
                                                                            # "mltost"
                                                                            )]]

plot(input_raster_stack)

v <- terra::vect(tracks_m)
v
bb <- terra::ext(v)
buff <- 10000 

bb <- sf::st_bbox(tracks)
# bb <- sf::st_bbox(tracks_trans)

bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 0.5),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 0.5),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 0.5),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 0.5)
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
names(input_raster)
str(input_occ)


tmap::tm_shape(input_raster[[names(input_raster) != "id"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = gray.colors(10)),
    col.legend = tm_legend_hide()
  ) +
  tmap::tm_shape(input_occ) +
  tmap::tm_dots(
    fill = "PA",
    fill.scale = tm_scale_categorical(),
    size = 0.5,
    fill_alpha = 0.5
  )



# test spatial autocreelation distance
set.seed(44)
sac1 <- blockCV::cv_spatial_autocor(r = input_raster, 
                                    num_sample = 10000,
                                    plot = TRUE)
sac1$range

summary(sac1)


sac2 <- blockCV::cv_spatial_autocor(r = input_raster,
                                    x = input_occ, 
                                    column =  "PA", 
                                    num_sample = 10000,
                                    plot = TRUE)
sac2$range

str(sac2)
summary(sac2)

sac3 <- blockCV::cv_spatial_autocor(x = input_occ,
                                    column = "PA",
                                    plot = TRUE)
sac3$range

library(automap)
sac2$variograms
plot(sac1$variograms[[1]])
plot(sac2$variograms[[1]])
plot(sac3$variograms[[1]]) 




graphics::plot(
  exp_var$dist,
  exp_var$gamma,
  xlim = c(0, 5e5),  # show only first 500 km
  xlab = "Distance (m)",
  ylab = "Semivariance"
)

exp_var |>
  dplyr::mutate(dist_km = dist / 1000) |>
  ggplot2::ggplot(ggplot2::aes(x = dist_km, y = gamma)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::coord_cartesian(xlim = c(0, 1000)) +  # e.g. 0–500 km window
  ggplot2::labs(
    x = "Distance (km)",
    y = "Semivariance"
  ) +
  ggplot2::theme_bw()





blockCV::cv_block_size(r = input_raster,
                       x = input_occ,
                       column = "PA",
                       min_size = 10000,
                       max_size = 1000000)








## Spatial Blocks
# random 
# range = 100000
# range = 350000 # size of the blocks in metres
range = 500000
pal <- cmocean::cmocean("balance")(100)

sb1 <- blockCV::cv_spatial(x = input_occ,
                           column = "PA", # the response column (binary or multi-class)
                           r = input_raster,
                           k = 5, # number of folds
                           size = range, 
                           hexagon = TRUE,
                           selection = "random", # random blocks-to-fold
                           iteration = 100, # find evenly dispersed folds
                           progress = TRUE,
                           seed = 666,
                           biomod2 = FALSE) # also create folds for biomod2
sb1$records

blockCV::cv_plot(cv = sb1, 
                 x = input_occ,
                 r = input_raster,
                 num_plots = 1:5,
                 nrow = 2, 
                 points_alpha = 0.5)


block <- blockCV::cv_plot(cv = sb1, 
                 x = input_occ,
                 r = input_raster,
                 num_plots = 5,
                 nrow = 2, 
                 points_alpha = 0.5) +
  theme(text = element_text(size =10))

block


sb1_sim <- blockCV::cv_similarity(cv = sb1,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:5,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  # points_colors = pal,
                                  progress = TRUE)


sb1_sim <- sb1_sim +
  ggplot2::scale_color_gradientn(colours = pal)
sb1_sim





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
                           hexagon = TRUE,
                           size = range, # size of the blocks in metres
                           selection = "systematic", 
                           iteration = 100, # find evenly dispersed folds
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
                                  num_plot = 1:5,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

sb3_sim


sb1_sim
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



# Nearest Neighbour Distance Matching (NNDM) LOO

nncv <- blockCV::cv_nndm(x = input_occ,
                column = "PA",
                r = input_raster,
                presence_bg = FALSE,   # we treat pseudo-absences as true absences as we carefully sampled them and can confdently say the animal has not been there 
                size = 500000,
                num_sample = 1000,
                sampling = "random",
                min_train = 0.25,
                plot = TRUE)


nncv$records

nncv$plot


blockCV::cv_plot(cv = nncv,
                 x = input_occ,
                 r = input_raster,
                 nrow = 2,
                 points_alpha = 0.5)




nncv_sim <- blockCV::cv_similarity(cv = nncv,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:100,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  points_colors = pal,
                                  progress = TRUE)

nncv_sim


# Save  -------------------------------------------------------------------

# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking.rds")
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_3days.rds")
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_3to7days.rds")
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_1to4days.rds")
# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_woPATS.rds")
# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final.rds")
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_test.rds")
# 
# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_2.rds")
# 
# saveRDS(sb3, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_3.rds")
# 
# 
saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_4.rds")
saveRDS(nncv, "data/processed/BlockCV_spatial_NNCV_folds_tracking_daily_mp_monthreduced_final_4.rds")

# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_5.rds")

# saveRDS(sb1, "data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_6.rds")
