

### Bluelink upward current BRAN2020



# Load required libraries
library(terra)
library(stringr)
library(lubridate)
library(ncdf4)







## load monthly composites downloaded from BRAN
nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Environmental_Varibales_Downloaded/Bluelink/Monthly"
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
nc_files

nc_bran <- terra::rast(nc_files)
print(nc_bran)
names(nc_bran)

# Extract the indices of the first layer of each month
first_layer_indices <- seq(1, nlyr(nc_bran), by = 51)

# Subset the SpatRaster to keep only the first layer of each month
Wz_monthly_comps_0.1 <- subset(nc_bran, first_layer_indices)
Wz_monthly_comps_0.1

n_layers <- terra::nlyr(Wz_monthly_comps_0.1)  # 228
n_layers
layer_dates <- seq.Date(
  from = as.Date("2005-01-16"),  # adjust start if needed
  by   = "1 month",
  length.out = n_layers
)

names(Wz_monthly_comps_0.1) <- layer_dates |>
  format("%Y-%m-%d")   # e.g. "2010_01", "2010_02", ...


names(Wz_monthly_comps_0.1)
print(Wz_monthly_comps_0.1)

plot(Wz_monthly_comps_0.1[[165]], range  = c(-0.00005, 0.00005))



# Define the extent of the region of interest
temp_raster <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_BRT_mean_rev_mp_crwPA.tif")


roi_extent <- terra::ext(temp_raster)
roi_extent


# Crop raster
Wz_monthly_comps_0.1_crop <- crop(Wz_monthly_comps_0.1, roi_extent)

bran_stack <- Wz_monthly_comps_0.1_crop


## Copernicus Marine 


reference_raster_10km <- bran_stack

output_path <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files"
nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/CMEMS_Global/Monthly/wo_after2024/"

nc_folder <- "/Volumes/Ingo_PhD_2/CMEMS_CoralSea/monthly/w/2025update"
nc_folder
depth_layer <- 5  # User-defined depth layer (e.g., 5m)


# create_monthly_rasters_from_files_CM(nc_folder, depth_layer, reference_raster = reference_raster_10km, output_path = output_path, output_filename = "Copernicus_W_5m_Monthly_Stack_2022to2025_10km")


create_monthly_rasters_from_files_CM(nc_folder, desired_depth = depth_layer,
                                     reference_raster = reference_raster_10km,
                                     output_path = output_path, output_filename = "COP_stack")


COP_stack

d <- as.Date(terra::time(COP_stack))
# Set each to the 16th of its month
t16 <- as.Date(strftime(d, format = "%Y-%m-16"))

# Update the stack's time and names
terra::time(COP_stack) <- as.POSIXct(t16, tz = "UTC")
names(COP_stack)       <- format(t16, "%Y-%m-%d")



# assure same ref grid
crs(COP_stack) <- "EPSG:4326"
crs(bran_stack) <- "EPSG:4326"

COP_stack
bran_stack


# --- Helper: parse dates from layer names -------------------------------------
get_dates <- function(x) {
  base::as.Date(base::trimws(terra::names(x)))
}

dates_bran <- get_dates(bran_stack)
dates_cop  <- get_dates(COP_stack)

# Keep only overlapping months (intersection of dates)
common_dates <- base::sort(base::intersect(dates_bran, dates_cop))
common_dates
common_dates <- as.Date(common_dates, origin = "1970-01-01")

idx_bran <- match(common_dates, dates_bran)
idx_cop  <- match(common_dates, dates_cop)

stopifnot(!anyNA(idx_bran), !anyNA(idx_cop))  # sanity check

bran_olap <- bran_stack[[idx_bran]]
cop_olap  <- COP_stack [[idx_cop ]]




# --- Align Copernicus to BRAN grid (CRS/res/extent) ---------------------------
# Reproject/resample ONLY if needed
if (!terra::compareGeom(bran_olap, cop_olap, stopOnError = FALSE)) {
  # Crop to overlapping extent first (cheap), then resample
  e_int      <- terra::intersect(terra::ext(bran_olap), terra::ext(cop_olap))
  bran_olap  <- terra::crop(bran_olap, e_int)
  cop_olap   <- terra::crop(cop_olap,  e_int)
  cop_olap   <- terra::resample(cop_olap, bran_olap, method = "bilinear")
}

# # Keep cells non-NA in BOTH stacks for ALL months
# keep <- !is.na(terra::app(bran_olap, fun = sum, na.rm = TRUE)) &
#   !is.na(terra::app(cop_olap,  fun = sum, na.rm = TRUE))
# 
# bran_olap <- terra::mask(bran_olap, keep)
# cop_olap  <- terra::mask(cop_olap,  keep)



# --- Area-weighting by latitude (for lon/lat EPSG:4326) -----------------------
lat_r   <- terra::init(bran_olap, fun = "y")                # per-cell latitude (degrees)
w_r     <- base::cos(lat_r * base::pi / 180)                # weights ∝ cell area
w_sum   <- terra::global(w_r, fun = "sum", na.rm = TRUE)[1,1]
w_mean  <- function(r) terra::global(r * w_r, "sum", na.rm = TRUE)[1,1] / w_sum
w_rmse  <- function(r) base::sqrt(w_mean(r^2))

# Weighted correlation helper
wcor <- function(x, y, w) {
  wx <- base::sum(w * x) / base::sum(w)
  wy <- base::sum(w * y) / base::sum(w)
  covw <- base::sum(w * (x - wx) * (y - wy)) / base::sum(w)
  sx <- base::sqrt(base::sum(w * (x - wx)^2) / base::sum(w))
  sy <- base::sqrt(base::sum(w * (y - wy)^2) / base::sum(w))
  covw / (sx * sy)
}

# --- Compute monthly metrics ---------------------------------------------------
# We’ll also sample points for correlation/regression (fast & memory-safe)
sample_size <- 80000L  # adjust to your machine; 20k–200k typical

monthly_stats <- base::lapply(seq_along(common_dates), function(i) {
  bran <- bran_olap[[i]]
  cop  <- cop_olap [[i]]
  d    <- cop - bran
  
  ## Unweighted stats
  bias_mean <- terra::global(d, fun = "mean", na.rm = TRUE)[1,1]
  bias_med  <- terra::global(d, fun = stats::median, na.rm = TRUE)[1,1]
  rmse      <- sqrt(terra::global(d^2, fun = "mean", na.rm = TRUE)[1,1])
  
  ## Area-weighted stats
  bias_wmean <- w_mean(d)
  rmse_w     <- w_rmse(d)
  
  ## MAE & robust stats (this month only)
  mae       <- terra::global(abs(d), fun = "mean",   na.rm = TRUE)[1,1]
  mae_w     <- terra::global(abs(d) * w_r, "sum",    na.rm = TRUE)[1,1] / w_sum
  medae     <- terra::global(abs(d), fun = stats::median, na.rm = TRUE)[1,1]
  med_e     <- terra::global(d,      fun = stats::median, na.rm = TRUE)[1,1]
  mad_sigma <- 1.4826 * terra::global(abs(d - med_e), fun = stats::median, na.rm = TRUE)[1,1]
  
  ## Correlation / regression using samples (avoid masking base::c)
  s  <- terra::spatSample(base::c(bran, cop, lat_r), size = sample_size,
                          method = "random", na.rm = TRUE, xy = FALSE)
  df <- base::as.data.frame(s)
  base::names(df) <- c("bran", "cop", "lat")
  df$w <- base::cos(df$lat * base::pi / 180)
  
  corr   <- stats::cor(df$bran, df$cop)
  corr_w <- {
    wx <- sum(df$w * df$bran) / sum(df$w)
    wy <- sum(df$w * df$cop)  / sum(df$w)
    covw <- sum(df$w * (df$bran - wx) * (df$cop - wy)) / sum(df$w)
    sx <- sqrt(sum(df$w * (df$bran - wx)^2) / sum(df$w))
    sy <- sqrt(sum(df$w * (df$cop  - wy)^2) / sum(df$w))
    covw / (sx * sy)
  }
  
  fit   <- stats::lm(cop ~ bran, data = df)
  a     <- stats::coef(fit)[1]
  bcoef <- stats::coef(fit)[2]
  
  fit_w <- stats::lm(cop ~ bran, data = df, weights = df$w)
  aw    <- stats::coef(fit_w)[1]
  bw    <- stats::coef(fit_w)[2]
  
  data.frame(
    date        = common_dates[i],
    year        = lubridate::year(common_dates[i]),
    month       = lubridate::month(common_dates[i]),
    bias_mean   = bias_mean,
    bias_median = bias_med,
    rmse        = rmse,
    bias_wmean  = bias_wmean,
    rmse_w      = rmse_w,
    mae         = mae,
    mae_w       = mae_w,
    medae       = medae,
    mad_sigma   = mad_sigma,
    corr        = corr,
    corr_w      = corr_w,
    slope       = bcoef,
    intercept   = a,
    slope_w     = bw,
    intercept_w = aw,
    n_sample    = nrow(df),
    stringsAsFactors = FALSE
  )
})


monthly_stats_w <- dplyr::bind_rows(monthly_stats)

# --- Optional: per-month quick plots ------------------------------------------
# 1) Time series of area-weighted mean bias (Cop - BRAN)
monthly_stats_w |>
  dplyr::mutate(ym = lubridate::ymd(sprintf("%04d-%02d-15", year, month))) |>
  (\(df) {
    graphics::plot(df$ym, df$bias_wmean, type = "b",
                   xlab = "Month", ylab = "Area-weighted mean bias (m/s)",
                   main = "Copernicus − BRAN (monthly, area-weighted)")
    df
  })()

# 2) Time series of area-weighted RMSE
monthly_stats_w |>
  dplyr::mutate(ym = lubridate::ymd(sprintf("%04d-%02d-15", year, month))) |>
  (\(df) {
    graphics::plot(df$ym, df$rmse_w, type = "b",
                   xlab = "Month", ylab = "Area-weighted RMSE (m/s)",
                   main = "Monthly RMSE (Cop vs BRAN, area-weighted)")
    df
  })()

# --- Optional: save results ----------------------------------------------------
# readr::write_csv(monthly_stats, "bran_vs_cop_monthly_stats_2022_2024.csv")
monthly_stats_w


range_bias <- range(monthly_stats_w$bias_wmean, na.rm = TRUE)
range_rmse <- range(monthly_stats_w$rmse_w,     na.rm = TRUE)
range_mae  <- range(monthly_stats_w$mae_w,      na.rm = TRUE)
sprintf("bias_wmean: % .2e to % .2e;  RMSE_w: %.2e–%.2e;  MAE_w: %.2e–%.2e",
        range_bias[1], range_bias[2], range_rmse[1], range_rmse[2], range_mae[1], range_mae[2])


op <- par(mfrow = c(2,1), mar = c(4,5,2,2))

# 1a) area-weighted mean bias
plot(monthly_stats_w$date, monthly_stats_w$bias_wmean, type = "b", pch = 16,
     xlab = "", ylab = expression("Bias"[w]*" (m s"^{-1}*")"),
     main = "Copernicus − BRAN: Area-weighted mean bias")
abline(h = 0, lty = 3)

# 1b) area-weighted errors
yr <- range(c(monthly_stats_w$rmse_w, monthly_stats_w$mae_w), na.rm = TRUE)

plot(monthly_stats_w$date, monthly_stats_w$rmse_w,
     type = "b", pch = 16, lwd = 2,
     ylim = yr,
     xlab = "Month",
     ylab = expression("Error (m s"^{-1}*")"),
     main = "Monthly errors")

lines(monthly_stats_w$date, monthly_stats_w$mae_w,
      type = "b", pch = 1, lwd = 2, col = 2)   # different style/color

legend("topleft", bty = "n",
       pch = c(16,1), lty = 1, lwd = 2, col = c(1,2),
       legend = c("RMSE","MAE"))

par(op)



err  <- cop_olap - bran_olap
absE <- abs(err)

## (1) Zero-mean offset across months (bias_wmean CI)
b <- monthly_stats_w$bias_wmean
ci <- mean(b) + c(-1,1) * qt(0.975, df=length(b)-1) * sd(b)/sqrt(length(b))
list(mean_bias_w = mean(b), CI95 = ci)

## (2) Errors relative to typical magnitude of the field
lat_r <- terra::init(bran_olap, fun = "y"); w_r <- cos(lat_r*pi/180)
w_sum <- terra::global(w_r, "sum", na.rm=TRUE)[1,1]

aw_mean_abs_bran <- sapply(1:terra::nlyr(bran_olap), function(i)
  terra::global(abs(bran_olap[[i]])*w_r, "sum", na.rm=TRUE)[1,1] / w_sum)

rel_MAE <- monthly_stats_w$mae_w / aw_mean_abs_bran     # unitless
rel_RMSE <- monthly_stats_w$rmse_w / aw_mean_abs_bran   # unitless
summary(rel_MAE); summary(rel_RMSE)

## (3) Errors relative to natural variability (z-RMSE and exceedance)
bran_sd <- terra::app(bran_olap, fun = stats::sd, na.rm = TRUE)  # per-pixel SD across months
sd_floor <- 1e-6
scaleA <- terra::ifel(bran_sd < sd_floor, sd_floor, bran_sd)

z2 <- (err/scaleA)^2
zRMSE_w <- sqrt( terra::global(z2 * w_r, "sum", na.rm=TRUE)[1,1] / w_sum )

frac_gt2 <- terra::global( abs(err/scaleA) > 2, fun = "mean", na.rm = TRUE)[1,1]  # per-layer fraction
summary(as.numeric(frac_gt2))
zRMSE_w







# combined plot -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)


# Build tidy data for both variables
ms_w <- monthly_stats_w %>%
  transmute(
    date = as.Date(date),
    variable = "W",
    Bias = bias_wmean, RMSE = rmse_w, MAE = mae_w
  )

ms_w



df <- ms_w |> 
  pivot_longer(c(Bias, RMSE, MAE), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Bias","RMSE","MAE")))

# Zero line only for Bias panels
bias_hline <- df |>
  filter(metric == "Bias") |>
  distinct(variable) |>
  mutate(yintercept = 0)

str(df)




scale_w <- 1  # display W in µm s^-1 (set to 1 to keep m s^-1)

w_df <- df |>
  dplyr::filter(variable == "W") |>
  dplyr::mutate(
    value = value * scale_w,                 # rescale for display
    panel = if_else(metric == "Bias", "Bias", "Errors")
  )

P.w <- w_df |>
  ggplot(aes(date, value)) +
  # dashed zero only on Bias panel
  geom_hline(
    data = w_df |> filter(panel == "Bias"),
    aes(yintercept = 0),
    linetype = 3, colour = "grey50", linewidth = 0.4, inherit.aes = FALSE
  ) +
  # Bias line
  geom_line(
    data = w_df |> filter(panel == "Bias"),
    linewidth = 0.8, colour = "black"
  ) +
  geom_point(
    data = w_df |> filter(panel == "Bias"),
    size = 1.6, colour = "black"
  ) +
  # Errors panel: both RMSE & MAE together
  geom_line(
    data = w_df |> filter(panel == "Errors"),
    aes(colour = metric, linetype = metric),
    linewidth = 0.9
  ) +
  geom_point(
    data = w_df |> filter(panel == "Errors"),
    aes(colour = metric, shape = metric),
    size = 1.8
  ) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(RMSE = "steelblue", MAE = "firebrick")) +
  scale_shape_manual(values  = c(RMSE = 16,     MAE = 1)) +
  scale_linetype_manual(values = c(RMSE = "solid", MAE = "solid")) +
  labs(x = "Month",
       # y = if (scale_w == 1) NULL else "Value (µm/s)",
       y = "Vertical Velocity (W) [m/s]",
       title = "",
       subtitle = "") +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position    = c(0.04, 0.47),    # inside bottom facet
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )


P.w



ggsave("BRAN_Copernicus_sitch_justification_plot.png", plot = P.w, path ="outputs/final_figures", scale =1, width = 18, height = 18, units = "cm", dpi = 300)



















# 
# 
# 
# # MLD ---------------------------------------------------------------------
# 
# 
# 
# ## load monthly composites downloaded from BRAN
# nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Environmental_Varibales_Downloaded/Bluelink/MLD/Monthly/"
# nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
# nc_files
# 
# nc_bran <- terra::rast(nc_files)
# print(nc_bran)
# names(nc_bran)
# 
# mld_monthly_comps_0.1 <- nc_bran
# 
# names_monthly_composites <- readRDS("names_monthly_comps.rds")
# names_monthly_composites_2023 <- names_monthly_composites[1:168]
# 
# # change layer names to dates 
# names(mld_monthly_comps_0.1) <- names_monthly_composites_2023
# 
# names(mld_monthly_comps_0.1)
# print(mld_monthly_comps_0.1)
# 
# plot(mld_monthly_comps_0.1[[165]])
# 
# 
# 
# # Define the extent of the region of interest
# SST_mean <- rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files/SST_means_2010_2023_MUR_erddap.tif")
# 
# roi_extent <- terra::ext(SST_mean)
# roi_extent
# 
# 
# 
# # Crop raster
# mld_monthly_comps_0.1_crop <- crop(mld_monthly_comps_0.1, roi_extent)
# 
# bran_stack <- mld_monthly_comps_0.1_crop
# 
# 
# ## Copernicus Marine 
# 
# 
# reference_raster_10km <- bran_stack
# 
# output_path <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Processed_Raster_Files"
# nc_folder <- "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Copernicus/MonthlyDownloads/MLD/Monthly/2025update/"
# nc_folder
# 
# nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
# 
# 
# nc_cop <- terra::rast(nc_files)
# print(nc_cop)
# names(nc_cop)
# 
# 
# 
# COP_stack <- nc_cop
# 
# d <- as.Date(terra::time(COP_stack))
# # Set each to the 16th of its month
# t16 <- as.Date(strftime(d, format = "%Y-%m-16"))
# 
# # Update the stack's time and names
# terra::time(COP_stack) <- as.POSIXct(t16, tz = "UTC")
# names(COP_stack)       <- format(t16, "%Y-%m-%d")
# 
# COP_stack
# 
# # assure same ref grid
# crs(COP_stack) <- "EPSG:4326"
# crs(bran_stack) <- "EPSG:4326"
# 
# COP_stack
# bran_stack
# 






# # --- Helper: parse dates from layer names -------------------------------------
# get_dates <- function(x) {
#   base::as.Date(base::trimws(terra::names(x)))
# }
# 
# dates_bran <- get_dates(bran_stack)
# dates_cop  <- get_dates(COP_stack)
# 
# # Keep only overlapping months (intersection of dates)
# common_dates <- base::sort(base::intersect(dates_bran, dates_cop))
# common_dates
# common_dates <- as.Date(common_dates, origin = "1970-01-01")
# 
# idx_bran <- match(common_dates, dates_bran)
# idx_cop  <- match(common_dates, dates_cop)
# 
# stopifnot(!anyNA(idx_bran), !anyNA(idx_cop))  # sanity check
# 
# bran_olap <- bran_stack[[idx_bran]]
# cop_olap  <- COP_stack [[idx_cop ]]
# 
# # --- Align Copernicus to BRAN grid (CRS/res/extent) ---------------------------
# 
# if (!terra::compareGeom(bran_olap, cop_olap, stopOnError = FALSE)) {
#   e_int      <- terra::intersect(terra::ext(bran_olap), terra::ext(cop_olap))
#   bran_olap  <- terra::crop(bran_olap, e_int)
#   cop_olap   <- terra::crop(cop_olap,  e_int)
#   cop_olap   <- terra::resample(cop_olap, bran_olap, method = "bilinear")
# }
# 
# 
# 
# # --- Area-weighting by latitude (for lon/lat EPSG:4326) -----------------------
# lat_r   <- terra::init(bran_olap, fun = "y")                # per-cell latitude (degrees)
# w_r     <- base::cos(lat_r * base::pi / 180)                # weights ∝ cell area
# w_sum   <- terra::global(w_r, fun = "sum", na.rm = TRUE)[1,1]
# w_mean  <- function(r) terra::global(r * w_r, "sum", na.rm = TRUE)[1,1] / w_sum
# w_rmse  <- function(r) base::sqrt(w_mean(r^2))
# 
# # Weighted correlation helper
# wcor <- function(x, y, w) {
#   wx <- base::sum(w * x) / base::sum(w)
#   wy <- base::sum(w * y) / base::sum(w)
#   covw <- base::sum(w * (x - wx) * (y - wy)) / base::sum(w)
#   sx <- base::sqrt(base::sum(w * (x - wx)^2) / base::sum(w))
#   sy <- base::sqrt(base::sum(w * (y - wy)^2) / base::sum(w))
#   covw / (sx * sy)
# }
# 
# # --- Compute monthly metrics ---------------------------------------------------
# # We’ll also sample points for correlation/regression (fast & memory-safe)
# sample_size <- 80000L  # adjust to your machine; 20k–200k typical
# 
# monthly_stats <- base::lapply(seq_along(common_dates), function(i) {
#   bran <- bran_olap[[i]]
#   cop  <- cop_olap [[i]]
#   d    <- cop - bran
#   
#   ## Unweighted stats
#   bias_mean <- terra::global(d, fun = "mean", na.rm = TRUE)[1,1]
#   bias_med  <- terra::global(d, fun = stats::median, na.rm = TRUE)[1,1]
#   rmse      <- sqrt(terra::global(d^2, fun = "mean", na.rm = TRUE)[1,1])
#   
#   ## Area-weighted stats
#   bias_wmean <- w_mean(d)
#   rmse_w     <- w_rmse(d)
#   
#   ## MAE & robust stats (this month only)
#   mae       <- terra::global(abs(d), fun = "mean",   na.rm = TRUE)[1,1]
#   mae_w     <- terra::global(abs(d) * w_r, "sum",    na.rm = TRUE)[1,1] / w_sum
#   medae     <- terra::global(abs(d), fun = stats::median, na.rm = TRUE)[1,1]
#   med_e     <- terra::global(d,      fun = stats::median, na.rm = TRUE)[1,1]
#   mad_sigma <- 1.4826 * terra::global(abs(d - med_e), fun = stats::median, na.rm = TRUE)[1,1]
#   
#   ## Correlation / regression using samples (avoid masking base::c)
#   s  <- terra::spatSample(base::c(bran, cop, lat_r), size = sample_size,
#                           method = "random", na.rm = TRUE, xy = FALSE)
#   df <- base::as.data.frame(s)
#   base::names(df) <- c("bran", "cop", "lat")
#   df$w <- base::cos(df$lat * base::pi / 180)
#   
#   corr   <- stats::cor(df$bran, df$cop)
#   corr_w <- {
#     wx <- sum(df$w * df$bran) / sum(df$w)
#     wy <- sum(df$w * df$cop)  / sum(df$w)
#     covw <- sum(df$w * (df$bran - wx) * (df$cop - wy)) / sum(df$w)
#     sx <- sqrt(sum(df$w * (df$bran - wx)^2) / sum(df$w))
#     sy <- sqrt(sum(df$w * (df$cop  - wy)^2) / sum(df$w))
#     covw / (sx * sy)
#   }
#   
#   fit   <- stats::lm(cop ~ bran, data = df)
#   a     <- stats::coef(fit)[1]
#   bcoef <- stats::coef(fit)[2]
#   
#   fit_w <- stats::lm(cop ~ bran, data = df, weights = df$w)
#   aw    <- stats::coef(fit_w)[1]
#   bw    <- stats::coef(fit_w)[2]
#   
#   data.frame(
#     date        = common_dates[i],
#     year        = lubridate::year(common_dates[i]),
#     month       = lubridate::month(common_dates[i]),
#     bias_mean   = bias_mean,
#     bias_median = bias_med,
#     rmse        = rmse,
#     bias_wmean  = bias_wmean,
#     rmse_w      = rmse_w,
#     mae         = mae,
#     mae_w       = mae_w,
#     medae       = medae,
#     mad_sigma   = mad_sigma,
#     corr        = corr,
#     corr_w      = corr_w,
#     slope       = bcoef,
#     intercept   = a,
#     slope_w     = bw,
#     intercept_w = aw,
#     n_sample    = nrow(df),
#     stringsAsFactors = FALSE
#   )
# })
# 
# 
# monthly_stats_mld <- dplyr::bind_rows(monthly_stats)
# dev.off()
# # --- Optional: per-month quick plots ------------------------------------------
# # 1) Time series of area-weighted mean bias (Cop - BRAN)
# monthly_stats_mld |>
#   dplyr::mutate(ym = lubridate::ymd(sprintf("%04d-%02d-15", year, month))) |>
#   (\(df) {
#     graphics::plot(df$ym, df$bias_wmean, type = "b",
#                    xlab = "Month", ylab = "Area-weighted mean bias (m)",
#                    main = "Copernicus − BRAN (monthly, area-weighted)")
#     df
#   })()
# 
# # 2) Time series of area-weighted RMSE
# monthly_stats_mld |>
#   dplyr::mutate(ym = lubridate::ymd(sprintf("%04d-%02d-15", year, month))) |>
#   (\(df) {
#     graphics::plot(df$ym, df$rmse_w, type = "b",
#                    xlab = "Month", ylab = "Area-weighted RMSE (m)",
#                    main = "Monthly RMSE (Cop vs BRAN, area-weighted)")
#     df
#   })()
# 
# # --- Optional: save results ----------------------------------------------------
# # readr::write_csv(monthly_stats, "bran_vs_cop_monthly_stats_2022_2024.csv")
# monthly_stats_mld
# 
# 
# range_bias <- range(monthly_stats_mld$bias_wmean, na.rm = TRUE)
# range_rmse <- range(monthly_stats_mld$rmse_w,     na.rm = TRUE)
# range_mae  <- range(monthly_stats_mld$mae_w,      na.rm = TRUE)
# sprintf("bias_wmean: % .2e to % .2e;  RMSE_w: %.2e–%.2e;  MAE_w: %.2e–%.2e",
#         range_bias[1], range_bias[2], range_rmse[1], range_rmse[2], range_mae[1], range_mae[2])
# 
# 
# op <- par(mfrow = c(2,1), mar = c(4,5,2,2))
# 
# # 1a) area-weighted mean bias
# plot(monthly_stats_mld$date, monthly_stats_mld$bias_wmean, type = "b", pch = 16,
#      xlab = "", ylab = expression("Bias"[w]*" (m)"),
#      main = "Copernicus − BRAN: Area-weighted mean bias")
# abline(h = 0, lty = 3)
# 
# # 1b) area-weighted errors
# yr <- range(c(monthly_stats_mld$rmse_w, monthly_stats_mld$mae_w), na.rm = TRUE)
# 
# plot(monthly_stats_mld$date, monthly_stats_mld$rmse_w,
#      type = "b", pch = 16, lwd = 2,
#      ylim = yr,
#      xlab = "Month",
#      ylab = expression("Error (m)"),
#      main = "Monthly errors")
# 
# lines(monthly_stats_mld$date, monthly_stats_mld$mae_w,
#       type = "b", pch = 1, lwd = 2, col = 2)   # different style/color
# 
# legend("topleft", bty = "n",
#        pch = c(16,1), lty = 1, lwd = 2, col = c(1,2),
#        legend = c("RMSE","MAE"))
# 
# par(op)
# 
# 
# err  <- cop_olap - bran_olap
# absE <- abs(err)
# 
# ## (1) Zero-mean offset across months (bias_wmean CI)
# b <- monthly_stats_mld$bias_wmean
# ci <- mean(b) + c(-1,1) * qt(0.975, df=length(b)-1) * sd(b)/sqrt(length(b))
# list(mean_bias_w = mean(b), CI95 = ci)
# 
# ## (2) Errors relative to typical magnitude of the field
# lat_r <- terra::init(bran_olap, fun = "y"); w_r <- cos(lat_r*pi/180)
# w_sum <- terra::global(w_r, "sum", na.rm=TRUE)[1,1]
# 
# aw_mean_abs_bran <- sapply(1:terra::nlyr(bran_olap), function(i)
#   terra::global(abs(bran_olap[[i]])*w_r, "sum", na.rm=TRUE)[1,1] / w_sum)
# 
# rel_MAE <- monthly_stats_mld$mae_w / aw_mean_abs_bran     # unitless
# rel_RMSE <- monthly_stats_mld$rmse_w / aw_mean_abs_bran   # unitless
# summary(rel_MAE); summary(rel_RMSE)
# 
# ## (3) Errors relative to natural variability (z-RMSE and exceedance)
# bran_sd <- terra::app(bran_olap, fun = stats::sd, na.rm = TRUE)  # per-pixel SD across months
# sd_floor <- 1e-6
# scaleA <- terra::ifel(bran_sd < sd_floor, sd_floor, bran_sd)
# 
# z2 <- (err/scaleA)^2
# zRMSE_w <- sqrt( terra::global(z2 * w_r, "sum", na.rm=TRUE)[1,1] / w_sum )
# 
# frac_gt2 <- terra::global( abs(err/scaleA) > 2, fun = "mean", na.rm = TRUE)[1,1]  # per-layer fraction
# summary(as.numeric(frac_gt2))
# zRMSE_w
# 




# combined plot -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)


# Build tidy data for both variables
ms_w <- monthly_stats_w %>%
  transmute(
    date = as.Date(date),
    variable = "W",
    Bias = bias_wmean, RMSE = rmse_w, MAE = mae_w
  )

ms_mld <- monthly_stats_mld %>%
  transmute(
    date = as.Date(date),
    variable = "MLD",
    Bias = bias_wmean, RMSE = rmse_w, MAE = mae_w
  )

df <- bind_rows(ms_w, ms_mld) |>
  pivot_longer(c(Bias, RMSE, MAE), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Bias","RMSE","MAE")))

# Zero line only for Bias panels
bias_hline <- df |>
  filter(metric == "Bias") |>
  distinct(variable) |>
  mutate(yintercept = 0)

str(df)




scale_w <- 1  # display W in µm s^-1 (set to 1 to keep m s^-1)

w_df <- df |>
  dplyr::filter(variable == "W") |>
  dplyr::mutate(
    value = value * scale_w,                 # rescale for display
    panel = if_else(metric == "Bias", "Bias", "Errors")
  )

P.w <- w_df |>
  ggplot(aes(date, value)) +
  # dashed zero only on Bias panel
  geom_hline(
    data = w_df |> filter(panel == "Bias"),
    aes(yintercept = 0),
    linetype = 3, colour = "grey50", linewidth = 0.4, inherit.aes = FALSE
  ) +
  # Bias line
  geom_line(
    data = w_df |> filter(panel == "Bias"),
    linewidth = 0.8, colour = "black"
  ) +
  geom_point(
    data = w_df |> filter(panel == "Bias"),
    size = 1.6, colour = "black"
  ) +
  # Errors panel: both RMSE & MAE together
  geom_line(
    data = w_df |> filter(panel == "Errors"),
    aes(colour = metric, linetype = metric),
    linewidth = 0.9
  ) +
  geom_point(
    data = w_df |> filter(panel == "Errors"),
    aes(colour = metric, shape = metric),
    size = 1.8
  ) +
  facet_grid(panel ~ ., scales = "free_y") +
  scale_colour_manual(values = c(RMSE = "steelblue", MAE = "firebrick")) +
  scale_shape_manual(values  = c(RMSE = 16,     MAE = 1)) +
  scale_linetype_manual(values = c(RMSE = "solid", MAE = "solid")) +
  labs(x = "Month",
       # y = if (scale_w == 1) NULL else "Value (µm/s)",
       y = "Vertical Velocity (W) [m/s]",
       title = "",
       subtitle = "") +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position    = c(0.02, 0.47),    # inside bottom facet
    legend.justification = c("left", "top"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )




mld_df <- df |>
  dplyr::filter(variable == "MLD") |>
  dplyr::mutate(panel = if_else(metric == "Bias", "Bias", "Errors"),
                date = as.Date(date))

mld_df





P.mld <- mld_df |>
  ggplot(aes(date, value)) +
  
  geom_hline(
    data = mld_df |> filter(panel == "Bias"),
    aes(yintercept = 0),
    linetype = 3, colour = "grey50", linewidth = 0.4
  ) +
  geom_line(data = mld_df |> filter(panel == "Bias"),
            linewidth = 0.8, colour = "black") +
  geom_point(data = mld_df |> filter(panel == "Bias"),
             size = 1.6, colour = "black") +
  geom_line(data = mld_df |> filter(panel == "Errors"),
            aes(colour = metric, linetype = metric),
            linewidth = 0.9) +
  geom_point(data = mld_df |> filter(panel == "Errors"),
             aes(colour = metric, shape = metric),
             size = 1.8) +
  facet_grid(panel ~ ., scales = "free_y") +
  scale_colour_manual(values = c(RMSE = "steelblue", MAE = "firebrick")) +
  scale_shape_manual(values  = c(RMSE = 16,     MAE = 1)) +
  scale_linetype_manual(values = c(RMSE = "solid", MAE = "solid")) +
  labs(x = "Month",
       y = "Mixed Layer Depth (MLD) [m]") +
  theme_bw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position   = "none",   
    legend.justification = c("left", "top"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

P.mld

library(patchwork)
PLOT <- P.w + P.mld +
  plot_annotation(tag_levels = 'A', tag_suffix = ")") &
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 2))

PLOT

ggsave("BRAN_Copernicus_siwth_justification_plot.png", plot = PLOT, path ="/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/Output_Plots", scale =1, width = 18, height = 12, units = "cm", dpi = 300)

