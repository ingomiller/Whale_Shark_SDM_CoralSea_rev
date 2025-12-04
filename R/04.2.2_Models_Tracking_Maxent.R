#_____________________________________________________________________________
#                        Models: Tracking - Maxent
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




# Import DATA -------------------------------------------------------------


# tracks <- readRDS("data/processed/Tracks_PA_w_3days_dynSDM_10_2018_2025_extract_processed.rds")
# # tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced_final.rds")

# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced.rds")

# tracks |>
#   dplyr::group_by(id) |>
#   dplyr::summarise(Pre = sum(PA==1), Abs = sum(PA==0), ratio = Abs/Pre, .groups="drop") |>
#   print(n = 200)
  

tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")

# random backgroudn sampled pseudo absences; rep=15 (~1:10 ratio PA); thetao outlier removed too
tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_5.rds")

# tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_6.rds")
str(tracks)

tracks |> dplyr::select(-Date) |>   mapview::mapview()



model_dt <- tracks |>
  
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year)
  ) |>
  sf::st_drop_geometry() |>
  as.data.frame()



# model_dt <- tracks |>
#   # dplyr::rename(depth = Depth,
#   #               slope = Slope,
#   #               roughness = Roughness,
#   #               wz = Wz) |> 
#   dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
#                 year  = factor(year),
#                 depth = case_when(depth >= 0 ~ -5, TRUE ~ depth),
#                 id = as.factor(id),
#                 dist2000 = dist2000/1000) |> 
#   sf::st_drop_geometry() |> 
#   as.data.frame()
# tracks_old <- readRDS("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM/Modelling_Simulations_PA_Predictors_Meta_AC_thinned_FINAL.rds")
# str(tracks_old)
# 
# model_dt <- tracks_old |>
#   dplyr::mutate(month = factor(Month, levels = sort(unique(Month))),  # 1..12 etc.
#                 year  = factor(lubridate::year(date)),
#                 #PA = as.numeric(PA)
#                 pres_abs = case_when(rep == 0 ~ 1, TRUE ~ 0)) |>
#   sf::st_drop_geometry() |>
#   as.data.frame()

str(model_dt)

model_dt |> dplyr::summarise(start = min(Date),
                             end = max(Date))

model_dt |> 
  rstatix::get_summary_stats(depth)


model_dt |> 
  as.data.frame() |>
  dplyr::group_by(PA) |> 
  rstatix::get_summary_stats(thetao, type = "common")


# Explore Data ------------------------------------------------------------


plots <- model_dt |> dplyr::mutate(PA = as.factor(PA))
plots

(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist200, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist1000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)



na_counts <- model_dt |> 
  dplyr::summarise(across(everything(), ~ sum(is.na(.))))
na_counts

str(model_dt)

# model_dt <- model_dt_trans

model_dt$month






swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = model_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = model_dt |>
    dplyr::select(thetao, 
                  mltost, 
                  chl, 
                  uv, 
                  wz, 
                  slope, 
                  # roughness,
                  depth, 
                  # dist200,
                  # dist1000,
                  dist2000, 
                  # dist_seamount,
                  # dist_knoll,
                  month) |> 
    as.data.frame(),
  pa = model_dt$PA)

swd

# 
# swd <- SDMtune::SWD(
#   species = "Rhincodon_typus",
#   coords = model_dt |> dplyr::select(lon, lat) |> as.data.frame(),
#   data = model_dt |>
#     dplyr::select(MUR_SST, MLD, CM_Chl, CM_uv, Wz, Depth, Slope, dist2000, month) |> as.data.frame(),
#   pa = model_dt$pres_abs)

swd

swd@data




# # lets use a test and train dataset (which is a bit tricky when having a small sample size): Method = "hold-out' validation 
# library(zeallot)
# output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
# colnames(output) <- c("seed", "trainAUC", "testAUC")
# 
# set.seed(25)
# seeds <- sample.int(1000, 10) # Create 10 different random seeds
# 
# for (i in seq_along(seeds)) { # Loop through the seeds
#   c(train, test) %<-% trainValTest(swd, 
#                                    test = 0.2, 
#                                    seed = seeds[i], 
#                                    only_presence = TRUE) # Make the train/test split
#   
#   m <- train("Maxent", 
#              data = train) # train the model
#   
#   # Populate the output data.frame
#   output[i, 1] <- seeds[i]
#   output[i, 2] <- auc(m)
#   output[i, 3] <- auc(m, test = test)
# }
# 
# output
# 
# m




# Using Cross Validation Method: ------------------------------------------

# get the blockCV spatial block folds
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_3days.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_3to7days.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_1to4days.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final.rds")

sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_4.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_5.rds")

# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_6.rds")






blockCV::cv_plot(cv = sp_blocks,
                 # x = input_occ,
                 # r = predictors_mean[[1]],
                 nrow = 2, 
                 points_alpha = 0.5)

set.seed(28)
cv_m0 <- SDMtune::train(method = "Maxent", 
                        data = swd,
                        folds = sp_blocks)



# cv_m0 <- SDMtune::train(method = "Maxent",
#                         data = swd,
#                         folds = sp_blocks,
#                         fc = "lqph",
#                         reg = 1.5,
#                         iter = 5000)


cv_m0

cat("Training AUC: ", SDMtune::auc(cv_m0))
cat("Testing AUC: ", SDMtune::auc(cv_m0, test = TRUE))
cat("Training TSS: ", SDMtune::tss(cv_m0))
cat("Testing TSS: ", SDMtune::tss(cv_m0, test = TRUE))






ROCplots <- lapply(seq_len(ncol(cv_m0@folds$test)), function(i){
  idx <- cv_m0@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = cv_m0@data@species,
    coords  = cv_m0@data@coords[idx, , drop = FALSE],
    data    = cv_m0@data@data  [idx, , drop = FALSE],
    pa      = cv_m0@data@pa    [idx]
  )
  SDMtune::plotROC(cv_m0@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})




ROCplots_combined <- wrap_plots(ROCplots) +
  plot_layout(ncol = 3, nrow = 2, axes = "collect") &
  ggplot2::theme(
    legend.position = c(0.95, 0.03),  # x, y coordinates (right–bottom)
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ROCplots_combined

## Variable Importance


vi_maxent <- SDMtune::varImp(cv_m0, 
                    permut = 10)

vi_maxent

SDMtune::plotVarImp(vi_maxent)


jk1 <- SDMtune::doJk(cv_m0, 
                     metric = "tss", 
                     test = TRUE)
jk1

SDMtune::plotJk(jk1, 
                type = "train", 
                ref = SDMtune::tss(cv_m0))

SDMtune::plotJk(jk1, 
                type = "test", 
                ref = SDMtune::tss(cv_m0, test = TRUE))


jk2 <- SDMtune::doJk(cv_m0, 
                     metric = "auc", 
                     test = TRUE)
jk2

SDMtune::plotJk(jk2, 
                type = "train", 
                ref = SDMtune::auc(cv_m0))

SDMtune::plotJk(jk2, 
                type = "test", 
                ref = SDMtune::auc(cv_m0, test = TRUE))


m <- cv_m0

P_sst <- SDMtune::plotResponse(m, 
                               var = "thetao", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("sst ("*degree*"C)")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_chl <- SDMtune::plotResponse(m, 
                               var = "chl", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log (chl)") +
  #scale_x_log10() +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_chl

P_uv <- SDMtune::plotResponse(m, 
                              var = "uv", 
                              type = "cloglog", 
                              only_presence = TRUE, 
                              marginal = TRUE, 
                              rug = TRUE,
                              fun = mean,
                              color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("uv (m s"^{-1}*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_wz <- SDMtune::plotResponse(m, 
                              var = "wz", 
                              type = "cloglog", 
                              only_presence = TRUE, 
                              marginal = TRUE, 
                              rug = TRUE,
                              fun = mean,
                              color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("wz (m s"^{-1}*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_mld <- SDMtune::plotResponse(m, 
                               var = "mltost", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = median,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "mld (m)") +
  # scale_x_log10() +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_mld

P_depth <- SDMtune::plotResponse(m, 
                                 var = "depth", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  labs(x = "depth (m)") +
  ggplot2::scale_y_continuous(limits = c(-0.05, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_slope <- SDMtune::plotResponse(m, 
                                 var = "slope", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("slope ("*degree*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_dist <- SDMtune::plotResponse(m, 
                                var = "dist2000", 
                                type = "cloglog", 
                                only_presence = TRUE, 
                                marginal = TRUE, 
                                rug = TRUE,
                                fun = mean,
                                color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "dist2000 (m)") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_month <- SDMtune::plotResponse(m, 
                                 var = "month", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

marginal_plots <- (P_depth + P_dist + P_sst + P_chl + P_uv + P_mld + P_wz + P_slope + P_month) +
  patchwork::plot_layout(ncol = 3, guides = "collect", axes = "collect")

marginal_plots

# Variable Correlations ---------------------------------------------------


SDMtune::plotCor(swd, 
                 method = "spearman", 
                 cor_th = NULL)

## let's draw random background points from the mean environment over study period to see which ones correlate 
predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log_ML.tif")
names(predictors_mean)
plot(predictors_mean)




## according to package vignette, we draw random samples from predictors, however, calibration is not equal prdiction area, thus this is not a good representation for our data, plus we use dynamic SDM, so mean environment is also not suitable torepresent the background. So, we will use our extacted pseudo-absences instead
# set.seed(66)
# bg <- terra::spatSample(predictors_mean,
#                         size = 10000,
#                         method = "random",
#                         na.rm = TRUE,
#                         xy = TRUE,
#                         values = FALSE)

# bg <- SDMtune::prepareSWD(species = "Bgs", 
#                  a = bg, 
#                  env = predictors_mean)


# bg_dt <- model_dt |> dplyr::filter(PA == 0)
bg_dt <- model_dt 
bg4cor <- SDMtune::SWD(
  species = "Bgs",
  coords  = bg_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data    = bg_dt |> dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, 
                                   roughness,
                                   dist2000, 
                                   # dist_seamount, dist_knoll, 
                                   month) |> as.data.frame(),
  pa      = bg_dt$PA)



SDMtune::plotCor(bg4cor, 
                 method = "spearman", 
                 cor_th = NULL)

SDMtune::corVar(bg4cor, 
                method = "spearman", 
                cor_th = 0.6)


cv_m1 <- SDMtune::varSel(cv_m0, 
                         metric = "auc", 
                         test = TRUE, 
                         bg4cor = bg4cor,
                         method = "spearman", 
                         cor_th = 0.6,
                         permut = 10)

cv_m1



# Fine tune model  --------------------------------------------------------


SDMtune::getTunableArgs(cv_m1)

h <- list(reg = seq(0.5, 5, by = 0.5), 
          fc = c("l", "lq", "lh", "lqp", "lqph"),
          # fc  = c("l", "lh", "lqph"),
          iter = c(2500, 5000, 7500, 10000)
)

h


expected_fits(pop = 20, 
              gen = 5, 
              keep_best = 0.4,
              keep_random = 0.2,
              rounding= "round")


# cv_m2 <- SDMtune::optimizeModel(cv_m1, 
#                                 hypers = h, 
#                                 metric = "auc",
#                                 test = TRUE,
#                                 pop = 20,
#                                 gen = 5,
#                                 keep_best = 0.4,
#                                 keep_random = 0.2,
#                                 mutation_chance = 0.4,
#                                 interactive = TRUE,
#                                 progress = TRUE,
#                                 seed = 694
# )

cv_m2 <- SDMtune::optimizeModel(cv_m1, 
                                hypers = h, 
                                metric = "tss",
                                test = TRUE,
                                pop = 20,
                                gen = 5,
                                keep_best = 0.4,
                                keep_random = 0.2,
                                mutation_chance = 0.4,
                                interactive = TRUE,
                                progress = TRUE,
                                seed = 694
)



cv_m2


cv_m2@results
#ordered 
# cv_m2@results[order(-cv_m2@results$test_AUC), ]
cv_m2@results[order(-cv_m2@results$test_TSS), ]

# Index of the best model in the experiment
# index <- terra::which.max(cv_m2@results$test_AUC)
index <- terra::which.max(cv_m2@results$test_TSS)
index <- 19
index


best_cv_m2 <- cv_m2@models[[index]]

best_cv_m2
cv_m2@results[index, ]



cat("Training AUC: ", SDMtune::auc(best_cv_m2))
cat("Testing AUC: ", SDMtune::auc(best_cv_m2, test = TRUE))
cat("Training TSS: ", SDMtune::tss(best_cv_m2))
cat("Testing TSS: ", SDMtune::tss(best_cv_m2, test = TRUE))


ROCplots <- lapply(seq_len(ncol(best_cv_m2@folds$test)), function(i){
  idx <- best_cv_m2@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = best_cv_m2@data@species,
    coords  = best_cv_m2@data@coords[idx, , drop = FALSE],
    data    = best_cv_m2@data@data  [idx, , drop = FALSE],
    pa      = best_cv_m2@data@pa    [idx]
  )
  SDMtune::plotROC(best_cv_m2@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})



ROCplots_combined <- wrap_plots(ROCplots) +
  plot_layout(ncol = 3, nrow = 2, axes = "collect") &
  ggplot2::theme(
    legend.position = c(0.95, 0.03),  # x, y coordinates (right–bottom)
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ROCplots_combined



# Refine using gridSearch -------------------------------------------------

h_refine <- list(
  reg = seq(0.5, 2.5, by = 0.5),
  fc   = c("lqp", "lqph"),
  iter = seq(4000, 10000, by =1000)
)




cv_m3 <- SDMtune::gridSearch(
  model  = cv_m1,
  hypers = h_refine,
  metric = "tss",
  interactive = TRUE,
  progress = TRUE
)



cv_m3

cv_m3@results
#ordered
cv_m3@results[order(-cv_m3@results$test_TSS), ]

# Index of the best model in the experiment
# index <- terra::which.max(cv_m3@results$test_TSS)
index <- 8


best_cv_m3 <- cv_m3@models[[index]]
best_cv_m3
cv_m3@results[index, ]


cat("Training TSS: ", SDMtune::tss(best_cv_m3))
cat("Testing TSS: ", SDMtune::tss(best_cv_m3, test = TRUE))
cat("Training AUC: ", SDMtune::auc(best_cv_m3))
cat("Testing AUC: ", SDMtune::auc(best_cv_m3, test = TRUE))


# check if refine made it better:
cat("Testing AUC: ", SDMtune::auc(cv_m1, test = TRUE))
cat("Testing AUC: ", SDMtune::auc(best_cv_m3, test = TRUE))
cat("Testing TSS: ", SDMtune::tss(cv_m1, test = TRUE))
cat("Testing TSS: ", SDMtune::tss(best_cv_m3, test = TRUE))

cat("Testing AUC before tuning: ", SDMtune::auc(best_cv_m2, test = TRUE))
cat("Testing AUC after tuning: ", SDMtune::auc(best_cv_m3, test = TRUE))
cat("Testing TSS before tuning: ", SDMtune::tss(best_cv_m2, test = TRUE))
cat("Testing TSS after tuning: ", SDMtune::tss(best_cv_m3, test = TRUE))

cat("Training AUC before tuning: ", SDMtune::auc(best_cv_m2))
cat("Training AUC after tuning: ", SDMtune::auc(best_cv_m3))
cat("Training TSS before tuning: ", SDMtune::tss(best_cv_m2))
cat("Training TSS after tuning: ", SDMtune::tss(best_cv_m3))


m <- best_cv_m3
ROCplots <- lapply(seq_len(ncol(m@folds$test)), function(i){
  idx <- m@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = m@data@species,
    coords  = m@data@coords[idx, , drop = FALSE],
    data    = m@data@data  [idx, , drop = FALSE],
    pa      = m@data@pa    [idx]
  )
  SDMtune::plotROC(m@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})

ROCplots








# Removing variables with low importance  --------------------------------

SDMtune::varImp(best_cv_m3,
                permut = 10)

# SDMtune::varImp(cv_m0,
#                 permut = 10)


best_cv_m3

cv_m4 <- SDMtune::reduceVar(best_cv_m3, 
                            th = 0.5, 
                            metric = "tss", 
                            test = TRUE, 
                            permut = 10, 
                            use_jk = TRUE,
                            interactive = TRUE,
                            verbose = TRUE)

# cv_m4 <- readRDS("models/objects/Tracks_Maxent_CV_Models_SDMtune_mp.rds")

cat("Testing TSS before: ", SDMtune::tss(best_cv_m2, test = TRUE))
cat("Testing TSS after: ", SDMtune::tss(cv_m4, test = TRUE))
cat("Testing AUC before: ", SDMtune::auc(best_cv_m2, test = TRUE))
cat("Testing AUC after: ", SDMtune::auc(cv_m4, test = TRUE))

cat("Traning TSS before: ", SDMtune::tss(best_cv_m2))
cat("Traning TSS after: ", SDMtune::tss(cv_m4))
cat("Traning AUC before: ", SDMtune::auc(best_cv_m2))
cat("Traning AUC after: ", SDMtune::auc(cv_m4))

# cv_m4 <- best_cv_m2

m <- cv_m4
cv_m4

ROCplots <- lapply(seq_len(ncol(m@folds$test)), function(i){
  idx <- m@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = m@data@species,
    coords  = m@data@coords[idx, , drop = FALSE],
    data    = m@data@data  [idx, , drop = FALSE],
    pa      = m@data@pa    [idx]
  )
  SDMtune::plotROC(m@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})

m@folds$train

ROCplots_combined <- wrap_plots(ROCplots) +
  plot_layout(ncol = 3, nrow = 2, axes = "collect") &
  ggplot2::theme(
    legend.position = c(0.95, 0.03),  # x, y coordinates (right–bottom)
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ROCplots_combined


# helper to build an SWD from row indices of the original SWD
.build_swd <- function(obj, idx) {
  SDMtune::SWD(
    species = obj@data@species,
    coords  = obj@data@coords[idx, , drop = FALSE],
    data    = obj@data@data  [idx, , drop = FALSE],
    pa      = obj@data@pa    [idx]
  )
}

k <- ncol(m@folds$test)

metrics_per_fold <- data.frame(
  fold      = seq_len(k),
  auc_train = vapply(seq_len(k), function(i) {
    idx_tr <- m@folds$train[, i]
    swd_tr <- .build_swd(m, idx_tr)
    SDMtune::auc(m@models[[i]], test = swd_tr)
  }, numeric(1)),
  auc_test  = vapply(seq_len(k), function(i) {
    idx_te <- m@folds$test[, i]
    swd_te <- .build_swd(m, idx_te)
    SDMtune::auc(m@models[[i]], test = swd_te)
  }, numeric(1)),
  tss_train = vapply(seq_len(k), function(i) {
    idx_tr <- m@folds$train[, i]
    swd_tr <- .build_swd(m, idx_tr)
    SDMtune::tss(m@models[[i]], test = swd_tr)
  }, numeric(1)),
  tss_test  = vapply(seq_len(k), function(i) {
    idx_te <- m@folds$test[, i]
    swd_te <- .build_swd(m, idx_te)
    SDMtune::tss(m@models[[i]], test = swd_te)
  }, numeric(1))
)

# rounded view and fold means (optional)
metrics_per_fold_rounded <- metrics_per_fold |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ base::round(.x, 2)))


cv_summary <- base::list(
  metrics_per_fold |>
    dplyr::summarise(
      dplyr::across(c(auc_train, auc_test, tss_train, tss_test), ~ base::mean(.x))
    ) |>
    dplyr::mutate(stat = "mean"),
  metrics_per_fold |>
    dplyr::summarise(
      dplyr::across(c(auc_train, auc_test, tss_train, tss_test), ~ stats::sd(.x))
    ) |>
    dplyr::mutate(stat = "sd")
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(dplyr::across(-stat, ~ base::round(.x, 3))) |>
  dplyr::relocate(stat)

metrics_per_fold_rounded
cv_summary



cv_m4


max_vi <- SDMtune::varImp(cv_m4, permut = 10)

# 2) Turn it into a tidy data.frame/tibble
max_vi_df <- max_vi |>
  tibble::as_tibble(rownames = "variable") |>
  dplyr::rename(importance = 2L) |>
  dplyr::mutate(algorithm = "MaxEnt")


max_vi_df

saveRDS(max_vi_df, "models/tables/MaxEnt_Var_Importance_CV.rds")


jk4 <- SDMtune::doJk(cv_m4, 
                     metric = "tss", 
                     test = TRUE)
jk4

SDMtune::plotJk(jk4, 
                type = "train", 
                ref = SDMtune::auc(cv_m4))

SDMtune::plotJk(jk4, 
                type = "test", 
                ref = SDMtune::auc(cv_m4, test = TRUE))


m <- cv_m4
# m <- best_cv_m2
P_sst <- SDMtune::plotResponse(m, 
                               var = "thetao", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("sst ("*degree*"C)")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_chl <- SDMtune::plotResponse(m, 
                               var = "chl", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("chl (log(mg m"^{-3}*"))")) +
  #scale_x_log10() +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_chl

P_uv <- SDMtune::plotResponse(m, 
                              var = "uv", 
                              type = "cloglog", 
                              only_presence = TRUE, 
                              marginal = TRUE, 
                              rug = TRUE,
                              fun = mean,
                              color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("uv (m s"^{-1}*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_wz <- SDMtune::plotResponse(m, 
                              var = "wz", 
                              type = "cloglog", 
                              only_presence = TRUE, 
                              marginal = TRUE, 
                              rug = TRUE,
                              fun = mean,
                              color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("wz (m s"^{-1}*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_mld <- SDMtune::plotResponse(m, 
                               var = "mltost", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = median,
                               color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "mld (m)") +
  # scale_x_log10() +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_mld

P_depth <- SDMtune::plotResponse(m, 
                                 var = "depth", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  labs(x = "depth (m)") +
  ggplot2::scale_y_continuous(limits = c(-0.05, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_slope <- SDMtune::plotResponse(m, 
                                 var = "slope", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("slope ("*degree*")")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_dist <- SDMtune::plotResponse(m, 
                                var = "dist2000", 
                                type = "cloglog", 
                                only_presence = TRUE, 
                                marginal = TRUE, 
                                rug = TRUE,
                                fun = mean,
                                color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "dist2000 (m)") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_month <- SDMtune::plotResponse(m, 
                                 var = "month", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "hotpink") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

marginal_plots <- (P_depth + P_dist + P_sst + P_chl + P_uv + P_mld + P_wz + P_slope + P_month) +
  patchwork::plot_layout(ncol = 3, guides = "collect", axes = "collect")

marginal_plots


ggsave("MaxEnt_RepsonsePlots_Marginal_Tracking_mp.png", plot = marginal_plots, path ="outputs/final_figures", scale =1, width = 18, height = 20, units = "cm", dpi = 300)




plot_brt_interaction_3d <- function(model,
                                    var1,
                                    var2,
                                    n = 50,
                                    fun = stats::median) {
  
  if (!methods::is(model, "SDMmodel")) {
    base::stop("`model` must be an SDMmodel object (e.g., a single BRT, not SDMmodelCV).")
  }
  
  # Training data
  x <- model@data@data |>
    base::as.data.frame()
  
  if (!(var1 %in% base::names(x))) {
    base::stop("Variable ", var1, " not found in model@data@data.")
  }
  if (!(var2 %in% base::names(x))) {
    base::stop("Variable ", var2, " not found in model@data@data.")
  }
  
  # Sequences over the two focal variables
  v1_seq <- base::seq(
    from = base::min(x[[var1]], na.rm = TRUE),
    to   = base::max(x[[var1]], na.rm = TRUE),
    length.out = n
  )
  
  v2_seq <- base::seq(
    from = base::min(x[[var2]], na.rm = TRUE),
    to   = base::max(x[[var2]], na.rm = TRUE),
    length.out = n
  )
  
  # Grid over the two variables
  grid <- base::expand.grid(
    var1 = v1_seq,
    var2 = v2_seq
  )
  base::names(grid)[1:2] <- c(var1, var2)
  
  # Hold other predictors constant
  other_vars <- base::setdiff(base::names(x), c(var1, var2))
  
  for (nm in other_vars) {
    if (base::is.numeric(x[[nm]])) {
      grid[[nm]] <- fun(x[[nm]], na.rm = TRUE)
    } else if (base::is.factor(x[[nm]]) || base::is.character(x[[nm]])) {
      tab <- base::table(x[[nm]])
      grid[[nm]] <- base::names(tab)[base::which.max(tab)]
    } else {
      grid[[nm]] <- fun(x[[nm]], na.rm = TRUE)
    }
  }
  
  # Predict
  grid$pred <- SDMtune::predict(model, data = grid)
  
  # Turn predictions into matrix for surface plot
  z_mat <- base::matrix(
    grid$pred,
    nrow = base::length(v1_seq),
    ncol = base::length(v2_seq),
    byrow = FALSE
  )
  
  # 3D surface with plotly
  p <- plotly::plot_ly(
    x = v1_seq,
    y = v2_seq,
    z = z_mat
  ) |>
    plotly::add_surface() |>
    plotly::layout(
      scene = list(
        xaxis = list(title = var1),
        yaxis = list(title = var2),
        zaxis = list(title = "BRT prediction")
      )
    )
  
  return(p)
}


m <- final_m
plot_brt_interaction_3d(
  model = final_brt,
  var1  = "thetao",
  var2  = "chl",
  fun = mean
)


plot_brt_interaction_3d(
  model = final_brt,
  var1  = "uv",
  var2  = "chl",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "wz",
  var2  = "chl",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "thetao",
  var2  = "wz",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "depth",
  var2  = "wz",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "slope",
  var2  = "wz",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "depth",
  var2  = "thetao",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "depth",
  var2  = "dist2000",
  fun = mean
)

plot_brt_interaction_3d(
  model = final_brt,
  var1  = "slope",
  var2  = "dist2000",
  fun = mean
)
























# get final model 
set.seed(25)
final_m <- SDMtune::combineCV(cv_m4)
# final_m <- SDMtune::combineCV(best_cv_m2)

SDMtune::plotROC(final_m)
final_m@model
final_m

cv_m0


relevant_vars <- names(final_m@data@data)
relevant_vars
predictors <- predictors_mean[[relevant_vars]]
names(predictors_mean)

e <- terra::ext(c(140, 170, -40, 0))


predictors <-  predictors |>
  terra::crop(e, snap = "out") 

plot(predictors)


map <- SDMtune::predict(final_m,
                        data = predictors,
                        type = "cloglog",
                        const = (data.frame(month = "0"))
                        )

SDMtune::plotPred(map,
                  lt = "Habitat\nsuitability",
                  colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

map

ths <- SDMtune::thresholds(final_m, 
                           type = "cloglog")

ths

SDMtune::plotPA(map, 
                th = ths[3, 2])

SDMtune::plotPA(map, 
                th = 0.7)







SDMtune::modelReport(final_m, 
                     type = "cloglog", 
                     folder = "mmodels/objects/Tracks_Maxent_tuned_mp_crwPA_final", 
                     test = NULL, 
                     response_curves = TRUE, 
                     only_presence = TRUE, 
                     jk = FALSE, 
                     env = predictors_mean)



# saveRDS(cv_m4, "models/objects/Tracks_Maxent_CV_Models_SDMtune_mp.rds")
# saveRDS(final_m, "models/objects/Tracks_Maxent_final_Model_SDMtune_mp.rds")

saveRDS(cv_m4, "models/objects/Tracks_Maxent_CV_Models_SDMtune_mp_crwPA_final.rds")
saveRDS(final_m, "models/objects/Tracks_Maxent_final_Model_SDMtune_mp_crwPA_final.rds")

cv_m4 <- readRDS("models/objects/Tracks_Maxent_CV_Models_SDMtune_mp_crwPA_final.rds")
final_m <- readRDS( "models/objects/Tracks_Maxent_final_Model_SDMtune_mp_crwPA_final.rds")
final_m
# final_m_good <- readRDS("models/objects/Tracks_Maxent_final_Model_SDMtune_mp.rds")
# str(final_m_good@data@data)
# final_m_good
# final_m

# SDMtune::varImp(final_m_good, 
#                 permut = 10)


# create extrnal validation data from sightings 
sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full.rds")
str(sight)
sight |> dplyr::group_by(PA) |>  dplyr::filter(lon > 140) |>  dplyr::summarise(N = n())


mapview::mapview(sight |> dplyr::select(-Date) |>  dplyr::filter(PA==1 & lon > 140)) 


val_dt <- sight |>
  # dplyr::filter(PA == 1) |>
  dplyr::filter(lon >150 & lon < 180) |> # making external dataset spatially independent of calibrationn data
  # dplyr::filter(lon >142 & lon < 160) |>
  # dplyr::filter(lat >-17 & lat < -5) |>
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year)) |> 
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |>
  sf::st_drop_geometry() |> 
  # dplyr::filter(!month %in% c(7, 8, 9, 10)) |>
  as.data.frame()

str(val_dt)


val_dt <- val_dt |> 
  dplyr::mutate(
    # slope = log1p(slope),
    # roughness = log1p(roughness),
    chl = log(chl),
    dist2000 = dist2000/1000)



plots <- val_dt |> dplyr::mutate(PA = as.factor(PA))

(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


val_dt |> dplyr::summarise(start = min(Date),
                           end = max(Date))


mapview::mapview(val_dt |> dplyr::select(-Date) |> dplyr::filter(PA == 1) |> sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE))

val.swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = val_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = val_dt |>
    dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, dist2000, month) |> as.data.frame(),
  pa = val_dt$PA)

val.swd


# quartz()
SDMtune::auc(final_m, 
             test = val.swd)
SDMtune::tss(final_m, 
             test = val.swd)
SDMtune::plotROC(final_m, test = val.swd)  +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "External validation (AUC = 0.77)")



final_m@model

pred <- SDMtune::predict(final_m, data = val.swd)
labels <- val.swd@pa

# Use pROC for reliable plotting
pROC::roc(labels, pred, plot = TRUE, col = "steelblue", lwd = 2, legacy.axes = TRUE)

roc_obj <- pROC::roc(labels, pred)
roc_df <- data.frame(
  tpr = rev(roc_obj$sensitivities),
  fpr = rev(1 - roc_obj$specificities)
)

ggplot2::ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = sprintf("ROC curve (AUC = %.3f)", roc_obj$auc))



ths_ext <- SDMtune::thresholds(
  model = final_m,
  test  = val.swd    # external validation SWD
  # type = "cloglog" # only needed for Maxent/Maxnet with specific output type
)

ths_ext


pred_val <- SDMtune::predict(final_m, data = val.swd)
labels   <- val.swd@pa

# Boyce index on external validation set
boyce_ext <- ecospat::ecospat.boyce(
  fit   = pred_val,                 # predictions for ALL external points
  obs   = pred_val[labels == 1],    # predictions at external presences
  PEplot = TRUE
)



boyce_ext$cor







# Residuals autocorrealtion -----------------------------------------------

m <- best_cv_m3

# 1) Extract SWD data used to fit the model
swd.r <- m@data   # SDMtune::SWD object

# 2) Fitted probabilities for those points
# For MaxEnt: use cloglog or logistic, NOT raw
fitted <- SDMtune::predict(m, data = swd.r, type = "cloglog")
# or type = "prob" depending on how you set it up; usually "raw" is fine for Bernoulli

# 3) Observed 0/1 (presence/absence or presence/background)
obs <- swd.r@pa

# 4) Simple residuals (you can also use Pearson residuals if you want)
resid <- obs - fitted

# Put everything into a data.frame with coordinates
resid_df <- data.frame(
  x      = swd.r@coords[, 1],
  y      = swd.r@coords[, 2],
  obs    = obs,
  fitted = fitted,
  resid  = resid
)


resid_df
resid_sf <- resid_df |>
  sf::st_as_sf(coords = c("x", "y"), crs = 4326) |>
  sf::st_transform(crs = 3577)  # or another equal-area / metric CRS you’re using

resid_sf
str(resid_sf)
predictors
ac_resid <- blockCV::cv_spatial_autocor(
  # r      = predictors,        # SpatRaster / Raster* is fine
  x      = resid_sf,     
  column = "resid"       # use residuals, not the raw PA
  # other args: size, k, etc. you can leave default to let it explore lags
)

summary(ac_resid)






# Boyce index -------------------------------------------------------------



boyce_from_cv_by_fold <- function(cv_model,
                                  type = "cloglog",
                                  use_background = TRUE,
                                  boyce_res = 100L,
                                  boyce_nclass = 0L,
                                  plot = FALSE) {
  k <- base::ncol(cv_model@folds$test)
  
  boyce_train  <- base::numeric(k)
  boyce_test   <- base::numeric(k)
  n_train_bins <- base::integer(k)
  n_test_bins  <- base::integer(k)
  
  # helper to compute Boyce (Spearman rho) from preds + pa
  compute_boyce <- function(preds, pa) {
    obs <- preds[pa == 1L]
    fit <- if (use_background) preds[pa == 0L] else preds
    
    obs <- obs[stats::complete.cases(obs)]
    fit <- fit[stats::complete.cases(fit)]
    
    boyce <- ecospat::ecospat.boyce(
      fit      = fit,
      obs      = obs,
      window.w = "default",
      nclass   = boyce_nclass,
      res      = boyce_res,
      PEplot   = plot
    )
    
    df <- base::data.frame(HS = boyce$HS, PE = boyce$F.ratio) |>
      dplyr::filter(stats::complete.cases(HS, PE))
    
    ct <- stats::cor.test(df$HS, df$PE, method = "spearman", exact = FALSE)
    
    list(
      rho     = unname(ct$estimate),
      p_value = ct$p.value,
      n_bins  = base::nrow(df)
    )
  }
  
  for (i in base::seq_len(k)) {
    # --- indices ---
    test_idx <- cv_model@folds$test[, i]
    
    if ("train" %in% base::names(cv_model@folds)) {
      train_idx <- cv_model@folds$train[, i]
    } else {
      train_idx <- !test_idx
    }
    
    # --- build SWDs ---
    train_swd <- SDMtune::SWD(
      species = cv_model@data@species[train_idx],
      coords  = cv_model@data@coords[train_idx, , drop = FALSE],
      data    = cv_model@data@data  [train_idx, , drop = FALSE],
      pa      = cv_model@data@pa    [train_idx]
    )
    
    test_swd <- SDMtune::SWD(
      species = cv_model@data@species[test_idx],
      coords  = cv_model@data@coords[test_idx, , drop = FALSE],
      data    = cv_model@data@data  [test_idx, , drop = FALSE],
      pa      = cv_model@data@pa    [test_idx]
    )
    
    # --- predictions ---
    preds_train <- SDMtune::predict(cv_model@models[[i]],
                                    data = train_swd,
                                    type = type)
    
    preds_test  <- SDMtune::predict(cv_model@models[[i]],
                                    data = test_swd,
                                    type = type)
    
    # --- Boyce train ---
    bt <- compute_boyce(preds_train, train_swd@pa)
    boyce_train[i]  <- bt$rho
    n_train_bins[i] <- bt$n_bins
    
    # --- Boyce test ---
    bte <- compute_boyce(preds_test, test_swd@pa)
    boyce_test[i]  <- bte$rho
    n_test_bins[i] <- bte$n_bins
  }
  
  fold_results <- base::data.frame(
    fold         = base::seq_len(k),
    boyce_train  = boyce_train,
    boyce_test   = boyce_test,
    n_train_bins = n_train_bins,
    n_test_bins  = n_test_bins
  )
  
  summary <- list(
    train_mean = mean(boyce_train, na.rm = TRUE),
    train_sd   = sd(boyce_train,   na.rm = TRUE),
    test_mean  = mean(boyce_test,  na.rm = TRUE),
    test_sd    = sd(boyce_test,    na.rm = TRUE)
  )
  
  list(
    folds   = fold_results,
    summary = summary
  )
}



boyce_cv_m <- boyce_from_cv_by_fold(
  cv_model       = cv_m4,
  type           = "cloglog",
  use_background = TRUE,
  boyce_res      = 100L,
  boyce_nclass   = 0L,
  plot           = FALSE
)

boyce_cv_m$folds    # per-fold Boyce for train + test
boyce_cv_m$summary  # mean ± sd (train/test)





set.seed(1)

# 1) data
x  <- res_boyce$boyce
df <- data.frame(HS = x$HS, PE = x$F.ratio)
df <- df[is.finite(df$HS) & is.finite(df$PE), ]

ct <- suppressWarnings(cor.test(df$HS, df$PE,
                                method = "spearman", exact = FALSE))

# 2) grid + bootstrap loess
hs_grid <- seq(min(df$HS), max(df$HS), length.out = 200)
B <- 500L      # increase to 1000+ for publication

pred_mat <- matrix(NA_real_, nrow = length(hs_grid), ncol = B)
for (b in seq_len(B)) {
  i <- sample.int(nrow(df), replace = TRUE)
  fit <- stats::loess(PE ~ HS, data = df[i, ], span = 0.5, degree = 2, surface = "direct")
  pred_mat[, b] <- stats::predict(fit, newdata = data.frame(HS = hs_grid))
}

band <- data.frame(
  HS  = hs_grid,
  lo  = apply(pred_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
  med = apply(pred_mat, 1, stats::quantile, probs = 0.500, na.rm = TRUE),
  hi  = apply(pred_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
)

# 3) (optional) bootstrap CI for Spearman ρ too
boot_rho <- replicate(B, {
  i <- sample.int(nrow(df), replace = TRUE)
  suppressWarnings(stats::cor(df$HS[i], df$PE[i], method = "spearman"))
})
rho_hat <- unname(stats::cor(df$HS, df$PE, method = "spearman"))
rho_ci  <- stats::quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)


df_plot   <- subset(df, is.finite(HS) & is.finite(PE))
df_plot
band_plot <- subset(band, is.finite(HS) & is.finite(lo) & is.finite(hi) & is.finite(med))

ggplot(df_plot, aes(x = HS, y = PE)) +
  geom_ribbon(data = band_plot, aes(x = HS, ymin = lo, ymax = hi),
              inherit.aes = FALSE, alpha = 0.15, fill = "blue") +
  geom_line(data = band_plot, aes(x = HS, y = med), linewidth = 1, color = "blue") +
  geom_point(shape = 21, fill = NA, stroke = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x = "Habitat suitability (HS)", y = "P/E ratio",
       subtitle = sprintf("Boyce (Spearman ρ) = %.3f,  p = %.3g, (95%% CI: %.3f–%.3f)",
                          rho_hat, ct$p.value, rho_ci[1], rho_ci[2])) +
  coord_cartesian(ylim = c(0, max(df_plot$PE, na.rm = TRUE))) +
  theme_bw()





# Kappa  ------------------------------------------------------------------

preds <- SDMtune::predict(final_m, data = swd, type = "cloglog")
obs   <- swd@pa
obs



library(flexsdm)
preds_pres <- preds[swd@pa == 1]
preds_abs  <- preds[swd@pa == 0]
eval_tbl   <- sdm_eval(p = preds_pres,
                       a = preds_abs,
                       thr = c("max_sens_spec","max_jaccard"))
print(eval_tbl)


# Extract numeric threshold value
thr <- eval_tbl$thr_value[eval_tbl$threshold == "max_sens_spec"]

# 3. Create a dataframe with observed and predicted classes
df <- tibble(
  obs = factor(obs, levels = c(0, 1)),
  pred = preds,
  pred_bin = factor(ifelse(preds >= thr, 1, 0), levels = c(0, 1))
)

# 4. Compute AUC, Kappa, TSS, Boyce
auc_val   <- SDMtune::auc(final_m, test = swd)
tss_val   <- SDMtune::tss(final_m, test = swd)
kappa_val <- yardstick::kap(df, truth = obs, estimate = pred_bin)$.estimate


performance_tbl <- tibble::tibble(
  Metric = c("AUC", "TSS", "Kappa", "Boyce"),
  Value  = round(c(auc_val, tss_val, kappa_val, rho_hat), 3)
)

performance_tbl


# Omission rate -----------------------------------------------------------


omission_cv_maxTSS_train_test <- function(cv_model) {
  
  k <- base::ncol(cv_model@folds$test)
  
  # Determine prediction type once
  type <- if (inherits(cv_model@models[[1]]@model, "Maxent")) "raw" else "link"
  
  get_fold_omission <- function(i) {
    test_idx <- cv_model@folds$test[, i]
    
    # If train indices exist, use them; otherwise use complement of test
    if ("train" %in% base::names(cv_model@folds)) {
      train_idx <- cv_model@folds$train[, i]
    } else {
      train_idx <- !test_idx
    }
    
    # --- build SWD objects ---
    train_swd <- SDMtune::SWD(
      species = cv_model@data@species[train_idx],
      coords  = cv_model@data@coords[train_idx, , drop = FALSE],
      data    = cv_model@data@data  [train_idx, , drop = FALSE],
      pa      = cv_model@data@pa    [train_idx]
    )
    
    test_swd <- SDMtune::SWD(
      species = cv_model@data@species[test_idx],
      coords  = cv_model@data@coords[test_idx, , drop = FALSE],
      data    = cv_model@data@data  [test_idx, , drop = FALSE],
      pa      = cv_model@data@pa    [test_idx]
    )
    
    # --- TRAIN confusion matrix + maxTSS ---
    cm_train <- SDMtune::confMatrix(
      model = cv_model@models[[i]],
      test  = train_swd,
      type  = type
    )
    
    tpr_train <- cm_train$tp / (cm_train$tp + cm_train$fn)
    tnr_train <- cm_train$tn / (cm_train$fp + cm_train$tn)
    tss_train <- tpr_train + tnr_train - 1
    
    i_best_train   <- which.max(tss_train)
    omission_train <- 1 - tpr_train[i_best_train]
    
    # --- TEST confusion matrix + maxTSS ---
    cm_test <- SDMtune::confMatrix(
      model = cv_model@models[[i]],
      test  = test_swd,
      type  = type
    )
    
    tpr_test <- cm_test$tp / (cm_test$tp + cm_test$fn)
    tnr_test <- cm_test$tn / (cm_test$fp + cm_test$tn)
    tss_test <- tpr_test + tnr_test - 1
    
    i_best_test   <- which.max(tss_test)
    omission_test <- 1 - tpr_test[i_best_test]
    
    tibble::tibble(
      fold              = i,
      threshold_train   = cm_train$th[i_best_train],
      omission_train    = omission_train,
      sensitivity_train = tpr_train[i_best_train],
      specificity_train = tnr_train[i_best_train],
      tss_train         = tss_train[i_best_train],
      threshold_test    = cm_test$th[i_best_test],
      omission_test     = omission_test,
      sensitivity_test  = tpr_test[i_best_test],
      specificity_test  = tnr_test[i_best_test],
      tss_test          = tss_test[i_best_test]
    )
  }
  
  folds_df <- purrr::map_dfr(base::seq_len(k), get_fold_omission)
  
  summary_df <- folds_df |>
    dplyr::summarise(
      mean_omission_train = base::mean(omission_train),
      sd_omission_train   = stats::sd(omission_train),
      mean_tss_train      = base::mean(tss_train),
      sd_tss_train        = stats::sd(tss_train),
      mean_omission_test  = base::mean(omission_test),
      sd_omission_test    = stats::sd(omission_test),
      mean_tss_test       = base::mean(tss_test),
      sd_tss_test         = stats::sd(tss_test)
    )
  
  list(
    folds   = folds_df,   # per-fold omission/TSS at each set's own maxTSS
    summary = summary_df  # mean ± sd across folds for train and test
  )
}

om_cv <- omission_cv_maxTSS_train_test(cv_m4)

om_cv$folds    # per-fold train + test omission/TSS/threshold
om_cv$summary  # mean ± sd train vs test



# MOnthly predictions -----------------------------------------------------


relevant_vars <- names(final_m@data@data)
relevant_vars

input_list <- monthly_stacks_lst_0.1_trans


input_list <- lapply(input_list, function(stack) {
  stack[[relevant_vars]] # Extract layers that match the model variables
})

plot(input_list[[1]])

# 1) Get start/end from your model data
rng <- range(as.Date(model_dt$Date), na.rm = TRUE)
start_mon <- as.Date(format(rng[1], "%Y-%m-01"))
end_mon   <- as.Date(format(rng[2], "%Y-%m-01"))
want <- seq(start_mon, end_mon, by = "1 month")

# 2) Align your monthly rasters to that range
avail <- as.Date(names(monthly_stacks_lst_0.1_trans))  # "YYYY-MM-01"
idx <- match(want, avail)

if (anyNA(idx)) {
  warning("Missing months in raster list: ",
          paste(as.character(want[is.na(idx)]), collapse = ", "))
}

monthly_predictors   <- monthly_stacks_lst_0.1_trans[idx[!is.na(idx)]]
monthly_dates  <- avail[idx[!is.na(idx)]]


e <- terra::ext(c(140, 170, -40, 0))
e
tictoc::tic("Monthly dynamic predictions took: " )
monthly_predictions_stack <- sdmtune_predict_monthly(
  model = final_m,
  monthly_list = monthly_predictors,
  dates = monthly_dates,  # already "YYYY-MM-01"
  type = "cloglog",
  extent = e,
  verbose = TRUE
)
tictoc::toc()


monthly_predictions_stack
names(monthly_predictions_stack)
plot(monthly_predictions_stack[[56:68]], range = c(0, 1))



SDMtune::plotPred(monthly_predictions_stack[[64]],
                  lt = "Habitat\nsuitability",
                  colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


ths <- SDMtune::thresholds(final_m, 
                           type = "cloglog")

ths

SDMtune::plotPA(monthly_predictions_stack[[64]],
                th = ths[3, 2])


writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_2018_2025_rev_mp_crwPA.tif", overwrite = TRUE)

# monthly_predictions_stack <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_2018_2025_rev_mp.tif")

### Monthly Means over tracking period

# Define list for monthly averages:
monthly_means <- list()

# Loop through each month:
for (i in 1:12){
  
  # Extract the layers corresponding to the current month
  month_layers <- monthly_predictions_stack[[grep(sprintf("%02d", i), names(monthly_predictions_stack))]]
  
  # Check if layers exist for the current month
  if (length(month_layers) > 0) {
    
    # Calculate the mean across all layers for the current month
    #monthly_mean <- calc(month_layers, mean, na.rm = TRUE)
    monthly_mean <- terra::app(month_layers, fun = mean, na.rm = TRUE)
    
    # Assign meaningful names:
    names(monthly_mean) <- paste(month.abb[i], "Mean")
    
    # Print current layer being processed
    print(paste("Calculating mean for", month.abb[i]))
    
    # Store in list
    monthly_means[[i]] <- monthly_mean
  }
}


monthly_means

monthly_means_stack <- terra::rast(monthly_means)
names(monthly_means_stack) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
monthly_means_stack


# plot(monthly_means_stack, range = c(0, 1))
plot(monthly_means_stack[[12]], range = c(0, 1))
plot(monthly_means_stack[[5]], range = c(0, 1))


# Save for ensemble modelling
writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_monthly_means_rev_crwPA.tif", overwrite = TRUE)


## Overall mean 


mean_climate <- terra::app(monthly_means_stack, fun = mean, na.rm = TRUE)
mean_climate
plot(mean_climate)

# Save for ensemble modelling
# writeRaster(mean_climate, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_Maxent_mean_rev_mp.tif", overwrite = TRUE)

writeRaster(mean_climate, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_Maxent_mean_rev_mp_crwPA.tif", overwrite = TRUE)





### Monsoon vs Trade Wind::



seasons_2 <- c("Monsoon", "Dry")

# Define the months corresponding to each season
season_months_2 <- list(
  Monsoon = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr"),
  Dry = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
)


# Initialize an empty stack to store seasonal means
seasonal_means_stack_2 <- terra::rast()

# Loop through each season
for (season in seasons_2) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months_2[[season]]
  season_stack <- monthly_means_stack[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"Mean", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_means_stack_2 <- c(seasonal_means_stack_2, seasonal_mean)
}


seasonal_means_stack_2


seasonal_means_stack_2 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp_crwPA.tif")

plot(seasonal_means_stack_2, range = c(0, 1))
plot(seasonal_means_stack_2, range = c(0, 1), xlim = c(140, 170), ylim = c(-40, 0), axes =FALSE, legend =FALSE)

plot(seasonal_means_stack_2, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))



max_values <- terra::app(seasonal_means_stack_2, fun = "max", na.rm = TRUE)
high <- global(max_values, fun = "max", na.rm = TRUE) |> as.numeric()
high = 1

Mean_Season_monsoon_SDM_plot <-  rasterVis::levelplot(seasonal_means_stack_2,
                                                      layout=c(2, 1),
                                                      main = "Seasonal whale shark habitat suitability",
                                                      names.attr = c("Monsoon (Nov - Apr)", 
                                                                     "Dry (May - Oct)"),
                                                      zlim = c(0, high),
                                                      at = seq(0, high, length.out = 100),
                                                      # col.regions =  rev(topo.colors(100)),
                                                      # col.regions = cm_ocean_palette,
                                                      # col.regions =  viridisLite::inferno(200),
                                                      col.regions =  colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100),
                                                      colorkey = list(space = "right", 
                                                                      length = 0.5, 
                                                                      height = 0.75,
                                                                      labels = list(
                                                                        at = c(0, high),  # Positions for "Low" and "High"
                                                                        labels = c("Low", "High"))))  # Labels to use

Mean_Season_monsoon_SDM_plot




# Save for ensemble modelling
# writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp.tif", overwrite = TRUE)

writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp_crwPA.tif", overwrite = TRUE)




seasons <- c("Q1", "Q2", "Q3", "Q4")
season_months <- list(
  Q1 = c("Jan", "Feb", "Mar"),
  Q2 = c("Apr", "May", "Jun"),
  Q3 = c("Jul", "Aug", "Sep"),
  Q4 = c("Oct", "Nov", "Dec"))


# Initialize an empty stack to store seasonal means
seasonal_means_stack_4 <- terra::rast()

# Loop through each season
for (season in seasons) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months[[season]]
  season_stack <- monthly_means_stack[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"Mean", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_means_stack_4 <- c(seasonal_means_stack_4, seasonal_mean)
}


seasonal_means_stack_4

plot(seasonal_means_stack_4, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))

Mean_Season_quartal_SDM_plot <-  rasterVis::levelplot(seasonal_means_stack_4,
                                                      layout=c(4, 1),
                                                      main = "Seasonal whale shark habitat suitability",
                                                      names.attr = c("Jan - Mar", 
                                                                     "Apr - Jun", 
                                                                     "Jul - Sep", 
                                                                     "Oct - Dec"),
                                                      zlim = c(0, high),
                                                      at = seq(0, high, length.out = 100),
                                                      #col.regions =  rev(topo.colors(200)),
                                                      # col.regions = cm_ocean_palette, 
                                                      col.regions =  colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100),
                                                      colorkey = list(space = "right", 
                                                                      length = 0.5, 
                                                                      height = 0.75,
                                                                      labels = list(
                                                                        at = c(0, high),  # Positions for "Low" and "High"
                                                                        labels = c("Low", "High"))))  # Labels to use

Mean_Season_quartal_SDM_plot


# writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons4_means_rev_mp.tif", overwrite = TRUE)

writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons4_means_rev_mp_crwPA.tif", overwrite = TRUE)




