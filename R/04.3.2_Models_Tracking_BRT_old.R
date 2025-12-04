#_____________________________________________________________________________
#                        Models: Tracking - BRT
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


tracks <- readRDS("data/processed/Tracks_PA_w_dynSDM_10_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_1to4days_dynSDM_30_2018_2025_extract.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_seamounts.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")
tracks


tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced_woPATS.rds")


tracks |>
  dplyr::group_by(id) |>
  dplyr::summarise(Pre = sum(PA==1), Abs = sum(PA==0), ratio = Abs/Pre, .groups="drop") |>
  print(n = 200)

str(tracks)

model_dt <- tracks |>
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |> 
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year),
                depth = case_when(depth >= 0 ~ -5, TRUE ~ depth),
                id = as.factor(id),
                dist2000 = dist2000/1000
                ) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()


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




# Explore Data ------------------------------------------------------------


plots <- model_dt  |> dplyr::mutate(PA = as.factor(PA))

(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots |>  ggplot(aes(x = dist200, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


# Transformations

# Transformations

model_dt_trans <- model_dt |>
  dplyr::mutate(
    # depth = sign(depth) * abs(depth)^(1/3), # cube root\
    # depth = -depth,
    # depth = log1p(depth),
    # depth = log1p(abs(depth)),
    # slope = (slope)^(1/3),
    # slope = log1p(slope),
    # roughness = (roughness)^(1/3),
    # roughness = log1p(roughness),
    # uv = (uv)^(1/3),
    # uv = log1p(uv),
    chl = log(chl),
    # chl = log(chl)
    # mltost = log(mltost)
    # mltost = log(mltost)
  )


model_dt_trans |> 
  rstatix::get_summary_stats() |> 
  print(n=50)
model_dt |> 
  rstatix::get_summary_stats() |> 
  print(n=50)



plots <- model_dt_trans |> dplyr::mutate(PA = as.factor(PA))
plots

(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots |>  ggplot(aes(x = dist200, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
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

model_dt <- model_dt_trans

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
                  month) |> 
    as.data.frame(),
  pa = model_dt$PA)

swd




## lets try to tweak the package brt function



w <- ifelse(swd@pa == 1L, 1,
            sum(swd@pa == 1L) / sum(swd@pa == 0L))

w

# 
# # Make sure every row has a stable name, and key weights by those names
# rownames(swd@data) <- sprintf("r%07d", base::seq_len(nrow(swd@data)))
# names(w) <- rownames(swd@data)
# 
# # Store in an attribute (not a column!)
# attr(swd@data, "weights_by_row") <- w
# 
# 
# str(swd@data)
# 




#  override SDMtune's internal binding for this session with our modified trainBRT functuon from the helpers
assignInNamespace("trainBRT", patched_trainBRT, ns = "SDMtune")




# Using Cross Validation Method: ------------------------------------------

# get the blockCV spatial block folds

sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_1to4days.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")

sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_woPATS.rds")
sp_blocks

set.seed(28)
cv_brt0 <- SDMtune::train(method = "BRT", 
                        data = swd,
                        folds = sp_blocks)


cv_brt0

model <- cv_brt0

cat("Training AUC: ", SDMtune::auc(model))
cat("Testing AUC: ", SDMtune::auc(model, test = TRUE))
cat("Training TSS: ", SDMtune::tss(model))
cat("Testing TSS: ", SDMtune::tss(model, test = TRUE))





ROCplots <- lapply(seq_len(ncol(model@folds$test)), function(i){
  idx <- model@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = model@data@species,
    coords  = model@data@coords[idx, , drop = FALSE],
    data    = model@data@data  [idx, , drop = FALSE],
    pa      = model@data@pa    [idx]
  )
  SDMtune::plotROC(model@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})

ROCplots

## Variable Importance


vi_brt <- varImp(model, 
                    permut = 10)

vi_brt

SDMtune::plotVarImp(vi_brt)


jk1 <- SDMtune::doJk(model, 
                     metric = "tss", 
                     test = TRUE)
jk1

SDMtune::plotJk(jk1, 
                type = "train", 
                ref = SDMtune::tss(model))

SDMtune::plotJk(jk1, 
                type = "test", 
                ref = SDMtune::tss(model, test = TRUE))


jk2 <- SDMtune::doJk(model, 
                     metric = "auc", 
                     test = TRUE)
jk2

SDMtune::plotJk(jk2, 
                type = "train", 
                ref = SDMtune::auc(model))

SDMtune::plotJk(jk2, 
                type = "test", 
                ref = SDMtune::auc(model, test = TRUE))


m <- model
SDMtune::plotResponse(m, 
                      var = "thetao", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "red") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "uv", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  #scale_x_log10() +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "wz", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "mltost", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = median,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  # scale_x_log10() +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "chl", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = median,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  # scale_x_log10() +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "depth", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "slope", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

SDMtune::plotResponse(m, 
                      var = "dist2000", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()


SDMtune::plotResponse(m,
                      var = "month",
                      type = "logistic",
                      only_presence = TRUE,
                      marginal = TRUE,
                      rug = TRUE,
                      fun = mean,
                      color = "orange") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()



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


bg_dt <- model_dt |> dplyr::filter(PA == 0)
bg4cor <- SDMtune::SWD(
  species = "Bgs",
  coords  = bg_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data    = bg_dt |> dplyr::select(thetao, mltost, chl, uv, wz, 
                                   depth, 
                                   slope,
                                   # roughness, 
                                   dist2000, 
                                   month) |> as.data.frame(),
  pa      = bg_dt$PA)



SDMtune::plotCor(bg4cor, 
                 method = "spearman", 
                 cor_th = NULL)

SDMtune::corVar(bg4cor, 
                method = "spearman", 
                cor_th = 0.6)


cv_brt1 <- SDMtune::varSel(cv_brt0, 
                         metric = "auc", 
                         test = TRUE, 
                         bg4cor = bg4cor,
                         method = "spearman", 
                         cor_th = 0.6,
                         permut = 10)

cv_brt1



# Fine tune model  --------------------------------------------------------


SDMtune::getTunableArgs(cv_brt1)

h_brt <- list(
  n.trees           = c(2500L, 5000L, 10000L, 20000L),
  interaction.depth = c(1L, 2L, 3L),
  shrinkage         = c(0.01, 0.005, 0.001),
  bag.fraction      = c(0.5, 0.75),
  distribution      = "bernoulli"
)



expected_fits(pop = 20, 
              gen = 5, 
              keep_best = 0.4,
              keep_random = 0.2,
              rounding= "round")


cv_brt2 <- SDMtune::optimizeModel(cv_brt1, 
                                hypers = h_brt, 
                                metric = "auc",
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

# cv_brt2.2 <- SDMtune::optimizeModel(cv_brt1, 
#                                   hypers = h_brt, 
#                                   metric = "tss",
#                                   test = TRUE,
#                                   pop = 20,
#                                   gen = 5,
#                                   keep_best = 0.4,
#                                   keep_random = 0.2,
#                                   mutation_chance = 0.4,
#                                   interactive = TRUE,
#                                   progress = TRUE,
#                                   seed = 694
# )
# 
# 
# cv_brt2.2

cv_brt2

plot(cv_brt2,
     title = "Model tuning",
     interactive = TRUE)

cv_brt2@results
#ordered 
cv_brt2@results[order(-cv_brt2@results$test_AUC), ]
# cv_brt2.2@results[order(-cv_brt2@results$test_AUC), ]

# Index of the best model in the experiment
index <- terra::which.max(cv_brt2@results$test_AUC)
index


best_cv_brt2 <- cv_brt2@models[[index]]
best_cv_brt2
cv_brt2@results[index, ]

m <- best_cv_brt2



cat("Training AUC: ", SDMtune::auc(m))
cat("Testing AUC: ", SDMtune::auc(m, test = TRUE))
cat("Training TSS: ", SDMtune::tss(m))
cat("Testing TSS: ", SDMtune::tss(m, test = TRUE))


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

# h_refine <- list(
#   n.trees           = c(5000L, 10000L, 20000L),
#   interaction.depth = c(2L, 3L, 4L, 5L),
#   shrinkage         = c(0.01, 0.001, 0.005),
#   bag.fraction      = c(0.5, 0.75),
#   distribution      = "bernoulli"
# )
# 
# 
# 
# 
# cv_brt3 <- SDMtune::gridSearch(
#   model  = cv_brt1,         
#   hypers = h_refine,
#   metric = "auc",
#   interactive = TRUE,
#   progress = TRUE
# )
# 
# 
# 
# cv_brt3
# 
# cv_brt3@results
# #ordered 
# cv_brt3@results[order(-cv_brt3@results$test_AUC), ]
# 
# # Index of the best model in the experiment
# index <- terra::which.max(cv_brt3@results$test_AUC)
# index
# 
# 
# best_cv_brt3 <- cv_brt3@models[[index]]
# best_cv_brt3
# cv_brt3@results[index, ]
# 
# 
# 
# m <- best_cv_brt2
# m
# 
# cat("Training TSS: ", SDMtune::tss(m))
# cat("Testing TSS: ", SDMtune::tss(m, test = TRUE))
# cat("Training AUC: ", SDMtune::auc(m))
# cat("Testing AUC: ", SDMtune::auc(m, test = TRUE))
# 
# 
# # check if refine made it better:
# cat("Testing AUC: ", SDMtune::auc(cv_brt1, test = TRUE))
# cat("Testing AUC: ", SDMtune::auc(best_cv_brt3, test = TRUE))
# cat("Testing TSS: ", SDMtune::tss(cv_brt1, test = TRUE))
# cat("Testing TSS: ", SDMtune::tss(best_cv_brt3, test = TRUE))
# 
# cat("Training AUC before tuning: ", SDMtune::auc(cv_brt1))
# cat("Training AUC after tuning: ", SDMtune::auc(best_cv_brt3))
# cat("Training TSS before tuning: ", SDMtune::tss(cv_brt1))
# cat("Training TSS after tuning: ", SDMtune::tss(best_cv_brt3))
# 
# 
# 
# ROCplots <- lapply(seq_len(ncol(m@folds$test)), function(i){
#   idx <- m@folds$test[, i]
#   test_i <- SDMtune::SWD(
#     species = m@data@species,
#     coords  = m@data@coords[idx, , drop = FALSE],
#     data    = m@data@data  [idx, , drop = FALSE],
#     pa      = m@data@pa    [idx]
#   )
#   SDMtune::plotROC(m@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
# })
# 
# ROCplots[[1]]
# 
# 
# ROCplots_combined <- wrap_plots(ROCplots) +
#   plot_layout(ncol = 3, nrow = 2, axes = "collect") &
#   ggplot2::theme(
#     legend.position = c(0.95, 0.03),  # x, y coordinates (right–bottom)
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
#     legend.key.size = unit(10, "pt"),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 10)
#   )
# 
# ROCplots_combined
# 




# Removing variables with low importance  --------------------------------

SDMtune::varImp(best_cv_brt2, 
                permut = 10)




cv_brt4 <- SDMtune::reduceVar(best_cv_brt2, 
                            th = 1, 
                            metric = "auc", 
                            test = TRUE, 
                            permut = 10, 
                            use_jk = TRUE)

# cv_brt4 <- SDMtune::reduceVar(best_cv_brt2,
#                               th = 2, 
#                               metric = "tss", 
#                               test = TRUE, 
#                               permut = 10, 
#                               use_jk = TRUE)

cat("Testing TSS before: ", SDMtune::tss(best_cv_brt2, test = TRUE))
cat("Testing TSS after: ", SDMtune::tss(cv_brt4, test = TRUE))
cat("Testing AUC before: ", SDMtune::auc(best_cv_brt2, test = TRUE))
cat("Testing AUC after: ", SDMtune::auc(cv_brt4, test = TRUE))

cv_brt4


vi_brt <- SDMtune::varImp(cv_brt4, 
                             permut = 10)
vi_brt
SDMtune::plotVarImp(vi_brt)


jk4 <- SDMtune::doJk(cv_brt4, 
                     metric = "tss", 
                     test = TRUE)
jk4

SDMtune::plotJk(jk4, 
                type = "train", 
                ref = SDMtune::auc(cv_brt4))

SDMtune::plotJk(jk4, 
                type = "test", 
                ref = SDMtune::auc(cv_brt4, test = TRUE))


m <- cv_brt4
P_sst <- SDMtune::plotResponse(m, 
                      var = "thetao", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

P_chl <- SDMtune::plotResponse(m, 
                              var = "chl", 
                              type = "cloglog", 
                              only_presence = TRUE, 
                              marginal = TRUE, 
                              rug = TRUE,
                              fun = mean,
                              color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  #scale_x_log10() +
  theme_bw()

P_uv <- SDMtune::plotResponse(m, 
                      var = "uv", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  #scale_x_log10() +
  theme_bw()

P_wz <- SDMtune::plotResponse(m, 
                      var = "wz", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

P_mld <- SDMtune::plotResponse(m, 
                      var = "mltost", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = median,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  # scale_x_log10() +
  theme_bw()

P_mld

P_depth <- SDMtune::plotResponse(m, 
                      var = "depth", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(-0.05, 1)) +
  theme_bw()

P_slope <- SDMtune::plotResponse(m, 
                      var = "slope", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

P_dist <- SDMtune::plotResponse(m, 
                      var = "dist2000", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

P_month <- SDMtune::plotResponse(m, 
                      var = "month", 
                      type = "cloglog", 
                      only_presence = TRUE, 
                      marginal = TRUE, 
                      rug = TRUE,
                      fun = mean,
                      color = "orange") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

(P_depth + P_dist + P_sst + P_chl + P_uv + P_mld + P_wz + P_slope + P_month) +
  patchwork::plot_layout(nrow = 3, guides = "collect", axes = "collect")


# get final model 
set.seed(25)
# final_m <- SDMtune::combineCV(best_cv_brt3)
final_m <- SDMtune::combineCV(cv_brt4)


# create extrnal validation data from sightings 
sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full.rds")
str(sight)



mapview::mapview(sight |> dplyr::select(-Date))


val_dt <- sight |>
  # dplyr::filter(PA == 1) |> 
  dplyr::filter(lon >140 & lon < 170) |>
  # dplyr::filter(lat >-35 & lat < -5) |>
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year)) |> 
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |>
  sf::st_drop_geometry() |> 
  as.data.frame()

str(val_dt)


val_dt <- val_dt |> 
  dplyr::mutate(
    # slope = log1p(slope),
    # roughness = log1p(roughness),
    chl = log(chl),
    dist2000 = dist2000/1000,
    PA = as.factor(PA))


plots <- val_dt

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
SDMtune::auc(final_m)
SDMtune::tss(final_m)

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






# Predictions -------------------------------------------------------------
relevant_vars <- names(final_m@data@data)
predictors <- predictors_mean[[relevant_vars]]


e <- terra::ext(c(140, 170, -40, 0))


map <- SDMtune::predict(final_m,
                        data = predictors,
                        type = "cloglog",
                        const = (data.frame(month = "0")),
                        extent = e)

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
                     folder = "models/objects/Tracks_BRT_tuned_mp", 
                     test = NULL, 
                     response_curves = TRUE, 
                     only_presence = TRUE, 
                     jk = FALSE, 
                     env = predictors_mean)



saveRDS(cv_brt4, "models/objects/Tracks_BRT_CV_Models_SDMtune_mp.rds")
saveRDS(final_m, "models/objects/Tracks_BRT_final_Model_SDMtune_mp.rds")


# Boyce index -------------------------------------------------------------

oof_boyce_from_cv <- function(cv_model, type = "cloglog", use_test_background = TRUE,
                              boyce_res = 100L, boyce_nclass = 0L, plot = TRUE) {
  k <- base::ncol(cv_model@folds$test)
  obs_all <- base::numeric(0)
  fit_all <- base::numeric(0)
  
  for (i in base::seq_len(k)) {
    idx <- cv_model@folds$test[, i]
    test_i <- SDMtune::SWD(
      species = cv_model@data@species,
      coords  = cv_model@data@coords[idx, , drop = FALSE],
      data    = cv_model@data@data  [idx, , drop = FALSE],
      pa      = cv_model@data@pa    [idx]
    )
    preds <- SDMtune::predict(cv_model@models[[i]], data = test_i, type = type)
    obs_all <- base::c(obs_all, preds[test_i@pa == 1L])
    fit_all <- base::c(fit_all, if (use_test_background) preds[test_i@pa == 0L] else preds)
  }
  
  obs_all <- obs_all[stats::complete.cases(obs_all)]
  fit_all <- fit_all[stats::complete.cases(fit_all)]
  
  boyce <- ecospat::ecospat.boyce(
    fit      = fit_all,
    obs      = obs_all,
    window.w = "default",
    nclass   = boyce_nclass,
    res      = boyce_res,
    PEplot   = plot
  )
  
  df <- base::data.frame(HS = boyce$HS, PE = boyce$F.ratio) |>
    dplyr::filter(stats::complete.cases(HS, PE))
  
  ct <- stats::cor.test(df$HS, df$PE, method = "spearman", exact = FALSE)
  
  list(
    boyce      = boyce,
    df         = df,                 # for ggplot
    rho        = unname(ct$estimate),
    p_value    = ct$p.value,
    n_points   = base::nrow(df)
  )
}



# --- run on best CV trackss model ---
res_boyce <- oof_boyce_from_cv(cv_brt4, type = "cloglog", use_test_background = TRUE)
res_boyce$p_value
res_boyce$boyce
res_boyce$rho
res_boyce$df



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
avail <- as.Date(names(input_list))  # "YYYY-MM-01"
idx <- match(want, avail)

if (anyNA(idx)) {
  warning("Missing months in raster list: ",
          paste(as.character(want[is.na(idx)]), collapse = ", "))
}

monthly_predictors   <- input_list[idx[!is.na(idx)]]
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



SDMtune::plotPred(monthly_predictions_stack[[68]],
                  lt = "Habitat\nsuitability",
                  colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


ths <- SDMtune::thresholds(final_m, 
                           type = "cloglog")

ths

SDMtune::plotPA(monthly_predictions_stack[[66]],
                th = ths[3, 2])


writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_2019_2025_rev_mp.tif", overwrite = TRUE)



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


plot(monthly_means_stack, range = c(0, 1))
plot(monthly_means_stack[[12]], range = c(0, 1))



# Save for ensemble modelling
writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_means_rev_mp.tif", overwrite = TRUE)


monthly_means_stack <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_monthly_means_rev_mp.tif")

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

plot(seasonal_means_stack_2, range = c(0, 1))

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
                                                      #col.regions =  rev(topo.colors(200)),
                                                      # col.regions = cm_ocean_palette, 
                                                      # col.regions =  viridisLite::inferno(200),
                                                      col.regions = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100),
                                                      colorkey = list(space = "right", 
                                                                      length = 0.5, 
                                                                      height = 0.75,
                                                                      labels = list(
                                                                        at = c(0, high),  # Positions for "Low" and "High"
                                                                        labels = c("Low", "High"))))  # Labels to use

Mean_Season_monsoon_SDM_plot




# Save for ensemble modelling
writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp.tif", overwrite = TRUE)




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



writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp.tif", overwrite = TRUE)
