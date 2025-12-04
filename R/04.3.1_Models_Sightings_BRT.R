#_____________________________________________________________________________
#                        Models: Sightings - Boosted Regression Tree
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
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------
# 
# 
# sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full_processed.rds")

sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed_SUPPS_MODEL.rds")



sight


model_dt <- sight |>
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year),
                dist2000 = dist2000/1000,
                chl = log(chl)) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()

str(model_dt)

model_dt |> dplyr::summarise(start = min(Date),
                             end = max(Date))




swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = model_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = model_dt |>
    dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, roughness, dist2000, month) |> as.data.frame(),
  pa = model_dt$PA)

swd

swd@data



## lets try to tweak the package brt function

# make class-balancing weights (pres = 1; abs scaled so totals match)
# swd@data$case_w <- ifelse(swd@pa == 1L, 1, sum(swd@pa == 1L) / sum(swd@pa == 0L))


# w <- ifelse(swd@pa == 1L,
#             1,
#             sum(swd@pa == 1L) / sum(swd@pa == 0L))
# 
# w
# 
# 
# # Make sure every row has a stable name, and key weights by those names
# base::rownames(swd@data) <- sprintf("r%07d", base::seq_len(base::nrow(swd@data)))
# base::names(w) <- base::rownames(swd@data)
# 
# # Store in an attribute (not a column!)
# base::attr(swd@data, "weights_by_row") <- w
# 
# 
# str(swd@data)





#  override SDMtune's internal binding for this session with our modified trainBRT functuon from the helpers
assignInNamespace("trainBRT", patched_trainBRT, ns = "SDMtune")








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
blocks <- readRDS("data/processed/BlockCV_spatial_folds_sightings_groups.rds")

blocks <- readRDS("data/processed/BlockCV_spatial_folds_sightings_SUPPS_MODEL.rds")

set.seed(25)
cv_brt_s0 <- SDMtune::train(method = "BRT", 
                        data = swd,
                        folds = blocks)


cv_brt_s0

cat("Training AUC: ", SDMtune::auc(cv_brt_s0))
cat("Testing AUC: ", SDMtune::auc(cv_brt_s0, test = TRUE))
cat("Training TSS: ", SDMtune::tss(cv_brt_s0))
cat("Testing TSS: ", SDMtune::tss(cv_brt_s0, test = TRUE))


m <- cv_brt_s0
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

## Variable Importance

vi_brt <- varImp(cv_brt_s0, 
                 permut = 10)

vi_brt

SDMtune::plotVarImp(vi_brt)


jk1 <- SDMtune::doJk(cv_brt0, 
                     metric = "tss", 
                     test = TRUE)
jk1

SDMtune::plotJk(jk1, 
                type = "train", 
                ref = auc(cv_brt0))

SDMtune::plotJk(jk1, 
                type = "test", 
                ref = auc(cv_brt0, test = TRUE))


jk2 <- SDMtune::doJk(cv_brt0, 
                     metric = "auc", 
                     test = TRUE)
jk2

SDMtune::plotJk(jk2, 
                type = "train", 
                ref = auc(cv_brt0))

SDMtune::plotJk(jk2, 
                type = "test", 
                ref = auc(cv_brt0, test = TRUE))



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
  scale_x_log10() +
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
  scale_x_log10() +
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
predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1.tif")
names(predictors_mean)





bg_dt <- model_dt |> dplyr::filter(PA == 0)
bg4cor <- SDMtune::SWD(
  species = "Bgs",
  coords  = bg_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data    = bg_dt |> dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, roughness, dist2000, month) |> as.data.frame(),
  pa      = bg_dt$PA)



SDMtune::plotCor(bg4cor, 
                 method = "spearman", 
                 cor_th = NULL)

SDMtune::corVar(bg4cor, 
                method = "spearman", 
                cor_th = 0.7)


cv_brt_s1 <- SDMtune::varSel(cv_brt_s0, 
                         metric = "tss", 
                         test = TRUE, 
                         bg4cor = bg4cor,
                         method = "spearman", 
                         cor_th = 0.7,
                         permut = 10)

cv_brt_s1


vi_brt <- SDMtune::varImp(cv_brt_s1, 
                             permut = 10)
vi_brt
SDMtune::plotVarImp(vi_brt)


jk1 <- SDMtune::doJk(cv_brt1, 
                     metric = "tss", 
                     test = TRUE)
jk1

SDMtune::plotJk(jk1, 
                type = "train", 
                ref = auc(cv_brt1))

SDMtune::plotJk(jk1, 
                type = "test", 
                ref = auc(cv_brt1, test = TRUE))


m <- cv_brt1
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
  scale_x_log10() +
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






# Fine tune model  --------------------------------------------------------


SDMtune::getTunableArgs(cv_brt_s1)

h_brt <- list(
  n.trees           = c(2500L, 5000L, 7500L, 10000L, 15000L),
  interaction.depth = c(1L, 2L, 3L, 5L),
  shrinkage         = c(0.01, 0.005, 0.001),
  bag.fraction      = c(0.5, 0.75),
  distribution      = "bernoulli"
)






cv_brt_s2 <- SDMtune::optimizeModel(cv_brt_s1, 
                                hypers = h_brt, 
                                metric = "tss",
                                pop = 20,
                                gen = 5,
                                keep_best = 0.4,
                                keep_random = 0.2,
                                mutation_chance = 0.4,
                                interactive = TRUE,
                                progress = TRUE,
                                seed = 799
                                )


cv_brt_s2

summary(cv_brt_s2)
slotNames(cv_brt_s2)

cv_brt_s2@results
#ordered 
cv_brt_s2@results[order(-cv_brt_s2@results$test_TSS), ]

# Index of the best model in the experiment
# index <- terra::which.max(cv_brt_s2@results$test_TSS)
index <- 3


best_cv_brt_s2 <- cv_brt_s2@models[[index]]
best_cv_brt_s2
cv_brt_s2@results[index, ]



best_cv_brt_s2
cat("Training AUC: ", SDMtune::auc(best_cv_brt_s2))
cat("Testing AUC: ", SDMtune::auc(best_cv_brt_s2, test = TRUE))
cat("Training TSS: ", SDMtune::tss(best_cv_brt_s2))
cat("Testing TSS: ", SDMtune::tss(best_cv_brt_s2, test = TRUE))

m <- best_cv_brt_s2
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




# Refine using gridSearch -------------------------------------------------

h_refine <- list(
  n.trees           = c(8000L, 9000L, 10000L, 12000L),  #around 7k
  interaction.depth = c(2L, 3L),                             # 2 was best; 3 as a safety check
  shrinkage         = c(0.0005, 0.00075, 0.001, 0.0015, 0.002), # 0.001 was best, check a little up and down 
  bag.fraction      = c(0.75, 0.5), # 0.75 most, but check 0.5 anyway
  distribution      = "bernoulli"
)


cv_brt3 <- SDMtune::gridSearch(
  model  = cv_brt1,         
  hypers = h_refine,
  metric = "tss",
  interactive = TRUE,
  progress = TRUE
)

cv_brt3

cv_brt3@results
#ordered 
cv_brt3@results[order(-cv_brt3@results$test_TSS), ]

# Index of the best model in the experiment
index <- terra::which.max(cv_brt3@results$test_TSS)
index


best_cv_brt3 <- cv_brt3@models[[index]]
best_cv_brt3
cv_brt3@results[index, ]


cat("Training TSS: ", SDMtune::tss(best_cv_brt3))
cat("Testing TSS: ", SDMtune::tss(best_cv_brt3, test = TRUE))
cat("Training AUC: ", SDMtune::auc(best_cv_brt3))
cat("Testing AUC: ", SDMtune::auc(best_cv_brt3, test = TRUE))


# check if refine made it better:
cat("Testing AUC before refining: ", SDMtune::auc(best_cv_brt2, test = TRUE))
cat("Testing AUC after refining: ", SDMtune::auc(best_cv_brt3, test = TRUE))
cat("Testing TSS before refining: ", SDMtune::tss(best_cv_brt2, test = TRUE))
cat("Testing TSS after refining: ", SDMtune::tss(best_cv_brt3, test = TRUE))

# --> based on this we continue with the best model before refining


# test the startin model against best model after tuning:
cat("Testing TSS before tuning: ", SDMtune::tss(cv_brt0, test = TRUE))
cat("Testing TSS after tuning: ", SDMtune::tss(best_cv_brt2, test = TRUE))

cat("Testing AUC before tuning: ", SDMtune::auc(cv_brt0, test = TRUE))
cat("Testing AUC after tuning: ", SDMtune::auc(best_cv_brt2, test = TRUE))



m <- best_cv_m2
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

SDMtune::varImp(best_cv_brt_s2, 
                permut = 10)

best_cv_brt_s2
cv_brt_s3 <- SDMtune::reduceVar(best_cv_brt_s2, 
                            th = 1, 
                            metric = "tss", 
                            test = TRUE, 
                            permut = 10, 
                            use_jk = TRUE)

cat("Testing TSS before: ", SDMtune::tss(best_cv_brt2_sight, test = TRUE))
cat("Testing TSS after: ", SDMtune::tss(cv_brt3_sight, test = TRUE))
cat("Testing AUC before: ", SDMtune::auc(best_cv_brt2_sight, test = TRUE))
cat("Testing AUC after: ", SDMtune::auc(cv_brt3_sight, test = TRUE))

cv_brt3

## --> nothign removed


# Get the FINAL MODEL -----------------------------------------------------

set.seed(25)
final_brt_sight <- SDMtune::combineCV(best_cv_brt_s2)
final_brt_sight

SDMtune::plotROC(final_brt_sight)
final_brt_sight@model




m <- best_cv_brt_s2
m
vi_brt <- SDMtune::varImp(m, 
                          permut = 10)
vi_brt
SDMtune::plotVarImp(vi_brt)


P_sst <- SDMtune::plotResponse(m, 
                               var = "thetao", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "steelblue") +
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
                               color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("chl (log(mg m"^{-3}*"))")) +
  # scale_x_log10() +
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
                              color = "steelblue") +
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
                              color = "steelblue") +
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
                               color = "steelblue") +
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
                                 color = "steelblue") +
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
                                 color = "steelblue") +
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
                                color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = "dist2000 (km)") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# P_dist200 <- SDMtune::plotResponse(m,
#                                 var = "dist200",
#                                 type = "cloglog",
#                                 only_presence = TRUE,
#                                 marginal = TRUE,
#                                 rug = TRUE,
#                                 fun = mean,
#                                 color = "steelblue") +
#   ggplot2::scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "dist200 (km)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 10),
#         axis.text = element_text(size = 8))

# P_rough <- SDMtune::plotResponse(m,
#                                 var = "roughness",
#                                 type = "cloglog",
#                                 only_presence = TRUE,
#                                 marginal = TRUE,
#                                 rug = TRUE,
#                                 fun = mean,
#                                 color = "steelblue") +
#   ggplot2::scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "roughness") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 10),
#         axis.text = element_text(size = 8))

# P_seam <- SDMtune::plotResponse(m,
#                                 var = "dist_seamount",
#                                 type = "cloglog",
#                                 only_presence = TRUE,
#                                 marginal = TRUE,
#                                 rug = TRUE,
#                                 fun = mean,
#                                 color = "steelblue") +
#   ggplot2::scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "seamount (m)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 10),
#         axis.text = element_text(size = 8))


P_month <- SDMtune::plotResponse(m, 
                                 var = "month", 
                                 type = "cloglog", 
                                 only_presence = TRUE, 
                                 marginal = TRUE, 
                                 rug = TRUE,
                                 fun = mean,
                                 color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

P_month



marginal_plots <- (
  P_depth +
    P_dist +
    # P_dist200 +
    # P_seam +
    P_sst + 
    P_chl + 
    P_uv +
    P_mld + 
    P_wz + 
    P_slope +
    # P_rough +
    P_month) +
  patchwork::plot_layout(ncol = 3, guides = "collect", axes = "collect") &
  ggplot2::labs(y = "Rel. Habitat Suitability")

marginal_plots


# get final model 


map <- SDMtune::predict(final_brt_sight,
                        data = predictors_mean,
                        type = "cloglog",
                        const = (data.frame(month = "0")))

SDMtune::plotPred(map,
                  lt = "Habitat\nsuitability",
                  colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

map

ths <- SDMtune::thresholds(final_brt, 
                           type = "cloglog")

ths

SDMtune::plotPA(map, 
                th = ths[3, 2])

SDMtune::plotPA(map, 
                th = 0.7)





SDMtune::modelReport(final_brt_sight, 
                     type = "cloglog", 
                     folder = "models/objects/Sightings_BRT_tuned", 
                     test = NULL, 
                     response_curves = TRUE, 
                     only_presence = TRUE, 
                     jk = FALSE, 
                     env = predictors_mean)



saveRDS(best_cv_brt_s2, "models/objects/Sightings_BRT_CV_Models_SDMtune_SUPPS_MODEL.rds")
saveRDS(final_brt_sight, "models/objects/Sightings_BRT_final_Model_SDMtune_SUPPS_MODEL.rds")





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



# --- run it on best CV sightings model ---
res_boyce <- oof_boyce_from_cv(best_cv_brt_s2, type = "cloglog", use_test_background = TRUE)
res_boyce$p_value
res_boyce$boyce
res_boyce$rho
res_boyce$df



set.seed(1)

# 1) data
x  <- res_boyce$boyce
df <- data.frame(HS = x$HS, PE = x$F.ratio)
df <- df[is.finite(df$HS) & is.finite(df$PE), ]

ct <- suppressWarnings(cor.test(df_plot$HS, df_plot$PE,
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





# MOnthly predictions -----------------------------------------------------


relevant_vars <- names(final_brt@data@data)
relevant_vars

input_list <- monthly_stacks_lst_0.1_trans



# 1) Get start/end from your model data
rng <- range(as.Date(model_dt$Date), na.rm = TRUE)
start_mon <- as.Date(format(rng[1], "%Y-%m-01"))
end_mon   <- as.Date(format(rng[2], "%Y-%m-01"))
want <- seq(start_mon, end_mon, by = "1 month")

# 2) Align your monthly rasters to that range
avail <- as.Date(names(monthly_stacks_lst_0.1))  # "YYYY-MM-01"
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
  model = final_brt_sight,
  monthly_list = monthly_predictors,
  dates = monthly_dates,  # already "YYYY-MM-01"
  type = "cloglog",
  extent = e,
  verbose = TRUE
)
tictoc::toc()


monthly_predictions_stack
names(monthly_predictions_stack)
plot(monthly_predictions_stack[[164:179]])



SDMtune::plotPred(monthly_predictions_stack[[179]],
                  lt = "Habitat\nsuitability",
                  colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


ths <- SDMtune::thresholds(final_brt_sight, 
                           type = "cloglog")

ths

SDMtune::plotPA(monthly_predictions_stack[[179]],
                th = ths[3, 2])


writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Sightings_BRT_monthly_2010_2025_rev.tif", overwrite = TRUE)



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
writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Sightings_BRT_monthly_means_rev.tif", overwrite = TRUE)




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

# seasonal_means_stack_2 <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA.tif")

plot(seasonal_means_stack_2, range = c(0, 1))
plot(seasonal_means_stack_2, range = c(0, 1), xlim = c(140, 170), ylim = c(-40, 0), axes =FALSE, legend =FALSE)

plot(seasonal_means_stack_2, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))





# writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA.tif", overwrite = TRUE)

# writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA_RANDOM_FOLDS_CV.tif", overwrite = TRUE)




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
plot(seasonal_means_stack_4[[4]], range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))



# writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp_crwPA.tif", overwrite = TRUE)




# writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons4_means_rev_mp_crwPA_RANDOM_FOLDS_CV.tif", overwrite = TRUE)



