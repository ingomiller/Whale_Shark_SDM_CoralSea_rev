#_____________________________________________________________________________
#                        Models: Sightings - Maxent
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


sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
sight


model_dt <- sight |>
  dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
                year  = factor(year)) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()

str(model_dt)

model_dt |> dplyr::summarise(start = min(Date),
                             end = max(Date))




swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = model_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = model_dt |>
    dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, dist2000, month) |> as.data.frame(),
  pa = model_dt$PA)

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
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds.rds")
set.seed(25)
cv_m0 <- SDMtune::train(method = "Maxent", 
                     data = swd,
                     folds = sp_blocks)


cv_m0

cat("Training AUC: ", auc(cv_m0))
cat("Testing AUC: ", auc(cv_m0, test = TRUE))
cat("Training TSS: ", tss(cv_m0))
cat("Testing TSS: ", tss(cv_m0, test = TRUE))


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

ROCplots

## Variable Importance


vi_maxent <- varImp(cv_m0, 
                    permut = 10)

vi_maxent

SDMtune::plotVarImp(vi_maxent)


jk1 <- SDMtune::doJk(cv_m0, 
           metric = "tss", 
           test = TRUE)
jk1

SDMtune::plotJk(jk1, 
       type = "train", 
       ref = auc(cv_m0))

SDMtune::plotJk(jk1, 
       type = "test", 
       ref = auc(cv_m0, test = TRUE))


jk2 <- SDMtune::doJk(cv_m0, 
            metric = "auc", 
            test = TRUE)
jk2

SDMtune::plotJk(jk2, 
       type = "train", 
       ref = auc(cv_m0))

SDMtune::plotJk(jk2, 
       type = "test", 
       ref = auc(cv_m0, test = TRUE))


m <- cv_m0
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
  data    = bg_dt |> dplyr::select(thetao, mltost, chl, uv, wz, depth, slope, dist2000, month) |> as.data.frame(),
  pa      = bg_dt$PA)



SDMtune::plotCor(bg4cor, 
        method = "spearman", 
        cor_th = NULL)

SDMtune::corVar(bg4cor, 
       method = "spearman", 
       cor_th = 0.7)


cv_m1 <- SDMtune::varSel(cv_m0, 
                         metric = "tss", 
                         test = TRUE, 
                         bg4cor = bg4cor,
                         method = "spearman", 
                         cor_th = 0.7,
                         permut = 10)

cv_m1


vi_maxent <- SDMtune::varImp(cv_m1, 
                    permut = 10)
vi_maxent
SDMtune::plotVarImp(vi_maxent)


jk1 <- SDMtune::doJk(cv_m1, 
            metric = "tss", 
            test = TRUE)
jk1

SDMtune::plotJk(jk1, 
       type = "train", 
       ref = auc(cv_m1))

SDMtune::plotJk(jk1, 
       type = "test", 
       ref = auc(cv_m1, test = TRUE))


m <- cv_m1
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


SDMtune::getTunableArgs(cv_m1)

h <- list(reg = seq(0.1, 3, 0.1), 
          fc = c("lq", "lh", "lqp", "lqph", "lqpht"),
          iter = c(500, 1000, 2500, 5000, 10000, 20000, 50000)
)

h


expected_fits(pop = 30, 
              gen = 10, 
              keep_best = 0.4,
              keep_random = 0.2,
              rounding="round")# also 60 with these params


cv_m2 <- SDMtune::optimizeModel(cv_m1, 
                    hypers = h, 
                    metric = "tss",
                    test = TRUE,
                    pop = 30,
                    gen = 10,
                    keep_best = 0.4,
                    keep_random = 0.2,
                    mutation_chance = 0.4,
                    interactive = TRUE,
                    progress = TRUE,
                    seed = 699
)


cv_m2

summary(cv_m2)
slotNames(cv_m2)

cv_m2@results
#ordered 
cv_m2@results[order(-cv_m2@results$test_TSS), ]

# Index of the best model in the experiment
index <- terra::which.max(cv_m2@results$test_TSS)
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

ROCplots



# Refine using gridSearch -------------------------------------------------

h_refine <- list(
  reg = seq(2.4, 3.2, by = 0.1), 
  fc   = c("lqh", "lqph", "lh"), 
  iter = c(1000, 2000, 5000)
)




cv_m3 <- SDMtune::gridSearch(
  model  = cv_m1,         # same SWD + same spatial folds
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
index <- terra::which.max(cv_m3@results$test_TSS)
index


best_cv_m3 <- cv_m3@models[[index]]
best_cv_m3
cv_m3@results[index, ]


cat("Training TSS: ", SDMtune::tss(best_cv_m3))
cat("Testing TSS: ", SDMtune::tss(best_cv_m3, test = TRUE))
cat("Training AUC: ", SDMtune::auc(best_cv_m3))
cat("Testing AUC: ", SDMtune::auc(best_cv_m3, test = TRUE))


# check if refine made it better:
cat("Testing AUC: ", SDMtune::auc(best_cv_m2, test = TRUE))
cat("Testing AUC: ", SDMtune::auc(best_cv_m3, test = TRUE))
cat("Testing TSS: ", SDMtune::tss(best_cv_m2, test = TRUE))
cat("Testing TSS: ", SDMtune::tss(best_cv_m3, test = TRUE))


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

SDMtune::varImp(best_cv_m3 , 
       permut = 10)


cv_m4 <- SDMtune::reduceVar(best_cv_m3, 
                   th = 1, 
                   metric = "tss", 
                   test = TRUE, 
                   permut = 10, 
                   use_jk = TRUE)

cat("Testing TSS before: ", SDMtune::tss(best_cv_m3, test = TRUE))
cat("Testing TSS after: ", SDMtune::tss(cv_m4, test = TRUE))
cat("Testing AUC before: ", SDMtune::auc(best_cv_m3, test = TRUE))
cat("Testing AUC after: ", SDMtune::auc(cv_m4, test = TRUE))

cv_m4


vi_maxent <- SDMtune::varImp(cv_m4, 
                             permut = 10)
vi_maxent
SDMtune::plotVarImp(vi_maxent)


jk4 <- SDMtune::doJk(cv_m4, 
                     metric = "tss", 
                     test = TRUE)
jk4

SDMtune::plotJk(jk4, 
                type = "train", 
                ref = auc(cv_m1))

SDMtune::plotJk(jk4, 
                type = "test", 
                ref = auc(cv_m1, test = TRUE))


m <- cv_m4
SDMtune::plotResponse(cv_m4, 
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
  ggplot2::scale_y_continuous(limits = c(-0.05, 1)) +
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





# get final model 
set.seed(25)
final_m <- SDMtune::combineCV(cv_m4)

SDMtune::plotROC(final_m)
final_m@model





map <- SDMtune::predict(final_m,
               data = predictors_mean,
               type = "cloglog",
               const = (data.frame(month = "0")))

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
            folder = "models/objects/Sightings_Maxent_tuned", 
            test = NULL, 
            response_curves = TRUE, 
            only_presence = TRUE, 
            jk = FALSE, 
            env = predictors_mean)



saveRDS(cv_m1, "models/objects/Sightings_Maxent_CV_Models_SDMtune.rds")
saveRDS(cv_m4, "models/objects/Sightings_Maxent_final_Model_SDMtune.rds")


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



# --- run it on your best CV sightings model ---
res_boyce <- oof_boyce_from_cv(cv_m4, type = "cloglog", use_test_background = TRUE)
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


relevant_vars <- names(final_m@data@data)
relevant_vars

input_list <- monthly_stacks_lst_0.1



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

monthly_predictors   <- monthly_stacks_lst_0.1[idx[!is.na(idx)]]
monthly_dates  <- avail[idx[!is.na(idx)]]


e <- terra::ext(c(136, 170, -40, -5))
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
plot(monthly_predictions_stack[[164:179]])



SDMtune::plotPred(monthly_predictions_stack[[179]],
         lt = "Habitat\nsuitability",
         colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


ths <- SDMtune::thresholds(final_m, 
                  type = "cloglog")

ths

SDMtune::plotPA(monthly_predictions_stack[[179]],
       th = ths[3, 2])


writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Sightings_Maxent_monthly_2010_2025_rev.tif", overwrite = TRUE)



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



# Save for ensemble modelling
writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Sightings_Maxent_monthly_means_rev.tif", overwrite = TRUE)






