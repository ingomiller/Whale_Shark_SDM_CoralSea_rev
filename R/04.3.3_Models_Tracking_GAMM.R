#_____________________________________________________________________________
#                        Models: Tracking - GAMM
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
library(mgcv)
library(patchwork)
library(MuMIn)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


# tracks <- readRDS("data/processed/Tracks_PA_w_dynSDM_10_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_3days_dynSDM_10_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")

# this si the one that worked well before 
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")
# 
# str(tracks)
# 
# mapview::mapview(tracks |>  dplyr::select(-Date))
# 
# 
# tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced.rds")

tracks <- readRDS( "data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2019_2025_extract_final_processed_monthsreduced_4.rds")


str(tracks)

mapview::mapview(tracks |>  dplyr::select(-Date))


tracks |> 
  dplyr::group_by(PA) |> 
  dplyr::summarise(N = n())

model_dt <- tracks |>
  dplyr::mutate(id = as.factor(id),
                month_fc = as.factor(month)) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()



model_dt |> rstatix::get_summary_stats(depth, slope, roughness, thetao, uv, chl, dist2000, 
                                       # dist200, 
                                       wz, mltost, dist_seamount, dist_knoll )

plots <- model_dt |> dplyr::mutate(PA = as.factor(PA))


(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    # plots |>  ggplot(aes(x = dist200, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_seamount, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = dist_knoll, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = sst_slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)



model_dt |>
  dplyr::group_by(month) |>
  dplyr::summarise(
    presences = sum(PA == 1, na.rm = TRUE),
    absences  = sum(PA == 0, na.rm = TRUE),
    total     = dplyr::n()
  )

model_dt |>
  dplyr::summarise(
    presences = sum(PA == 1, na.rm = TRUE),
    absences  = sum(PA == 0, na.rm = TRUE),
    total     = dplyr::n()
  )



# Check Colinearity of repdictor variables


model_dt |>
  dplyr::filter(PA == 0) |> 
  dplyr::select(depth, slope, roughness, dist2000, thetao, chl, uv, wz, mltost, dist_seamount, dist_knoll, sst_slope) |>
  psych::pairs.panels(method = "spearman", cor = TRUE, scale = FALSE, stars = FALSE)



# 1) Do presences and buffered absences differ at all (univariate)?
library(pROC)
vars <- c("depth","slope","dist2000","uv","thetao","wz","mltost","chl","dist_seamount","dist_knoll")
sapply(vars, \(v) roc(model_dt$PA, model_dt[[v]], quiet=TRUE)$auc)

# 2) Fit *single-term* GAMs to see if *any* variable has signal
library(mgcv)
#w <- with(model_dt, ifelse(PA==1, 1, sum(PA==1)/sum(PA==0)))  # global class-balance
sig_check <- lapply(vars, \(v){
  f <- as.formula(paste0("PA ~ s(", v, ", bs='cs', k=6)"))
  m <- gam(f, data=model_dt, family=binomial(), method="REML")
  c(var=v, edf=summary(m)$s.table[1,"edf"], p=summary(m)$s.table[1,"p-value"],
    dev=summary(m)$dev.expl)
})
do.call(rbind, sig_check)





# Build Model -------------------------------------------------------------






set.seed(6666)

str(model_dt)
k_general <- 5


# ths oen for old daat that worked fien!; don;t change this
form <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  # s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(mltost,       month, bs = c("cr","cr"), k = c(4,3)) +
  # PROXY for east-coast vs offshore seasonal shift
  # ti(depth, month, bs = c("cr","cr"), k = c(4,3)) +
  s(id, bs = "re")



form <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  s(sst_slope, bs = "cr", k = 5) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  # ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  # ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(mltost,       month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(sst_slope,       month, bs = c("cr","cr"), k = c(4,3)) +
  # PROXY for east-coast vs offshore seasonal shift
  # ti(depth, month, bs = c("cr","cr"), k = c(4,3)) +
  s(id, bs = "re")





form
model_data <- model_dt
str(model_data)

# saveRDS(model_data, "data/processed/model_input_data_processed_crwPA.rds")
# 
# model_data <- readRDS("data/processed/model_input_data_processed.rds") # this one worked before; do not overwrite


weights_PA <- ifelse(model_data$PA == 1, 1, 0.1)
weights_PA

# month_counts <- table(model_data$month)
# month_factor <- 1 / as.numeric(month_counts[as.character(model_data$month)])
# weights_combined <- ifelse(model_data$PA == 1, 1, 0.1) * month_factor
# weights_combined <- weights_combined / mean(weights_combined)  # normalize
# weights_combined
# 
# summary(weights_combined)
# range(weights_combined)
# 
# weights_combined <- pmin(weights_combined, quantile(weights_combined, 0.99))
# weights_combined <- pmax(weights_combined, quantile(weights_combined, 0.01))
# weights_combined <- weights_combined / mean(weights_combined)
# summary(weights_combined)
# range(weights_combined)
# 
# 
# ggplot(model_data, aes(x=month, y=PA)) +
#   stat_summary(fun=mean, geom="line") +
#   geom_point(alpha=0.2) +
#   labs(y="Raw presence rate")


# model_jan <- model_data |> 
#   dplyr::filter(month == "1")
# weights_PA <- ifelse(model_jan$PA == 1, 1, 0.5)


# n1 <- sum(model_data$PA == 1)
# n0 <- sum(model_data$PA == 0)
# target_pi <- 0.25   # choose your target prevalence
# w_abs <- (n1 * (1 - target_pi)) / (n0 * target_pi)
# w_abs



gam.full <- mgcv::bam(form,
                 family = binomial(link="cloglog"), 
                 method="fREML", 
                 control =  mgcv::gam.control(maxit = 500, epsilon = 1e-5, nthreads = 6),
                 discrete = TRUE,
                 weights = weights_PA,
                 na.action = na.fail,
                 data = model_data,
                 select = TRUE,
                 gamma = 1
)


summary(gam.full)
AIC(gam.full)

plot(gam.full, pages = 1, scheme = 1, shade = TRUE,
     seWithMean = TRUE, rug = TRUE, scale = 0, residuals = FALSE, se = TRUE,
     shade.col = "steelblue1", 
     shift = coef(gam.full)[1],
     trans = stats::family(gam.full)$linkinv)

visreg::visreg(gam.full, xvar = "thetao", by = "month",
       overlay = TRUE, partial = FALSE,
       scale = "response", type = "conditional",
       rug = FALSE, 
       gg = TRUE) +
  ggplot2::labs(y = "Predicted probability", x = "thetao (°C)") +
  ggplot2::theme_bw()



# gam.2 <- MuMIn::dredge(gam.full,
#                        evaluate = TRUE,
#                        rank = "AICc",
#                        beta = "sd",
#                        trace = 2)
# 
# 
# nrow(gam.2)



# # saveRDS(gam.2, "models/objects/gamm_dredge_results")
# gam.2 <- readRDS("models/objects/gamm_dredge_results")
# 
# print(subset(gam.2, delta <= 10), abbrev.names = TRUE)
# 
# gam.2[1]
# 
# # 1) See the real column names the dredge table is using
# cols <- names(gam.2)
# print(cols, max.levels = 0)
# 
# # 2) Convert to a data.frame without mangling names
# 
# tab <- as.data.frame(gam.2, stringsAsFactors = FALSE, check.names = FALSE)
# 
# # helper: included terms are recorded as "+", 1, or TRUE in dredge table
# to01 <- function(x) as.integer(x %in% c("+", 1, TRUE))
# 
# # helper: does a row include ANY column whose name matches regex `pattern`?
# has_any <- function(pattern) {
#   cols_match <- grep(pattern, names(tab), value = TRUE)
#   if (length(cols_match) == 0L) return(rep(FALSE, nrow(tab)))
#   mat <- sapply(cols_match, function(nm) to01(tab[[nm]]))
#   if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
#   rowSums(mat) > 0
# }
# 
# # 3) Build must-have masks using REGEX that tolerates abbreviated or full names
# must_chl  <- has_any("^s\\(chl")                      # s(chl, ...)
# must_tht  <- has_any("^s\\((tht|thetao)")             # s(tht, ...) or s(thetao, ...)
# must_dpt  <- has_any("^s\\((dpt|depth)")              # s(dpt|depth, ...)
# must_d20  <- has_any("^s\\((d20|dist2000)")           # s(d20|dist2000, ...)
# must_uv   <- has_any("^s\\(uv")                       # s(uv, ...)
# must_wz   <- has_any("^s\\(wz")                       # s(wz, ...)
# 
# # 4) Month interaction: at least one of these present
# has_m_chl <- has_any("^ti\\(chl.*month")
# has_m_tht <- has_any("^ti\\((tht|thetao).*month")
# has_m_uv  <- has_any("^ti\\(uv.*month")
# has_m_wz  <- has_any("^ti\\(wz.*month")
# has_any_month_inter <- has_m_chl | has_m_tht | has_m_uv | has_m_wz
# 
# # 5) Final filter: all must-haves AND ≥1 month interaction
# keep_idx <- must_chl & must_tht & must_dpt & must_d20 & must_uv & must_wz & has_any_month_inter
# 
# filtered <- tab[keep_idx, , drop = FALSE]
# if (nrow(filtered) == 0L) {
#   stop("No model matches those constraints. Call `names(gam.2)` and paste the exact headers here; I’ll adjust the regex.")
# }
# 
# # 6) Pick the best by AICc among the filtered rows
# ord <- order(filtered$AICc)
# ord
# filtered_best <- filtered[ord[1], , drop = FALSE]
# print(filtered_best)
# 
# # 7) Recover the corresponding model object
# row_id <- as.integer(rownames(filtered_best)[1])
# best_model <- MuMIn::get.models(gam.2, subset = row_id)[[1]]
# 
# summary(best_model)
# 
# best.models <- MuMIn::get.models(gam.2, delta < 10)
# best.models
# 
# AIC(best_model)
# 
# best.gam <- MuMIn::get.models(gam.2, delta==5.76156)[[1]]
# best.gam
# 
# best.gam <- best.models[[3]]
# 
# summary(best.gam)
# 
# AIC(best.gam)
# AIC(gam.full)
# 
# plot(subset(gam.2, delta <= 10))



## selecting best model based ecological grounds not just stats

# gam <- best.gam
gam <- gam.full

summary(gam)

mgcv::k.check(gam)
mgcv::gam.check(gam)
mgcv::concurvity(gam, full = TRUE)
mgcv::concurvity(gam, full = FALSE)

plot(gam, pages = 1, scheme = 1, shade = TRUE,
     seWithMean = TRUE, rug = TRUE, scale = 0, residuals = FALSE, se = TRUE,
     shade.col = "steelblue1", 
     shift = coef(gam)[1],
     trans = function(eta) 1 - exp(-exp(eta)))

visreg::visreg(gam, scale = "response", type = "conditional", gg = TRUE)
visreg::visreg(gam, scale = "response", type = "contrast", gg = TRUE)



dev.off()

visreg::visreg2d(gam, xvar = "thetao", yvar = "uv", scale = "response", plot.type = "persp",
                 xlab = "thetao", ylab = "uv", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "thetao", yvar = "month", scale = "response", plot.type = "persp",
                 xlab = "thetao", ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 30, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "chl", yvar = "month", scale = "response", plot.type = "persp",
                 xlab = "chl", ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 30, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "uv", yvar = "month", scale = "response", plot.type = "persp",
                 xlab = "uv", ylab = "month", zlab = "Probability of presence",
                 theta = 120, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "wz", yvar = "month", scale = "response", plot.type = "persp",
                 xlab = "wz", ylab = "month", zlab = "Probability of presence",
                 theta = 120, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "depth", yvar = "month", scale = "response", plot.type = "persp",
                 xlab = "depth", ylab = "month", zlab = "Probability of presence",
                 theta = 120, phi = 15, zlim = c(0,1))


visreg::visreg2d(gam, xvar = "chl", yvar = "uv", scale = "response", plot.type = "rgl",
                 xlab = "chl", ylab = "uv", zlab = "Probability of presence",
                 theta = 70, phi = 30, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "thetao", yvar = "uv", scale = "response", plot.type = "rgl",
                 xlab = "thetao", ylab = "uv", zlab = "Probability of presence",
                 theta = 130, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "thetao", yvar = "month", scale = "response", plot.type = "rgl",
                 xlab = "thetao", ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "chl", yvar = "month", scale = "linear", plot.type = "rgl",
                 xlab = "chl", ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 15)

visreg::visreg2d(gam, xvar = "thetao", yvar = "chl", scale = "response", plot.type = "rgl",
                 xlab = "thetao", ylab = "chl", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))


visreg::visreg2d(gam, xvar = "lon", yvar = "month", scale = "response", plot.type = "rgl",
                 xlab = "lon", ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "lon", yvar = "lat", scale = "response", plot.type = "rgl",
                 xlab = "lon", ylab = "lat", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "depth", yvar = "month", scale = "response", plot.type = "rgl", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))
visreg::visreg2d(gam, xvar = "chl", yvar = "dist200", scale = "response", plot.type = "rgl", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))

visreg::visreg2d(gam, xvar = "dist2000", yvar = "month", scale = "response", plot.type = "rgl", zlab = "Probability of presence",
                 theta = 145, phi = 15, zlim = c(0,1))



mgcv::vis.gam(
  gam,
  view = c("chl","month"),
  plot.type = "contour",
  too.far = 0.05,
  zlim = NULL,                    # link scale
  color = "terrain",
  exclude = "s(chl)"              # <-- show the interaction surface only
)

mgcv::vis.gam(
  gam,
  view = c("chl","month"),
  plot.type = "link",
  theta = 145, phi = 20,
  zlim = NULL,
  color = "terrain",
  exclude = "s(chl)"
)

mgcv::vis.gam(
  gam,
  view = c("thetao","month"),
  plot.type = "persp",
  theta = 145, phi = 20,
  zlim = NULL,
  exclude = "s(chl)"
)

mgcv::vis.gam(
  gam,
  view = c("uv","month"),
  plot.type = "persp",
  theta = 145, phi = 20,
  zlim = NULL,
  exclude = "s(chl)"
)

mgcv::vis.gam(
  gam,
  view = c("wz","month"),
  plot.type = "persp",
  theta = 145, phi = 20,
  zlim = NULL,
  exclude = "s(chl)"
)

mgcv::vis.gam(
  gam,
  view = c("depth","month"),
  plot.type = "persp",
  theta = 145, phi = 20,
  zlim = NULL,
  exclude = "s(depth)"
)

mgcv::vis.gam(
  gam,
  view = c("chl","month"),
  plot.type = "contour",
  too.far = 0.05,
  zlim = NULL                     # link scale (no trans)
  # no exclude -> full prediction on link scale
)

gam.vcomp(gam)

AIC(gam)




## 1) Per-term contributions on the link scale (exclude random effect)
terms_link <- predict(
  object  = gam,
  newdata = model_data,
  type    = "terms",
  exclude = "s(id)"        # don't let RE swamp importance
)

terms_link

## Term-level importance: mean absolute contribution
term_imp <- colMeans(abs(terms_link), na.rm = TRUE)
term_imp_df <- data.frame(
  term = names(term_imp),
  mean_abs = as.numeric(term_imp)
) |>
  dplyr::arrange(dplyr::desc(mean_abs))

# quick look
print(term_imp_df)

base_preds <- c("thetao","uv","wz","chl","mltost", "depth","month", "slope","dist2000")

## Helper: which predictors appear in a term name?
pattern <- paste0("\\b(", paste(base_preds, collapse = "|"), ")\\b")

pred_imp <- stats::setNames(numeric(length(base_preds)), base_preds)


for (j in seq_along(term_imp)) {
  term_name <- names(term_imp)[j]
  vars_in_term <- stringr::str_extract_all(term_name, pattern)[[1]] |> unique()
  if (length(vars_in_term) > 0L) {
    share <- term_imp[[j]] / length(vars_in_term)       # split interaction fairly
    for (v in vars_in_term) pred_imp[[v]] <- pred_imp[[v]] + share
  }
}

## Sorted per-predictor importance + % of total
pred_imp_df <- data.frame(
  predictor = names(pred_imp),
  mean_abs  = as.numeric(pred_imp)
) |>
  dplyr::mutate(pct = 100 * mean_abs / sum(mean_abs)) |>
  dplyr::arrange(dplyr::desc(mean_abs))

print(pred_imp_df)

## Optional: bar plots
ggplot2::ggplot(term_imp_df, ggplot2::aes(x = reorder(term, mean_abs), y = mean_abs)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean |contribution| (link scale)", title = "GAM term importance")

ggplot2::ggplot(pred_imp_df, ggplot2::aes(x = reorder(predictor, mean_abs), y = mean_abs)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean |contribution| (link scale)", title = "Per-predictor importance")






predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log.tif")
# predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_gam.tif")
names(predictors_mean)
plot(predictors_mean)


gam <- gam.full

relevant_vars <- all.vars(gam$formula)[-1] # Excluding the intercept
gam$formula
relevant_vars
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general", "year", "month")) ]


predictors <- predictors_mean[[relevant_vars]]
predictors
names(predictors)
plot(predictors)


map <- terra::predict(predictors,
                      model = gam,
                      type = "response",
                      const = (data.frame( month = 12)),
                      exclude = c("s(id)", "s(year)", "s(month)"),
                      na.rm = TRUE
)


map


plot(map, range = c(0, 1))
terra::plot(map, range = c(0, 1), xlim = c(140, 170), ylim = c(-35, 0))

dev.off()


# mean contribution by month for each term (exclude RE)
terms_link <- predict(gam, newdata = model_data, type = "terms", exclude = "s(id)")
contrib <- cbind(month = model_dt_trans$month, as.data.frame(terms_link))

monthly_term_means <- contrib |>
  dplyr::group_by(month) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)))

print(monthly_term_means |> dplyr::select(month, dplyr::starts_with("ti(")))



model_data |>
  dplyr::group_by(month) |>
  dplyr::summarise(thetao_mean = mean(thetao, na.rm = TRUE),
                   chl_mean    = mean(chl,    na.rm = TRUE)) |>
  print()






# sample covariates to check concurvity/correlation
dat_samp <- terra::spatSample(predictors, size = 100000, method = "regular",
                              as.df = TRUE, na.rm = TRUE)
str(dat_samp)
# rank-based correlations (robust)
cor_mat <- stats::cor(dat_samp[, c("depth","slope","uv","thetao","wz","chl","dist2000")],
                      method = "spearman", use = "pairwise")
print(cor_mat)

cc <- mgcv::concurvity(gam, full = TRUE)

cc_est <- if (base::is.list(cc)) cc$estimate else cc
cc_wst <- if (base::is.list(cc)) cc$worst     else NULL
cc_obs <- if (base::is.list(cc)) cc$observed  else NULL

base::print(cc_est)  # main concurvity matrix (cols explain rows)

cc_simple <- mgcv::concurvity(gam, full = FALSE)
base::print(cc_simple)  # often a named list of vectors




k_general <- 5
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_3to7days.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
# sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")

sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced_final_4.rds")

sp_blocks$records

# GAMM/GAM formula (shrinkage to let terms drop if redundant)
fml <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(mltost,   month, bs = c("cr","cr"), k = c(4,3)) 
  # PROXY for east-coast vs offshore seasonal shift
  # ti(depth, month, bs = c("cr","cr"), k = c(4,3)) +
  # s(id, bs = "re")








sp_blocks <- sb1
gamm_cv_out <- gam_cv(
  dat   = model_data,             
  fid   = sp_blocks$folds_ids,                   
  fml   = fml,          
  family = binomial(link = "cloglog"),
  base_preds  = c("thetao", "mltost", "chl", "uv", "wz", "slope", "depth", "dist2000", "month"),
  re_term    = "s(id)"
)


gamm_cv_out


gamm_cv_out$summary

gamm_vi_df <- gamm_cv_out$varimp |>
  tibble::as_tibble() |>
  dplyr::mutate(algorithm = "GAMM")

gamm_vi_df


saveRDS(gamm_vi_df, "models/tables/GAMM_Var_Importance_CV.rds")



# Miller MCS

# gcv is what gam_cv() returned
models <- gamm_cv_out$models

# container for out-of-fold predictions
pred_oof <- numeric(nrow(dat))

kvals <- sort(unique(fid))
kvals
for (i in seq_along(kvals)) {
  k     <- kvals[i]
  te_id <- fid == k
  m_i   <- models[[i]]
  
  pred_oof[te_id] <- stats::predict(
    object  = m_i,
    newdata = dat[te_id, , drop = FALSE],
    type    = "response"
  )
}

mc_cv <- modEvA::MillerCalib(
  obs        = dat$PA,
  pred       = pred_oof,
  plot       = TRUE,
  main       = "Miller calibration – GAMM (CV out-of-fold predictions)",
  xlab       = "Predicted probability",
  ylab       = "Observed presence (logit scale)",
  line.col   = "darkblue",
  diag       = TRUE,
  diag.col   = "grey80",
  plot.values = TRUE
)

mc_cv  # list with intercept, slope, slopeDiff





# Model Residuals ---------------------------------------------------------



resid <- mgcv::residuals.gam(gam, type = "pearson")



# 3) Add coords + residuals
resid_df <- model_data |>
  dplyr::mutate(resid = resid) |>
  dplyr::select(lon, lat, PA, resid)   # adjust lon/lat names if needed

str(resid_df)

# 4) sf + projection
resid_sf <- resid_df |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  sf::st_transform(crs = 3577)

# 5) Spatial AC of residuals
ac_resid <- blockCV::cv_spatial_autocor(
  x      = resid_sf,
  column = "resid"
)

summary(ac_resid)



#  fold’s GAM
m <- cv$models[[5]]
summary(m)
plot(m, pages=1)


model <- m
class(model)
summary(model)

par(mfrow = c(3,3))
plot(model, pages = 1, scheme = 1, shade = TRUE,
     seWithMean = TRUE, rug = TRUE, scale = 0)


map <- terra::predict(predictors,
                      model = model,
                      type = "response",
                      const = (data.frame(month = 12)),
                      exclude = c( "s(id)"),
                      na.rm = TRUE
)


map

par(mfrow = c(1,1))
plot(map, range = c(0, 1))
terra::plot(map, range = c(0, 1), xlim = c(140, 170), ylim = c(-35, 0))


# saveRDS(gam, "models/objects/Tracks_GAMM_final_Model_mp.rds")
# saveRDS(cv, "models/objects/Tracks_GAMM_CV_Models_mp.rds")

saveRDS(gam, "models/objects/Tracks_GAMM_final_Model_mp_crwPA.rds")
saveRDS(cv, "models/objects/Tracks_GAMM_CV_Models_mp_crwPA.rds")


gam <- readRDS("models/objects/Tracks_GAMM_final_Model_mp_crwPA.rds")
cv <- readRDS("models/objects/Tracks_GAMM_CV_Models_mp_crwPA.rds")
# Response curves ---------------------------------------------------------
model <- gam
rug = 4

p.depth <- sjPlot::plot_model(model, type = "pred", terms = "depth", colors = c("turquoise1"),
                              pred.type = "fe",
                              axis.lim =  c(0, 1),
                              show.data = FALSE,
                              show.intercept = FALSE) +
  labs(x = expression("depth (m)"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = depth),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = depth),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.depth



p.slope <- sjPlot::plot_model(model, type = "pred", terms = "slope", colors = c("turquoise1"),
                              #pred.type = "re",
                              axis.lim =  c(0, 1),
                              show.data = FALSE,
                              show.intercept = FALSE) +
  # geom_jitter(data = model_dt, aes(x = Slope, y = PA), 
  #             color = "grey20", width = 0, height = 0, alpha = 0.5, size = 1, shape = 1) + 
  labs(x = expression("slope ("*degree*")"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = slope),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = slope),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.slope



p.dist <- sjPlot::plot_model(model, type = "pred", terms = "dist2000", colors = c("turquoise1"),
                             #pred.type = "re",
                             axis.lim =  c(0, 1),
                             show.data = FALSE,
                             show.intercept = FALSE) +
  labs(x = "dist2000 (km)", y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = dist2000),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = dist2000),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.dist



p.sst <- sjPlot::plot_model(model, type = "pred", terms = "thetao", colors = c("turquoise1"),
                            #pred.type = "re",
                            axis.lim =  c(0, 1),
                            show.data = FALSE,
                            show.intercept = FALSE) +
  labs(x = expression("sst ("*degree*"C)"),  y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = thetao),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = thetao),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.sst



p.uv <- sjPlot::plot_model(model, type = "pred", terms = "uv", colors = c("turquoise1"),
                           #pred.type = "re",
                           axis.lim =  c(0, 1),
                           show.data = FALSE,
                           show.intercept = FALSE,
                           show.p = FALSE) +
  # geom_jitter(data = model_dt, aes(x = uv, y = PA),
  #             color = "grey20", width = 0, height = 0, alpha = 0.5, size = 1, shape = 1) +
  labs(x = expression("uv (m s"^{-1}*")"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = uv),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = uv),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
p.uv


p.wz <- sjPlot::plot_model(model, type = "pred", terms = "wz", colors = c("turquoise1"),
                           #pred.type = "re",
                           axis.lim =  c(0, 1),
                           show.data = FALSE,
                           show.intercept = FALSE) +
  labs(x = expression("wz (m s"^{-1}*")"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = wz),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = wz),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.wz


p.chl <- sjPlot::plot_model(model, type = "pred", terms = "chl", colors = c("turquoise1"),
                           #pred.type = "re",
                           axis.lim =  c(0, 1),
                           show.data = FALSE,
                           show.intercept = FALSE) +
  labs(x = expression("chl (log(mg m"^{-3}*"))"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = chl),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = chl),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.chl



p.mld <- sjPlot::plot_model(model, type = "pred", terms = "mltost", colors = c("turquoise1"),
                            #pred.type = "re",
                            axis.lim =  c(0, 1),
                            show.data = FALSE,
                            show.intercept = FALSE) +
  labs(x = expression("mld (m)"), y = "Probability of Presence", title = NULL) +
  # presences (top)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = mltost),
    inherit.aes = FALSE,
    sides = "t",                    # top
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  # absences (bottom)
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = mltost),
    inherit.aes = FALSE,
    sides = "b",                    # bottom
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.mld



P_effects <- 
  (p.depth +  p.slope + p.dist + p.sst + p.chl + p.uv + p.wz + p.mld) +
  plot_layout(ncol = 3, axis_titles  = "collect_y") 

P_effects

P_effects2 <- cowplot::ggdraw(P_effects) + cowplot::draw_plot_label("A)", x = 0, y = 1, hjust = -.2, vjust = 1.5, size = 12)

P_effects2




# ggsave("GAMM_RepsonsePlots_Marginal_Tracking_mp.png", plot = P_effects2, path ="outputs/final_figures", scale =1, width = 18, height = 20, units = "cm", dpi = 300)




pal_bal <- cmocean::cmocean(
  name = "balance", alpha = 1,
  start = 0, end = 1, direction = -1
)(256)


model
p.sst_m <- visreg::visreg2d(model, xvar = "thetao", yvar = "month", scale = "response", plot.type = "gg",
                 xlab = expression("sst ("*degree*"C)"), ylab = "month", zlab = "Probability of presence",
                 theta = 145, phi = 30, zlim = c(0,1)) +
  ggplot2::scale_fill_gradientn(
    colours = cmocean::cmocean(
      name = "thermal", alpha = 1,
      start = 0.1, end = 0.9, direction = 1
    )(256),
    # limits = c(0, 1),
    name   = "Probability of\npresence"
  ) +
  # ggplot2::scale_fill_viridis_c(name = "Probability of\npresence", direction = 1) +
  # ggplot2::scale_y_continuous(breaks = 1:12) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        legend.title        = ggplot2::element_text(size = 6, face = "bold", vjust = 1),
        legend.text         = ggplot2::element_text(size = 4),
        legend.key.height   = grid::unit(2,  "mm"), 
        legend.key.width    = grid::unit(3,  "mm"),
        legend.margin       = ggplot2::margin(2, 2, 2, 2),
        legend.box.margin   = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing  = grid::unit(1, "mm"))

p.sst_m


p.chl_m <- visreg::visreg2d(model, xvar = "chl", yvar = "month", scale = "response", plot.type = "gg",
                            xlab = expression("chl (log(mg m"^{-3}*"))"), ylab = "month") +
  ggplot2::scale_fill_gradientn(
    colours = cmocean::cmocean(
      name = "algae", alpha = 1,
      start = 0, end = 1, direction = 1
    )(256),
    # limits = c(0, 1),
    name   = "Probability of\npresence"
  ) +
  # ggplot2::scale_fill_viridis_c(name = "Probability of\npresence", direction = 1) +
  # ggplot2::scale_y_continuous(breaks = 1:12) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        legend.title        = ggplot2::element_text(size = 6, face = "bold", vjust = 1),
        legend.text         = ggplot2::element_text(size = 4),
        legend.key.height   = grid::unit(2,  "mm"), 
        legend.key.width    = grid::unit(3,  "mm"),
        legend.margin       = ggplot2::margin(2, 2, 2, 2),
        legend.box.margin   = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing  = grid::unit(1, "mm"))


p.chl_m


p.uv_m <- visreg::visreg2d(model, xvar = "uv", yvar = "month", scale = "response", plot.type = "gg",
                            xlab = expression("uv (m s"^{-1}*")"), ylab = "month") +
  ggplot2::scale_fill_gradientn(
    colours = cmocean::cmocean(
      name = "speed", alpha = 1,
      start = 0, end = 1, direction = -1
    )(256),
    # limits = c(0, 1),
    name   = "Probability of\npresence"
  ) +
  # ggplot2::scale_fill_viridis_c(name = "Probability of\npresence", direction = 1) +
  # ggplot2::scale_y_continuous(breaks = 1:12) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        legend.title        = ggplot2::element_text(size = 6, face = "bold", vjust = 1),
        legend.text         = ggplot2::element_text(size = 4),
        legend.key.height   = grid::unit(2,  "mm"), 
        legend.key.width    = grid::unit(3,  "mm"),
        legend.margin       = ggplot2::margin(2, 2, 2, 2),
        legend.box.margin   = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing  = grid::unit(1, "mm"))


p.uv_m

p.wz_m <- visreg::visreg2d(model, xvar = "wz", yvar = "month", scale = "response", plot.type = "gg",
                           xlab = expression("wz (m s"^{-1}*")"), ylab = "month") +
  ggplot2::scale_fill_gradientn(
    colours = cmocean::cmocean(
      name = "haline", alpha = 1,
      start = 0, end = 1, direction = 1
    )(256),
    # limits = c(0, 1),
    name   = "Probability of\npresence"
  ) +
  # ggplot2::scale_fill_viridis_c(name = "Probability of\npresence", direction = 1) +
  # ggplot2::scale_y_continuous(breaks = 1:12) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        legend.title        = ggplot2::element_text(size = 6, face = "bold", vjust = 1),
        legend.text         = ggplot2::element_text(size = 4),
        legend.key.height   = grid::unit(2,  "mm"), 
        legend.key.width    = grid::unit(3,  "mm"),
        legend.margin       = ggplot2::margin(2, 2, 2, 2),
        legend.box.margin   = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing  = grid::unit(1, "mm"))


p.wz_m

p.mld_m <- visreg::visreg2d(model, xvar = "mltost", yvar = "month", scale = "response", plot.type = "gg",
                           xlab = expression("mld (m)"), ylab = "month") +
  ggplot2::scale_fill_gradientn(
    colours = cmocean::cmocean(
      name = "ice", alpha = 1,
      start = 0, end = 1, direction = 1
    )(256),
    # limits = c(0, 1),
    name   = "Probability of\npresence"
  ) +
  # ggplot2::scale_fill_viridis_c(name = "Probability of\npresence", direction = 1) +
  # ggplot2::scale_y_continuous(breaks = 1:12) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.position = "bottom",
        legend.title        = ggplot2::element_text(size = 6, face = "bold", vjust = 1),
        legend.text         = ggplot2::element_text(size = 4),
        legend.key.height   = grid::unit(2,  "mm"), 
        legend.key.width    = grid::unit(3,  "mm"),
        legend.margin       = ggplot2::margin(2, 2, 2, 2),
        legend.box.margin   = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing  = grid::unit(1, "mm"))


p.mld_m



P_int <- (p.sst_m + p.chl_m + p.uv_m + p.wz_m + p.mld_m) +
  patchwork::plot_layout(ncol = 3, 
                         # guides = "collect",
                         axes = "collect"
                         ) &
  ggplot2::theme(aspect.ratio = 1) &
  ggplot2::coord_cartesian(expand = FALSE) &
  ggplot2::theme(plot.margin = grid::unit(c(5,5,5,5), "pt"))

P_int


P_int2 <- cowplot::ggdraw(P_int) + cowplot::draw_plot_label("B)", x = 0, y = 0.9, hjust = -.2, vjust = 1.5, size = 12)

P_int2



effects_final <- (P_effects2 / P_int2) +
  patchwork::plot_layout(widths = c(1,1), heights = c(1, 1)) &
  ggplot2::theme(aspect.ratio = 1) &
  # ggplot2::coord_cartesian(expand = FALSE) &
  ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "pt"))

effects_final


ggsave("GAMM_RepsonsePlots_Marginal_Tracking_mp_crwPA.png", plot = effects_final, path ="outputs/final_figures", scale =1, width = 18, height = 30, units = "cm", dpi = 300)


# External Vaildation -----------------------------------------------------



# create extrnal validation data from sightings 

sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full.rds")
str(sight)



mapview::mapview(sight |> dplyr::select(-Date))


val_dt <- sight |>
  # dplyr::filter(PA == 1) |>
  # dplyr::filter(lon >150 & lon < 180) |> # making external dataset spatially independent of calibrationn data
  # dplyr::filter(lon >142 & lon < 160) |>
  # dplyr::filter(lat >-17 & lat < -5) |>
  # dplyr::mutate(month = factor(month, levels = sort(unique(month))),  # 1..12 etc.
  #               year  = factor(year)) |> 
  # dplyr::rename(depth = Depth,
  #               slope = Slope,
  #               roughness = Roughness,
  #               wz = Wz) |>
  sf::st_drop_geometry() |> 
  dplyr::filter(!month %in% c(7, 8, 9, 10)) |>
  as.data.frame()

str(val_dt)



val_dt <- val_dt |> 
  dplyr::mutate(
    # slope = log1p(slope),
    # roughness = log1p(roughness),
    chl = log(chl),
    dist2000 = dist2000/1000)




plots <- val_dt

(plots |>  ggplot(aes(x = depth, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = PA)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)

val_dt |>
  dplyr::select(-Date) |>
  dplyr::filter(PA == 1) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  mapview::mapview(col.regions = "white", color = "white")


# ## save for repo
# val_dt_repo <- val_dt |> 
#     dplyr::select(id, lon, lat, PA, rep, month, depth, slope, roughness, dist2000, thetao,uv, mltost, chl, wz)
# write_csv(val_dt_repo, "models/repo/ext_val_data_repo.csv")


pred_prob <- predict(
  gam,
  newdata = val_dt,
  type    = "response",
  exclude = "s(id)"
)

p <- pred_prob[val_dt$PA == 1]
a <- pred_prob[val_dt$PA == 0]

ev <- dismo::evaluate(p = p, a = a)

# AUC
auc_ext <- ev@auc
auc_ext

roc_obj <- pROC::roc(response = val_dt$PA, predictor = pred_prob, quiet = TRUE, plot = TRUE, col = "steelblue", lwd = 2, legacy.axes = TRUE)
roc_obj
auc_ext <- as.numeric(pROC::auc(roc_obj))
auc_ext

# TSS at the max-sum threshold (sens + spec)
thr     <- dismo::threshold(ev, stat = "spec_sens")  # same operating point SDMtune uses by default
tss_ext <- ev@TPR[which.max(ev@TPR + ev@TNR)] + ev@TNR[which.max(ev@TPR + ev@TNR)] - 1
tss_ext

roc_df <- data.frame(FPR = 1 - ev@TNR, TPR = ev@TPR)
gg_ext <- ggplot2::ggplot(roc_df, ggplot2::aes(FPR, TPR)) +
  ggplot2::geom_line() +
  ggplot2::geom_abline(linetype = 2) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = sprintf("External validation (AUC = %.2f)", auc_ext),
                x = "False positive rate", y = "True positive rate")

print(gg_ext)


cat(sprintf("External AUC = %.3f | TSS = %.3f | Threshold = %.3f\n", auc_ext, tss_ext, thr))


pred_bin <- ifelse(pred_prob >= thr, 1, 0)

y <- val_dt$PA  # 0/1

tp <- sum(pred_bin == 1 & y == 1)
fn <- sum(pred_bin == 0 & y == 1)

omission_ext <- fn / (tp + fn)
omission_ext


boyce_ext <- ecospat::ecospat.boyce(
  fit    = pred_prob,               # predictions for ALL external points (pres + abs)
  obs    = pred_prob[y == 1],       # predictions at external presences
  PEplot = TRUE                     # or FALSE if you don't want the plot
)

boyce_ext$cor 


# Monthlky predictions ----------------------------------------------------


model <- gam
model
relevant_vars <- all.vars(model$formula)[-1] # Excluding the intercept
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general")) ]
relevant_vars

# Extract relevant layers from each raster stack in the list

min(tracks$date) 
max(tracks$date) 

input_list <- monthly_stacks_lst_0.1_trans
names(input_list)

input_list <- lapply(input_list, function(stack) {
  stack[[relevant_vars]] # Extract layers that match the model variables
})

input_list <- input_list[119:186]

names(input_list)
print(input_list[[65]])
names(input_list[[65]])

plot(input_list[[65]])



# input_list <- input_list[175:186]
dates <- as.Date(names(input_list))

monthly_predictions <- list()
N <- length(input_list)

tictoc::tic("Loop run took: " )
for (i in seq_along(input_list)) {
  # Extract the current monthly raster stack
  monthly_stack <- input_list[[i]]
  
  # Extract the current raster's date"
  Month_Date <- names(input_list[i])
  
  # Print information about the current prediction
  cat("Prediction", i, "of", N, "\n")
  
  # Generate predictions for the current monthly raster stack
  predictions <- terra::predict(monthly_stack,
                                model,
                                # const = (data.frame(Year = "0")),
                                #cores = 6,
                                type = "response"
                                )
  
  #predictions <- dismo::predict(maxent_model, monthly_stack)
  
  # Assign Date as layer name:
  names(predictions) <- Month_Date
  
  # Store the predictions in the list
  monthly_predictions[[i]] <- predictions
}
tictoc::toc()

monthly_predictions

monthly_predictions_stack <- terra::rast(monthly_predictions)
monthly_predictions_stack
names(monthly_predictions_stack)


# writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_2018_2025_rev_mp.tif", overwrite = TRUE)

writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_2018_2025_rev_mp_crwPA.tif", overwrite = TRUE)





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

plot(monthly_means_stack, range = c(0,1), xlim = c(140, 170), ylim = c(-35, 0))
plot(monthly_means_stack[[12]], range = c(0,1), xlim = c(140, 170), ylim = c(-40, 0))


# Save for ensemble modelling
# writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_rev_mp.tif", overwrite = TRUE)

writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_rev_mp_crwPA.tif", overwrite = TRUE)



## Overall mean 

mean_climate <- terra::app(monthly_means_stack, fun = mean, na.rm = TRUE)
mean_climate
plot(mean_climate)

# Save for ensemble modelling
# writeRaster(mean_climate, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_mean_rev_mp.tif", overwrite = TRUE)

writeRaster(mean_climate, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_mean_rev_mp_crw_PA.tif", overwrite = TRUE)





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

plot(seasonal_means_stack_2, range = c(0, 1), xlim = c(140, 170), ylim = c(-40, 0))
plot(seasonal_means_stack_2, range = c(0, 1), xlim = c(140, 170), ylim = c(-40, 0), axes =FALSE, legend =FALSE)

plot(seasonal_means_stack_2, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100), xlim = c(140, 170), ylim = c(-40, 0))



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
                                                      col.regions = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100),
                                                      colorkey = list(space = "right", 
                                                                      length = 0.5, 
                                                                      height = 0.75,
                                                                      labels = list(
                                                                        at = c(0, high),  # Positions for "Low" and "High"
                                                                        labels = c("Low", "High"))))  # Labels to use

Mean_Season_monsoon_SDM_plot




# Save for ensemble modelling
# writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp.tif", overwrite = TRUE)

writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp_crwPA.tif", overwrite = TRUE)




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

plot(seasonal_means_stack_4, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100),  xlim = c(140, 170), ylim = c(-40, 0))

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



# writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons4_means_rev_mp.tif", overwrite = TRUE)

writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons4_means_rev_mp_crwPA.tif", overwrite = TRUE)



