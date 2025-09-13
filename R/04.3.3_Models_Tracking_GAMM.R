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


tracks <- readRDS("data/processed/Tracks_PA_w_dynSDM_10_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_3days_dynSDM_10_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_3to7days_dynSDM_10_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed.rds")
tracks <- readRDS("data/processed/Tracks_PA_w_daily_mp_dynSDM_30_2018_2025_extract_processed_monthsreduced.rds")

tracks


tracks |> 
  dplyr::group_by(PA) |> 
  dplyr::summarise(N = n())

model_dt <- tracks |>
  dplyr::mutate(year  = factor(year),
                depth = case_when(depth >= 0 ~ -5, TRUE ~ depth),
                id = as.factor(id),
                dist2000 = dist2000/1000) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()

tracks |> 
  group_by(PA) |>
  summarise(N = n())





model_dt |> 
  group_by(PA) |>
  summarise(N = n())


model_dt |> rstatix::get_summary_stats(depth, slope, roughness, thetao, uv, chl, dist2000, dist200, wz, mltost, dist_seamount, dist_knoll )

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
    plot_layout(guides = "collect")
)


# Transformations

model_dt_trans <- model_dt |>
  dplyr::mutate(
    chl = log(chl),
    month_fc = as.factor(month)
  )

model_dt_trans


model_dt_trans |> 
  rstatix::get_summary_stats(depth, slope, roughness, thetao, uv, chl, dist2000, wz, mltost, dist_seamount, dist_knoll )
model_dt |> 
  rstatix::get_summary_stats(depth, slope, roughness, thetao, uv, chl,  dist2000, wz, mltost, dist_seamount, dist_knoll )



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
  dplyr::select(depth, slope, roughness, dist2000, thetao, chl, uv, wz, mltost, dist_seamount, dist_knoll) |>
  psych::pairs.panels(method = "spearman", cor = TRUE, scale = FALSE, stars = FALSE)

model_dt_trans |>
  # dplyr::filter(PA == 0) |> 
  dplyr::select(depth, slope, roughness, dist2000, thetao, chl, uv, wz, mltost, dist_seamount, dist_knoll) |>
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

str(model_dt_trans)
k_general <- 5



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
  # PROXY for east-coast vs offshore seasonal shift
  # ti(depth, month, bs = c("cr","cr"), k = c(4,3)) +
  s(id, bs = "re")





form
model_data <- model_dt_trans

saveRDS(model_dt, "data/processed/model_input_data_processed.rds")

weights_PA <- ifelse(model_data$PA == 1, 1, 0.2)
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


n1 <- sum(model_data$PA == 1)
n0 <- sum(model_data$PA == 0)
target_pi <- 0.25   # choose your target prevalence
w_abs <- (n1 * (1 - target_pi)) / (n0 * target_pi)
w_abs



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
     shift = coef(gam)[1],
     trans = stats::family(gam.full)$linkinv)

visreg::visreg(gam.full, xvar = "thetao", by = "month",
       overlay = TRUE, partial = FALSE,
       scale = "response", type = "conditional",
       rug = FALSE, 
       gg = TRUE) +
  ggplot2::labs(y = "Predicted probability", x = "thetao (°C)") +
  ggplot2::theme_bw()

visreg::visreg(gam.full, xvar = "wz", by = "month_fc",
               overlay = TRUE, partial = FALSE,
               scale = "response", type = "conditional",
               rug = FALSE, 
               gg = TRUE) +
  ggplot2::labs(y = "Predicted probability", x = "thetao (°C)") +
  ggplot2::theme_bw()

visreg::visreg(gam.full, xvar = "chl", by = "month_fc",
               overlay = TRUE, partial = FALSE,
               scale = "response", type = "conditional",
               rug = FALSE, 
               gg = TRUE) +
  ggplot2::labs(y = "Predicted probability", x = "thetao (°C)") +
  ggplot2::theme_bw()

visreg::visreg(gam.full, xvar = "dist2000", by = "month_fc",
               overlay = TRUE, partial = FALSE,
               scale = "response", type = "conditional",
               rug = FALSE, 
               gg = TRUE) +
  ggplot2::labs(y = "Predicted probability", x = "thetao (°C)") +
  ggplot2::theme_bw()


gam.2 <- MuMIn::dredge(gam.full,
                       evaluate = TRUE,
                       rank = "AICc",
                       beta = "sd",
                       trace = 2)


nrow(gam.2)






# saveRDS(gam.2, "models/objects/gamm_dredge_results")
gam.2 <- readRDS("models/objects/gamm_dredge_results")

print(subset(gam.2, delta <= 5), abbrev.names = TRUE)

gam.2[1]

# 1) See the real column names the dredge table is using
cols <- names(gam.2)
print(cols, max.levels = 0)

# 2) Convert to a data.frame without mangling names

tab <- as.data.frame(gam.2, stringsAsFactors = FALSE, check.names = FALSE)

# helper: included terms are recorded as "+", 1, or TRUE in dredge table
to01 <- function(x) as.integer(x %in% c("+", 1, TRUE))

# helper: does a row include ANY column whose name matches regex `pattern`?
has_any <- function(pattern) {
  cols_match <- grep(pattern, names(tab), value = TRUE)
  if (length(cols_match) == 0L) return(rep(FALSE, nrow(tab)))
  mat <- sapply(cols_match, function(nm) to01(tab[[nm]]))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
  rowSums(mat) > 0
}

# 3) Build must-have masks using REGEX that tolerates abbreviated or full names
must_chl  <- has_any("^s\\(chl")                      # s(chl, ...)
must_tht  <- has_any("^s\\((tht|thetao)")             # s(tht, ...) or s(thetao, ...)
must_dpt  <- has_any("^s\\((dpt|depth)")              # s(dpt|depth, ...)
must_d20  <- has_any("^s\\((d20|dist2000)")           # s(d20|dist2000, ...)
must_uv   <- has_any("^s\\(uv")                       # s(uv, ...)
must_wz   <- has_any("^s\\(wz")                       # s(wz, ...)

# 4) Month interaction: at least one of these present
has_m_chl <- has_any("^ti\\(chl.*month")
has_m_tht <- has_any("^ti\\((tht|thetao).*month")
has_m_uv  <- has_any("^ti\\(uv.*month")
has_m_wz  <- has_any("^ti\\(wz.*month")
has_any_month_inter <- has_m_chl | has_m_tht | has_m_uv | has_m_wz

# 5) Final filter: all must-haves AND ≥1 month interaction
keep_idx <- must_chl & must_tht & must_dpt & must_d20 & must_uv & must_wz & has_any_month_inter

filtered <- tab[keep_idx, , drop = FALSE]
if (nrow(filtered) == 0L) {
  stop("No model matches those constraints. Call `names(gam.2)` and paste the exact headers here; I’ll adjust the regex.")
}

# 6) Pick the best by AICc among the filtered rows
ord <- order(filtered$AICc)
ord
filtered_best <- filtered[ord[1], , drop = FALSE]
print(filtered_best)

# 7) Recover the corresponding model object
row_id <- as.integer(rownames(filtered_best)[1])
best_model <- MuMIn::get.models(gam.2, subset = row_id)[[1]]

summary(best_model)

best.models <- MuMIn::get.models(gam.2, delta < 10)
best.models

AIC(best_model)

best.gam <- MuMIn::get.models(gam.2, delta==5.76156)[[1]]
best.gam

best.gam <- best.models[[3]]

summary(best.gam)

AIC(best.gam)
AIC(gam.full)

plot(subset(gam.2, delta <= 10))



## selecting best model based ecological grounds not just stats

gam <- best.gam
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
  plot.type = "persp",
  theta = 145, phi = 20,
  zlim = NULL,
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
  newdata = model_dt_trans,
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

base_preds <- c("thetao","uv","wz","chl","depth","month", "lon", "lat", "slope","dist200", "dist1000", "dist2000", "dist_seamount", "dist_knoll")

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
                      const = (data.frame( month = 2)),
                      exclude = c("s(id)", "s(year)", "s(month)"),
                      na.rm = TRUE
)


map


plot(map, range = c(0, 1))
terra::plot(map, range = c(0, 1), xlim = c(140, 170), ylim = c(-35, 0))

dev.off()


# mean contribution by month for each term (exclude RE)
terms_link <- predict(gam, newdata = model_dt_trans, type = "terms", exclude = "s(id)")
contrib <- cbind(month = model_dt_trans$month, as.data.frame(terms_link))

monthly_term_means <- contrib |>
  dplyr::group_by(month) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)))

print(monthly_term_means |> dplyr::select(month, dplyr::starts_with("ti(")))



model_dt_trans |>
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
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_3to7days.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp.rds")
sp_blocks <- readRDS("data/processed/BlockCV_spatial_folds_tracking_daily_mp_monthreduced.rds")



# GAMM/GAM formula (shrinkage to let terms drop if redundant)
fml <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  # s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  ti(thetao, chl, bs = c("tp","cr"), k = c(4,4)) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  # ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  # PROXY for east-coast vs offshore seasonal shift
  ti(depth, month, bs = c("cr","cr"), k = c(4,3)) +
  s(id, bs = "re")






gam_cv <- function(dat, fid, fml, family = binomial(link="cloglog")){
  kvals <- sort(unique(fid))
  
  tss_at <- function(y, p, thr){
    y <- as.integer(y == 1)
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & y == 1); tn <- sum(pred == 0 & y == 0)
    fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
    if ((tp+fn)==0 || (tn+fp)==0) return(NA_real_)
    sens <- tp/(tp+fn); spec <- tn/(tn+fp)
    sens + spec - 1
  }
  
  models  <- vector("list", length(kvals))
  metrics <- vector("list", length(kvals))
  
  for (i in seq_along(kvals)){
    k <- kvals[i]
    tr_idx <- fid != k; te_idx <- fid == k
    train  <- dat[tr_idx, , drop=FALSE]
    test   <- dat[te_idx, , drop=FALSE]
    
    # class-balanced weights on TRAIN subset
    # pres <- sum(train$PA == 1L); absn <- sum(train$PA == 0L)
    # train$w <- if (absn > 0) ifelse(train$PA == 1L, 2, 0.1) else 1
    train$w <- ifelse(train$PA == 1L, 1, 0.1)
    
    m <- mgcv::bam(fml, 
                   data = train, 
                   family = family, 
                   method="fREML", 
                   control =  mgcv::gam.control(maxit = 500, epsilon = 1e-5, nthreads = 6),
                   discrete = FALSE,
                   weights = w,
                   na.action = na.fail,
                   select = TRUE,
                   gamma = 1
    )
    
    # predictions
    p_tr <- predict(m, newdata = train, type = "response")
    p_te <- predict(m, newdata = test,  type = "response")
    
    # TRAIN metrics (AUC + best threshold by Youden)
    tr_auc <- tryCatch(as.numeric(pROC::roc(train$PA, p_tr, quiet=TRUE)$auc),
                       error = function(e) NA_real_)
    co_tr  <- tryCatch(pROC::coords(pROC::roc(train$PA, p_tr, quiet=TRUE),
                                    x="best", best.method="youden",
                                    ret=c("threshold","sensitivity","specificity")),
                       error = function(e) c(threshold=NA, sensitivity=NA, specificity=NA))
    thr_tr <- as.numeric(co_tr["threshold"])
    tr_tss <- if (is.na(thr_tr)) NA_real_ else tss_at(train$PA, p_tr, thr_tr)
    
    # TEST metrics
    te_auc <- tryCatch(as.numeric(pROC::roc(test$PA, p_te, quiet=TRUE)$auc),
                       error = function(e) NA_real_)
    # (i) Test TSS at the TRAIN-chosen threshold (proper CV)
    te_tss_at_tr <- if (is.na(thr_tr)) NA_real_ else tss_at(test$PA, p_te, thr_tr)
    # (ii) (optional) Test-optimal TSS (not used for selection typically)
    co_te  <- tryCatch(pROC::coords(pROC::roc(test$PA, p_te, quiet=TRUE),
                                    x="best", best.method="youden",
                                    ret=c("threshold","sensitivity","specificity")),
                       error = function(e) c(threshold=NA, sensitivity=NA, specificity=NA))
    te_tss_opt <- if (any(is.na(co_te))) NA_real_
    else as.numeric(co_te["sensitivity"] + co_te["specificity"] - 1)
    
    models[[i]] <- m
    metrics[[i]] <- data.frame(
      fold = k,
      train_auc = tr_auc,
      test_auc  = te_auc,
      train_tss = tr_tss,
      test_tss  = te_tss_at_tr,     # <-- this mirrors rigorous CV (threshold set on train)
      test_tss_opt = te_tss_opt,    # <-- optional: test-optimised TSS (for info)
      thr_train = thr_tr
    )
  }
  
  metrics <- do.call(rbind, metrics)
  summary <- within(data.frame(
    train_auc = mean(metrics$train_auc, na.rm=TRUE),
    test_auc  = mean(metrics$test_auc,  na.rm=TRUE),
    train_tss = mean(metrics$train_tss, na.rm=TRUE),
    test_tss  = mean(metrics$test_tss,  na.rm=TRUE)
  ), diff_tss <- train_tss - test_tss)
  
  list(models = models, metrics = metrics, summary = summary)
}

dat <- model_dt_trans         
fid <- sp_blocks$folds_ids     

cv <- gam_cv(dat, fid, fml)
cv

cv$metrics   
cv$summary   


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
                      const = (data.frame(month = 7)),
                      exclude = c( "s(id)"),
                      na.rm = TRUE
)


map

par(mfrow = c(1,1))
plot(map, range = c(0, 1))
terra::plot(map, range = c(0, 1), xlim = c(140, 170), ylim = c(-35, 0))





# Monthlky predictions ----------------------------------------------------

saveRDS(gam, "models/objects/Tracks_GAMM_final_Model_mp.rds")


model <- gam
relevant_vars <- all.vars(model$formula)[-1] # Excluding the intercept
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general")) ]
relevant_vars
#for maxent
# relevant_vars <- c("Depth", "Slope", "dist2000", "CM_uv", "MUR_SST", "Wz", "MLD", "Month", "Year")

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


writeRaster(monthly_predictions_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_2018_2025_rev_mp.tif", overwrite = TRUE)




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
writeRaster(monthly_means_stack, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_rev_mp.tif", overwrite = TRUE)
# monthly_means_stack <- terra::rast("//Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracking_GAMM_monthly_rev_mp.tif")

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
writeRaster(seasonal_means_stack_2, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp.tif", overwrite = TRUE)




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



writeRaster(seasonal_means_stack_4, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons4_means_rev_mp.tif", overwrite = TRUE)



# MESS values  ------------------------------------------------------------

library(dismo)

# Use the same predictors as in your model
env_train <- model_data[, c("thetao", "uv", "wz", "chl", "depth", "slope", "dist2000")]

predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_ext_log.tif")

names(predictors_mean)
plot(predictors_mean)


relevant_vars <- all.vars(model$formula)[-1] # Excluding the intercept
model$formula
relevant_vars
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general", "year", "month", "id")) ]


predictors <- predictors_mean[[relevant_vars]]
predictors
names(predictors)
plot(predictors)
predictors_r <- raster::stack(predictors)
predictors_r

# Compute MESS
mess_map <- dismo::mess(predictors_r, env_train)
mess_map <- terra::rast(mess_map)
plot(mess_map)

extrapolated_mask <- mess_map < 0
plot(extrapolated_mask)
plot(map)

print(extrapolated_mask)
print(map)

# Convert to data frames
map_df <- terra::as.data.frame(map, xy = TRUE)
map_df <- terra::as.data.frame(seasonal_means_stack_2, xy = TRUE)
names(map_df)[3] <- "pred"

extrap_df <- terra::as.data.frame(extrapolated_mask, xy = TRUE)
names(extrap_df)[3] <- "mess"

# keep only extrapolated (TRUE) cells
extrap_df <- subset(extrap_df, mess)

ggplot() +
  geom_raster(data = map_df, aes(x = x, y = y, fill = pred)) +
  geom_raster(data = extrap_df, aes(x = x, y = y), fill = "red", alpha = 0.1) +
  scale_fill_viridis_c(name = "Suitability") +
  labs(title = "Prediction with MESS extrapolation mask") +
  coord_sf(expand = FALSE) +
  theme_bw()


writeRaster(mess_map, filename = "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp.tif", overwrite = TRUE)

saveRDS(extrap_df, "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/MESS_map_df.RDS")

