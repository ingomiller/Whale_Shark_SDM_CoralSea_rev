#_____________________________________________________________________________
#                        Models: Sightings - GAMM
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
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
sight




model_dt <- sight |>
  dplyr::mutate(year = format(date, "%Y"),
                month = format(date, "%m"),
                month = as.integer(month),
                #Month = as.factor(Month),
                season = case_when(
                  month %in% 10:12 | month %in% 1:3 ~ "Monsoon",
                  month %in% 4:9 ~ "Trade winds"),
                PA = as.integer(as.character(PA)),
                id = as.factor(id),
                year = as.factor(year),
                season = as.factor(season),
                rep = as.integer(rep), 
                depth = ifelse(depth >=0, -10, depth),
                # CHnage distance to km:
                dist2000 = dist2000) |> 
  sf::st_drop_geometry() |> 
  as.data.frame()

str(model_dt)

model_dt |> dplyr::summarise(start = min(Date),
                             end = max(Date))


model_dt |> 
  rstatix::get_summary_stats(depth)




# Explore Data ------------------------------------------------------------


plots <- model_dt

(plots |>  ggplot(aes(x = depth, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)



# Transformations

model_dt_trans <- model_dt |>
  dplyr::mutate(
    depth = sign(depth) * abs(depth)^(1/3), # cube root
    slope = (slope)^(1/3),
    roughness = (roughness)^(1/3),
    uv = (uv)^(1/3),
    chl = log(chl),
    dist2000 = (dist2000)^(1/3),
    mltost = log(mltost))


model_dt_trans |> 
  rstatix::get_summary_stats(dist2000) |> 
  print(n=50)
model_dt |> 
  rstatix::get_summary_stats(dist2000) |> 
  print(n=50)

plots <- model_dt_trans

(plots |>  ggplot(aes(x = depth, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = slope, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = roughness, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = thetao, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = uv, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = chl, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = dist2000, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots |>  ggplot(aes(x = wz, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plots|>  ggplot(aes(x = mltost, colour = Occ_Source)) + geom_histogram(bins = 50) + theme_bw() +
    plot_layout(guides = "collect")
)


na_counts <- model_dt |> 
  dplyr::summarise(across(everything(), ~ sum(is.na(.))))



# Check Colinearity of repdictor variables


model_dt |>
  dplyr::filter(PA == 0) |> 
  dplyr::select(depth, slope, roughness, dist2000, thetao, chl, uv, wz, mltost) |>
  psych::pairs.panels(method = "spearman", cor = TRUE, scale = FALSE, stars = FALSE)

model_dt_trans |>
  # dplyr::filter(PA == 0) |> 
  dplyr::select(depth, slope, roughness, dist2000, thetao, chl, uv, wz, mltost) |>
  psych::pairs.panels(method = "spearman", cor = TRUE, scale = FALSE, stars = FALSE)



# 1) Do presences and buffered absences differ at all (univariate)?
library(pROC)
vars <- c("depth","slope","dist2000","uv","thetao","wz","mltost","chl")
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
k_general = 5

form = PA ~ 
  s(thetao, bs = "cr",  k = k_general) + 
  s(uv, bs = "cr", k = k_general) + 
  s(wz, bs = "cr", k = k_general) +
  s(mltost, bs = "cr", k = 5) +
  s(chl, bs = "cr", k=k_general) +
  #s(depth, bs = "ts", k = k_general) + 
  s(slope, bs = "cr", k = 5) + 
  s(dist2000, bs = "cr", k = 5) +
  #ti(thetao, depth, bs=c('tp','tp'), k=c(6,6)) +   # key interaction
  #ti(chl, depth, bs = c("tp","tp"), k = c(5,5)) +  # <- key
  ti(chl, mltost, bs = c("tp","tp"), k = c(10,10)) +  # <- key
  ti(thetao, mltost, bs = c("tp","tp"), k = c(10,10)) +  # <- key
  ti(chl, thetao, bs = c("tp","tp"), k = c(10,10)) +  # <- key
  #ti(thetao, mltost, bs=c('tp','tp'), k=c(6,6)) +   # joint surface
  #ti(thetao, month, bs=c('cr','cc'), k=c(5,4)) +
  #ti(chl, month, bs=c('cr','cc'), k=c(5,4)) +
  #ti(uv, month, bs=c('cr','cc'), k=c(5,10)) +
  #ti(wz, dist2000, bs=c('cr','cr'), k=c(5,5)) +
  #ti(mltost, month, bs=c('cr','cc'), k=c(6,10)) +
  #ti(lon, lat, bs = c("cr","cr"), k = c(10,10)) +
  #s(year, bs = "re"),
  s(month, bs = "cc", k = 10) 
  


weights_PA <- ifelse(model_dt_trans$PA == 1, 1, 0.5)
# weights_PA

gam <- mgcv::bam(form,
                 family = binomial(link="cloglog"), 
                 method="REML", 
                 control = mgcv::gam.control(maxit = 1000, epsilon = 1e-5),
                 discrete = FALSE,
                 weights = weights_PA,
                 na.action = na.omit,
                 data = model_dt_trans,
                 select = TRUE,
                 gamma = 1.4
)


summary(gam)

visreg::visreg(gam, scale = "response", type = "conditional", gg = TRUE)

plot(gam, pages = 1, scheme = 1, shade = TRUE,
     seWithMean = TRUE, rug = TRUE, scale = 0)

AIC(gam)


predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_transf_month_predictor_stack_0.1_gam.tif")
# predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1_gam.tif")
names(predictors_mean)
plot(predictors_mean)


relevant_vars <- all.vars(gam$formula)[-1] # Excluding the intercept
gam$formula
relevant_vars
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general", "year")) ]


predictors <- predictors_mean[[relevant_vars]]
predictors
names(predictors)
plot(predictors)


map <- terra::predict(predictors,
                      model = gam,
                      type = "response",
                      const = (data.frame(year = 2024)),
                      exclude = "s(year)",
                      na.rm = TRUE
                      )


map


plot(map, range = c(0, 1))

dev.off()





# sample covariates to check concurvity/correlation
dat_samp <- terra::spatSample(predictors, size = 100000, method = "regular",
                              as.df = TRUE, na.rm = TRUE)
# rank-based correlations (robust)
cor_mat <- stats::cor(dat_samp[, c("depth","slope","dist2000","uv","thetao","wz","mltost", "chl")],
                      method = "spearman", use = "pairwise")
print(cor_mat)

cc <- mgcv::concurvity(gam, full = TRUE)

cc_est <- if (base::is.list(cc)) cc$estimate else cc
cc_wst <- if (base::is.list(cc)) cc$worst     else NULL
cc_obs <- if (base::is.list(cc)) cc$observed  else NULL

base::print(cc_est)  # main concurvity matrix (cols explain rows)

cc_simple <- mgcv::concurvity(gam, full = FALSE)
base::print(cc_simple)  # often a named list of vectors







# GAMM/GAM formula (shrinkage to let terms drop if redundant)
fml = PA ~ 
  s(thetao, bs = "cr",  k = k_general) + 
  #s(uv, bs = "cr", k = k_general) + 
  s(wz, bs = "cr", k = k_general) +
  s(mltost, bs = "cr", k = 3) +
  s(chl, bs = "cr", k=k_general) +
  #s(depth, bs = "cr", k = k_general) + 
  s(slope, bs = "cr", k = 5) + 
  s(dist2000, bs = "cr", k = 5) +
  #ti(thetao, depth, bs=c('tp','tp'), k=c(6,6)) +   # key interaction
  #ti(chl, depth, bs = c("cr","cr"), k = c(5,5)) +  # <- key
  #ti(thetao, mltost, bs=c('tp','tp'), k=c(6,6)) +   # joint surface
  ti(thetao, month, bs=c('cr','cc'), k=c(5,4)) +
  ti(chl, month, bs=c('cr','cc'), k=c(5,4)) +
  #ti(uv, month, bs=c('cr','cc'), k=c(5,10)) +
  #ti(mltost, month, bs=c('cr','cc'), k=c(6,10)) +
  #s(c(lon, lat), bs = "cr", k = 10) +
  #s(year, bs = "re"),
  s(month, bs = "cc", k = 4) 



gam_cv <- function(dat, fid, fml, family = binomial()){
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
    pres <- sum(train$PA == 1L); absn <- sum(train$PA == 0L)
    #train$w <- if (absn > 0) ifelse(train$PA == 1L, 1, pres/absn) else 1
    
    m <- mgcv::bam(fml, 
             data = train, 
             family = family, 
             method = "REML",
             select = TRUE,
             #weights = w,
             discrete = FALSE
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

dat <- model_dt_trans         # your transformed data
fid <- sp_blocks$folds_ids      # from blockCV

cv <- gam_cv(dat, fid, fml)

cv$metrics   
cv$summary   


# third fold’s GAM
m3 <- cv$models[[2]]
summary(m3)
plot(m3, pages=1)


model <- m3
class(model)
summary(model)

par(mfrow = c(3,3))
plot(model, pages = 1, scheme = 1, shade = TRUE,
     seWithMean = TRUE, rug = TRUE, scale = 0)








