# Ensemble Models ---------------------------------------------------------


library(terra)
library(tidyverse)
library(rasterVis)




## Themes/Palletes

#Generate a color palette using cmocean
cm_ocean_palette <- cmocean::cmocean(name = "balance", alpha = 1, start = 0.2, end = 0.8, direction = 1)

#for levelplot background
myTheme <- rasterVis::BTCTheme()
myTheme$panel.background$col = 'gray' 





# Tracking Data -----------------------------------------------------------




### Maxent Models 



months_all_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_monthly_2018_2024.tif")
print(months_all_max)


mean_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_mean_climate.tif")
print(mean_max)


annual_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_annual.tif")

monthly_max <-  terra::rast("SDM_whalesharks_Tracking_MAXENT_monthly.tif")


season_4_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_4seasons.tif")

season_2_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_2seasons.tif")
season_2_max


quartals_max <- terra::rast("SDM_whalesharks_Tracking_MAXENT_quartals.tif")

### GAMM Models




months_all_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_monthly_2018_2024.tif")
print(months_all_gam)


mean_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_mean_climate.tif")
print(mean_gam)

annual_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_annual.tif")

monthly_gam <-  terra::rast("SDM_whalesharks_Tracking_GAMM_monthly.tif")


season_4_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_4seasons.tif")

season_2_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_2seasons.tif")

quartals_gam <- terra::rast("SDM_whalesharks_Tracking_GAMM_quartals.tif")

### BRT




months_all_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_monthly_2018_2024.tif")
print(months_all_brt)


mean_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_mean_climate.tif")
print(mean_brt)
plot(mean_brt, range = c(0,1))

annual_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_annual.tif")

monthly_brt <-  terra::rast("SDM_whalesharks_Tracking_BRT_monthly.tif")

season_4_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_4seasons.tif")

season_2_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_2seasons.tif")

quartals_brt <- terra::rast("SDM_whalesharks_Tracking_BRT_quartals.tif")



setwd("~/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM")


### Ensemble 





models_mean <- c(mean_gam, mean_max, mean_brt)
names(models_mean) <- c("GAM", "MaxEnt", "BRT")

plot(models_mean, range = c(0,1))

## Now with a weighted mean model (weighted using model AUCs)
mean_climate_ensemble <- terra::weighted.mean(models_mean, w = c(0.7460, 0.7461, 0.7755))
mean_climate_ensemble_thres <- mean_climate_ensemble > 0.5

plot(mean_climate_ensemble)
plot(mean_climate_ensemble_thres)



season2 <- c(season_2_gam, season_2_max, season_2_brt)
season2
names(season2) <- c("GAM: Monsoon Season", "GAM: Dry Season ", "MaxEnt: Monsoon Season", "MaxEnt: Dry Season", "BRT: Monsoon Season", "BRT: Dry Season ")


plot(season2)

monsoon <- season2[[c(1, 3, 5)]]
dry <- season2[[c(2, 4, 6)]]
monsoon

plot(monsoon)
plot(dry)

weights_t = c(0.7460, 0.7461, 0.7755)

monsoon_ensemble <- terra::weighted.mean(monsoon, w = weights_t)
dry_ensemble <- terra::weighted.mean(dry, w = weights_t)



season_2_ensemble <- c(monsoon_ensemble, dry_ensemble)
names(season_2_ensemble) <- c("Monsoon_Season_Ensemble_Tracking", "Dry_Season_Ensemble_Tracking")
season_2_ensemble


season_2_ensemble_thres <- season_2_ensemble > 0.5

plot(season_2_ensemble, range = c(0,1))
plot(season_2_ensemble_thres)




max_values <- terra::app(season_2_ensemble, fun = "max", na.rm = TRUE)
high <- global(max_values, fun = "max", na.rm = TRUE) |> as.numeric()




season_2_ensemble_t_map <-  rasterVis::levelplot(season_2_ensemble,
                                                 layout=c(2, 1),
                                                 par.settings = myTheme,
                                                 main = "Seasonal whale shark habitat suitability (Tracking Ensemble)",
                                                 names.attr = c("Monsoon (Nov - Apr)", 
                                                                "Dry (May - Oct)"),
                                                 zlim = c(0, high),
                                                 at = seq(0, high, length.out = 200),
                                                 #col.regions =  rev(topo.colors(200)),
                                                 col.regions = cm_ocean_palette, 
                                                 #col.regions =  viridisLite::inferno(200),
                                                 colorkey = list(space = "right", 
                                                                 length = 0.5, 
                                                                 height = 0.75,
                                                                 labels = list(
                                                                   at = c(0, high),  # Positions for "Low" and "High"
                                                                   labels = c("Low", "High"))))  # Labels to use

season_2_ensemble_t_map





## Quartals 

quartals_t <- c(quartals_gam, quartals_max, quartals_brt)
quartals_t
names(quartals_t) <- c("GAM: Q1", "GAM: Q2", "GAM: Q3", "GAM: Q4",
                       "MaxEnt: Q1", "MaxEnt: Q2", "MaxEnt: Q3", "MaxEnt: Q4", 
                       "BRT: Q1", "BRT: Q2", "BRT: Q3", "BRT: Q4")

plot(quartals_t)

q1_t <- quartals_t[[c(1, 5, 9)]]
q2_t <- quartals_t[[c(2, 6, 10)]]
q3_t <- quartals_t[[c(3, 7, 11)]]
q4_t <- quartals_t[[c(4, 8, 12)]]
q1_t

plot(q1_t)
plot(q4_t)

weights_t = c(0.7460, 0.7461, 0.7755)

q1_ensemble_t <- terra::weighted.mean(q1_t, w = weights_t)
q2_ensemble_t <- terra::weighted.mean(q2_t, w = weights_t)
q3_ensemble_t <- terra::weighted.mean(q3_t, w = weights_t)
q4_ensemble_t <- terra::weighted.mean(q4_t, w = weights_t)



quartals_ensemble_t <- c(q1_ensemble_t, q2_ensemble_t, q3_ensemble_t, q4_ensemble_t)
quartals_ensemble_t

names(quartals_ensemble_t) <- c("Q1_Ensemble", "Q2_Ensemble", "Q3_Ensemble", "Q4_Ensemble")
quartals_ensemble_t
quartals_ensemble_thres_t <- quartals_ensemble_t > 0.5

plot(quartals_ensemble_t)
plot(quartals_ensemble_thres_t)




