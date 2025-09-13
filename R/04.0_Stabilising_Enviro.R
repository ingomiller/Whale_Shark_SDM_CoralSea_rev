 


library(tidyverse)
library(rstatix)
library(ggplot2)
library(sf)

sims_env <- readRDS("~/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM/Simulations_50_Stabilising_Test_bathy_dist_SST_CM_uv_Chl_MLD_Wz_FINAL.rds")
sims_env <- readRDS("C:/Users/jc563815/OneDrive - James Cook University/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM/Simulations_50_Stabilising_Test_bathy_dist_SST_CM_uv_Chl_MLD_Wz_FINAL.rds")


sims_env <- readRDS("data/work_files/Tracks_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
sims_env <- readRDS("data/work_files/Tracks_mp_sims_50_raw_2010_2025_bathy_dist_sst_uv.curr_mld_chl.rds")
str(sims_env)


sims_env |> 
  dplyr::mutate(Month = lubridate::month(date)) |> 
  group_by(Month) |> 
  summarise(reps = max(rep))

summary(sims_env$Depth)
table(is.na(sims_env$Depth))
unique(head(sims_env$Depth, 20))




# sight_env <- readRDS("Sightings_PA_50_Bathy_MUR_SST_CM_uv_CM_Chl_dist_MLD_Wz_2024.rds")



sims_env |> 
  dplyr::group_by(PA) |> 
  dplyr::summarise(N = sum(n()))


# sight_env <- sight_env |> 
#   #dplyr::rename(CM_chl = CM_Chl) |> 
#   dplyr::arrange(date, id, rep) |> 
#   dplyr::mutate(dist2000 = dist2000/1000)

sims_env <- sims_env |> 
  tidyr::drop_na(Depth, Slope, Roughness, thetao, chl, mltost, uv) |>
  dplyr::filter(Depth < 0) |>
  dplyr::group_by(id_split) |> 
  dplyr::mutate(
    REP = case_when(
      PA == 0 ~ as.integer(factor(rep, levels = sort(unique(rep[PA == 0])))),
      TRUE    ~ rep
    )
  ) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(dist2000 = dist2000/1000,
                Depth = -Depth) |> 
  as.data.frame()

str(sims_env)

sims_env |> 
  group_by(id) |> 
  summarise(n_abs_rep = n_distinct(REP[PA == 0])) |> 
  print(n = 50)




sims_env |> 
  rstatix::get_summary_stats(chl, thetao, Depth, mltost, Slope, dist2000, uv, Roughness)


# Function to calculate the mean and sd incrementally
calculate_incremental_stats <- function(data) {
  max_reps <- max(data$REP)
  result <- data.frame(Simulation = integer(),
                       mean_uv = numeric(),
                       sd_uv = numeric(),
                       # mean_wz = numeric(),
                       # sd_wz = numeric(),
                       mean_chl = numeric(),
                       sd_chl = numeric(),
                       mean_thetao = numeric(),
                       sd_thetao = numeric(),
                       mean_depth = numeric(),
                       sd_depth = numeric(),
                       mean_dlope = numeric(),
                       sd_dlope = numeric(),
                       mean_roughness = numeric(),
                       sd_roughness = numeric(),
                       mean_dist2000 = numeric(),
                       sd_dist2000 = numeric(),
                       mean_mltost = numeric(),
                       sd_mltost = numeric()
  )
  
  for (i in 1:max_reps) {
    sims_env_filtered <- data |> 
      dplyr::filter(rep <= (i)) # Adjust to include rep 0
    
    stats <- sims_env_filtered |> 
      dplyr::summarise(Simulation = i, # Adjust numbering to start from 2
                mean_uv = mean(uv, na.rm = TRUE),
                sd_uv = sd(uv, na.rm = TRUE),
                # mean_wz = mean(wz, na.rm = TRUE),
                # sd_wz = sd(wz, na.rm = TRUE),
                mean_chl = mean(chl, na.rm = TRUE),
                sd_chl = sd(chl, na.rm = TRUE),
                mean_thetao = mean(thetao, na.rm = TRUE),
                sd_thetao = sd(thetao, na.rm = TRUE),
                mean_depth = mean(Depth, na.rm = TRUE),
                sd_depth = sd(Depth, na.rm = TRUE),
                mean_slope = mean(Slope, na.rm = TRUE),
                sd_slope = sd(Slope, na.rm = TRUE),
                mean_roughness = mean(Roughness, na.rm = TRUE),
                sd_roughness = sd(Roughness, na.rm = TRUE),
                mean_dist2000 = mean(dist2000, na.rm = TRUE),
                sd_dist2000 = sd(dist2000, na.rm = TRUE),
                mean_mltost = mean(mltost, na.rm = TRUE),
                sd_mltost = sd(mltost, na.rm = TRUE)
      )
    
    result <- rbind(result, stats)
  }
  
  return(result)
}

# Calculate the stats incrementally for up to 100 replicates
incremental_stats <- calculate_incremental_stats(sims_env)

str(incremental_stats)
head(incremental_stats)

# Transform the wide table into a long table
sims_plot_dt <- incremental_stats |> 
  pivot_longer(cols = -Simulation,
               names_to = c("stat", "variable"),
               names_pattern = "(mean|sd)_(.*)",
               values_to = "value")


# Print the resulting long data frame
print(sims_plot_dt)
head(sims_plot_dt)



ggplot(sims_plot_dt |> dplyr::filter(stat == "mean"), aes(x = Simulation, y = value, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "blue") + # Add vertical line
  facet_wrap(~variable, scales = "free_y", ncol =2) +
  # Adding secondary y-axis for standard deviation (sd)
  scale_y_continuous(sec.axis = sec_axis(~ . * max(sims_plot_dt$value[sims_plot_dt$stat == "sd"]) / max(sims_plot_dt$value[sims_plot_dt$stat != "sd"]), name = "Standard Deviation")) +
  labs(title = "Tracking data",
       x = "Number of Simulations",
       y = "Value",
       color = "Statistic") +
  theme_bw()

ggplot(sims_plot_dt |> dplyr::filter(stat == "sd"), aes(x = Simulation, y = value, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "blue") + # Add vertical line
  facet_wrap(~variable, scales = "free_y", ncol =2) +
  # Adding secondary y-axis for standard deviation (sd)
  scale_y_continuous(sec.axis = sec_axis(~ . * max(sims_plot_dt$value[sims_plot_dt$stat == "sd"]) / max(sims_plot_dt$value[sims_plot_dt$stat != "sd"]), name = "Standard Deviation")) +
  labs(title = "Tracking data",
       x = "Number of Simulations",
       y = "Value",
       color = "Statistic") +
  theme_bw()


## for sightings 

# Calculate the stats incrementally for up to 100 replicates
incremental_stats <- calculate_incremental_stats(sight_env)

str(incremental_stats)
head(incremental_stats)

# Transform the wide table into a long table
sight_plot_dt <- incremental_stats %>%
  pivot_longer(cols = -Simulation,
               names_to = c("stat", "variable"),
               names_pattern = "(mean|sd)_(.*)",
               values_to = "value")


# Print the resulting long data frame
print(sight_plot_dt)
head(sight_plot_dt)

ggplot(sight_plot_dt, aes(x = Simulation, y = value, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 10, linetype = "dashed", color = "blue") + # Add vertical line
  facet_wrap(~variable, scales = "free_y", ncol =2) +
  labs(title = "Sightings data",
       x = "Number of Simulations",
       y = "Value",
       color = "Statistic") +
  theme_minimal()





### calculate rate of chnage 


# Calculate the range of each variable
range <- sims_plot_dt |> 
  dplyr::group_by(stat, variable) |> 
  dplyr::summarise(range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE), .groups = 'drop')

# Join the ranges with the long data frame
sims_plot_dt2 <- sims_plot_dt |> 
  dplyr::left_join(range, by = c("stat", "variable"))

# Calculate the rate of change and normalize it by the range
sims_plot_dt2 <- sims_plot_dt2 |> 
  dplyr::group_by(stat, variable) |> 
  dplyr::arrange(Simulation) |> 
  dplyr::mutate(rate_of_change = c(NA, diff(value)),
                normalised_rate_of_change = rate_of_change / range) |> 
  dplyr::ungroup()



sims_plot_dt2


# calculate Convergence 

# Find the convergence point for each variable based on the normalized rate of change
# Define a threshold for normalized convergence
normalised_threshold <- 0.05

convergence_points <- sims_plot_dt2 |> 
  dplyr::group_by(stat, variable) |> 
  dplyr::filter(!is.na(normalised_rate_of_change) & abs(normalised_rate_of_change) < normalised_threshold) |> 
  dplyr::summarise(convergence_point = min(Simulation, na.rm = TRUE), .groups = 'drop')

# Find the maximum convergence point across all variables
overall_convergence_point <- max(convergence_points$convergence_point, na.rm = TRUE)
overall_convergence_point <- mean(convergence_points$convergence_point, na.rm = TRUE)

# Print the convergence point
print(convergence_points)

convergence_points |>
  group_by(stat) |>
  dplyr::summarise(max = max(convergence_point),
                   mean = mean(convergence_point))

print(overall_convergence_point)


absences_Plot_track <- ggplot(sims_plot_dt2, aes(x = Simulation, y = normalised_rate_of_change, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "blue") + # Add vertical line for visualization
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  labs(title = "Tracking data: Normalised Rate of Change",
       x = "Number of Simulations",
       y = "Normalised Rate of Change",
       color = "Statistic") +
  # Set custom colors for each stat value
  scale_color_manual(values = c("mean" = "tomato2", "sd" = "skyblue3")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 12))


absences_Plot_track

ggsave("Tracking_Absence_Simualtions_50.png", plot = absences_Plot_track, path ="C:/Users/jc563815/OneDrive - James Cook University/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/Output_Plots", scale =1, width = 17, height = 12, units = "cm", dpi = 300)


## for sightings 

# Calculate the range of each variable
range <- sight_plot_dt |>
  dplyr::group_by(stat, variable) |>
  dplyr::summarise(range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE), .groups = 'drop')

# Join the ranges with the long data frame
sight_plot_dt2 <- sight_plot_dt |>
  dplyr::left_join(range, by = c("stat", "variable"))

# Calculate the rate of change and normalize it by the range
sight_plot_dt2 <- sight_plot_dt2 |>
  dplyr::group_by(stat, variable) |>
  dplyr::arrange(Simulation) |>
  dplyr::mutate(rate_of_change = c(NA, diff(value)),
                normalised_rate_of_change = rate_of_change / range) |>
  dplyr::ungroup()



sight_plot_dt2


# calculate Convergence 

# Find the convergence point for each variable based on the normalized rate of change
# Define a threshold for normalized convergence
normalised_threshold <- 0.05

convergence_points <- sight_plot_dt2 |> 
  dplyr::group_by(stat, variable) |> 
  dplyr::filter(!is.na(normalised_rate_of_change) & abs(normalised_rate_of_change) < normalised_threshold) |> 
  summarize(convergence_point = min(Simulation, na.rm = TRUE), .groups = 'drop')

# Find the maximum convergence point across all variables
overall_convergence_point <- max(convergence_points$convergence_point, na.rm = TRUE)
overall_convergence_point <- mean(convergence_points$convergence_point, na.rm = TRUE)


convergence_points |>
  group_by(stat) |>
  dplyr::summarise(max = max(convergence_point),
                   mean = mean(convergence_point))

# Print the convergence point
print(convergence_points)
print(overall_convergence_point)




absences_Plot_sight <- ggplot(sight_plot_dt2, aes(x = Simulation, y = normalised_rate_of_change, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") + # Add vertical line for visualization
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  labs(title = "Sightings data: Normalised Rate of Change",
       x = "Number of Simulations",
       y = "Normalised Rate of Change",
       color = "Statistic") +
  theme_bw() +
  scale_color_manual(values = c("mean" = "tomato2", "sd" = "skyblue3")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 10))


absences_Plot_sight

ggsave("Sightings_Absence_Simualtions_50.png", plot = absences_Plot_sight, path ="C:/Users/jc563815/OneDrive - James Cook University/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/Output_Plots", scale =1, width = 17, height = 12, units = "cm", dpi = 300)



### in one plot:

sims_plot_dt2$Source <- "Tracking"
sight_plot_dt2$Source <- "Sightings"


plot_dt <- bind_rows(sims_plot_dt2, sight_plot_dt2)
glimpse(plot_dt)

plot_dt$Source <- factor(plot_dt$Source, levels = c("Tracking", "Sightings"))

absence_sims_plot <- ggplot(plot_dt, aes(x = Simulation, y = normalised_rate_of_change, color = stat)) +
  geom_line() +
  geom_vline(xintercept = 10, linetype = "dotted", color = "grey20") + # Add vertical line for visualization
  facet_grid(variable ~ Source, scales = "free_y") +
  labs(x = "Number of Simulations",
       y = "Normalised Rate of Change",
       color = "Statistic") +
  theme_bw() +
  scale_color_manual(values = c("mean" = "tomato2", "sd" = "skyblue3")) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none")

absence_sims_plot

ggsave("Absence_Simualtions_50.png", plot = absence_sims_plot, path ="C:/Users/jc563815/OneDrive - James Cook University/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/Output_Plots", scale =1, width = 17, height = 22, units = "cm", dpi = 300)
