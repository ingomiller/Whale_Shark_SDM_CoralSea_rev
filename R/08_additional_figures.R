
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
library(rasterVis)
library(tidyterra)
library(ggspatial)
library(sf)
library(basemaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(grid)
library(patchwork)


predictors_mean <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/Predictor_Rasters_rev/mean_month_predictor_stack_0.1.tif")
names(predictors_mean)

plot(predictors_mean)

predictors_mean


predictors_mean[["chl"]] <- predictors_mean[["chl"]] |>
  log()















# Model Data plot  --------------------------------------------------------




locs <- readRDS( "Model_Data_final_R1.rds")


locs_PA <- locs |> 
  dplyr::mutate(PA = as.factor(PA)) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 



PA_plot <- ggplot() +
  geom_sf(data = world, fill = "grey50", color = "grey50") +
  # 1) Absences first (drawn underneath)
  ggplot2::geom_sf(
    data  = dplyr::filter(locs_PA, PA == "0"),
    ggplot2::aes(color = PA),
    size  = 0.01,
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  
  # 2) Presences last (drawn on top)
  ggplot2::geom_sf(
    data  = dplyr::filter(locs_PA, PA == "1"),
    ggplot2::aes(color = PA),
    size  = 0.03,
    shape = 21,
    alpha = 0.9,
    show.legend = FALSE   # legend comes from first layer
  ) +
  ggplot2::scale_color_manual(
    values = c("0" = "dodgerblue", "1" = "firebrick"),
    labels = c("0" = "Pseudo-absence", "1" = "Presence"),
    breaks = c("0", "1"),
    drop   = FALSE,
    name   = ""
  ) +
  coord_sf(
    xlim = range(sf::st_coordinates(locs_PA)[, 1]),
    ylim = range(sf::st_coordinates(locs_PA)[, 2])
  ) +
  theme_bw() +
  theme(
    # Legend inside, compact
    legend.position      = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 5),
    legend.key.size      = unit(1, "mm"),
    legend.margin        = margin(0, 0, 0, 0, "mm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.border      = ggplot2::element_blank(),
    
    
    # Tight plot margins so the 5 cm aren’t wasted
    plot.margin = margin(1, 1, 1, 1, "mm")
  )

PA_plot


# Save as ~5 cm square, high resolution
ggplot2::ggsave(
  filename = "PA_plot.png",
  path ="outputs/figures",
  plot     = PA_plot,
  width    = 5,       # cm
  height   = 5,       # cm
  units    = "cm",
  dpi      = 300,
  bg = "transparent",
  scale =1
)





# Sightings Data plot  ----------------------------------------------------

sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")


sight_sf <-sight |> 
  dplyr::filter(PA == 1) |>  
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 



sight_plot <- ggplot() +
  geom_sf(data = world, fill = "grey50", color = "grey50") +
  # 1) Absences first (drawn underneath)
  ggplot2::geom_sf(
    data  = sight_sf,
    size  = 0.5,
    # shape = 21,
    alpha = 0.5,
    color = "red",
    show.legend = FALSE
  ) +
  
  coord_sf(
    xlim = c(142, 170),
    ylim = c(-40, -5)
  ) +
  theme_bw() +
  theme(
    # Legend inside, compact
    legend.position      = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 5),
    legend.key.size      = unit(1, "mm"),
    legend.margin        = margin(0, 0, 0, 0, "mm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.border      = ggplot2::element_blank(),
    
    
    # Tight plot margins so the 5 cm aren’t wasted
    plot.margin = margin(1, 1, 1, 1, "mm")
  )

sight_plot


# Save as ~5 cm square, high resolution
ggplot2::ggsave(
  filename = "Sight_plot.png",
  path ="outputs/figures",
  plot     = sight_plot,
  width    = 5,       # cm
  height   = 5,       # cm
  units    = "cm",
  dpi      = 300,
  bg = "transparent",
  scale =1
)






# Predictor maps  ---------------------------------------------------------





locs <- locs |> dplyr::filter(PA==1) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

world <- ne_countries(scale = 10, returnclass = "sf")

ext <- predictors_mean |>
  terra::ext()

sst_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["thetao"]]
  ) +
  ggplot2::scale_fill_viridis_c(
    name  = "SST (°C)"
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  # Locations on top
  ggplot2::geom_sf(
    data        = locs,
    # mapping     = ggplot2::aes(color = "white"),
    size        = 0.7,
    color = "firebrick",
    alpha       = 0.9,
    shape = 21,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

sst_plot


ggplot2::ggsave(
  filename = "outputs/figures/sst_locations_transparent.png",
  plot     = sst_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)





library(rayshader)


sst_rast <- predictors_mean[["thetao"]] |>
  setNames("thetao")
elmat <- sst_rast |>
  terra::as.matrix(wide = TRUE)

# rayshader expects [1,1] at bottom-left, so flip rows:
elmat <- elmat[base::nrow(elmat):1, ]

# Scale SST to 0–1 for heightmap (optional but nice)
elmat_scaled <- (elmat - base::min(elmat, na.rm = TRUE)) /
  (base::max(elmat, na.rm = TRUE) - base::min(elmat, na.rm = TRUE))

elmat_scaled[base::is.na(elmat_scaled)] <- 0


# 2) Build hillshade/texture (0–1 RGB array) ------------------------------

tex <- elmat_scaled |>
  rayshader::height_shade(
    texture = viridis::viridis(256)
  )


# 3) Plot in 3D -----------------------------------------------------------

rayshader::plot_3d(
  hillshade  = tex,          
  heightmap  = elmat_scaled, 
  zscale     = 50,
  fov        = 0,
  theta      = 315,
  phi        = 35,
  zoom       = 0.7,
  windowsize = c(800, 800)
)

# Save a snapshot
rgl::rgl.snapshot("sst_3d.png")




chl_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["chl"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

chl_plot


ggplot2::ggsave(
  filename = "outputs/figures/chl_map.png",
  plot     = chl_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)





uv_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["uv"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

uv_plot


ggplot2::ggsave(
  filename = "outputs/figures/uv_map.png",
  plot     = uv_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)



mld_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["mltost"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

mld_plot


ggplot2::ggsave(
  filename = "outputs/figures/mld_map.png",
  plot     = mld_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)





wz_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["wz"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

wz_plot


ggplot2::ggsave(
  filename = "outputs/figures/wz_map.png",
  plot     = wz_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)




depth_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["depth"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

depth_plot


ggplot2::ggsave(
  filename = "outputs/figures/bathy_map.png",
  plot     = depth_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)



dist_plot <- ggplot2::ggplot() +
  # SST raster
  tidyterra::geom_spatraster(
    data = predictors_mean[["dist2000"]]
  ) +
  ggplot2::scale_fill_viridis_c(
  ) +
  # Coastline
  ggplot2::geom_sf(
    data        = world,
    fill        = "grey40",
    colour      = "grey20",
    linewidth   = 0.2,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    xlim   = c(142, 170),
    ylim   = c(-30, 0),
    expand = FALSE
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid      = ggplot2::element_blank(),
    axis.title      = ggplot2::element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    legend.key.size = grid::unit(3, "mm"),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    panel.border      = ggplot2::element_blank(),
    plot.margin     = ggplot2::margin(2, 2, 2, 2, "mm")
  )

dist_plot


ggplot2::ggsave(
  filename = "outputs/figures/dist_map.png",
  plot     = dist_plot,
  width    = 10,
  height   = 10,
  units    = "cm",
  dpi      = 300,
  bg       = "transparent"
)




# Overview Map ------------------------------------------------------------



cities <- data.frame(Loc = c("Cooktown", "Cairns",  "Mackay", "Brisbane", "Sydney"),
                     Group = c("Town", "Town",  "Town", "Town", "Town"),
                     # Season = c("Monsoon Season (Nov - Apr)"),
                     # lyr = 1,
                     lat = c(-15.4758, -16.918246, -21.1434, -27.4705, -33.911609),
                     lon = c(145.2471, 145.771359, 149.1868, 153.026, 151.186715)) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

cities2 <- data.frame(Loc = c( "Port\nMoresby"),
                      lat = c(  -9.4790),
                      lon = c( 147.1494)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

WB <- data.frame(Loc = c( "Wreck\nBay"),
                 lat = c(  -12.132504),
                 lon = c( 143.893818)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 

str(cities)





australia_map_data <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(name == "Australia" & name == "Papua")
print(australia_map_data)

australia_map_data

world <- ne_countries(scale = 10, returnclass = "sf")

AUS_PNG <- world |> filter(name == "Australia" | name == "Papua New Guinea" | name == "Indonesia")
plot(AUS_PNG)

# Define the extent of the study area 
ext <- c(xmin = 140, xmax = 170, ymin = -40, ymax = 0)

# Create the inset map of Australia + Papua New Guinea and add the study area rectangle
AUS_map <- ggplot() +
  geom_sf(data = world, fill = "grey25", colour = NA, size = 0.1) +
  coord_sf(xlim = c(112, 170), ylim = c(-44, 1), expand = FALSE) +
  geom_rect(aes(xmin = ext["xmin"], xmax = ext["xmax"], ymin = ext["ymin"], ymax = ext["ymax"]), 
            fill = NA, color = "#0A3A5A", 
            linetype = "solid", size = 0.5) +  # Study area rectangle
  theme_void() +  
  theme(panel.border = element_blank(),
        panel.background = element_blank())  

AUS_map



GBR_zone <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/GBR_Zones/Great_Barrier_Reef_Marine_Park_Boundary.shp")


GBR_boundary_dt <- st_boundary(GBR_zone) %>%
  st_cast("LINESTRING") %>%
  st_coordinates() %>%
  as.data.frame()

GBR_boundary_dt



GBR_features <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/GBR_Zones/Great_Barrier_Reef_Features.shp")

GBR_features


iucn_range <- sf::st_read("/Volumes/Ingo_PhD/PhD_Data_Analysis/IUCN_whaleshark_range/data_0.shp")

sz <- 2.5

Fig1 <- ggplot() +
  
  geom_sf(data = iucn_range, fill = "steelblue1", color = NA, size = 0.1, alpha = 0.25) +  
  
  geom_sf(data = GBR_features |> dplyr::filter(FEAT_NAME == "Reef"), fill = "gray25", color = NA, size = 0.1) +  
  # add GBRMPA area
  geom_sf(data = GBR_zone, fill = NA, color = "gray25", size = 0.1) +  
  
  
  
  
  geom_sf(data = world, fill = "lightgray", colour = "gray25", size = 0.1) +
  
  
  
  geom_sf(data = cities,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 1, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = -2.5, 
                                nudge_y = 0.25, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 10) +
  
  
  geom_sf(data = cities2,
          shape = 21,
          colour = "black", 
          fill = "yellow", 
          alpha = 1, 
          size = 1.5,
          show.legend = FALSE) +
  coord_sf(xlim = c(140, 170), ylim = c(-40, 0), expand = FALSE) + 
  
  ggsflabel::geom_sf_text_repel(data = cities2,
                                colour = "black", 
                                aes(label = Loc), 
                                nudge_x = 1, 
                                nudge_y = 0, 
                                size = 2, 
                                force = 1,
                                force_pull = 10,
                                seed = 12) +
  
  ggsflabel::geom_sf_text_repel(data = WB,
                                colour = "black",
                                aes(label = Loc),
                                nudge_x = 2.5,
                                nudge_y = -0.5,
                                size = 2,
                                #fontface = "bold",
                                force = 1,
                                force_pull = 10,
                                seed = 15) +
  
  
  
  geom_sf(data = tag.locs,
          shape = 25,
          colour = "black", 
          fill = "red", 
          alpha = 1, 
          size = 1.5,
          show.legend = FALSE) +
  
  ggsflabel::geom_sf_text_repel(data = tag.locs,
                                colour = "red", 
                                aes(label = Loc_short), 
                                nudge_x = .1, 
                                nudge_y = 0, 
                                size = 2.5,
                                force = 1,
                                seed = 10) +
  
  ## Adding Area of Intrest for predictions 
  geom_rect(aes(xmin = ext["xmin"], xmax = ext["xmax"], ymin = ext["ymin"], ymax = ext["ymax"]), 
            fill = NA, color = "#0A3A5A", 
            linetype = "solid", linewidth = 0.5) +  # Study area rectangle
  
  
  coord_sf(xlim = c(138, 172), ylim = c(-42, 2), expand = FALSE) + 
  labs( x = "Longitude", y = "Latitude", title = "") +
  
  scale_y_continuous(limits = c(-41, 1), breaks = seq(-40, 0, by = 10), expand = c(0,0)) +
  scale_x_continuous(limits = c(139, 171), breaks = seq(140, 170, by = 10),expand = c(0,0)) +
  
  # Add the expanded inset map as an annotation
  annotation_custom(
    grob = ggplotGrob(AUS_map),  # Convert the inset map to a grob
    xmin = 141, xmax = 150, ymin = -33, ymax = -25  # Adjust position & size of inset map
  ) +
  
  # Add manual text at specific coordinates
  annotate("text", x = 144, y = -29.5, label = "AUS", color = "black", size = sz, fontface = "bold") +
  annotate("text", x = 144, y = -6, label = "PNG", color = "black", size = sz, fontface = "bold") +
  annotate("text", x = 165, y = -21, label = "NC", color = "black", size = sz, fontface = "bold") +
  annotate("text", x = 158.5, y = -8, label = "SI", color = "black", size = sz, fontface = "bold") +
  annotate("text", x = 166.75, y = -15.5, label = "VU", color = "black", size = sz, fontface = "bold") +
  annotate("text", x = 152, y = -15, label = "Coral Sea", color = "blue4", size = sz, fontface = "italic") +
  annotate("text", x = 154, y = -8.5, label = "Solomon\nSea", color = "blue4", size = sz, fontface = "italic")+
  # annotate("text", x = 139, y = -14, label = "Gulf of\nCarpentaria", color = "blue4", size = 3, fontface = "italic")+
  annotate("text", x = 160, y = -35, label = "Tasman Sea", color = "blue4", size = sz, fontface = "italic")+
  annotate("text", x = 145, y = -9, label = "Gulf of\nPapua", color = "blue4", size = sz, fontface = "italic")+
  annotate("text", x = 142.5, y = -10, label = "Torres Str.", color = "blue4", size = 2, fontface = "italic")+
  
  ggspatial::annotation_north_arrow(location = "br",
                                    which_north = "true", 
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"),
                                    pad_x = unit(1, "cm"),
                                    pad_y = unit(1.3, "cm"),
                                    style =  north_arrow_fancy_orienteering) +
  
  guides(alpha = "none") +
  
 
  
  # Adding sclase bar
  ggspatial::annotation_scale(location = "br", 
                              width_hint = 0.15, 
                              tick_height = 0.2,
                              text_cex = 0.7, 
                              style = "ticks", 
                              pad_x = unit(2.1, "cm"),
                              pad_y = unit(1.3, "cm"),
                              line_width = 1,
  ) +
  
  guides(alpha = "none") +
  
  
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "lightblue"),  # ocean color
    panel.grid.major = element_line(color = "gray", linetype = "dotted"), 
    panel.grid.minor = element_line(color = "gray", linetype = "dotted"),
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.line = element_blank(),
    #panel.background = element_blank(),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 



Fig1


ggsave("Study_Area_map.png", plot = Fig1, path ="/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/Output_Plots", scale =2, width = 14, height = 17, units = "cm", dpi = 300)


ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_1.png", plot = Fig1, scale = 1, width = 14, height = 17, units = "cm", dpi = 600)

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_1.pdf", plot = Fig1, scale =1, width = 14, height = 17, units = "cm", dpi = 600, device = "pdf")

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_1.eps", plot = Fig1, scale =1, width = 14, height = 17, units = "cm", dpi = 600)






# Predicotr Importance ----------------------------------------------------






# Variable Importance  ----------------------------------------------------

BRT_Var_Importance_CV <- readRDS("~/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM_rev/models/tables/BRT_Var_Importance_CV.rds")

MaxEnt_Var_Importance_CV <- readRDS("~/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM_rev/models/tables/MaxEnt_Var_Importance_CV.rds")

GAMM_Var_Importance_CV <- readRDS("~/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/Data_Analysis/R_workfolder/WhaleSharks_SDM_rev/models/tables/GAMM_Var_Importance_CV.rds")


str(BRT_Var_Importance_CV)

BRT_Var_Importance_CV
str(MaxEnt_Var_Importance_CV)
str(GAMM_Var_Importance_CV)

brt_imp <- BRT_Var_Importance_CV |> 
  dplyr::select(-variable) |> 
  dplyr::rename(variable = importance,
                importance = Permutation_importance)

max_imp <- MaxEnt_Var_Importance_CV |> 
  dplyr::select(-variable) |> 
  dplyr::rename(variable = importance,
                importance = Permutation_importance)

unique(brt_imp$variable)

sum_var_imp <- dplyr::bind_rows(brt_imp, max_imp, GAMM_Var_Importance_CV) |> 
  dplyr::mutate(
    variable = case_when(
      variable == "mltost" ~ "mld",
      variable == "thetao" ~ "sst",
      TRUE ~ variable  # Keep all other values unchanged
    ),
    Group = dplyr::case_when(
      variable == "month" ~ "temporal",
      variable == "depth" ~ "static",
      variable == "slope" ~ "static",
      variable == "dist2000" ~ "static",
      TRUE ~ "dynamic"
    )
  )

str(sum_var_imp)

unique(sum_var_imp$variable)

sum_var <- sum_var_imp |> 
  dplyr::group_by(algorithm) |> 
  dplyr::mutate(panel_fill = importance / max(importance),
                ymin = pmax(0, importance - sd),   # don’t go below 0
                ymax = importance + sd) |>   # Scale within each facet
  dplyr::ungroup() |> 
  dplyr::mutate(algorithm = factor(algorithm, levels = c("BRT", "MaxEnt", "GAMM")),
                variable = factor(variable, levels = c("depth", "slope", "dist2000", "sst", "chl", "uv", "wz", "mld", "month")))

sum_var

rel_inf_plot <- ggplot(sum_var, aes(x = variable, y = importance, fill = panel_fill)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,                       # whisker width
    linewidth = 0.3
  ) +
  coord_flip() +
  labs(x = "Predictor", y = "Relative Influence (%)") +
  scale_y_continuous(expand = c(0, 1), limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  scale_fill_gradientn(colors = scales::brewer_pal(palette = "Blues")(9), 
                       name = "Relative Influence",
                       limits = c(0, 1)) +
  # geom_text(aes(label = round(importance, 1)), hjust = -0.1, size = 2.5, fontface = "italic") +  # Add labels slightly outside bars
  facet_grid(cols = vars(algorithm), scales = "free") +
  theme_bw() +
  theme(
    #lot.margin = unit(c(0.5, 0.5, 0.5, 0.1), "cm"),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8))


rel_inf_plot



ggsave("RelInfluence_final_models.png", plot = rel_inf_plot, path ="outputs/final_figures/", scale =1, width = 17, height = 12, units = "cm", dpi = 300)


ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models.png", plot = rel_inf_plot, scale =1, width = 17, height = 12, units = "cm", dpi = 600)

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models.pdf", plot = rel_inf_plot, scale =1, width = 17, height = 12, units = "cm", dpi = 600, device = "pdf")

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models.eps", plot = rel_inf_plot, scale =1, width = 17, height = 12, units = "cm", dpi = 600)






sum_group <- sum_var_imp |> 
  dplyr::group_by(algorithm, Group) |> 
  dplyr::summarise(importance_g = sum(importance)) |> 
  dplyr::mutate(
    panel_fill = importance_g / max(importance_g)) |>   # Scale within each facet
  dplyr::ungroup() |> 
  dplyr::mutate(algorithm = factor(algorithm, levels = c("BRT", "MaxEnt", "GAMM")),
                Group = factor(Group, levels = c("static", "dynamic", "temporal")))

sum_group

rel_inf_groups <- ggplot(sum_group, aes(x = Group, y = importance_g, fill = panel_fill)) +
  geom_bar(stat = "identity", width = 0.4) +
  coord_flip() +
  labs(x = "Predictor Group", y = "Relative Influence (%)") +
  scale_y_continuous(expand = c(0, 1.05), limits = c(0, 105), breaks = seq(0, 90, by = 20)) +
  scale_fill_gradientn(colors = scales::brewer_pal(palette = "Greens")(3),
                       name = "Relative Influence",
                       limits = c(0, 1)) +
  # geom_text(aes(label = round(rel_inf_g, 1)), hjust = -0.05, size = 2.5, fontface = "italic") +  # Add labels slightly outside bars
  facet_grid(cols = vars(algorithm), scales = "free") +
  theme_bw() +
  theme(
    #lot.margin = unit(c(0.5, 0.5, 0.5, 0.1), "cm"),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8))


rel_inf_groups


ggsave("RelInfluence_final_models_groups.png", plot = rel_inf_groups, path ="outputs/final_figures/", scale =1, width = 15, height = 10,  units = "cm", dpi = 300)


ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models_groups.png", plot = rel_inf_groups, scale =1, width = 15, height = 10, units = "cm", dpi = 600)

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models_groups.pdf", plot = rel_inf_groups, scale =1, width = 15, height = 10, units = "cm", dpi = 600, device = "pdf")

ggsave("/Users/ingo/Library/CloudStorage/OneDrive-JamesCookUniversity/02_PhD/06_Chapters/DataChapters/Chapter2_WhaleSharks_Mantas/00_Final_Manuscript_Files/Revision_1/Figure_S__RelImp_models_groups.eps", plot = rel_inf_groups, scale =1, width = 15, height = 10, units = "cm", dpi = 600)




#  external sightings map  ------------------------------------------------



sight <- readRDS("data/processed/Sightings_PA_w_dynSDM_10_2010_2025_extract_processed.rds")
# sight <- readRDS("data/work_files/Sightings_Validation_data_extract_full.rds")
str(sight)
sight <- sight  |>  dplyr::filter(lon > 140) 


mapview::mapview(sight |> dplyr::select(-Date))



locs_PA <- sight |> 
  dplyr::mutate(PA = as.factor(PA)) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) 



PA_plot <- ggplot() +
  geom_sf(data = world, fill = "grey50", color = "grey50") +
  # 1) Absences first (drawn underneath)
  ggplot2::geom_sf(
    data  = dplyr::filter(locs_PA, PA == "0"),
    ggplot2::aes(color = PA),
    size  = 1,
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  
  # 2) Presences last (drawn on top)
  ggplot2::geom_sf(
    data  = dplyr::filter(locs_PA, PA == "1"),
    ggplot2::aes(color = PA),
    size  = 1,
    shape = 21,
    alpha = 0.9,
    show.legend = FALSE   # legend comes from first layer
  ) +
  ggplot2::scale_color_manual(
    values = c("0" = "dodgerblue", "1" = "firebrick"),
    labels = c("0" = "Pseudo-absence", "1" = "Presence"),
    breaks = c("0", "1"),
    drop   = FALSE,
    name   = ""
  ) +
  coord_sf(
    xlim = range(sf::st_coordinates(locs_PA)[, 1]),
    ylim = range(sf::st_coordinates(locs_PA)[, 2])
  ) +
  theme_bw() +
  theme(
    # Legend inside, compact
    legend.position      = c(0.7, 0.95),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 5),
    legend.key.size      = unit(1, "mm"),
    legend.margin        = margin(0, 0, 0, 0, "mm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # Transparent backgrounds, no border
    plot.background   = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.background  = ggplot2::element_rect(fill = NA, colour = NA),
    # panel.border      = ggplot2::element_blank(),
    
    
    # Tight plot margins so the 5 cm aren’t wasted
    plot.margin = margin(1, 1, 1, 1, "mm")
  )

PA_plot


# Save as ~5 cm square, high resolution
ggplot2::ggsave(
  filename = "PA_plot_sightings.png",
  path ="outputs/final_figures",
  plot     = PA_plot,
  width    = 10,       # cm
  height   = 10,       # cm
  units    = "cm",
  dpi      = 300,
  bg = "transparent",
  scale =1
)


