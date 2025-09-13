
#_____________________________________________________________________________
#                          HELPER FUNCTIONS 
#_____________________________________________________________________________




# fix to dynamicSDM::spatiotemp_pseudoabs()  ------------------------------

spatiotemp_pseudoabs_fix <- function(
    spatial.method,
    temporal.method,
    occ.data,
    spatial.ext,
    temporal.ext,
    spatial.buffer,
    temporal.buffer,
    n.pseudoabs = 100,
    prj = "EPSG:4326",      # input CRS (lon/lat by default)
    work_crs = "EPSG:3577"  # metric CRS for buffering/sampling
) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")
  
  temporal.method <- base::match.arg(temporal.method, c("buffer", "random"))
  spatial.method  <- base::match.arg(spatial.method,  c("buffer", "random"))
  if (!is.numeric(n.pseudoabs)) stop("n.pseudoabs should be numeric")
  
  # coerce occ.data to needed cols
  if (!("x" %in% names(occ.data)) && "lon" %in% names(occ.data)) occ.data$x <- as.numeric(occ.data$lon)
  if (!("y" %in% names(occ.data)) && "lat" %in% names(occ.data)) occ.data$y <- as.numeric(occ.data$lat)
  if (!all(c("year","month","day") %in% names(occ.data)) && "date" %in% names(occ.data)) {
    occ.data$year  <- as.integer(format(as.Date(occ.data$date), "%Y"))
    occ.data$month <- as.integer(format(as.Date(occ.data$date), "%m"))
    occ.data$day   <- as.integer(format(as.Date(occ.data$date), "%d"))
  }
  needed <- c("x","y","year","month","day")
  if (!all(needed %in% names(occ.data))) {
    stop("occ.data must contain x, y, year, month, day (or date to derive them).")
  }
  occ.data <- occ.data |>
    dplyr::mutate(
      x = as.numeric(.data$x),
      y = as.numeric(.data$y),
      year = as.integer(.data$year),
      month = as.integer(.data$month),
      day = as.integer(.data$day)
    ) |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$y),
                  !is.na(.data$year), !is.na(.data$month), !is.na(.data$day))
  
  n_occ <- nrow(occ.data)
  if (n_occ == 0) stop("occ.data has no valid rows after cleaning.")
  value <- ceiling(n.pseudoabs / n_occ)  # draw per occurrence, then trim to n.pseudoabs at end
  
  # helper: extent to sf polygon in prj
  as_sf_extent <- function(ext, crs_in) {
    if (missing(ext)) return(NULL)
    if (inherits(ext, "sf") || inherits(ext, "sfc")) return(sf::st_transform(ext, crs = crs_in))
    if (is.numeric(ext) && length(ext) == 4) {
      m <- matrix(c(ext[1], ext[2], ext[3], ext[2], ext[3], ext[4], ext[1], ext[4], ext[1], ext[2]),
                  ncol = 2, byrow = TRUE)
      ply <- sf::st_sfc(sf::st_polygon(list(m)), crs = crs_in)
      return(sf::st_as_sf(sf::st_sf(geometry = ply)))
    }
    stop("spatial.ext must be an sf/sfc polygon or numeric bbox c(xmin,ymin,xmax,ymax) in 'prj' CRS.")
  }
  
  spatial.ext.sf <- try(as_sf_extent(spatial.ext, prj), silent = TRUE)
  if (inherits(spatial.ext.sf, "try-error")) spatial.ext.sf <- NULL
  
  # build occ sf in metric CRS
  occ_sf <- sf::st_as_sf(occ.data, coords = c("x","y"), crs = prj) |>
    sf::st_transform(work_crs)
  if (!is.null(spatial.ext.sf)) spatial.ext.sf <- sf::st_transform(spatial.ext.sf, work_crs)
  
  # SPATIAL
  PA_coords_list <- vector("list", length = n_occ)
  if (spatial.method == "random") {
    if (is.null(spatial.ext.sf)) stop("spatial.ext is required for spatial.method = 'random'.")
    for (i in seq_len(n_occ)) {
      samp <- sf::st_sample(spatial.ext.sf, type = "random", size = value)
      coords <- sf::st_coordinates(samp)
      PA_coords_list[[i]] <- data.frame(x = coords[,1], y = coords[,2])
    }
  } else { # buffer
    if (missing(spatial.buffer)) stop("No spatial.buffer specified for 'buffer' method.")
    if (!is.numeric(spatial.buffer)) stop("spatial.buffer must be numeric (metres).")
    if (length(spatial.buffer) > 2) stop("spatial.buffer must be length 1 or 2.")
    if (length(spatial.buffer) == 1) spatial.buffer <- c(1e-8, spatial.buffer)
    if (spatial.buffer[1] > spatial.buffer[2]) stop("Second spatial.buffer must be greater than first.")
    
    for (i in seq_len(n_occ)) {
      outer <- sf::st_buffer(occ_sf[i,], dist = spatial.buffer[2])
      inner <- sf::st_buffer(occ_sf[i,], dist = spatial.buffer[1])
      ring  <- suppressWarnings(sf::st_difference(outer, inner)) |> sf::st_make_valid()
      if (!is.null(spatial.ext.sf)) {
        ring <- suppressWarnings(sf::st_intersection(ring, spatial.ext.sf))
        if (length(ring) == 0) stop("For at least one occurrence, buffer does not intersect spatial.ext.")
      }
      samp <- sf::st_sample(ring, type = "random", size = value)
      coords <- sf::st_coordinates(samp)
      PA_coords_list[[i]] <- data.frame(x = coords[,1], y = coords[,2])
    }
  }
  PA_coords <- dplyr::bind_rows(PA_coords_list)
  
  # TEMPORAL
  PA_dates_list <- vector("list", length = n_occ)
  if (temporal.method == "random") {
    if (missing(temporal.ext)) stop("temporal.ext is required for temporal.method = 'random'.")
    if (!is.character(temporal.ext) || length(temporal.ext) != 2)
      stop("temporal.ext must be c('YYYY-MM-DD','YYYY-MM-DD').")
    d1 <- as.Date(temporal.ext[1]); d2 <- as.Date(temporal.ext[2])
    if (is.na(d1) || is.na(d2) || d1 > d2) stop("Invalid temporal.ext.")
    for (i in seq_len(n_occ)) {
      rdates <- d1 + sample.int(as.integer(d2 - d1), size = value, replace = TRUE)
      PA_dates_list[[i]] <- data.frame(
        year  = as.integer(format(rdates, "%Y")),
        month = as.integer(format(rdates, "%m")),
        day   = as.integer(format(rdates, "%d"))
      )
    }
  } else { # buffer
    if (missing(temporal.buffer)) stop("temporal.buffer is required for 'buffer' method.")
    if (!is.numeric(temporal.buffer)) stop("temporal.buffer must be numeric.")
    if (length(temporal.buffer) > 2) stop("temporal.buffer must be length 1 or 2.")
    if (length(temporal.buffer) == 1) temporal.buffer <- c(0, temporal.buffer)
    if (temporal.buffer[1] > temporal.buffer[2]) stop("Second temporal.buffer must be >= first.")
    base_dates <- as.Date(with(occ.data, paste(year, month, day, sep = "-")), "%Y-%m-%d")
    if (any(is.na(base_dates))) stop("occ.data year/month/day produced invalid dates.")
    for (i in seq_len(n_occ)) {
      center <- base_dates[i]
      left  <- seq(center - temporal.buffer[2], center - temporal.buffer[1], by = "day")
      right <- seq(center + temporal.buffer[1], center + temporal.buffer[2], by = "day")
      if (!missing(temporal.ext) && !is.null(temporal.ext)) {
        dmin <- as.Date(temporal.ext[1]); dmax <- as.Date(temporal.ext[2])
        if (!is.na(dmin)) { left  <- left[left  >= dmin]; right <- right[right >= dmin] }
        if (!is.na(dmax)) { left  <- left[left  <= dmax]; right <- right[right <= dmax] }
      }
      pool <- unique(c(left, right))
      if (length(pool) == 0) stop("Temporal buffer produced an empty date set.")
      rdates <- sample(pool, size = value, replace = TRUE)
      PA_dates_list[[i]] <- data.frame(
        year  = as.integer(format(rdates, "%Y")),
        month = as.integer(format(rdates, "%m")),
        day   = as.integer(format(rdates, "%d"))
      )
    }
  }
  PA_dates <- dplyr::bind_rows(PA_dates_list)
  
  # COMBINE + TRIM to n.pseudoabs (mimic original behaviour)
  pseudo.df <- dplyr::bind_cols(PA_coords, PA_dates)
  if (nrow(pseudo.df) > n.pseudoabs) {
    pseudo.df <- dplyr::slice_sample(pseudo.df, n = n.pseudoabs)
  }
  
  # coords are in work_crs → convert back to prj
  pts   <- sf::st_as_sf(pseudo.df, coords = c("x","y"), crs = work_crs)
  ptsll <- sf::st_transform(pts, crs = prj)
  cxy   <- sf::st_coordinates(ptsll)
  pseudo.df$x <- cxy[,1]; pseudo.df$y <- cxy[,2]
  
  pseudo.df[, c("x","y","year","month","day")]
}



# # Spatiotemporal thinning -------------------------------------------------
# #### REPLACED BY dynamicSDM fucntion 
# 
# sp_thinning <- function(input_df, distance_degrees = 0.25, verbose = FALSE) {
#   library(geosphere)
#   library(future.apply)
#   library(progressr)
#   
#   sharks <- unique(input_df$id)
#   df.save <- NULL
#   df.deleted <- NULL  # DataFrame to save the deleted rows
#   total_thinned <- 0  # Counter for the number of thinned points
#   
#   # Initialize progress bar
#   p <- progressor(along = sharks)
#   
#   # Define the function to process each shark in parallel
#   process_shark <- function(shark_id) {
#     aux.shark <- subset(input_df, id == shark_id)
#     if (verbose) cat("Processing shark ID:", shark_id, "with", nrow(aux.shark), "rows\n")
#     
#     reps <- unique(aux.shark$rep)  # Get unique reps for each shark
#     shark_save <- NULL
#     shark_deleted <- NULL
#     shark_thinned <- 0
#     
#     for (j in seq_along(reps)) {
#       aux.track <- subset(aux.shark, rep == reps[j])
#       if (verbose) cat("  Processing rep:", reps[j], "with", nrow(aux.track), "rows\n")
#       
#       if (nrow(aux.track) <= 1) {
#         if (verbose) cat("  Skipping track with less than 2 rows for id:", shark_id, "rep:", reps[j], "\n")
#         next
#       }
#       
#       index.save <- 1  # Always include the first point
#       index.base <- 1
#       
#       for (v in 2:nrow(aux.track)) {
#         point1 <- c(aux.track$lon[index.base], aux.track$lat[index.base])
#         point2 <- c(aux.track$lon[v], aux.track$lat[v])
#         dist <- distHaversine(point1, point2)
#         
#         if (verbose) cat("    Distance between points", index.base, "and", v, ":", dist, "meters\n")
#         
#         if (dist > distance_degrees * 111320) {  # Convert degrees to meters (~111.32 km per degree)
#           index.save <- c(index.save, v)
#           index.base <- v
#           if (verbose) cat("    Point", v, "retained with distance", dist, "meters\n")
#         }
#       }
#       
#       thinned_count <- nrow(aux.track) - length(index.save)  # Calculate the number of thinned points
#       shark_thinned <- shark_thinned + thinned_count  # Update the total thinned points
#       
#       aux.stand <- aux.track[index.save, ]
#       aux.deleted <- aux.track[-index.save, ]  # Identify the deleted rows
#       
#       if (verbose) cat("  Thinned data for rep:", reps[j], "- Retained:", nrow(aux.stand), "Deleted:", nrow(aux.deleted), "\n")
#       
#       shark_save <- rbind(shark_save, aux.stand)
#       shark_deleted <- rbind(shark_deleted, aux.deleted)
#     }
#     
#     p()  # Update progress bar
#     return(list(thinned_data = shark_save, deleted_data = shark_deleted, thinned_count = shark_thinned))
#   }
#   
#   # Run the shark processing function in parallel using future.apply
#   results <- future_lapply(sharks, process_shark)
#   
#   # Combine the results
#   df.save <- do.call(rbind, lapply(results, function(res) res$thinned_data))
#   df.deleted <- do.call(rbind, lapply(results, function(res) res$deleted_data))
#   total_thinned <- sum(sapply(results, function(res) res$thinned_count))
#   
#   if (is.null(df.save)) {
#     cat("No valid data after thinning.\n")
#     return(NULL)
#   }
#   
#   cat("Thinning complete. Number of spatially correlated points deleted:", total_thinned,
#       "| Final locations retained:", nrow(df.save), "\n")
#   
#   return(list(thinned_data = df.save, deleted_data = df.deleted))
# }










# CMEMS extractions  ------------------------------------------------------

## Function to extract surface layer values from CMEMS files

CMEMS_extract_by_date_surf <- function(rstack, sf_points, date_col = "Date", varname = "value", buffer_m = 10000) {
  stopifnot(inherits(rstack, "SpatRaster"))
  stopifnot(inherits(sf_points, "sf"))
  stopifnot(date_col %in% names(sf_points))
  
  # 1) normalise dates to match layer names ("YYYY-MM-DD")
  loc_dates_chr <- format(as.POSIXct(sf_points[[date_col]], tz = "UTC"), "%Y-%m-%d")
  lyr_names     <- names(rstack)
  
  # 2) find layer index for each row
  lyr_idx <- match(loc_dates_chr, lyr_names)
  
  if (all(is.na(lyr_idx))) {
    stop("No dates in locations matched any raster layer names", call. = FALSE)
  }
  if (any(is.na(lyr_idx))) {
    warning("Some locations had no matching date in the raster stack; returning NA for those rows.", call. = FALSE)
  }
  
  # 3) ensure geometries align
  v <- terra::vect(sf_points)
  if (!terra::same.crs(v, rstack)) {
    v <- terra::project(v, rstack)
  }
  
  # 4) extract per-layer
  vals <- rep(NA_real_, nrow(sf_points))
  idx_split <- split(seq_len(nrow(sf_points)), lyr_idx, drop = TRUE)
  
  buffer_m <- as.numeric(buffer_m)
  
  for (k in names(idx_split)) {
    li  <- as.integer(k)
    ids <- idx_split[[k]]
    vi  <- v[ids]
    val <- terra::extract(rstack[[li]], vi, method = "bilinear", search_radius = buffer_m)[, 2]
    vals[ids] <- val
  }
  
  # 5) append with dynamic column name
  sf_points[[varname]] <- vals
  sf_points
}




# Bluelink BRAN2020 -------------------------------------------------------

windDir <-function(u, v){
  (270 - atan2(u, v) * 180 / pi) %% 360 
}

windSpd <- function(u, v){
  sqrt(u^2 + v^2)
}

dataDownload <- function(type, year, month = NULL, dir, varname, quiet = TRUE) {
  if (type %in% c("day", "month") == FALSE) {
    stop("Please provide type as 'day' or 'month'.")
  }
  if (is.null(year)) {
    stop("Please provide a year of interest to download data.")
  }  
  if(!varname %in% c('ocean_temp', 'ocean_salt', 'ocean_u', 'ocean_v', 'ocean_w', 'ocean_eta_t', 'ocean_mld', 'atm_flux_diag')){
    stop("Environmental variable not recognised, options include:\n'ocean_temp', 'ocean_salt', 'ocean_u', 'ocean_v', 'ocean_w', 'ocean_eta_t', 'ocean_mld', 'atm_flux_diag'")}
  # Change variable name to download correct netCDF file from BRAN
  if (varname %in% c('air_wind')) {
    varname <- "atm_flux_diag"
  }
  # Increase timeout for large file downloads and slow internet connexions
  options(timeout = 1000000000) 
  # Daily resolution
  if (type == "day") {
    if (year < 2023) {
      if (is.null(month)) {
        month <- 1:12
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      } else {
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      }      
    } else {
      if (is.null(month)) {
        month <- 1:12 # 2023 data now goes to december!
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      } else {
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      }
    }
  }
  # Monthly resolution
  if (type == "month") {
    if (is.null(month)) {
      month <- 1:12
      for (i in 1:length(month)) {
        if (month[i] < 10) {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_0", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        } else {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        }
      }
    } else {
      for (i in 1:length(month)) {
        if (month[i] < 10) {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_0", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        } else {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        }
      }
    }
  }   
}   




extractWz <- function(df,
                      X,
                      Y,
                      datetime,
                      folder_name,                 # directory with files like ocean_w_YYYY_MM.nc
                      export_path = NULL,
                      max_depth = -200,
                      fill_gaps = TRUE,
                      buffer = 10000,             # metres (geodesic with s2)
                      verbose = TRUE,
                      export_step = TRUE) {
  
  # ---- checks ---------------------------------------------------------------
  if (!dir.exists(folder_name)) {
    base::stop(base::sprintf("NetCDF directory not found: %s", folder_name), call. = FALSE)
  }
  if (!("Depth" %in% base::names(df))) {
    base::warning("Depth column not found; returning input df unchanged.")
    return(df)
  }
  
  # clamp depths and prepare aux.date (YYYY-MM)
  df <- df |>
    dplyr::mutate(
      Depth_Wz = base::ifelse(Depth > -5, -6, base::ifelse(Depth < max_depth, max_depth, Depth)),
      Depth_Wz = base::round(Depth_Wz, 1),
      aux.date = base::substr(base::as.character(.data[[datetime]]), 1, 7)
    )
  
  if ("Wz" %in% base::names(df)) {
    if (verbose) base::message("Previous Wz values found. Continuing only for rows with NA Wz...")
    dates <- base::unique(df$aux.date[base::is.na(df$Wz)])
  } else {
    df$Wz <- base::as.numeric(NA)
    dates <- base::unique(df$aux.date)
  }
  if (base::length(dates) == 0L) {
    if (verbose) base::message("Nothing to do (no dates needing computation).")
    return(df |> dplyr::select(-.data$aux.date))
  }
  
  # helper: read time origin/units from NC and convert offsets to POSIXct
  .read_nc_time_info <- function(nc_path) {
    nc <- ncdf4::nc_open(nc_path)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    
    # pick a plausible time dim name
    time_name <- base::intersect(base::names(nc$dim), c("Time", "time", "TIME"))
    if (base::length(time_name) == 0L) time_name <- "time"
    
    units_att <- try(ncdf4::ncatt_get(nc, time_name, "units")$value, silent = TRUE)
    if (base::inherits(units_att, "try-error") || base::is.null(units_att)) units_att <- ""
    
    list(units = units_att)
  }
  
  .convert_offsets_to_datetime <- function(offset_numeric, units_att) {
    # units like "days since 1990-01-01 00:00:00" or "hours since 1979-01-01"
    if (!base::grepl("since", units_att, ignore.case = TRUE)) {
      # fallback: treat as days since 1979-01-01 (common for BRAN/Bluelink)
      origin <- base::as.POSIXct("1979-01-01 00:00:00", tz = "UTC")
      return(origin + base::as.difftime(offset_numeric, units = "days"))
    }
    unit_part <- base::tolower(base::sub(" since .*", "", units_att))
    origin_str <- base::sub(".*since\\s+", "", units_att)
    origin     <- base::as.POSIXct(origin_str, tz = "UTC")
    if (base::is.na(origin)) origin <- base::as.POSIXct(base::paste0(origin_str, " 00:00:00"), tz = "UTC")
    
    if (unit_part %in% c("day", "days")) {
      origin + base::as.difftime(offset_numeric, units = "days")
    } else if (unit_part %in% c("hour", "hours")) {
      origin + base::as.difftime(offset_numeric, units = "hours")
    } else if (unit_part %in% c("sec", "second", "seconds")) {
      origin + base::as.difftime(offset_numeric, units = "secs")
    } else {
      # unknown unit → assume days
      origin + base::as.difftime(offset_numeric, units = "days")
    }
  }
  
  # progress bar
  if (verbose) {
    base::message("Processing Wz from existing NetCDF files...")
    pb <- utils::txtProgressBar(min = 0, max = base::length(dates), initial = 0, style = 3, width = 60)
  }
  
  # main loop over months present in df
  for (i in base::seq_along(dates)) {
    aux.date <- dates[i]
    yyyy <- base::substr(aux.date, 1, 4)
    mm   <- base::substr(aux.date, 6, 7)
    
    nc_file <- base::file.path(folder_name, base::sprintf("ocean_w_%s_%s.nc", yyyy, mm))
    if (!base::file.exists(nc_file)) {
      base::stop(base::sprintf("Missing NetCDF for %s: %s", aux.date, nc_file), call. = FALSE)
    }
    
    # read NC as SpatRaster
    nc.bran <- terra::rast(nc_file)
    
    # parse depth & raw time offsets from layer names: "w_sw_ocean=X_Time=Y"
    nm <- terra::names(nc.bran) |>
      stringr::str_remove("w_sw_ocean=") |>
      stringr::str_remove("Time=") |>
      stringr::str_split("_")
    
    aux.depth <- base::numeric(0)
    aux.time  <- base::numeric(0)
    for (k in base::seq_along(nm)) {
      aux.depth <- base::c(aux.depth, base::as.numeric(nm[[k]][1]))
      aux.time  <- base::c(aux.time,  base::as.numeric(nm[[k]][2]))
    }
    
    # read time units/origin from the file and convert offsets → POSIXct
    time_info <- .read_nc_time_info(nc_file)
    meta_time <- .convert_offsets_to_datetime(aux.time, time_info$units)
    
    # meta table: depth (negative down), time (POSIXct), plus layer index
    meta <- dplyr::tibble(
      Depth = -1 * aux.depth,
      Time  = meta_time,
      lyr   = base::seq_along(aux.depth)
    )
    
    # rows for this month
    idx_month <- base::which(df$aux.date == aux.date)
    
    for (ii in base::seq_along(idx_month)) {
      r <- idx_month[ii]
      dt_i <- base::as.Date(df[[datetime]][r])
      
      # match by calendar day (ignore time-of-day in NC)
      sel <- meta |>
        dplyr::filter(base::as.Date(.data$Time) == dt_i)
      
      if (nrow(sel) == 0L) {
        # no matching time slice → 0
        df$Wz[r] <- 0
        next
      }
      
      # keep layers down to local bottom (Depth_Wz is negative)
      sel <- sel |>
        dplyr::filter(.data$Depth >= df$Depth_Wz[r])
      
      if (nrow(sel) == 0L) {
        df$Wz[r] <- 0
        next
      }
      
      # sort shallow (near 0) → deeper (more negative)
      sel <- sel |> dplyr::arrange(dplyr::desc(.data$Depth))
      
      # extract W at each retained layer
      sel$W <- NA_real_
      for (jj in base::seq_len(nrow(sel))) {
        aux.layer <- nc.bran[[sel$lyr[jj]]]
        
        # point extract
        val <- terra::extract(x = aux.layer, y = df[r, c(X, Y)])[1, 2]
        
        # gap-fill via buffer (geodesic metres with s2)
        if (base::is.na(val) && fill_gaps) {
          pos_sf <- df[r, c(X, Y)] |>
            sf::st_as_sf(coords = base::c(1, 2), crs = 4326, remove = FALSE) 
          
          # Silence "dist is assumed to be in decimal degrees (arc_degrees)." message
          pos_sf <- suppressMessages(sf::st_buffer(pos_sf, buffer))
          
          buf_vals <- terra::extract(x = aux.layer, y = pos_sf)
          if (base::is.data.frame(buf_vals) && base::ncol(buf_vals) >= 2) {
            base::names(buf_vals)[2] <- "Buffer"
            val <- buf_vals |>
              dplyr::group_by(.data$ID) |>
              dplyr::summarise(Buffer_mean = base::mean(.data$Buffer, na.rm = TRUE)) |>
              dplyr::pull(.data$Buffer_mean)
            if (base::is.nan(val)) val <- NA_real_
          }
        }
        
        sel$W[jj] <- val
      }
      
      # layer thickness (positive metres)
      sel$Height <- NA_real_
      for (jj in base::seq_len(nrow(sel))) {
        if (jj == 1L) {
          sel$Height[jj] <- 0 - sel$Depth[jj]            # from surface (0) to first layer
        } else {
          sel$Height[jj] <- sel$Depth[jj - 1] - sel$Depth[jj]  # delta between layers (Depth is negative)
        }
      }
      
      # weighted mean upward velocity
      if (base::all(base::is.na(sel$W)) || base::sum(sel$Height, na.rm = TRUE) == 0) {
        aux.wz <- 0
      } else {
        sel$Each <- sel$W * sel$Height
        aux.wz <- base::sum(sel$Each, na.rm = TRUE) / base::sum(sel$Height, na.rm = TRUE)
        if (base::is.na(aux.wz)) aux.wz <- 0
      }
      
      df$Wz[r] <- aux.wz
      base::gc()
    }
    
    if (base::isTRUE(export_step) && !base::is.null(export_path)) {
      utils::write.csv(df |> dplyr::select(-.data$aux.date),
                       base::paste0(export_path, ".csv"),
                       row.names = FALSE)
    }
    
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  
  if (verbose) base::close(pb)
  df |> dplyr::select(-.data$aux.date)
}




### This one for Copernicus Marine vertical current data to complement 2024 data to Bluelink BRAN 2020 data which will not be updated anymore until new product released

## Note: The NetCDF files need to be downloaded from Copernicus Mairne sepearately in Python using the Copernicus Marine Toolbox

# Everything else in this function is in accordance with the Bluelink function 

extractWz_CM <- function(df, X, Y, datetime, bathy_path, folder_name, export_path, max_depth = -200, 
                         fill_gaps = TRUE, buffer = 10000, verbose = TRUE, export_step = TRUE, keep_nc_files = TRUE) {
  
  if (!dir.exists(folder_name)) {
    base::stop(base::sprintf("NetCDF directory not found: %s", folder_name), call. = FALSE)
  }
  if (!("Depth" %in% base::names(df))) {
    base::warning("Depth column not found; returning input df unchanged.")
    return(df)
  }
      
  df <- df |>
    dplyr::mutate(
      Depth_Wz = base::ifelse(Depth > -5, -6, base::ifelse(Depth < max_depth, max_depth, Depth)),
      Depth_Wz = base::round(Depth_Wz, 1),
      aux.date = base::substr(base::as.character(.data[[datetime]]), 1, 7)
    )
  
  if (dir.exists(folder_name) == FALSE) 
    dir.create(folder_name)
  df$aux.date <- substr(df[,which(names(df) == datetime)], 1, 7)
  
  if ("Wz" %in% names(df)) {
    message("Previous Wz calculations found. Continuing...")
    dates <- unique(df$aux.date[is.na(df$Wz)])  
  } else {
    dates <- unique(df$aux.date)
    df$Wz <- NA
  }
  
  nc_files <- list.files(folder_name, pattern = "\\.nc$", full.names = TRUE)
  
  # cache all rasters once
  nc_idx <- lapply(nc_files, function(nc_file) {
    r <- terra::rast(nc_file, subds = "ocean_w")
    if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- "OGC:CRS84"
    list(
      file = nc_file,
      rast = r,
      tvec = as.Date(terra::time(r), tz = "UTC"),  # per-layer timestamps (length = nlyr)
      dvec = as.numeric(terra::depth(r))           # matching depth per layer (length = nlyr)
    )
  })
  
  if (verbose) {
    message("Processing Wz data")
    pb <- txtProgressBar(min = 0, max = nrow(df), initial = 0, style = 3, width = 60)
  }
  
  for (ri in seq_len(nrow(df))) {
    
    # skip if already computed
    if (!is.na(df$Wz[ri])) { if (verbose) setTxtProgressBar(pb, ri); next }
    
    dt_i     <- as.Date(df[[datetime]][ri])
    bottom_m <- abs(df$Depth_Wz[ri])
    found    <- FALSE
    
    # try each file until we find the date
    for (k in seq_along(nc_idx)) {
      idx <- nc_idx[[k]]
      
      lyr_day <- which(idx$tvec == dt_i)
      if (!length(lyr_day)) next
      
      # order depths shallow->deep for that day and keep to bottom_m
      d_day <- idx$dvec[lyr_day]
      o     <- order(d_day, na.last = NA)
      lyr_day <- lyr_day[o]
      d_day   <- d_day[o]
      
      keep <- which(d_day <= bottom_m)
      if (!length(keep)) next
      
      L      <- lyr_day[keep]
      d_keep <- d_day[keep]
      
      # extract W at each kept layer (with optional buffer fill)
      Wvals <- rep(NA_real_, length(L))
      for (j in seq_along(L)) {
        rlyr <- idx$rast[[L[j]]]
        val  <- terra::extract(rlyr, df[ri, c(X, Y)])[1, 2]
        
        if (is.na(val) && fill_gaps) {
          pt  <- sf::st_as_sf(df[ri, c(X, Y)], coords = c(1, 2), crs = 4326, remove = FALSE)
          buf <- suppressMessages(suppressWarnings(sf::st_buffer(pt, buffer)))  # degrees by default
          ex  <- terra::extract(rlyr, buf)
          if (is.data.frame(ex) && ncol(ex) >= 2) {
            val <- mean(ex[[2]], na.rm = TRUE)
            if (is.nan(val)) val <- NA_real_
          }
        }
        Wvals[j] <- val
      }
      
      # positive layer thicknesses (m), assuming d_keep is ascending
      Height <- c(d_keep[1], diff(d_keep))
      
      if (all(is.na(Wvals)) || sum(Height, na.rm = TRUE) == 0) {
        df$Wz[ri] <- NA_real_
      } else {
        df$Wz[ri] <- sum(Wvals * Height, na.rm = TRUE) / sum(Height, na.rm = TRUE)
      }
      
      found <- TRUE
      break  # handled this row; go to next row
    }
    
    if (!found && verbose) {
      message("No matching date in any file for row ", ri, " (", dt_i, ").")
    }
    
    if (export_step && (ri %% 50 == 0 || ri == nrow(df))) {
      utils::write.csv(df[, -which(names(df) == "aux.date")], paste0(export_path, ".csv"), row.names = FALSE)
    }
    
    if (verbose) setTxtProgressBar(pb, ri)
  }
  
  if (verbose) close(pb)
  
  if (!keep_nc_files)
    invisible(file.remove(folder_name))
  return(df)
}




# Downlaod MUR SST monthly global data   ----------------------------------

download_mursst41_monthly <- function(from_year = 2005L,
                                     out_dir,
                                     until_yearmon = NULL,   # e.g. "2025-08"; NULL = last full month
                                     retries = 4L,
                                     backoff_sec = 2) {
  if (base::missing(out_dir) || base::is.null(out_dir)) {
    base::stop("Please provide 'out_dir'.", call. = FALSE)
  }
  if (!base::requireNamespace("curl", quietly = TRUE)) {
    base::stop("Please install.packages('curl') first.", call. = FALSE)
  }
  
  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)
  
  # Determine last full month (default) or use user-supplied "YYYY-MM"
  last_day <- if (base::is.null(until_yearmon)) {
    first_this_month <- base::as.Date(base::format(base::Sys.Date(), "%Y-%m-01"))
    first_this_month - 1L
  } else {
    ym <- base::as.Date(paste0(until_yearmon, "-01"))
    seq_to_next <- base::seq.Date(ym, by = "month", length.out = 2L)
    seq_to_next[2L] - 1L
  }
  
  start_date <- base::as.Date(sprintf("%d-01-01", from_year))
  first_month <- base::as.Date(base::format(start_date, "%Y-%m-01"))
  last_month  <- base::as.Date(base::format(last_day,   "%Y-%m-01"))
  
  months_seq <- base::seq.Date(from = first_month, to = last_month, by = "month")
  n <- base::length(months_seq)
  if (n == 0L) {
    base::message("No months to download for the requested range.")
    return(invisible(base::data.frame(file = character(), status = character())))
  }
  
  # Start and end date per month
  next_months <- base::seq.Date(from = months_seq[1L], by = "month", length.out = n + 1L)
  month_start <- months_seq
  month_end   <- next_months[-1L] - 1L
  
  # Build expected filenames and URLs
  fnames <- base::sprintf("%s%s-GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc",
                          base::format(month_start, "%Y%m%d"),
                          base::format(month_end,   "%Y%m%d"))
  base_url <- "https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/"
  urls <- base::paste0(base_url, fnames)
  dest <- base::file.path(out_dir, fnames)
  
  # Skip what you already have
  need <- !base::file.exists(dest)
  urls <- urls[need]; dest <- dest[need]
  month_start <- month_start[need]; month_end <- month_end[need]
  
  if (base::length(urls) == 0L) {
    base::message("All files already present.")
    return(invisible(base::data.frame(file = character(), status = character())))
  }
  
  # Sequential downloads with retry/backoff
  status <- base::character(base::length(urls))
  for (i in base::seq_along(urls)) {
    ok <- FALSE
    for (attempt in 0L:retries) {
      res <- base::tryCatch({
        curl::curl_download(url = urls[i], destfile = dest[i], mode = "wb")
        TRUE
      }, error = function(e) {
        FALSE
      })
      if (res) { ok <- TRUE; break }
      # backoff
      slp <- backoff_sec * (2 ^ attempt)
      base::Sys.sleep(slp)
    }
    status[i] <- if (ok) "ok" else "failed"
    base::message(base::sprintf("[%d/%d] %s : %s",
                                i, base::length(urls), base::basename(dest[i]), status[i]))
  }
  
  base::data.frame(file = dest,
                   month_start = month_start,
                   month_end = month_end,
                   status = status,
                   row.names = NULL)
}




# Downloads MUR monthly SST files from 2005-01 through the last full month
# Status lines + optional cli progress bar + final summary + optional CSV log
download_mursst41_monthly <- function(from_year = 2005L,
                                          out_dir,
                                          until_yearmon = NULL,   # e.g., "2025-08"; NULL = last full month
                                          retries = 4L,
                                          backoff_sec = 2,
                                          progress = c("cli","text","none"),
                                          write_log = TRUE) {
  if (base::missing(out_dir) || base::is.null(out_dir)) {
    base::stop("Please provide 'out_dir'.", call. = FALSE)
  }
  if (!base::requireNamespace("curl", quietly = TRUE)) {
    base::stop("Please install.packages('curl') first.", call. = FALSE)
  }
  
  progress <- base::match.arg(progress)
  
  has_cli <- base::requireNamespace("cli", quietly = TRUE)
  if (progress == "cli" && !has_cli) progress <- "text"
  
  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)
  
  # Determine last full month
  last_day <- if (base::is.null(until_yearmon)) {
    first_this_month <- base::as.Date(base::format(base::Sys.Date(), "%Y-%m-01"))
    first_this_month - 1L
  } else {
    ym <- base::as.Date(paste0(until_yearmon, "-01"))
    seq_to_next <- base::seq.Date(ym, by = "month", length.out = 2L)
    seq_to_next[2L] - 1L
  }
  
  # Month sequence
  start_date  <- base::as.Date(sprintf("%d-01-01", from_year))
  first_month <- base::as.Date(base::format(start_date, "%Y-%m-01"))
  last_month  <- base::as.Date(base::format(last_day,   "%Y-%m-01"))
  months_seq  <- base::seq.Date(from = first_month, to = last_month, by = "month")
  if (base::length(months_seq) == 0L) {
    base::message("No months to download for the requested range.")
    return(invisible(base::data.frame()))
  }
  
  # Start/end per month
  next_months <- base::seq.Date(from = months_seq[1L], by = "month", length.out = base::length(months_seq) + 1L)
  month_start <- months_seq
  month_end   <- next_months[-1L] - 1L
  
  # Filenames + URLs
  make_fname <- function(s, e) {
    base::sprintf("%s%s-GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc",
                  base::format(s, "%Y%m%d"), base::format(e, "%Y%m%d"))
  }
  fnames  <- base::mapply(make_fname, month_start, month_end, USE.NAMES = FALSE)
  base_url <- "https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/"
  urls    <- base::paste0(base_url, fnames)
  dest    <- base::file.path(out_dir, fnames)
  
  # Skip existing files
  need <- !base::file.exists(dest)
  urls <- urls[need]; dest <- dest[need]
  month_start <- month_start[need]; month_end <- month_end[need]
  fnames <- fnames[need]
  
  if (base::length(urls) == 0L) {
    base::message("All files already present.")
    return(invisible(base::data.frame()))
  }
  
  # Helpers
  fmt_bytes <- function(n) {
    if (!base::is.finite(n) || base::is.na(n)) return("NA")
    units <- c("B","kB","MB","GB","TB")
    i <- base::floor(base::log(max(n, 1), 1024))
    i <- max(0, min(i, base::length(units)-1))
    base::sprintf("%.1f %s", n / (1024^i), units[i+1])
  }
  
  n <- base::length(urls)
  results <- base::data.frame(
    file = dest,
    month_start = month_start,
    month_end = month_end,
    status = base::rep(NA_character_, n),
    size_bytes = base::rep(NA_real_, n),
    attempt = base::rep(NA_integer_, n),
    stringsAsFactors = FALSE
  )
  
  # Progress
  if (progress == "cli") {
    cli::cli_progress_bar("Downloading MUR monthly files",
                          total = n, clear = FALSE, format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} {cli::pb_eta}")
  }
  
  # Download loop with retries + friendly status lines
  for (i in base::seq_len(n)) {
    ok <- FALSE
    last_err <- ""
    for (attempt in 0L:retries) {
      res <- base::tryCatch({
        curl::curl_download(url = urls[i], destfile = dest[i], mode = "wb")
        TRUE
      }, error = function(e) {
        last_err <<- base::conditionMessage(e)
        FALSE
      })
      if (res) { ok <- TRUE; results$attempt[i] <- attempt + 1L; break }
      base::Sys.sleep(backoff_sec * (2 ^ attempt))
    }
    
    # size and status
    if (ok) {
      sz <- base::file.info(dest[i])$size
      results$status[i] <- "ok"
      results$size_bytes[i] <- if (base::is.na(sz)) NA_real_ else sz
      line <- base::sprintf("[%d/%d] %s : ok (%s)",
                            i, n, base::basename(dest[i]), fmt_bytes(results$size_bytes[i]))
    } else {
      results$status[i] <- paste0("error: ", last_err)
      line <- base::sprintf("[%d/%d] %s : FAILED (%s)",
                            i, n, base::basename(dest[i]), last_err)
    }
    
    if (progress == "cli") {
      cli::cli_progress_update()
      # also echo the line so you get the nice per-file status you liked
      cli::cli_inform(line)
    } else if (progress == "text") {
      base::message(line)
    }
  }
  
  if (progress == "cli") cli::cli_progress_done()
  
  # Final summary
  n_ok <- base::sum(results$status == "ok", na.rm = TRUE)
  n_fail <- base::sum(!base::startsWith(results$status, "ok"), na.rm = TRUE)
  total_sz <- fmt_bytes(base::sum(results$size_bytes[results$status == "ok"], na.rm = TRUE))
  
  base::message("")
  base::message("===== MUR monthly download summary =====")
  base::message(base::sprintf("Saved to: %s", out_dir))
  base::message(base::sprintf("OK: %d | Failed: %d | Total size: %s", n_ok, n_fail, total_sz))
  if (n_fail > 0) {
    base::message("Failed files:")
    base::message(paste0("  - ", base::basename(results$file[!base::startsWith(results$status, "ok")])), sep = "\n")
  }
  
  # Optional CSV log
  if (write_log) {
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M%S")
    log_path <- base::file.path(out_dir, paste0("mursst_download_log_", ts, ".csv"))
    utils::write.csv(results, log_path, row.names = FALSE)
    base::message(base::sprintf("Log written: %s", log_path))
  }
  
  results
}







# Predictor Raster Functions  ---------------------------------------------



# Extract the month from the layer names
extract_month_from_names <- function(spat_raster) {
  # Extract the names of the layers
  layer_names <- names(spat_raster)
  
  # Extract the month part assuming the format "YYYY-MM-DD"
  months <- as.integer(substring(layer_names, 6, 7))
  
  # Handle potential NAs and return the result
  if (any(is.na(months))) {
    warning("Some month values could not be coerced to integers. Check the layer names format.")
  }
  
  return(months)
}




## make monthly raster stacks from CMEMS Wo dowanloads for Wo variable 

create_monthly_rasters_from_files_CM <- function(nc_folder, desired_depth, reference_raster = NULL, output_path, output_filename, start_date = NULL) {
  
  # Function to process auxiliary data from nc file
  process_aux_data_CM <- function(nc_file) {
    # Load the NetCDF file
    nc_CM <- terra::rast(nc_file)
    gc()  # Garbage collection to free up memory
    
    # Extract and process variable names
    aux_names <- stringr::str_remove(names(nc_CM), pattern = "wo_depth=")
    aux_names <- stringr::str_split(aux_names, pattern = "_")
    
    aux_depth <- NULL
    aux_time <- NULL
    for (i in 1:length(aux_names)) {
      aux_depth <- c(aux_depth, aux_names[[i]][1])
      aux_time <- c(aux_time, aux_names[[i]][2])
    }
    
    aux_names_df <- data.frame(Depth = as.numeric(aux_depth), Time = as.numeric(aux_time))
    
    # ---- start_date handling (no structural changes elsewhere) ----
    # If start_date is NULL, keep original origin of "2024-01-16"
    origin <- if (is.null(start_date)) as.Date("2024-01-16", tz = "UTC") else as.Date(start_date, tz = "UTC")
    idx <- as.integer(aux_names_df$Time)  # 1,2,3,... month indices
    date_seq <- seq.Date(from = origin, by = "1 month", length.out = max(idx, na.rm = TRUE))
    aux_names_df$Time <- date_seq[idx]
    # ----------------------------------------------------------------
    
    gc()  # Garbage collection to free up memory before returning
    return(aux_names_df)
  }
  
  # List all NetCDF files in the directory
  nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    stop("No NetCDF files found in the specified folder.")
  }
  
  all_layers <- list()
  available_depths <- NULL
  
  for (nc_file in nc_files) {
    message(paste("Processing file:", nc_file))
    aux_data <- process_aux_data_CM(nc_file)
    available_depths <- unique(c(available_depths, aux_data$Depth))
  }
  
  # Find the nearest available depth
  differences <- abs(available_depths - desired_depth)
  nearest_depth <- available_depths[which.min(differences)]
  message(paste("Nearest available depth layer selected:", nearest_depth))
  
  for (nc_file in nc_files) {
    aux_data <- process_aux_data_CM(nc_file)
    
    # Filter by the nearest depth layer
    depth_filter <- aux_data$Depth == nearest_depth
    if (!any(depth_filter)) {
      message(paste("No matching depth layer found in file:", nc_file))
      next  # Skip this file if no layers match the specified depth
    }
    
    filtered_data <- aux_data[depth_filter, ]
    
    # Load the NetCDF file and filter by the selected layers
    nc_CM <- terra::rast(nc_file)
    filtered_layers <- nc_CM[[which(depth_filter)]]
    gc()  # Garbage collection to free up memory
    
    for (i in 1:nrow(filtered_data)) {
      current_layer <- filtered_layers[[i]]
      names(current_layer) <- as.character(filtered_data$Time[i])
      all_layers <- c(all_layers, list(current_layer))
      gc()  # Garbage collection to free up memory
    }
  }
  
  if (length(all_layers) == 0) {
    stop("No layers were extracted. Please check the depth layer and NetCDF files.")
  }
  
  combined_raster <- rast(all_layers)
  
  # Resample and crop the combined raster to match the reference raster if provided
  if (!is.null(reference_raster)) {
    combined_raster <- terra::resample(combined_raster, reference_raster)
    combined_raster <- terra::crop(combined_raster, reference_raster)
  }
  
  # Construct the full output file path with .tif extension
  output_file <- file.path(output_path, paste0(output_filename, ".tif"))
  
  # Save the combined raster stack to the specified output file
  terra::writeRaster(combined_raster, output_file, overwrite = TRUE)
  message(paste("Saved combined raster stack to:", output_file))
  
  assign(output_filename, terra::rast(output_file), envir = .GlobalEnv)
  
  gc()
}



## Calculate mean velocity:

calculate_uv <- function(stack) {
  layers <- nlyr(stack) / 2
  result_stack <- terra::rast(ncol=ncol(stack), nrow=nrow(stack), extent=ext(stack), crs=crs(stack))
  
  for (i in 1:layers) {
    uo <- stack[[i]]
    vo <- stack[[i + layers]]
    uv_layer <- sqrt(uo^2 + vo^2)
    names(uv_layer) <- names(stack)[i]  # Preserve the original date name
    result_stack <- c(result_stack, uv_layer)
    gc()  # Perform garbage collection
  }
  return(result_stack)
}






# SDMtune Predictions -----------------------------------------------------

sdmtune_predict_monthly <- function(model,
                                    monthly_list,       # list of SpatRaster (one per month)
                                    dates = NULL,       # Date vector; if NULL we try names(monthly_list)
                                    type = "cloglog",
                                    extent = NULL,      # crop to this extent before predicting
                                    verbose = TRUE) {
  
  stopifnot(inherits(model, "SDMmodel"))
  
  # Required variables (used to subset each stack)
  req_vars <- colnames(model@data@data)
  if (is.null(req_vars) || !length(req_vars)) {
    stop("Couldn't detect model variables from model@data@data.")
  }
  
  # Dates
  if (is.null(dates)) {
    nm <- names(monthly_list)
    if (is.null(nm)) stop("Provide 'dates' or set names(monthly_list) to ISO dates.")
    dates <- as.Date(nm)
    if (any(is.na(dates))) stop("names(monthly_list) must be parseable as Date; otherwise pass 'dates='.")
  }
  if (length(dates) != length(monthly_list))
    stop("Length of 'dates' must match length of 'monthly_list'.")
  
  # Coerce 'extent' to a SpatExtent if provided
  ex <- NULL
  if (!is.null(extent)) {
    if (inherits(extent, "SpatExtent")) {
      ex <- extent
    } else if (inherits(extent, "SpatRaster") || inherits(extent, "SpatVector")) {
      ex <- terra::ext(extent)
    } else if (is.numeric(extent) && length(extent) == 4) {
      ex <- terra::ext(extent)
    } else {
      stop("`extent` must be a SpatExtent, SpatRaster/SpatVector, or numeric c(xmin, xmax, ymin, ymax).")
    }
  }
  
  N <- length(monthly_list)
  preds <- vector("list", N)
  
  for (i in seq_len(N)) {
    r <- monthly_list[[i]]
    if (!inherits(r, "SpatRaster")) stop("monthly_list[[", i, "]] is not a SpatRaster.")
    
    # Crop to extent if requested
    if (!is.null(ex)) {
      r <- terra::crop(r, ex)
      if (terra::ncell(r) == 0) {
        stop(sprintf("Cropped raster for %s has zero cells; check the extent/CRS.",
                     format(dates[i], "%Y-%m-%d")))
      }
    }
    
    # Make sure all required variables are present
    missing <- setdiff(req_vars, names(r))
    if (length(missing)) {
      stop(sprintf("Month %s is missing layers: %s",
                   format(dates[i], "%Y-%m-%d"), paste(missing, collapse = ", ")))
    }
    
    # Subset and keep model order
    r <- r[[req_vars]]
    
    if (verbose) cat("Predicting", format(dates[i], "%Y-%m"), sprintf("(%d/%d)\n", i, N))
    
    pr <- SDMtune::predict(model, data = r, type = type)
    names(pr) <- as.character(dates[i])   # e.g., "2010-01-01"
    preds[[i]] <- pr
  }
  
  out <- terra::rast(preds)
  names(out) <- as.character(dates)
  out
}



expected_fits <- function(pop = 20, gen = 5, keep_best = 0.4, keep_random = 0.2,
                          rounding = c("floor", "round", "ceiling")) {
  rounding <- match.arg(rounding)
  
  # basic checks
  stopifnot(pop > 0, gen >= 0,
            keep_best >= 0, keep_random >= 0,
            keep_best + keep_random <= 1)
  
  rfun <- switch(rounding,
                 floor   = base::floor,
                 round   = base::round,
                 ceiling = base::ceiling)
  
  kept_per_gen <- rfun(pop * keep_best) + rfun(pop * keep_random)
  new_per_gen  <- pop - kept_per_gen
  total_fits   <- pop + gen * new_per_gen
  
  as.integer(total_fits)
}





# define a patched trainBRT that reads weights from data@data$case_w if present ---
patched_trainBRT <- function(data,
                             distribution = "bernoulli",
                             n.trees = 100,
                             interaction.depth = 1,
                             shrinkage = 0.1,
                             bag.fraction = 0.5) {
  # Build the SDMmodel shell
  result <- SDMtune::SDMmodel(data = data)
  
  # Predictors (no weights column here)
  df <- base::cbind(pa = data@pa, data@data)
  
  # ---- compute weights PER SUBSET (robust to CV subsetting) ----
  w <- NULL
  
  if ("month" %in% base::names(data@data)) {
    # Month-wise balancing: in each month, sum(w_pres) ≈ sum(w_abs)
    w <- base::numeric(nrow(df))
    mlev <- base::unique(data@data$month)
    for (m in mlev) {
      idx  <- which(data@data$month == m)
      if (length(idx)) {
        pres <- base::sum(data@pa[idx] == 1L)
        absn <- base::sum(data@pa[idx] == 0L)
        w[idx] <- if (absn > 0) ifelse(data@pa[idx] == 1L, 1, pres/absn) else 1
      }
    }
  } else {
    # Global balancing: sum(w_pres) ≈ sum(w_abs)
    pres <- base::sum(data@pa == 1L)
    absn <- base::sum(data@pa == 0L)
    w <- if (absn > 0) ifelse(data@pa == 1L, 1, pres/absn) else base::rep(1, nrow(df))
  }
  
  # Final safety: no NAs/negatives/zeros
  if (base::anyNA(w)) {
    w[base::is.na(w)] <- 1
  }
  w[w < 0] <- 0
  
  
  # ---- fit gbm with weights ----
  gbm_fit <- gbm::gbm(
    pa ~ .,
    data = df,
    distribution = distribution,
    n.trees = n.trees,
    interaction.depth = interaction.depth,
    shrinkage = shrinkage,
    bag.fraction = bag.fraction,
    weights = w
  )
  
  # Wrap back into SDMtune classes
  brt_obj <- SDMtune::BRT(
    n.trees = n.trees,
    distribution = distribution,
    interaction.depth = interaction.depth,
    shrinkage = shrinkage,
    bag.fraction = bag.fraction,
    model = gbm_fit
  )
  result@model <- brt_obj
  result
}




patched_trainANN <- function(data,
                             size,
                             decay = 0,
                             rang  = 0.7,
                             maxit = 100) {
  result <- SDMtune::SDMmodel(data = data)
  
  # One-hot encode factors (month etc.)
  X <- stats::model.matrix(~ . - 1, data = data@data)
  y <- data@pa
  
  # ---- GLOBAL class-balanced weights per subset ----
  pres <- sum(y == 1L); absn <- sum(y == 0L)
  w <- if (absn > 0) ifelse(y == 1L, 1, pres/absn) else rep(1, length(y))
  w[is.na(w)] <- 1
  w[w < 0] <- 0
  
  # ---- Compute required network size and set MaxNWts locally ----
  p <- ncol(X)
  needed  <- size * (p + 2) + 1L         # total params for 1-hidden-layer net
  MaxNWts <- max(1000L, as.integer(needed * 1.25))  # small safety margin
  MaxNWts <- min(MaxNWts, 5e6)           # hard cap to avoid silly values
  
  # Temporarily raise the global option, then restore on exit
  old_opts <- options(nnet.MaxNWts = MaxNWts)
  on.exit(options(old_opts), add = TRUE)
  
  utils::capture.output(
    mdl <- nnet::nnet(
      x = X, y = y, weights = w,
      size = size, decay = decay, rang = rang, maxit = maxit,
      entropy = TRUE, linout = FALSE, trace = FALSE,
      MaxNWts = MaxNWts   # pass explicitly too
    )
  )
  
  result@model <- SDMtune::ANN(
    size = size, decay = decay, rang = rang, maxit = maxit, model = mdl
  )
  result
}
