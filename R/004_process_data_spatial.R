#########################
#########################
#### 004_process_data_streams.R

#### Aims
# 1) Process stream/section data 
# ... * Calculate distances from tagging sections to the lake

#### Prerequisites
# 1) Build stream/section KML files for necessary streams


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls()) 
try(pacman::p_unload("all"), silent = TRUE) 
dv::clear() 

#### Essential packages
library(dv)
library(sf)

#### Load data
# Define connection to KML files
con <- here_data_raw("spatial", "streams")
# Define stream names
con |> 
  list.files(pattern = ".kml") |> 
  basename() |>
  substr(1, 3) |> 
  unique()

#########################
#########################
#### Calculate distances to the lake

#### Define stream & stream sections
# Define layer of interest
name     <- "GLU"
# Load streams (linestring) and stream sections (points)
stream   <- st_read(file.path(con, paste0(name, "-linestring.kml")))
sections <- st_read(file.path(con, paste0(name, "-sections.kml")))
# Guarantee correct section order
sections <- 
  sections |>
  janitor::clean_names() |> 
  dplyr::mutate(name = as.integer(as.character(name))) |> 
  dplyr::arrange(name)
# Convert to UTM coordinates (m)
stream      <- st_transform(stream, 2056)
sections    <- st_transform(sections, 2056)
stream_sp   <- as(stream, "Spatial")
sections_sp <- as(sections, "Spatial")

#### Discretise streams into a series of points 
# Interpolate points along the stream every (e.g.) 0.5 m
s <- seq(0, as.numeric(st_length(stream)), by = 0.1)
stream_int <- rgeos::gInterpolate(stream_sp, d = s)
if (FALSE) {
  raster::plot(stream_sp)
  raster::plot(stream_int, add = TRUE, pch = ".", col = "red")
}
# Calculate distances between sequential points (0.5 m)
dists <- c(0, sp::spDists(stream_int, segments = TRUE))

#### Approximate section midpoints 
mids_pos <- 
  lapply(seq_len(nrow(sections) - 1), function(i) {
    # Define river section pair
    print(i)
    sect_start  <- sections_sp[i, ]
    sect_end    <- sections_sp[i + 1, ]
    # Find the nearest locations on the river grid
    nearest_start <- which.min(raster::pointDistance(sect_start, stream_int)) |> as.integer()
    nearest_end   <- which.min(raster::pointDistance(sect_end, stream_int)) |> as.integer()
    # Find the position of the mid point
    pos <- floor(median(nearest_start:nearest_end))
    if (FALSE) {
      raster::plot(stream_sp)
      points(sections_sp, col = "dimgrey", lwd = 0.75)
      points(sect_start, pch = 21, bg = "green", col = "green")
      points(sect_end, pch = 21, bg = "red", col = "red")
      points(stream_int[nearest_start, ], pch = ".")
      points(stream_int[nearest_end, ], pch = ".")
      points(stream_int[pos, ], col = "purple", pch = "_", lwd = 3)
      readline("Press [Enter] to continue...")
    }
    pos
  }) 

#### Visualise the positions of section midpoints 
# Get coordinates 
mids <- lapply(mids_pos, function(pos) {
  stream_int[pos, ]
})
mids <- do.call(rbind, mids)
# Plot 
if (FALSE) {
  raster::plot(stream_sp)
  points(sections_sp[1, ]) # mouth of river
  points(sections_sp, col = "blue")
  points(mids, pch = "_", col = "red")
}

#### Get the cumulative distances from the mouth of the river to section midpoints
dist_to_lake <- 
  sapply(mids_pos, function(pos) {
  sum(dists[seq_len(pos)])
})
# Check distances manually (e.g., using TEST file)
if (FALSE) {
  dists_on_page <- flapper::dist_btw_clicks(lonlat = FALSE)
  dists_on_page$dist |> sum()
}

#### Define data frame with distances to the lake
data.frame(stream = name, 
           section = sections$name, 
           dist = c(0, dist_to_lake)) |>
  # Drop the first point (the lake) 
  dplyr::filter(section != 0)


#### End of code.
#########################
#########################