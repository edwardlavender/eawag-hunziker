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
library(dplyr)
library(sf)
library(raster)
library(riverdist)
library(ggplot2) 
library(tictoc)

#### Load data
# Load altitude data
alt <- terra::rast(here_data_raw("spatial", "altitude.tif"))
# alt <- terra::project(alt, "epsg:4326") 
# Define connection to stream KML files
con <- here_data_raw("spatial", "streams", "geometry", "linestrings")
# Define stream names
stream_abbr <- 
  con |> 
  list.files(pattern = ".kml") |> 
  basename() |>
  substr(1, 3) |> 
  unique() |>
  sort()

#### Global parameters
manual <- FALSE


#########################
#########################
#### Identify coordinates of the mouth 

#### Extract coordinates
mouth_xy <- 
  lapply(stream_abbr, function(name) {
    # Define stream
    stream <- 
      con |>
      file.path(paste0(name, "-linestrings.kml")) |> 
      st_read(quiet = TRUE) |> 
      janitor::clean_names() |>
      filter(description %in% c("0", "1")) 
    # Extract coordinates
    xy <- 
      stream |>
      filter(description == "0") |>
      st_coordinates() |>
      as.data.frame() |>
      slice(1L) |>
      mutate(stream = name) |> 
      dplyr::select(stream, lon = X, lat = Y)
    # Visualise coordinates
    if (manual) {
      p <- 
        ggplot(stream) + 
        geom_sf(aes(colour = factor(description))) + 
        ggsflabel::geom_sf_label(aes(label = description)) + 
        geom_point(data = xy, aes(x = lon, y = lat), shape = 1, size = 4) +
        labs(title = xy$stream)
      print(p)
    } 
    if (manual) {
      readline("Press [Enter] to continue...")
    }
    # Return coordinates
    xy
}) |> 
  bind_rows()


#########################
#########################
#### Calculate altitude of each stream section

altitudes <- 
  lapply(stream_abbr, function(name) {
    
    #### Define stream
    # name <- "SGN"
    stream <- 
      con |>
      file.path(paste0(name, "-linestrings.kml")) |> 
      st_read(quiet = TRUE) |> 
      janitor::clean_names()  |>
      # Fix stream order
      mutate(description = as.integer(as.character(description))) |> 
      arrange(description)
    
    #### Define coordinates of stream sections
    xy <- 
      lapply(split(stream, stream$description), function(d) {
      xy <- st_coordinates(d)[1, c("X", "Y")] 
      xy <- data.frame(x = xy[1], y = xy[2])
      xy$stream  <- name
      xy$section <- d$description
      xy
    }) |> bind_rows()
    rownames(xy) <- NULL
    # Convert to LV95 (NB: this does not work via sf)
    sp::coordinates(xy) <- ~x + y
    xy@proj4string      <- sp::CRS("+init=EPSG:4326")
    xy <- sp::spTransform(xy, sp::CRS("+init=epsg:2056"))
    xy <- sp::coordinates(xy)
    
    #### Extract altitudes
    altm <- terra::extract(alt, xy)[, 1]
    data.frame(stream = name, section = stream$description, x = xy[, 1], y = xy[, 2], alt = altm)
    
  }) |> bind_rows()


#########################
#########################
#### Calculate distances to the lake (~6 s)

tic()
distances <- 
  lapply(stream_abbr, function(name) {
    
    #### Define stream & stream sections
    # name     <- "SGN"
    message(name)
    stream   <- 
      file.path(con, paste0(name, "-linestrings.kml")) |> 
      st_read(quiet = TRUE) |> 
      janitor::clean_names() |> 
      # Fix stream order
      # This is essential for correct functioning of {riverdist}
      mutate(description = as.integer(as.character(description))) |> 
      arrange(description)
    # Stream sections are named as follows:
    # * 0 (the lake to antenna/start of first tagging section)
    # * 1 (antenna/first tagging section to start of section section etc.)
    stream$description
    stopifnot(!is.unsorted(stream$description))
    # {riverdist} renames sections from 1 to the number of sections
    # * 1 = the lake to the antenna
    # * 2 = section 1 etc.
    length(stream$description)
    if (manual) {
      # plot(stream["description"])
      ggplot(stream) + 
        geom_sf(aes(colour = factor(description))) + 
        ggsflabel::geom_sf_label(aes(label = description))
      }
    
    #### Format stream for distance calculations with {riverdist}
    # Define stream in UTM
    stream    <- st_transform(stream, 2056)
    # Pull out lengths (m)
    section_lengths <- st_length(stream) |> as.numeric()
    # Convert to spatial 
    stream    <- as(stream, "Spatial")
    if (manual) plot(stream)
    # Define network
    # ... It is important that the tolerance is less than the minimum section length
    tol <- ifelse(min(section_lengths) < 30, 10, 30)
    stopifnot(tol < min(section_lengths))
    net_file <- here_data_raw("spatial", "streams", "geometry", "networks", paste0(name, ".rds"))
    if (!file.exists(net_file)) {
      # Define network 
      net <- line2network(stream, tolerance = tol)
      # (optional) clean up network
      # This is necessary for: 
      # a) SGN
      # ... The connection between 12 and 8 (as defined by {riverdist})
      # ... is not properly recognised because connections need to 
      # ... start/end at a specific segment.
      # ... This is fixed by moving the start of section 12 a few m 
      # ... so it starts at the end of section 7. 
      # net <- cleanup(net)
      saveRDS(net, net_file) 
    } else {
      net <- readRDS(net_file)
    }
    if (manual) plot(net)
    
    #### Calculate distances
    # Calculate the distances from river mouth to antenna
    dist_lake_to_antenna <- 
      riverdistance(startseg = 1, endseg = 2,
                    startvert = 1, endvert = 1, 
                    rivers = net, 
                    map = TRUE)
    # Calculate distances from the river mouth to the start of each section
    png(here_fig("checks", "riverdist", paste0(name, ".png")), 
        height = 10, width = 10, units = "in", res = 300)
    pp <- par(mfrow = prettyGraphics::par_mf(max(net$lineID$rivID)), 
              mar = c(1, 1, 1, 1))
    ids <- net$lineID$rivID
    sections <- ids[2:length(ids)]
    dist_lake_to_section <- 
      sapply(sections, function(end) {
        print(paste0("... ", end))
        riverdistance(startseg = 1, endseg = end,
                      startvert = 1, endvert = 1, 
                      rivers = net, 
                      map = TRUE)
      })
    par(pp)
    dev.off()
    
    #### Checks (passed) 
    # * In figures, all river sections should be correctly placed & identified 
    # * There should be no missing sections
    # * Distance arrows should correctly highlight from the lake to the 
    # * ... start of each section. 
    
    #### Define data frame with distances to the lake
    # Note that the first distance is also the distance from lake to antenna
    data.frame(stream = name, 
               section = sections - 1,
               dist_to_lake = dist_lake_to_section)
  
}) |> dplyr::bind_rows()
toc()


#########################
#########################
#### Save outputs

#### Collate stream characteristics (distances, altitudes)
stations              <- altitudes
stations$key          <- paste(stations$stream, stations$section)
distances$key         <- paste(distances$stream, distances$section)
stations$dist_to_lake <- distances$dist_to_lake[match(stations$key, distances$key)]
stations$dist_to_lake[stations$section == 0] <- 0

#### Save file
saveRDS(stations, here_data_raw("spatial", "streams", "characteristics", "stations.rds"))


#### End of code.
#########################
#########################