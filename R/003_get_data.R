#########################
#########################
#### get_data.R

#### Aims
# 1) Acquires necessary data for this project
# * The fish datasets have been collated by YD and are in data-raw/
# * We also need spatial datasets to map the study area
# * This code simply copies those spatial datasets into RStudio Project for QGIS 
# * This code is only designed to run on SIA-LAVENDED-M

#### Prerequisites
# 1) Obtain open-source spatial datasets from online sources


#########################
#########################
#### Set up

#### Wipe workspace
rm(list = ls()) 
try(pacman::p_unload("all"), silent = TRUE) 
dv::clear() 


#########################
#########################
#### Copy data

#### Define file paths 
con <- file.path("/Users", "lavended", "Documents", 
                 "eawag", "research", "datasets")
sink <- here::here("QGIS", "extdata")

#### Obtain GADM (country boundary) datasets
files <- list.files(file.path(con, "data-raw", "spatial", "gadm"), 
                    full.names = TRUE, recursive = TRUE)
files <- files[stringr::str_detect(files, "_0")]
file.copy(files, file.path(sink, basename(files)))

#### Obtain water feature (rivers, lakes) shapefiles
files <- list.files(file.path(con, "data", "water-features"), 
                    full.names = TRUE, recursive = TRUE)
file.copy(files, file.path(sink, basename(files)))


#### End of code.
#########################
#########################