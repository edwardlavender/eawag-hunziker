#########################
#########################
#### process_data.R

#### Aims
# 1) Process raw (fish) data for analysis

#### Prerequisites
# 1) Acquire raw (fish) data (assembled by YD) 


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

#### Load data
fish <- readxl::read_excel(here_data_raw("5. Real_residents_and_migs_below_SL243mm_without_N2&LRN.xlsx"))
# View(fish)


#########################
#########################
#### Process data

#### Select and process relevant columns
fish$stream <- as.character(fish$stream)
fish$stream[fish$stream == "Klosterbach UR"] <- "Klosterbach (UR)"
fish$stream[fish$stream == "Klosterbach SZ"] <- "Klosterbach (SZ)"
fish <- 
  fish |> 
  mutate(id = row_number()) |>
  select(id, 
         date = tag.date,
         stream,
         sex = sex, 
         length = SL,
         migration = downmig.bin, 
         migration_date = outmig.date, 
         ) |>
  mutate(id = factor(id), 
         sex = factor(sex),
         stream = factor(stream),
         migration = factor(migration),
         date = as.Date(date),
         yday = lubridate::yday(date),
         period = as.integer(as.Date("2015-06-30") - date), 
         length = length/10, 
         migration_date = as.Date(migration_date),
         migration_yday = lubridate::yday(migration_date)
         ) |>
  select(id, date, yday, period, stream, 
         sex, length, migration, migration_date, migration_yday)

#### Define migrants
migrants <- 
  fish |> 
  filter(migration == 1L) |> 
  select(length, sex, migration_date, migration_yday, yday, stream) 

#### Save data
saveRDS(fish, here_data("fish.rds"))
saveRDS(migrants, here_data("migrants.rds"))


#### End of code.
#########################
#########################