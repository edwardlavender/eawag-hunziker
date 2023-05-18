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
fish <- readxl::read_excel(here_data_raw("fish", "5. Real_residents_and_migs_below_SL243mm_without_N2&LRN.xlsx"))
pit  <- readxl::read_excel(here_data_raw("fish", "pit.xlsx"))
tag_sites_abb <- readxl::read_excel(here_data_raw("spatial", "streams", "stream-abbreviations.xlsx"), 
                                    sheet = "data")
tag_sites_dist <- readRDS(here_data_raw("spatial", "streams", "characteristics", "stations.rds"))
# View(fish)


#########################
#########################
#### Process migration data

#### Select and process relevant columns
# Process streams 
fish$stream <- as.character(fish$stream)
fish$stream[fish$stream == "Klosterbach UR"] <- "Klosterbach (UR)"
fish$stream[fish$stream == "Klosterbach SZ"] <- "Klosterbach (SZ)"
# Check size (length/weight) columns
table(is.na(fish$SL))
table(fish$mass == "na" | is.na(fish$mass))
plot(as.numeric(fish$mass), fish$SL)
# Define data frame
fish <- 
  fish |> 
  mutate(id = row_number()) |>
  select(id, 
         pit_id = pit.tag.no,
         fec_id = fishec.no,
         date = tag.date,
         stream,
         sex = sex, 
         length = SL, mass,
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
         mass = round(as.numeric(mass), digits = 1)/10,
         migration_date = as.Date(migration_date),
         migration_yday = lubridate::yday(migration_date)
         ) |>
  select(id, pit_id, fec_id, date, yday, period, stream, 
         sex, length, mass, migration, migration_date, migration_yday) |>
  as.data.frame()

#### Check individuals are defined by Fish Ec/PIT ID
# Fish Ec IDs are all defined and not duplicated
class(fish$fec_id)
any(is.na(fish$fec_id))
any(duplicated(fish$fec_id))
# One PIT is 'na'
table(fish$pit_id == "na")
any(is.na(fish$pit_id))
# No PITs in fish are duplicated
any(duplicated(fish$pit_id))


#########################
#########################
#### Pull out individual characteristics

#### Clean PIT data 
pit <- 
  pit |> 
  janitor::clean_names() |> 
  janitor::remove_empty(which = "cols") |>
  mutate(date = as.Date(date)) |>
  rename(fec_id = fishec) |>
  as.data.frame()

#### Check for duplicate names
# There are duplicate fish Ec IDs and PIT numbers
any(duplicated(pit$fec_id))
any(duplicated(pit$pit_id))
# These are the duplicated IDs
dup_fec <- pit$fec_id[duplicated(pit$fec_id)]
dup_fec <- dup_fec[!is.na(dup_fec)]
dup_fec
dup_pit <- pit$pit_id[duplicated(pit$pit_id)]
dup_pit <- dup_pit[!is.na(dup_pit)]
dup_pit

#### Confirm sampled fish are in the full dataset and not duplicates 
# All fish in the Hunziker data are in the overall dataset 
fish$fec_id[!(fish$fec_id %in% pit$fec_id)]
if (FALSE) {
  fish[fish$fec_id == 112964, ]
  pit |> 
    filter(species == "Trout") |> 
    filter(date == "2015-03-05") |>
    filter(river == "Schützenbrunnen") |> 
    arrange(fec_id) |> 
    janitor::remove_empty(which = "cols")
}
# Each fish should be recorded once in the overall database
# This is correct according to fishEc IDs (but not PIT IDs)
# I.e., None of the fish in the Hunziker data have duplicated IDs
table(fish$fec_id %in% dup_fec)
col_id <- "fec_id"
recorded <- 
  sapply(split(fish, fish$id), function(d) {
  nrow(pit[which(pit$date == d$date & pit[, col_id] == d[, col_id]), ])
})
table(recorded)

#### Pull information from pit into Hunziker data
# Define relevant data
pit <- pit |> filter(fec_id %in% fish$fec_id)
# Pull information
ind <- match(fish$fec_id, pit$fec_id)
fish$section <- factor(pit$section[ind])
# Final checks
# * All wild fish
table(pit$hatchery)

#### Recode section labels from 1 to the number of sections sampled in each stream
# (For modelling purposes)
fish <- 
  lapply(split(fish, fish$stream), function(d){
    # d <- split(fish, fish$stream)[[1]]
    lookup <- data.frame(old = sort(unique(d$section)))
    lookup$new <- seq_len(nrow(lookup))
    d$rc_section <- lookup$new[match(d$section, lookup$old)]
    d
  }) |> 
  bind_rows() |>
  mutate(rc_section = factor(rc_section))
str(fish)
# Check the number of sections per stream
fish |> 
  group_by(stream) |> 
  summarise(n = length(unique(rc_section)))


#########################
#########################
#### Define tagging site locations properties

#### Clean tag sites data.frame
tag_sites_abb <-
  tag_sites_abb |>
  janitor::clean_names()

#### Link fish dataframe with tagging sites dataframe
tag_sites_dist <- 
  tag_sites_dist |> 
  mutate(stream = tag_sites_abb$location[match(stream, tag_sites_abb$code)],
         ss = paste(stream, section))
fish <- 
  fish |> 
  mutate(ss = paste(stream, section)) 
# Checks 
stopifnot(!any(is.na(tag_sites_dist$stream)))
head(tag_sites_dist$ss)
head(fish$ss)
unique(fish$ss) |> sort()
# fish$ss[!(fish$ss %in% tag_sites_dist$ss)] |> unique() |> sort()
stopifnot(all(fish$ss %in% tag_sites_dist$ss))

#### Get altitudes & distances to lake/antenna
# Altitudes
fish$altitude <- tag_sites_dist$alt[match(fish$stream, tag_sites_dist$stream)]
stopifnot(!any(is.na(fish$altitude)))
# Distances of antenna to lake
dist_ant <- 
  tag_sites_dist |> 
  filter(stream %in% fish$stream) |> 
  filter(section == "1") |>
  select(stream, dist_to_lake)
fish$dist_to_lake_from_ant <- dist_ant$dist_to_lake[match(fish$stream, dist_ant$stream)]
stopifnot(!any(is.na(fish$dist_to_lake_from_ant)))
# Distances of tagging sites to lake
fish$dist_to_lake_from_tag <- tag_sites_dist$dist[match(fish$ss, tag_sites_dist$ss)]
stopifnot(!any(is.na(fish$dist_to_lake)))
# Distances of tagging sections to antenna
fish$dist_to_ant_from_tag <- fish$dist_to_lake_from_tag - fish$dist_to_lake_from_ant
stopifnot(!any(is.na(fish$dist_to_ant_from_tag)))


#########################
#########################
#### Define migrants

#### Define migrants
migrants <- 
  fish |> 
  filter(migration == 1L) |> 
  select(length, mass, sex, migration_date, migration_yday, yday, stream, section, rc_section,
         altitude, dist_to_lake_from_ant, dist_to_lake_from_tag, dist_to_ant_from_tag) 


#########################
#########################
#### Save data

#### Save data
saveRDS(tag_sites_dist, here_data("distances-to-lake.rds"))
saveRDS(fish, here_data("fish.rds"))
saveRDS(migrants, here_data("migrants.rds"))


#### End of code.
#########################
#########################