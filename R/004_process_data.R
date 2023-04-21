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


#########################
#########################
#### Process data

#### Select and process relevant columns
fish <- 
  fish |> 
  mutate(id = row_number()) |>
  select(id, 
         stream,
         sex = sex, 
         length = SL,
         day = tag.date.jul, 
         migration = downmig.bin
         ) |>
  mutate(id = factor(id), 
         sex = factor(sex),
         stream = factor(stream),
         day = as.integer(day), 
         migration = factor(migration), 
         )

### Save data
saveRDS(fish, here_data("fish.rds"))


#### End of code.
#########################
#########################