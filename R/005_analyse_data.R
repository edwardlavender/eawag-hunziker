#########################
#########################
#### analyse_data.R

#### Aims
# 1) Analyse sample data 
# ... * Produce summary tables
# ... * Produce summary statistics

#### Prerequisites
# 1) Process data via process_data.R


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
library(prettyGraphics)
library(ggplot2)

#### Load data
source(here_r("001_define_global_param.R"))
source(here_r("002_define_helpers.R"))
fish <- readRDS(here_data("fish.rds"))


#########################
#########################
#### Summary tables

#### Summary table of number of sampled fish per stream
fish |> 
  group_by(stream) |> 
  summarise(visits = length(unique(date)), 
            date = str_range(format(sort(date), "%d-%b")), 
            n = n(), 
            length = str_range(length), 
            n_migrant = length(which(migration == 1L))) |>
  mutate(pr_migrant = add_lagging_point_zero(round(n_migrant/n, 2), 2)) |>
  arrange(stream) |>
  tidy_write(here_fig("tables", "tag_summary.txt"))
  
#### Summary table of tagged fish 
fish |> 
  select(id, date, stream, sex, length, migration) |> 
  arrange(date, stream, sex, length) |> 
  mutate(id = row_number()) |>
  tidy_write(here_fig("tables", "tag_data.txt"))
  

#########################
#########################
#### Summary statistics

#### Size distribution of fish 
png(here_fig("fish-size-distribution.png"), 
    height = 4, width = 5, units = "in", res = 600)
bw <- 0.2
bks <- seq(plyr::round_any(min(fish$length), bw, floor), 
           plyr::round_any(max(fish$length), bw, ceiling), 
           by = bw)
h <- hist(fish$length, breaks = bks, 
     xlim = range(bks), axes = FALSE,
     xlab = "", ylab = "", main = "",
     col = scales::alpha("lightgrey", 0.75))
axis(side = 1, 
     at = prettyGraphics::pretty_seq(bks, lim = range(bks), pretty = list(n = 10))$at, 
     pos = 0)
axis(side = 2, 
     at = prettyGraphics::pretty_seq(h$counts, pretty = list(n = 10))$at,
     pos = min(bks), 
     las = TRUE)
mtext(side = 1, "Standard length (cm)", line = 2)
mtext(side = 2, "Frequency", line = 2)
HDInterval::hdi(fish$length, allowSplit = TRUE)
dev.off()

#### Count the number of streams/captures per stream
length(unique(fish$stream))
table(fish$stream)

#### Check the correlation between tagging date and length
# (Are fish thaty were tagged earlier smaller?)
plot(fish$yday, fish$length)
cor(fish$yday, fish$length)

#### Calculate the observed proportion of migrants by sex and size (ss)
size_class_cm <- 1.5
prop_ss <- 
  fish |>
  group_by(sex) |>
  mutate(bin = cut(length, seq(min(length), max(length), by = size_class_cm))) |>
  ungroup() |>
  group_by(sex, bin) |>
  summarise(pr = length(which(migration == "1"))/n(), 
            n = n()) |>
  mutate(length = parse_cut(bin), 
         col = as.character(scales::alpha(cols, alpha_pt)[as.character(sex)])) |> 
  select(sex, n, bin, length, pr, col) |>
  filter(!is.na(bin)) |> 
  as.data.frame()
saveRDS(prop_ss, here_data("prop_ss.rds"))

#### Calculate the observed proportion of migrants by sex, size and stream (sss)
# This code is the same as above but stream is included in the grouping structures
prop_sss <- 
  fish |>
  group_by(stream, sex) |>
  mutate(bin = cut(length, seq(min(length), max(length), by = size_class_cm))) |>
  ungroup() |>
  group_by(stream, sex, bin) |>
  summarise(pr = length(which(migration == "1"))/n(), 
            n = n()) |>
  mutate(length = parse_cut(bin), 
         col = as.character(scales::alpha(cols, alpha_pt)[as.character(sex)])) |> 
  select(stream, sex, n, bin, length, pr, col) |>
  filter(!is.na(bin)) |> 
  as.data.frame()
saveRDS(prop_sss, here_data("prop_sss.rds"))

#### For migrant individuals, check sizes and the timing of migration
# Size distribution
ggplot(fish) + 
  geom_histogram(aes(length)) + 
  facet_wrap(~sex + migration)
# Timing 
ggplot(fish) + 
  geom_histogram(aes(date)) + 
  facet_wrap(~sex + migration)
# Timing statistics
lapply(list(min, median, mean,  max), function(f) {
  f(fish$migration_date[fish$migration == 1L])
})



#### End of code.
#########################
#########################