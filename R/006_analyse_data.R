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
tag_sites_dist <- readRDS(here_data("distances-to-lake.rds"))


#########################
#########################
#### Data processing

#### Define altitude (a) of antenna and (b) other sections
tag_sites_dist_1 <- tag_sites_dist[tag_sites_dist$section == 1, ]
fish$altitude_1  <- tag_sites_dist_1$alt[match(fish$stream, tag_sites_dist_1$stream)] |> round()
fish$altitude    <- tag_sites_dist$alt[match(fish$ss, tag_sites_dist$ss)] |> round()


#########################
#########################
#### Summary tables

#### Summary table of stream-level data 
fish |> 
  group_by(stream) |> 
  summarise(altitude_1 = altitude_1[1], 
            visits = length(unique(date)), 
            date = str_range(format(sort(date), "%d-%b")), 
            n = n(), 
            length = str_range(length), 
            n_migrant = length(which(migration == 1L))) |>
  mutate(pr_migrant = add_lagging_point_zero(round(n_migrant/n, 2), 2)) |>
  arrange(stream) |>
  tidy_write(here_fig("tables", "stream_summary.txt"))

#### Summary of section-level data 
# Distances of antenna from lake
dist_ant <- 
  tag_sites_dist |> 
  filter(stream %in% fish$stream) |> 
  filter(section == "1") |>
  select(stream, dist_to_lake)
utils.add::basic_stats(dist_ant$dist_to_lake)
# Distances of tagging sections to the lake and other information (as above)
dist_tag <- 
  fish |> 
  group_by(stream, section) |>
  mutate(section = as.numeric(as.character(section)), 
         visits = length(unique(date)), 
         date = str_range(format(sort(date), "%d-%b")), 
         n = n(), 
         length = str_range(length), 
         n_migrant = length(which(migration == 1L)), 
         pr_migrant = add_lagging_point_zero(round(n_migrant/n, 2), 2)) |>
  slice(1L) |>
  select(stream, section, 
         altitude, dist_to_ant_from_tag, dist_to_lake_from_tag, 
         date,
         visits,
         n, 
         length, 
         n_migrant, 
         pr_migrant
  ) |> 
  arrange(stream, section)
# Do distances uniquely define each section? Essentially, yes. 
nrow(dist_tag)
length(unique(dist_tag$dist_to_ant_from_tag))
length(unique(dist_tag$dist_to_lake_from_tag))
# Do tagging dates uniquely define each stream or section? 
# ... Streams were fished on 2-3 occasions
# ... On any given date, 1-2 streams were fished
length(unique(fish$stream))
fish |> 
  group_by(stream, date) |> 
  slice(1L) |>
  ungroup() |>
  dplyr::select(stream, date) |> 
  arrange(date) 
# Check, for each stream, the number of unique tagging dates shared with other streams
sapply(split(fish, fish$stream), function(d) {
  other <- fish[!(fish$stream %in% d$stream), ]
  any(d$date %in% other$date)
  table(unique(d$date) %in% other$date)["TRUE"]
})
# Number of captures per section
utils.add::basic_stats(dist_tag$n)
# Distances of tagging sections to lake/antenna
utils.add::basic_stats(dist_tag$dist_to_lake_from_tag)
utils.add::basic_stats(dist_tag$dist_to_ant_from_tag)
# Tidy table
dist_tag$dist_to_ant_from_tag  <- round(dist_tag$dist_to_ant_from_tag)
dist_tag$dist_to_lake_from_tag <- round(dist_tag$dist_to_lake_from_tag)
dist_tag |> tidy_write(here_fig("tables", "stream_and_section_summary.txt"))
  
#### Summary table of tagged fish 
fish |> 
  select(id, date, stream, section, sex, length, migration) |> 
  arrange(date, stream, section, sex, length, migration) |> 
  mutate(id = row_number()) |>
  tidy_write(here_fig("tables", "tag_data.txt"))
  

#########################
#########################
#### Summary statistics

#### Stream profiles
tag_sites_dist |>
  group_by(stream) |> 
  arrange(section, .by_group = TRUE) |>
  mutate(alt = alt - alt[1]) |>
  ggplot() +
  geom_point(aes(dist_to_lake, alt, colour = section)) + 
  # scale_x_continuous(expand = c(0, 0)) + 
  facet_wrap(~stream, scales = "free")

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

#### For migrant individuals, check sex, sizes and the timing of migration
# Sex distribution (calculate approximately comparable statistics to Aarestrup et al., 2018)
aare <- fish[fish$migration == 1 & fish$length >= 9 & fish$length <= 21, ]
tbl <- table(aare$sex)
tbl/sum(tbl) * 100
mod <- glm(length ~ sex, data = aare)
pretty_predictions_1d(mod)
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

#### For migrants, compare the proportion of males versus females
# Calculate proportions
prop_ms <- 
  fish |>
  filter(migration == 1L) |>
  mutate(bin = cut(length, seq(min(length), max(length), by = size_class_cm))) |>
  group_by(bin) |> 
  summarise(n = n(), 
            nm = length(which(sex == "M")), 
            nf = length(which(sex == "F")), 
            pm = nm/n, 
            pf = nf/n,) |> 
  mutate(bin = parse_cut(bin))
# Define plotting properties
squash_param <- 20
prop_ms$cex <- prop_ms$n/squash_param
tiny <- prop_ms[prop_ms$cex <= 0.5, ]
# Make plot
# ... This code is adapted from the plotting routines in analyse_p1.R
png(here_fig("migrant-prop.png"), 
    height = 5, width = 7, units = "in", res = 600)
pp <- par(oma = c(0, 0, 0, 4))
paa <- list(lim = list(x = lim_length, 
                       y = c(-0.1, 1.1)),
            axis = list(x = list(at = at_length), 
                        y = list(at = seq(0, 1, 0.2))), 
            control_axis = list(las = TRUE, cex.axis = 1.25)
            )
pretty_plot(prop_ms$bin, prop_ms$pm,
            pretty_axis_args = paa,
            xlab = "", ylab = "",
            type = "n")
lines(prop_ms$bin, prop_ms$pf, 
      col = scales::alpha(cols["F"], alpha_fit))
points(prop_ms$bin, prop_ms$pf, 
       pch = 21, 
       bg = scales::alpha(cols["F"], alpha_pt), 
       col = scales::alpha(cols["F"], alpha_pt),
       cex = prop_ms$n/squash_param)
points(tiny$bin, tiny$pf, col = "red", lwd = 0.5)
lines(prop_ms$bin, prop_ms$pm)
points(prop_ms$bin, prop_ms$pm, 
       pch = 21, bg = cols["M"], col = cols["M"], 
       cex = prop_ms$n/squash_param)
points(tiny$bin, tiny$pm, col = "red", lwd = 0.5)
pos_f <- which(fish$migration == 1L & fish$sex == "F")
points(fish$length[pos_f], rep(1, length(pos_f)), 
    col = scales::alpha(cols["F"], 0.5), cex = 0.5)
pos_m <- which(fish$migration == 1L & fish$sex == "M")
points(fish$length[pos_m], rep(0, length(pos_m)), 
       col = scales::alpha(cols["M"], 0.5), cex = 0.5)
add_axes_labels <- function(cex = 1.25, line = 2, ...) {
  mtext(side = 1, "Standard length (cm)", cex = cex, line = line, ...)
  mtext(side = 2, "Proportion", cex = cex, line = line, ...)
}
add_axes_labels()
px <- par(xpd = NA)
cex.leg <- 1.1
legend(25, 1,
       title = "Sex", title.font = 2,
       legend = c("F", "M"), 
       lty = 1, pch = 21, 
       col = scales::alpha(cols, alpha_pt),
       pt.bg = scales::alpha(cols, alpha_pt), 
       bty = "n",
       cex = cex.leg)
ns <- c(10, 25, 50)
legend(25, 0.6,
       title = "Sample (N)", title.font = 2,
       legend = ns, 
       pch = 21, 
       col = scales::alpha(cols[2], alpha_pt),
       pt.bg = scales::alpha(cols[2], alpha_pt), 
       pt.cex = ns/squash_param,
       y.intersp = 1.4,
       cex = cex.leg, 
       bty = "n")
par(px)
par(pp)
dev.off()



#### End of code.
#########################
#########################