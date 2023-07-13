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
# Define summary data 
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
  tidy_write(here_fig("summaries", "stream_summary.txt"))
# Check the number of sections per stream with observations 
fish |> 
  group_by(stream) |> 
  summarise(n_sections = length(unique(section)))

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
dist_tag |> tidy_write(here_fig("summaries", "stream_and_section_summary.txt"))
  
#### Summary table of tagged fish 
fish |> 
  select(id, date, stream, section, sex, length, migration) |> 
  arrange(date, stream, section, sex, length, migration) |> 
  mutate(id = row_number()) |>
  tidy_write(here_fig("summaries", "tag_data.txt"))
  

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

#### Relationship between size metrics
plot(fish$length, fish$mass)
sl <- seq(min(fish$length), max(fish$length), length.out = 100)
lm(mass ~ poly(length, 3), data = fish)
ggplot(fish, aes(x = length, y = mass)) + 
  geom_point() + 
  geom_smooth(method = lm, formula = y ~ poly(x, 3), data = fish)

#### Size distribution of fish 
png(here_fig("fish-length-distribution.png"), 
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

#### Relationship between SL and TL 
# Define SL and TL
fish$SL <- fish$length
fish$TL <- fish$total_length
# A minimum TL of 105 mm was used
# This corresponds to a smallest size of 89 mm
range(fish$SL, na.rm = TRUE)
range(fish$TL, na.rm = TRUE)
sort(fish$SL[fish$TL == min(fish$TL)])
# ... Which is very similar to the average (expected) standard length
mod <- lm(SL ~ TL, data = fish)
ptl <- 105/10
psl <- predict(mod, newdata = data.frame(TL = ptl))
if (FALSE) {
  pp <- par(mfrow = c(2, 2))
  plot(mod)
  par(pp)
}
# Visualise SL vs TL 
png(here_fig("sl-vs-tl.png"), 
    height = 5, width = 5, units = "in", res = 600)
xlim <- c(7, 30)
ylim <- c(7, 30)
axis_ls <- 
  prettyGraphics::pretty_predictions_1d(mod, 
                                        xlim = xlim, ylim = ylim, 
                                        add_points = list(cex = 0.5, lwd = 0.5, col = "grey20"),
                                        add_xlab = list(text = ""), add_ylab = list(text = ""), 
                                        add_main = NULL)
# lines(rep(ptl, 2), c(ylim[1], psl), lty = 3, col = "red", lwd = 2)
# lines(c(xlim[1], ptl), rep(psl, 2), lty = 3, col = "red", lwd = 2)
len <- 0.05
col <- "red"
lwd <- 1.5
arrows(x0 = ptl, x1 = ptl, y0 = ylim[1], y1 = psl, 
       length = 0, col = col, lwd = lwd)
arrows(x0 = ptl, x1 = xlim[1], y0 = psl, y1 = psl, 
       length = len, col = col, lwd = lwd)
lines(xlim, xlim, lty = 3)
cex.mtext <- 1.25
mtext(side = 1, "Total length (cm)", cex = cex.mtext, line = 2)
mtext(side = 2, "Standard length (cm)", cex = cex.mtext, line = 2.5)
dev.off()

#### Count the number of streams/captures per stream
length(unique(fish$stream))
table(fish$stream)

#### Check the correlation between tagging date and length
# (Are fish thaty were tagged earlier smaller?)
plot(fish$yday, fish$length)
cor(fish$yday, fish$length)

#### Calculate the observed proportion of migrants by sex and size (ss)
size_class <- list(length = 1.5,
                   mass = 1.5)
size_vars <- c("length", "mass")
prop_ss <- 
  lapply(size_vars, function(size_var) {
    fish$size <- fish[, size_var]
    out <- 
      fish |>
      # group_by(sex) |>
      filter(!is.na(size)) |>
      mutate(bin = cut(size, seq(min(size), max(size) + size_class[[size_var]], 
                                 by = size_class[[size_var]]), 
                       include.lowest = TRUE, right = FALSE)
             ) |>
      # ungroup() |>
      group_by(sex, bin) |>
      summarise(pr = length(which(migration == "1"))/n(), 
                n = n()) |>
      mutate(size = parse_cut(bin), 
             col = as.character(scales::alpha(cols, alpha_pt)[as.character(sex)])) |> 
      select(sex, n, bin, size, pr, col) |>
      as.data.frame()
    stopifnot(!any(is.na(out$bin)))
    out
  })
names(prop_ss) <- size_vars
saveRDS(prop_ss, here_data("prop_ss.rds"))

#### Calculate the observed proportion of migrants by sex, size and stream (sss)
# This code is the same as above but stream is included in the grouping structures
prop_sss <- 
  lapply(size_vars, function(size_var) {
    fish$size <- fish[, size_var]
    out <- 
      fish |>
      filter(!is.na(size)) |>
      group_by(stream, sex) |>
      mutate(bin = cut(size, seq(min(size), max(size) + size_class[[size_var]], 
                                 by = size_class[[size_var]]), 
                       include.lowest = TRUE, right = FALSE)
             ) |>
      ungroup() |>
      group_by(stream, sex, bin) |>
      summarise(pr = length(which(migration == "1"))/n(), 
                n = n()) |>
      mutate(size = parse_cut(bin), 
             col = as.character(scales::alpha(cols, alpha_pt)[as.character(sex)])) |> 
      select(stream, sex, n, bin, size, pr, col) |>
      as.data.frame()
    stopifnot(!any(is.na(out$bin)))
    out
  })
names(prop_sss) <- size_vars
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
lapply(c("length", "mass"), function(size_var) {
  fish$size <- fish[, size_var]
  prop_ms <- 
    fish |>
    filter(migration == 1L) |>
    filter(!is.na(size)) |>
    mutate(bin = cut(size, seq(min(size), max(size) + size_class[[size_var]], 
                               by = size_class[[size_var]]), 
                     include.lowest = TRUE, right = FALSE)
    ) |>
    group_by(bin) |> 
    summarise(n = n(), 
              nm = length(which(sex == "M")), 
              nf = length(which(sex == "F")), 
              pm = nm/n, 
              pf = nf/n) |> 
    mutate(bin = parse_cut(bin))
  stopifnot(!any(is.na(prop_ms$bin)))
  # Define plotting properties
  squash_param <- 20
  prop_ms$cex <- prop_ms$n/squash_param
  tiny <- prop_ms[prop_ms$cex <= 0.5, ]
  # Make plot
  # ... This code is adapted from the plotting routines in analyse_p1.R
  png(here_fig("relationships", "fig", glue::glue("by-{size_var}"), "migrant-prop.png"), 
      height = 5, width = 7, units = "in", res = 600)
  pp <- par(oma = c(0, 0, 0, 4))
  paa <- list(lim = list(x = lim_size[[size_var]], 
                         y = c(-0.1, 1.1)),
              axis = list(x = list(at = at_size[[size_var]]), 
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
    if (size_var == "length") {
      xlab <- "Standard length (cm)"
    } else if (size_var == "mass") {
      xlab <- "Mass (10 g)"
    }
    mtext(side = 1, xlab, cex = cex, line = line, ...)
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
  
})


#### End of code.
#########################
#########################