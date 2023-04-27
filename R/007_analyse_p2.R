#########################
#########################
#### analyse_p2.R

#### Aims
# 1) Analyse P2 (the timing of migration)

#### Prerequisites
# 1) Process data for analysis (004_process_data.R)


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
library(mgcv)
library(DHARMa)
library(ggeffects)
library(prettyGraphics)
library(ggplot2)

#### Load data
source(here_r("001_define_global_param.R"))
source(here_r("002_define_helpers.R"))
fish     <- readRDS(here_data("fish.rds"))
migrants <- readRDS(here_data("migrants.rds"))


#########################
#########################
#### Data exploration

#### Examine relationships among variables
# Full matrix
pairs(migrants)
# Migration day ~ tagging date
plot(migrants$yday, migrants$migration_date)
cor(migrants$yday, migrants$migration_yday) # ~0.2
# Migration day ~ sex and length
plot(migrants$sex, migrants$migration_yday)
plot(migrants$length, migrants$migration_yday)

#### Check numbers of M/F migrants by stream
migrants_by_stream <- 
  migrants |> 
  group_by(stream, sex) |> 
  summarise(n = n())
# Simple statistics
utils.add::basic_stats(migrants_by_stream$n[migrants_by_stream$sex == "F"])
utils.add::basic_stats(migrants_by_stream$n[migrants_by_stream$sex == "M"])
# Barplot of counts
pp <- par(oma = c(4, 2, 2, 2))
barplot(migrants_by_stream$n, 
        names.arg = paste(migrants_by_stream$stream, migrants_by_stream$sex), 
        las = 2)
par(pp)

#### Check size distributions of migrants by stream
migrants |> 
  group_by(stream, sex) |> 
  summarise(min(length), median(length), mean(length), max(length))
ggplot(migrants) + 
  geom_histogram(aes(length)) + 
  facet_wrap(~stream)
boxplot(length ~ stream, data = migrants, las = 2)

#### Check distribution of migration dates by stream
# There is evidence of earlier migration in some streams
ggplot(migrants) + 
  geom_histogram(aes(migration_yday)) + 
  facet_wrap(~stream)


#########################
#########################
#### Implement modelling 

#### Fit models
# We will fit a model of day ~ sex * length
# * We expect migration timing to be delayed (larger) for smaller individuals 
# ... earlier migrants should be larger
# ... later migrants may be smaller (minimise predation risk)
# * We expect this decline to be faster for females
# ... I.e., small females will risk migration sooner than males
# ... of the equivalent size 

# mod_1: Initial model: sex-specific smoothers of length (with varying wiggliness)
# ... This model is closest to our hypothesis in that it allows M/F to vary in timing by size
mod_1 <- gam(migration_yday ~ 
               sex +
               s(length, by = sex, bs = "tp", m = 2) + 
               s(yday, bs = "cc") + 
               s(stream, bs = "re"), 
             knots = list(yday = c(0, 365)),
             family = gaussian(link = "identity"), data = migrants, 
             method = "REML")
plot(mod_1, pages = 1, scheme = 1, all.terms = TRUE)

# mod_2: mod_1 but with log(length) as predictor
mod_2 <- gam(migration_yday ~
               sex +
               s(log(length), by = sex, bs = "tp", m = 2) + 
               s(yday, bs = "cc") + 
               s(stream, bs = "re"), 
             knots = list(yday = c(0, 365)),
             family = gaussian(link = "identity"), data = migrants, 
             method = "REML")
plot(mod_2, pages = 1, scheme = 1, all.terms = TRUE)

# mod_3: individual sex AND stream specific smoothers (for effect of length)
# The initial models (e.g., mod_1) describe the data reasonably well 
# ... on average across streams but the predictions for individual streams
# ... are poor. There is evidence that the shape of the relationship between migration
# ... timing by body size varies substantially by stream; while normally negative
# ... the steepness of the shape varies. This suggests a model with effects that
# ... vary by stream would be appropriate. 
# Thus we amend the previous model to permit sex/stream specific smoothers
mod_3 <- gam(migration_yday ~ 
               sex +
               s(stream, bs = "re") + 
               # Group-level smoothers with different wiggliness:
               s(length, by = interaction(sex, stream), bs = "tp", m = 2) + 
               s(yday, bs = "cc"),
             knots = list(yday = c(0, 365)),
             family = gaussian(link = "identity"), data = migrants, 
             method = "REML")
plot(mod_3, pages = 1, scheme = 1, all.terms = TRUE)

# mod_4: mod_3 but with stream-specific smoothers only 
# According to AIC, this similar model is preferable
# I.e., the effect of size is similar between the sexes
# Putting P1 and P2 together, we can say that 
# * Females, and especially small females, are more likely to migrate than similarly sized males
# * Of the individuals that do migrate, larger individuals migrate sooner than smaller ones
# * Of the migrants, females do not appear to migrate sooner than males of equivalent size
mod_4 <- gam(migration_yday ~ 
               sex +
               s(stream, bs = "re") + 
               # Group level smoothers with different wiggliness:
               s(length, by = stream, bs = "tp", m = 2) + 
               s(yday, bs = "cc"),
             knots = list(yday = c(0, 365)),
             family = gaussian(link = "identity"), data = migrants, 
             method = "REML")
plot(mod_4, pages = 1, scheme = 1, all.terms = TRUE)

# mod_5: mod_4 but with Gamma likelihood
mod_5 <- gam(migration_yday ~ 
               sex +
               s(stream, bs = "re") + 
               s(length, by = stream, bs = "tp", m = 2) + 
               s(yday, bs = "cc"),
             knots = list(yday = c(0, 365)),
             family = Gamma(link = "log"), data = migrants, 
             method = "REML")

#### Compare models (with the same likelihood)
# mod_4 is preferred according to AIC
data.frame(mod = c(1, 2, 3, 4),
           aic = c(AIC(mod_1), AIC(mod_2), AIC(mod_3), AIC(mod_4))) |> 
  arrange(aic) |>
  mutate(delta_aic = aic - aic[1])
mod <- mod_4

##### Check model summary 
summary(mod)
AIC(mod)
plot(mod, pages = 1, scheme = 1, all.terms = TRUE)
sink(here_fig("tables", "migration-timing-mod.txt"))
print(summary.gam(mod, digits = 1))
sink()
  

#########################
#########################
#### Visualise observations and predictions

#### Quick visuals
if (FALSE) {
  pretty_predictions_1d(mod)
  do.call(gridExtra::grid.arrange, plot(ggpredict(mod), add.data = TRUE))
}

#### Visualise predictions across all streams
# This is appropriate for models 1--2

response  <- "migration_yday"
predictor <- "length"
mframe <- model.frame(mod)
mo <- seq(as.Date("2015-03-01"), as.Date("2015-06-01"), by = "months")
jd <- lubridate::yday(mo)
mo <- format(mo, "%b")
ylim <- c(min(migrants$migration_yday) - 2, max(migrants$migration_yday))
paa <- list(lim = list(x = lim_length, y = ylim), 
            axis = list(x = list(at = at_length), 
                        y = list(at = jd, labels = mo)))
add_axes_labels <- function(cex = 1.25, line = 2, ...) {
  mtext(side = 1, "Standard length (cm)", cex = cex, line = line, ...)
  mtext(side = 2, "Time of migration (months)", cex = cex, line = line, ...)
}
if (FALSE) {
  png(here_fig("migration-timing.png"), 
      height = 5, width = 6, units = "in", res = 600)
  pp <- par(oma = c(1, 1, 1, 2))
  pretty_blank(mframe, predictor, response, pretty_axis_args = paa)
  pred <- gen_pred(mod, predictor = predictor,
                   exclude = "s(stream)", newdata.guaranteed = TRUE)
  add_error_envelopes_by_sex(pred, predictor)
  pt.cex <- 0.5
  add_obs_by_sex(mframe, predictor, response, cex = 0.5)
  px <- par(xpd = NA)
  cex.leg <- 1.1
  legend(25, 180,
         title = "Sex", title.font = 2,
         legend = c("F", "M"), 
         lty = 1, pch = 21, 
         col = scales::alpha(cols, alpha_pt),
         pt.bg = scales::alpha(cols, alpha_pt), 
         pt.cex = pt.cex,
         box.lty = 3,
         cex = cex.leg)
  par(px)
  add_axes_labels()
  dev.off() 
}

#### Zoom into stream-specific relationships
# These plots are motivated by the diagnostic plots (see below)
# which suggest that Schutzenburnnen is an outlier
# It is difficult to see here that Schutzenburnnen is an outlier
# It is true that most migrants are a bit larger than expected
# across the range of migration days
# This is shown more clearly below where we plot the effects for each stream
# The intercept/average value of the smooth is higher in Schützenbrunnen
# than in the other streams.
plot(ggpredict(mod, terms = c("sex", "stream")), add.data = TRUE)
plot(ggpredict(mod, terms = c("length", "stream", "sex")), add.data = TRUE)
plot(ggpredict(mod, terms = c("yday", "stream")), add.data = TRUE)

#### Examine the partial effect of stream
# Define streams 
streams <- sort(unique(migrants$stream))
streams <- factor(streams, levels = levels(migrants$stream))
# Get partial effects 
# ... since type = "terms" in the code below
# ... the values of the other variables in newdata don't matter
nd <- 
  lapply(streams, function(stream) {
    d <- model.frame(mod)[1, ]
    d$stream <- stream
    d
  }) |> dplyr::bind_rows()
nd
nd$partial <- predict(mod, newdata = nd, type = "terms")[, "s(stream)"]
# Visualise partial effects
qn <- qqnorm(nd$partial)
basicPlotteR::addTextLabels(qn$x, qn$y, labels = nd$stream)
qqline(nd$partial)
shapiro.test(nd$partial)
# We can show these have been correctly computed
# e.g., via comparison to gratia::draw() output
# (Hover over the values and compare them to those computed in nd$partial, above)
if (FALSE) gratia::draw(mod) |> plotly::ggplotly()

#### Examine stream-specific predictions 
png(here_fig("migration-timing-by-stream.png"), 
    height = 5, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(2, 4), oma = c(3, 3, 1, 1), mar = c(2, 2, 2, 2))
lapply(seq_len(length(unique(fish$stream))), function(i) {
  
  #### Define variables
  response  <- "migration_yday"
  predictor <- "length"
  mframe <- model.frame(mod)
  stream <- sort(unique(fish$stream))[i]
  mframe_for_stream <- mframe[mframe$stream == stream, ]
  
  #### Create plot
  pretty_blank(mframe, predictor, response, pretty_axis_args = paa)
  pred <- gen_pred(mod, stream, predictor, exclude = "s(stream)")
  add_error_envelopes_by_sex(pred, predictor)
  legend("bottomright", legend = paste0("n = ", nrow(mframe_for_stream)), bty = "n")
  # legend("bottomright", legend = bquote(italic(n) * " = " * .(nrow(mframe_for_stream))), bty = "n")
  add_obs_by_sex(mframe_for_stream, predictor, response)
  mtext(side = 3, stream, font = 2)
}) |> invisible()

#### Add axes
add_axes_labels(outer = TRUE, line = 1)
par(pp)
dev.off()


#########################
#########################
#### Model diagnostics 

#### GAM checks
# * mod_4: shows some evidence of overdispersion
# * mod_5: broadly similar
pp <- par(mfrow = c(2, 2))
gam.check(mod, rep = 1e3)
par(pp)
k.check(mod)

#### Simulate DHARMa residuals
# mod_4: ok. 
# mod_5: ok.
mframe <- model.frame(mod)
res    <- simulateResiduals(mod, refit = FALSE, plot = TRUE, re.form = NULL)

#### Check residuals versus predictors
plotResiduals(res, form = mframe$sex)
plotResiduals(res, form = mframe[, "length"])
plotResiduals(res, form = mframe[, "yday"])
plotResiduals(res, form = mframe[, "stream"])

#### Run additional DHARMa checks
testDispersion(res)


#### End of code. 
#########################
#########################