#########################
#########################
#### analyse_h2.R

#### Aims
# 1) Analyse H2 (the size of migrants)

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


#########################
#########################
#### Data exploration

#### Define migrants
migrants <- 
  fish |> 
  filter(migration == 1L) |> 
  select(length, sex, migration_date, migration_yday, yday, stream) 
# migrants$length <- migrants$length * 10 

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

### Check size distributions of migrants by stream
migrants |> 
  group_by(stream, sex) |> 
  summarise(min(length), median(length), mean(length), max(length))
ggplot(migrants) + 
  geom_histogram(aes(length)) + 
  facet_wrap(~stream)
boxplot(length ~ stream, data = migrants, las = 2)

#### Check migration dates of individuals by stream
# There is evidence of earlier migration in some streams
ggplot(migrants) + 
  geom_histogram(aes(migration_yday)) + 
  facet_wrap(~stream)

#### Check relationships between tagging and migration dates
plot(migrants$yday, migrants$migration_date)
cor(migrants$yday, migrants$migration_yday) # ~0.2


#########################
#########################
#### Implement modelling 

#### Fit models
# We will fit a model of length ~ sex * time of migration: 
# * We expect larger individuals to migrate earlier
# * Smaller individuals may migrate later because predation risk is higher
# * This relationship may be steeper for females, i.e., 
# ... while large individuals of both sexes may migrate early
# ... as time passes, females of a given size may choose to migrate before
# ... males of the equivalent size 
# ... i.e., small females risk migration earlier

# mod_1: Initial model
mod_1 <- gam(length ~ 
               sex +
               s(migration_yday, by = sex, bs = c("cc", "fs")) + 
               s(yday, bs = "cc") + 
               s(stream, bs = "re"), 
             knots = list(yday = c(0, 365), 
                          migration_yday = c(0, 365)),
             family = gaussian, data = migrants, 
             method = "REML")
plot(mod_1, pages = 1, scheme = 1, all.terms = TRUE)

# mod_2: mod_1 without the factor (sex) - smooth (migration day) interaction
mod_2 <- gam(length ~ 
               sex +
               s(migration_yday, bs = "cc") + 
               s(yday, bs = "cc") + 
               s(stream, bs = "re"), 
             knots = list(yday = c(0, 365), 
                          migration_yday = c(0, 365)),
             family = gaussian, data = migrants, 
             method = "REML")
plot(mod_2, pages = 1, scheme = 1, all.terms = TRUE)

# mod_3: smoothers for each stream
mod_3 <- gam(length ~ 
               sex +
               s(migration_yday, by = stream, bs = c("cc", "re")) + 
               s(yday, bs = "cc"),
             knots = list(yday = c(0, 365), 
                          migration_yday = c(0, 365)),
             family = gaussian, data = migrants, 
             method = "REML")
plot(mod_3, pages = 1, scheme = 1, all.terms = TRUE)

# mod_4: mod_2 with Gamma likelihood 
mod_4 <- gam(length ~ 
               sex +
               s(migration_yday, bs = "cc") + 
               s(yday, bs = "cc") + 
               s(stream, bs = "re"), 
             knots = list(yday = c(0, 365), 
                          migration_yday = c(0, 365)),
             family = Gamma(link = "log"), data = migrants, 
             method = "REML")
plot(mod_4, pages = 1, scheme = 1, all.terms = TRUE)

#### Compare models 
# Check AIC (for models that share the same likelihood)
# * mod_3 has lowest AIC, but predictions look too wiggly (overfitting)
# * mod_2 without the interaction is preferable to mod_1 (by ~2.2 AIC)
# * Residual diagnostics are problematic for mod_2, 
# ... e.g., evidence of overdispersion in gam.check() and issues in DHARMa resids
data.frame(mod = c(1, 2, 3),
           aic = c(AIC(mod_1), AIC(mod_2), AIC(mod_3))) |>
  arrange(aic)
AIC(mod_1) - AIC(mod_2)

##### Check model summary 
mod <- mod_4
summary(mod)
AIC(mod)
plot(mod, pages = 1, scheme = 1, all.terms = TRUE)
sink(here_fig("tables", "migrant-length-mod.txt"))
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
png(here_fig("migrant-lengths.png"), 
    height = 5, width = 6, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 1))
response  <- "length"
predictor <- "migration_yday"
mframe <- model.frame(mod)
mo <- seq(as.Date("2015-03-01"), as.Date("2015-06-01"), by = "months")
jd <- lubridate::yday(mo)
mo <- format(mo, "%b")
xlim <- c(min(migrants$migration_yday) - 2, max(migrants$migration_yday))
paa <- list(lim = list(x = xlim, y = NULL), 
            axis = list(x = list(at = jd, labels = mo), 
                        y = list(NULL)))
pretty_blank(mframe, predictor, response, pretty_axis_args = paa)
pred <- gen_pred(mod, predictor = predictor,
                 exclude = "s(stream)", newdata.guaranteed = TRUE)
add_error_envelopes_by_sex(pred, predictor)
pt.cex <- 0.5
add_obs_by_sex(mframe, predictor, response, cex = 0.5)
add_axes_labels <- function(cex = 1.25, line = 2, ...) {
  mtext(side = 1, "Time of migration (months)", cex = cex, line = line, ...)
  mtext(side = 2, "Standard length (cm)", cex = cex, line = line, ...)
}
px <- par(xpd = NA)
cex.leg <- 1.1
legend(max(migrants$migration_yday) - 10, 25,
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
plot(ggpredict(mod, terms = c("migration_yday", "stream", "sex")), add.data = TRUE)
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
# We can show these have been correctly computed
# e.g., via comparison to gratia::draw() output
# (Hover over the values and compare them to those computed in nd$partial, above)
if (FALSE) gratia::draw(mod) |> plotly::ggplotly()

#### Examine stream-specific predictions 
png(here_fig("migrant-lengths-by-stream.png"), 
    height = 5, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(2, 4), oma = c(3, 3, 1, 1), mar = c(2, 2, 2, 2))
lapply(seq_len(length(unique(fish$stream))), function(i) {
  
  #### Define variables
  response  <- "length"
  predictor <- "migration_yday"
  mframe <- model.frame(mod)
  stream <- sort(unique(fish$stream))[i]
  
  #### Create plot
  pretty_blank(mframe, predictor, response, pretty_axis_args = paa)
  pred <- gen_pred(mod, stream, predictor)
  add_error_envelopes_by_sex(pred, predictor)
  add_obs_by_sex(mframe[mframe$stream == stream, ], predictor, response)
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
# * mod_2: shows evidence of overdispersion
# * mod_4: also overdispersion, but less strong
pp <- par(mfrow = c(2, 2))
gam.check(mod, rep = 1e3)
par(pp)
k.check(mod)

#### Simulate DHARMa residuals
# mod_2: evidence of overdispersion and other issues
# mod_4: residual diagnostics appear better, 
# ... although there is still some issue e.g.,
# ... pattern in residual versus predicted vals
mframe <- model.frame(mod)
res    <- simulateResiduals(mod, refit = FALSE, plot = TRUE, re.form = NULL)

#### Check residuals versus predictors
plotResiduals(res, form = mframe$sex)
plotResiduals(res, form = mframe[, "migration_yday"])
plotResiduals(res, form = mframe[, "yday"])

#### Run additional DHARMa checks
testDispersion(res)


#### End of code. 
#########################
#########################