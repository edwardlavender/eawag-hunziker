#########################
#########################
#### analyse_p1.R

#### Aims
# 1) Analyse P1 (migration Pr ~ sex * size)

#### Prerequisites
# 1) Process data for analysis (004_process_data.R)
# 2) Calculate observed proportions of migrants (005_analyse_data.R)


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
library(lme4)
library(DHARMa)
library(ggeffects)
library(prettyGraphics)
# {mgcViz} is required for {DHARMa} with GAMs
if (!requireNamespace("mgcViz", quietly = TRUE)) {
  install.packages("mgcViz")
}

#### Load data
source(here_r("001_define_global_param.R"))
source(here_r("002_define_helpers.R"))
fish     <- readRDS(here_data("fish.rds"))
prop_ss  <- readRDS(here_data("prop_ss.rds"))
prop_sss <- readRDS(here_data("prop_sss.rds"))


#########################
#########################
#### Implement modelling 

#### GLMM: model migration (0, 1) as a function of size and sex and other factors
# Controls
# * We will re-scale length to improve model identifiability 
mod_1 <- glmer(migration ~ 
                 log(length) * sex + 
                 yday + (1|stream) + (1|stream:rc_section), 
             data = fish, family = binomial(link = "logit"))

#### GAMM: model migration (0, 1) as a function of size and sex and other factors
# Motivation 
# * Examination of the model structure for mod_1 suggests 
# ... a non linear relationship for examples at small sizes 
# ... suggesting a GAM model would be more appropriate here. 
# Structure
# * Length is scaled via log()
# ... Scaling length is not necessary 
# ... and makes little visible difference to predictions
# ... but marginally improves model fit (AIC)
# ... scale() produces poor fit on predictive plot 
# * REML or ML recommended for fitting GAMs

# GAM parameters 
gamma <- 1

# mod_2: sex-specific smoothers of the length effect with their own degree of wiggliness
# ... This model is closest to our expectation that sexes may respond differently 
# ... We also include a nested random effect of section within stream.
# ... To be nested means an interaction term. 
# ... i.e., We can’t predict the results for section 'A' in any other stream. 
# ... hence, there is no overall 'section' effect.
mod_2 <- gam(migration ~ 
               sex + 
               s(log(length), by = sex, bs = "tp", m = 2) + 
               s(yday, k = 5, bs = "cc") + 
               s(stream, bs = "re") + 
               s(stream, rc_section, bs = "re"), 
             knots = list(yday = c(0, 365)),
             family = binomial, data = fish, 
             gamma = gamma,
             method = "REML")

# mod_3: sex-specific and stream-specific smoothers 
mod_3 <- gam(migration ~ 
               sex + 
               s(log(length), by = interaction(sex, stream), bs = "tp", m = 2) + 
               s(yday, k = 5, bs = "cc") + 
               s(stream, bs = "re") + 
               s(stream, rc_section, bs = "re"), 
             knots = list(yday = c(0, 365)),
             family = binomial, data = fish, 
             gamma = gamma,
             method = "REML")

#### Compare model AICs
data.frame(mod = c(1, 2, 3),
           aic = c(AIC(mod_1), AIC(mod_2), AIC(mod_3))) |>
  arrange(aic) |> 
  mutate(delta = aic - aic[1])
# Choose model
mod <- mod_2
is_glmer <- inherits(mod, "merMod")

##### Check model summary 
summary(mod)
AIC(mod)
if (is_glmer) {
  equatiomatic::extract_eq(mod)
  MuMIn::r.squaredGLMM(mod)
} else {
  fit <- plot(mod, pages = 1, scheme = 1, all.terms = TRUE, 
              residuals = TRUE, seWithMean = FALSE)
  if (FALSE) {
    pretty_smooth_1d(fit, select = seq_len(length(fit) - 1), 
                     add_resid = list(pch = "."), add_rug = list())
  }
}
sink(here_fig("tables", "migration-prob-mod.txt"))
print(summary.gam(mod, digits = 1))
sink()


#########################
#########################
#### Visualise observations and predictions

#########################
#### Visualise predictions/observations across all streams

#### Set up plot
png(here_fig("migration-prob.png"), 
    height = 5, width = 7, units = "in", res = 600)
pp <- par(oma = c(0, 0, 0, 4))
predictor <- "length"
response  <- "pr"

#### Create blank plot
# We define the vertical axis to range slightly beyond c(0, 1)
# This enables us to add rugs at the top for males/females (see add_outcomes())
paa <- list(lim = list(x = lim_length, 
                       y = c(-0.1, 1.1)),
            axis = list(x = list(at = at_length), 
                        y = list(at = seq(0, 1, 0.2))), 
            control_axis = list(las = TRUE, cex.axis = 1.25)
            )
pretty_blank(prop_ss, predictor, response, pretty_axis_args = paa)

#### Add model predictions with SEs
# Define sequence of body sizes for prediction
n <- 100
ms <- fish$length[fish$sex == "M"]
fs <- fish$length[fish$sex == "F"]
ms <- seq(min(ms), max(ms), length = n)
fs <- seq(min(fs), max(fs), length = n)
# Add CIs
if (is_glmer) {
  
  #### Generate model predictions across all body sizes (via ggeffects::ggpredict())
  # We can conveniently generate predictions via ggpredict()
  # ... (See https://strengejacke.github.io/ggeffects/articles/introduction_randomeffects.html)
  # ... For predictions, type = "fixed" is the default
  # ... I.e., predictions are on the population-level & do not account for the random effect vars 
  # ... Intervals are confidence intervals for the predicted values
  # ... CIs are calculated as described here:
  # ... https://strengejacke.github.io/ggeffects/articles/ggeffects.html
  # Generate predictions across the range of body sizes for males/females
  pr <- ggpredict(mod, terms = c("length [all]", "sex"), type = "fixed")
  pm <- ggpredict(mod, terms = "length [ms]", condition = c(sex = "M"))
  pf <- ggpredict(mod, terms = "length [fs]", condition = c(sex = "F"))
  if (FALSE) {
    plot(pr, add.data = TRUE)
    plot(pm, add.data = TRUE)
    plot(pf, add.data = TRUE)
  }
  pm <- as.data.frame(pm)
  pf <- as.data.frame(pf)
  
  #### Add error envelopes
  add_error_envelope(pm$x, ci = list(fit = pm$predicted, lowerCI = pm$conf.low, upperCI = pm$conf.high), 
                     add_fit = list(col = scales::alpha(cols["M"], alpha_fit)), 
                     add_ci = list(col = scales::alpha(cols["M"], alpha_ci), border = FALSE))
  add_error_envelope(pf$x, ci = list(fit = pf$predicted, lowerCI = pf$conf.low, upperCI = pf$conf.high), 
                     add_fit = list(col = scales::alpha(cols["F"], alpha_fit)), 
                     add_ci = list(col = scales::alpha(cols["F"], alpha_ci), border = FALSE))
  
} else {
  
  #### Generate model predictions via predict.gam()
  # We can exclude the effect of stream via the exclude argument
  # predict(mod, type = "terms") |> head()
  nd <- data.frame(sex = factor(c(rep("F", n), rep("M", n))),
                   yday = median(fish$yday),
                   length = c(fs, ms))
  p <- predict(mod, newdata = nd, 
               exclude = c("s(stream)", "s(stream,rc_section)"),
               se.fit = TRUE, 
               newdata.guaranteed = TRUE, type = "link")
  pred <- nd
  pred$fit <- as.numeric(p$fit)
  pred$se.fit <- as.numeric(p$se.fit)
  pred$lowerCI <- mod$family$linkinv(pred$fit - 1.96 * p$se.fit)
  pred$upperCI <- mod$family$linkinv(pred$fit + 1.96 * p$se.fit)
  pred$fit     <- mod$family$linkinv(pred$fit)
  pred$col <- cols[pred$sex]
  
  #### Add error envelopes for males/females
  lapply(split(pred, pred$sex), function(d) {
    add_error_envelope(d$length, 
                       ci = list(fit = d$fit, lowerCI = d$lowerCI, upperCI = d$upperCI), 
                       add_fit = list(col = scales::alpha(cols[d$sex[1]], alpha_fit)), 
                       add_ci = list(col = scales::alpha(cols[d$sex[1]], alpha_ci), border = FALSE))
  })
  
}

#### Add observed proportions 
squash_param <- 27
add_proportions <- function(data, add_lines = FALSE, squash = squash_param) {
  if (add_lines) {
    # Add lines for the observed proportion of migrants with body size 
    lapply(split(data, data$sex), function(d) {
      # d <- split(data, data$sex)[["M"]]
      lines(d$length, d$pr, col = d$col[1], lty = 3)
    })
  }
  # Add observed proportions for each size class
  data$cex <- data$n/squash
  points(data$length, data$pr,  
         pch = 21, bg = data$col, col = data$col, 
         cex = data$cex)
  # Highlight small points in red
  tiny <- data[data$cex < 0.5, ]
  points(tiny$length, tiny$pr, col = "red", lwd = 0.5)
}
add_proportions(prop_ss)

#### Add observed outcomes (migration/no migration)
# ... We shift the points for males up/down slightly to facilitate visualisation
add_outcomes <- function(data) {
  mi <- which(data$sex == "M")
  fi <- which(data$sex == "F")
  data$migration_int <- as.integer(as.character(data$migration))
  data$migration_p <- data$migration_int 
  data$migration_p[data$sex == "M" & data$migration == "0"] <- -0.1
  data$migration_p[data$sex == "M" & data$migration == "1"] <- 1.1
  data$migration_p[data$sex == "F" & data$migration == "0"] <- -0.05
  data$migration_p[data$sex == "F" & data$migration == "1"] <- 1.05
  points(data$length[mi], data$migration_p[mi], 
         col = scales::alpha(cols["M"], 0.5), cex = 0.5)
  points(data$length[fi], data$migration_p[fi], 
         col = scales::alpha(cols["F"], 0.5), cex = 0.5)
}
add_outcomes(fish)

#### Add axis labels
add_axes_labels <- function(cex = 1.25, line = 2, ...) {
  mtext(side = 1, "Standard length (cm)", cex = cex, line = line, ...)
  mtext(side = 2, "Probability of out-migration", cex = cex, line = line, ...)
}
add_axes_labels()

#### Write legend(s)
# (It is easier to manipulate these manually after the fact)
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

dev.off()


#########################
#### Visualise observations/predictions by stream

# This code is currently only implemented for the GAM
if (!is_glmer) {
  
  png(here_fig("migration-prob-by-stream.png"), 
      height = 5, width = 10, units = "in", res = 600)
  pp <- par(mfrow = c(2, 4), oma = c(3, 3, 1, 1), mar = c(2, 2, 2, 2))
  lapply(seq_len(length(unique(fish$stream))), function(i) {

    #### Define variables
    response  <- "pr"
    predictor <- "length"
    mframe <- model.frame(mod)
    stream <- sort(unique(fish$stream))[i]
    # Define fish_for_stream (rather than model.frame(mod)) because
    # ... `fish_for_stream` contains 'length' 
    # ... whereas model.frame(mod) contains log(length)
    fish_for_stream <- fish[fish$stream == stream, ]
    # Check size ranges
    # Size ranges are shown via the rug (add_outcomes())
    # But points mark proportions in size classes
    range(fish_for_stream$length[fish_for_stream$sex == "M"])
    range(fish_for_stream$length[fish_for_stream$sex == "F"])

    #### Create plot
    pretty_blank(prop_ss, predictor, response, pretty_axis_args = paa)
    pred <- gen_pred(mod, stream, predictor, mframe = fish_for_stream, 
                     exclude = "s(stream,rc_section)", newdata.guaranteed = TRUE)
    add_error_envelopes_by_sex(pred, predictor)
    add_proportions(prop_sss[prop_sss$stream == stream, ], squash = 5)
    add_outcomes(fish_for_stream)
    mtext(side = 3, stream, font = 2)
  }) |> invisible()
  
  #### Add axes
  add_axes_labels(outer = TRUE, line = 1)
  par(pp)
  dev.off()
  
}


#########################
#########################
#### Compare predictions for selected sizes

#### Compare predicted Pr migration for small/large males/females

if (is_glmer) {
  
  ggpredict(mod, terms = "sex", condition = c(length = small))
  ggpredict(mod, terms = "sex", condition = c(length = large))
  
} else {
  
  wrap_compare_gam <- function(model, newdata) {
    newdata$stream     <- fish$stream[1]
    newdata$rc_section <- fish$rc_section[1]
    compare_gam(model, newdata, 
                exclude = c("s(stream)", "s(rc_section,stream)"), 
                newdata.guaranteed = TRUE)
  }
  
  # Predicted migration pr for small/large males/females
  rbind(
    # Small females/males 
    wrap_compare_gam(mod, data.frame(sex = "F", length = small, yday = median(fish$yday))),
    wrap_compare_gam(mod, data.frame(sex = "M", length = small, yday = median(fish$yday))),
    # Large females/males
    wrap_compare_gam(mod, data.frame(sex = "F", length = large, yday = median(fish$yday))),
    wrap_compare_gam(mod, data.frame(sex = "M", length = large, yday = median(fish$yday)))
  ) |> round(digits = 2)
    
  # How much more likely are small females likely to migrate than smaller males? 
  # ... Difference in fitted Pr is 0.30 (females) versus 0.15 (males)
  # ... Difference in CIs 
  
}


#########################
#########################
#### Model diagnostics 

#### Simulate residuals
# ... By default, DHARMa re-simulated all levels in the model hierarchy 
# ... This tests the model structure as a whole (see also ?simulate.merMod)
# ... but we can check the diagnostics are stable with:
# ... * re.form = NA (default) -- condition on no random effects
# ... ... i.e., simulate new values for the random effects
# ... * re.form = NULL -- condition on all random effects
# ... ... i.e., random effects are set on fitted values
mframe <- model.frame(mod)
res    <- simulateResiduals(mod, refit = FALSE, plot = TRUE, re.form = NULL)

#### Check residuals versus predictors
plotResiduals(res, form = mframe$sex)
plotResiduals(res, form = mframe[, "log(length)"])

#### Run additional DHARMa checks
testDispersion(res)
testZeroInflation(res)
str(mod)

#### Check the distribution of the random effects
if (is_glmer) {
  
  re <- ranef(mod)[[1]][, 1]
  qqnorm(re)
  qqline(re)
  shapiro.test(re)

} else {
  
  pp <- par(mfrow = c(2, 2))
  gam.check(mod, rep = 500)
  par(pp)
  k.check(mod)

}


#### End of code. 
#########################
#########################