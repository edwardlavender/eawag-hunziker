#########################
#########################
#### analyse_h1.R

#### Aims
# 1) Analyse H1 (migration Pr ~ sex * size)

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
# * NAGQ controls the accuracy of the log likelihood evaluation
# * The default option is 1 but we can increase this (at the expense of speed) 

mod_1 <- glmer(migration ~ log(length) * sex + yday + (1|stream), 
             data = fish, family = binomial(link = "logit"), 
             nAGQ = 25)

#### GAMM: model migration (0, 1) as a function of size and sex and other factors
# Motivation 
# * Examination of the model structure for mod_1 suggests 
# ... a non linear relationship for examples at small sizes 
# ... suggesting a GAM model would be more appropriate here. 
# Structure
# * We fit a model with group-level smoothness & different trends
# ... i.e., 'Model I' in Pedersen et al. (2018)
# * Length is scaled via log()
# ... Scaling length is not necessary 
# ... and makes little visible difference to predictions
# ... but marginally improves model fit (AIC)
# ... scale() produces poor fit on predictive plot 
# * REML or ML recommended for fitting GAMs
mod_2 <- gam(migration ~ sex + s(log(length), by = sex) + s(stream, bs = "re"), 
             family = binomial, data = fish, gamma = 1.4, 
             method = "REML")

#### Compare model AICs
# For the glmer model
AIC(mod_1)
# For the GAM model
AIC(mod_2)
# * log(length): 769.2234 ("ML"), 767.7716 ("REML"), little visual difference 
# * length:      773.2707 ("ML"), 769.8456 ("REML"), little visual difference 
# * scale(length): poorer fit visually apparent
# Choose model
mod <- mod_2
is_glmer <- inherits(mod, "merMod")

##### Check model summary 
summary(mod)
AIC(mod)
if (is_glmer) {
  equatiomatic::extract_eq(mod)
  MuMIn::r.squaredGLMM(mod)
}


#########################
#########################
#### Visualise observations and predictions

#### Create blank plot 
pretty_plot(prop_ss$length, prop_ss$pr,
            type = "n", 
            xlab = "", ylab = "")

#### Add model predictions with SEs
n <- 100
ms <- fish$length[fish$sex == "M"]
fs <- fish$length[fish$sex == "F"]
ms <- seq(min(ms), max(ms), length = n)
fs <- seq(min(fs), max(fs), length = n)

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
  nd <- data.frame(sex = factor(c(rep("F", n), rep("M", n))),
                   length = c(fs, ms))
  p <- predict(mod, newdata = nd, 
               exclude = "s(stream)", 
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
# Add lines for the observed proportion of migrants with body size 
lapply(split(prop_ss, prop_ss$sex), function(d) {
  # d <- split(prop_ss, prop_ss$sex)[["M"]]
  lines(d$length, d$pr, col = d$col[1], lty = 3)
})
# Add observed proportions for each size class
points(prop_ss$length, prop_ss$pr,  
       pch = 21, bg = prop_ss$col, col = prop_ss$col, 
       cex = prop_ss$n/20)

#### Add observed outcomes (migration/no migration)
# ... We shift the points for males up/down slightly to facilitate visualisation
mi <- which(fish$sex == "M")
fi <- which(fish$sex == "F")
fish$migration_int <- as.integer(as.character(fish$migration))
fish$migration_p <- fish$migration_int 
fish$migration_p[fish$sex == "M" & fish$migration == "0"] <- 0.05
fish$migration_p[fish$sex == "M" & fish$migration == "1"] <- 0.95
points(fish$length[mi], fish$migration_p[mi], 
       col = scales::alpha(cols["M"], 0.5), cex = 0.5)
points(fish$length[fi], fish$migration_p[fi], 
       col = scales::alpha(cols["F"], 0.5), cex = 0.5)

#### Add axis labels
mtext(side = 1, "Total length (cm)", line = 2)
mtext(side = 2, "Probability of out-migration", line = 2)


#########################
#########################
#### Compare predictions for selected sizes

#### Compare predicted Pr migration for small/large males/females
small <- 100
large <- 150

if (is_glmer) {
  
  ggpredict(mod, terms = "sex", condition = c(length = small))
  ggpredict(mod, terms = "sex", condition = c(length = large))
  
} else {
  
  rbind(
    compare_gam(mod, data.frame(sex = "F", length = small)),
    compare_gam(mod, data.frame(sex = "F", length = large)),
    compare_gam(mod, data.frame(sex = "M", length = small)),
    compare_gam(mod, data.frame(sex = "M", length = large))
  )
    
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
res <- simulateResiduals(mod, refit = FALSE, plot = TRUE, re.form = NULL)

#### Check residuals versus predictors
plotResiduals(res, form = fish$sex)
plotResiduals(res, form = fish$length)

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