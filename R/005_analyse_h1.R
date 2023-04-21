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
library(lme4)
library(DHARMa)
library(ggeffects)
library(prettyGraphics)

#### Load data
source(here_r("001_define_global_param.R"))
source(here_r("002_define_helpers.R"))
fish <- readRDS(here_data("fish.rds"))


#########################
#########################
#### Summarise data

#### Count the number of streams/captures per stream
length(unique(fish$stream))
table(fish$stream)

#### Calculate the observed proportion of migrant males/females by size class
prs <- 
  fish |>
  group_by(sex) |>
  mutate(bin = cut(length, seq(min(length), max(length), length = 10))) |>
  ungroup() |>
  group_by(sex, bin) |>
  summarise(pr = length(which(migration == "1"))/n(), 
            n = n()) |>
  mutate(length = parse_cut(bin), 
         col = as.character(scales::alpha(cols, alpha_pt)[as.character(sex)])) |> 
  select(sex, n, bin, length, pr, col) |>
  filter(!is.na(bin))


#########################
#########################
#### Implement modelling 

#### Model migration (0, 1) as a function of size and sex 
# Controls
# * We will re-scale length to improve model identifiability 
# * NAGQ controls the accuracy of the log likelihood evaluation
# * The default option is 1 but we can increase this (at the expense of speed) 
mod <- glmer(migration ~ log10(length) * sex + (1|stream), 
             data = fish, family = binomial(link = "logit"), 
             nAGQ = 25)
# Check model summary 
summary(mod)

#### Generate model predictions (via predict())
# We can generate predictions for Pr migration via ?predict.merMod
# With this function, you can't compute standard errors though
# see ?predict.merMod and ?bootMer
# Define dataframe for predictions 
n <- 100
ms <- fish$length[fish$sex == "M"]
fs <- fish$length[fish$sex == "F"]
ms <- seq(min(ms), max(ms), length = n)
fs <- seq(min(fs), max(fs), length = n)
nd <- data.frame(sex = c(rep("M", n), rep("F", n)), 
                 length = c(ms, fs))
nd <- lapply(unique(fish$stream), function(stream) {
  nd$stream <- stream
  nd
}) |> dplyr::bind_rows()
nd$sex <- factor(nd$sex)
nd$stream <- factor(nd$stream)
head(nd)
# Make predictions with explicit specification of random effect
predict(mod, newdata = nd, re.form = ~(1|stream))

#### Generate model predictions (via ggeffects::ggpredict())
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


#########################
#########################
#### Visualise observations and predictions

#### Create blank plot 
pretty_plot(prs$length, prs$pr,
            type = "n", 
            xlab = "", ylab = "")

#### Add model predictions with SEs
add_error_envelope(pm$x, ci = list(fit = pm$predicted, lowerCI = pm$conf.low, upperCI = pm$conf.high), 
                   add_fit = list(col = scales::alpha(cols["M"], alpha_fit)), 
                   add_ci = list(col = scales::alpha(cols["M"], alpha_ci), border = FALSE))
add_error_envelope(pf$x, ci = list(fit = pf$predicted, lowerCI = pf$conf.low, upperCI = pf$conf.high), 
                   add_fit = list(col = scales::alpha(cols["F"], alpha_fit)), 
                   add_ci = list(col = scales::alpha(cols["F"], alpha_ci), border = FALSE))

#### Add observed proportions 
# Add lines for the observed proportion of migrants with body size 
lapply(split(prs, prs$sex), function(d) {
  # d <- split(prs, prs$sex)[["M"]]
  lines(d$length, d$pr, col = d$col[1], lty = 3)
})
# Add observed proportions for each size class
points(prs$length, prs$pr,  
       pch = 21, bg = prs$col, col = prs$col, 
       cex = prs$n/20)

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
#### Model diagnostics 

#### Simulate residuals
# ... By default, DHARMa re-simulated all levels in the model hierarchy 
# ... This tests the model structure as a whole (see also ?simulate.merMod)
res <- simulateResiduals(mod, refit = FALSE, plot = TRUE)

#### Check residuals versus predictors
plotResiduals(res, form = fish$sex)
plotResiduals(res, form = fish$length)

#### Run additional DHARMa checks
testDispersion(res)
testZeroInflation(res)
str(mod)

#### Check the distribution of the random effects
re <- ranef(mod)[[1]][, 1]
qqnorm(re)
qqline(re)
shapiro.test(re)


#### End of code. 
#########################
#########################