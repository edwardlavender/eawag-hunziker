#########################
#########################
#### define_global_param.R

#### Aims
# 1) Define global parameters

#### Prerequisites
# 1) NA


#########################
#########################
#### Parameter definitions

#### Define colours for males/females
# Define colours
cols <- setNames(c("royalblue", "black"), c("F", "M"))
# Define colour transparency for points/fitted lines/CIs
alpha_pt  <- 0.75
alpha_fit <- 0.50
alpha_ci  <- 0.25

#### Define length limits for plots
# range(fish$length)
lim_length <- c(8.5, 25)
at_length  <- seq(10, 24, by = 4)
# range(fish$mass, na.rm = TRUE)
lim_mass <- c(1, 30)
at_mass  <- seq(5, 30, by = 5)
# Collate information
lim_size <- list(length = lim_length, 
                 mass = lim_mass)
at_size <- list(length = at_length, 
                mass = at_mass)

#### Define 'small' and 'large' sizes for comparisons
# The units here need to be the same is the data used to fit models
small <- 10 # cm or kg
large <- 15 # cm or kg


#### End of code.
#########################
#########################