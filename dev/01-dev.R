###########################
###########################
#### 02-dev.R

#### Aims
# 1) Record project development

#### Prerequisites
# 1) NA


###########################
###########################
#### Implement dev

#### Use git
# usethis::use_git()
# usethis::use_github()
usethis::git_vaccinate()

#### Use dependency management
renv::init()

#### Install package(s)
renv::install("edwardlavender/dv")
# commonmark/(r)markdown packages (for README documentation)
if (!requireNamespace("commonmark", quietly = TRUE))
  renv::install("commonmark")
if (!requireNamespace("markdown", quietly = TRUE))
  renv::install("markdown")
if (!requireNamespace("rmarkdown", quietly = TRUE))
  renv::install("rmarkdown")

#### Use dv templates

# Set up template project structure
dv::use_template_proj()

# Update .gitignore
dv::use_template_gitignore()

# Add a README and associated files
# usethis::use_code_of_conduct("insert_email_here")
dv::use_template_readme(title = "Migratory propensity in juvenile brown trout depends on the interaction between sex and body size",
                        author = "Edward Lavender",
                        email = "edward.lavender@eawag.ch")

# Add template scripts
if (!requireNamespace("pacman", quietly = TRUE)) renv::install("pacman")
dv::use_template_script(here_r("insert_script_name_1.R"))
dv::use_template_script(here_r("insert_script_name_2.R"))

# Save the project directory 'tree'
# ... This enables the project directory tree to be rebuilt on another machine
# ... This function should be re-run when the directory tree is updated
dv::use_template_tree(save = dv::here_data("inst", "tree.rds"))

#### Update renv
## Take snapshot
renv::snapshot()
## Clean snapshot
# Note that this may attempt to drop 'suggested packages'
# ... that are required by (some) functions from other packages
# ... but which are not used directly. To guard against this,
# ... make an arbitrary call to the required packages
# ... where they are needed.
renv::clean()


#### End of code.
###########################
###########################
