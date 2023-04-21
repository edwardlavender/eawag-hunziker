###########################
###########################
#### define_helpers.R

#### Aims
# 1) Define helper functions

#### Prerequisites
# 1) NA


###########################
###########################
#### Define functions

#' @title Parse bins from `cut()`
#' @description This function parses bins from `cut()`. Bins are defined as the midpoint between each pair of values in the bin. 

parse_cut <- function(x) {
  x <- as.character(x)
  x <- substr(x, 2, nchar(x) - 1)
  xs <- stringr::str_split_fixed(x, ",", 2)
  xs <- apply(xs, 2, as.numeric)
  apply(xs, 1, mean)
}


#### End of code.
###########################
###########################