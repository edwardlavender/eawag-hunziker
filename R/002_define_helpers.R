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

#' @title Define ranges as a string ('a–b')

str_range <- function(x) {
  if (inherits(x, "character")) {
    paste0(c(x[1], x[length(x)]), collapse = "–")
  } else {
    paste0(c(min(x), max(x)), collapse = "–")
  }
}


#' @title Parse bins from `cut()`
#' @description This function parses bins from `cut()`. Bins are defined as the midpoint between each pair of values in the bin. 

parse_cut <- function(x) {
  x <- as.character(x)
  x <- substr(x, 2, nchar(x) - 1)
  xs <- stringr::str_split_fixed(x, ",", 2)
  xs <- apply(xs, 2, as.numeric)
  apply(xs, 1, mean)
}

#' @title Compare gam
#' @description This function compares the predictions of a GAM for a hypothetical individual. The `s(stream)` effect is excluded. 
#' @details This function is designed for the model(s) in analyse_h1.R

compare_gam <- function(model, newdata) {
  newdata$sex <- factor(newdata$sex, levels = c("F", "M"))
  p <- 
    predict(model, 
            newdata = newdata,
            exclude = "s(stream)", 
            newdata.guaranteed = TRUE, se.fit = TRUE) |> 
    list_CIs(inv_link = mod$family$linkinv, plot_suggestions = FALSE)
  do.call(cbind, p)
}

#### End of code.
###########################
###########################