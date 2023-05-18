###########################
###########################
#### define_helpers.R

#### Aims
# 1) Define helper functions

#### Prerequisites
# 1) NA


###########################
###########################
#### Define utils

#' @title Define ranges as a string ('a–b')
#' 
str_range <- function(x) {
  if (inherits(x, "character")) {
    x1 <- x[1]
    x2 <- x[length(x)]
  } else {
    x1 <- min(x)
    x2 <- max(x)
  }
  if (x1 == x2) {
    rng <- x1
  } else {
    rng <- paste0(c(x1, x2), collapse = "–")
  }
  rng
}


#' @title Parse bins from `cut()`
#' @description This function parses bins from `cut()`. Bins are defined as the midpoint between each pair of values in the bin. 
#' 
parse_cut <- function(x) {
  x <- as.character(x)
  x <- substr(x, 2, nchar(x) - 1)
  xs <- stringr::str_split_fixed(x, ",", 2)
  xs <- apply(xs, 2, as.numeric)
  apply(xs, 1, mean)
}


###########################
###########################
#### Plotting helpers

#' @title Create a blank plot
#' 
pretty_blank <- function(mframe, predictor, response, ...) {
  pretty_plot(mframe[, predictor], mframe[, response],
              ...,
              type = "n", 
              xlab = "", ylab = "")
}

#' @title Generate model predictions 
#' @details This function assumes the following objects exist in the workspace:
#' * `fish` data.frame
#' * `cols`
#' 
gen_pred <- function(mod, stream = NULL, predictor, mframe = NULL, n = 100, 
                     exclude = NULL, newdata.guaranteed = FALSE) {
  # Define data for stream 
  if (is.null(mframe)) mframe <- model.frame(mod)
  if (!is.null(stream)) {
    mframe <- mframe[mframe$stream == stream, ]
  }
  # Define a sequence of values of the predictor, separately for M/F, for prediction
  ms <- mframe[mframe$sex == "M", predictor]
  fs <- mframe[mframe$sex == "F", predictor]
  ms <- seq(min(ms, na.rm = TRUE), max(ms, na.rm = TRUE), length = n)
  fs <- seq(min(fs, na.rm = TRUE), max(fs, na.rm = TRUE), length = n)
  # Generate predictions for stream
  nd <- data.frame(sex = factor(c(rep("F", n), rep("M", n))),
                   predictor = c(fs, ms), 
                   yday = median(fish$yday))
  if (!is.null(stream)) nd$stream <- stream
  colnames(nd)[colnames(nd) %in% "predictor"] <- predictor
  p <- predict(mod, newdata = nd, 
               se.fit = TRUE, 
               exclude = exclude, newdata.guaranteed = newdata.guaranteed,
               type = "link")
  # Return dataframe with CIs
  pred <- nd
  pred$fit <- as.numeric(p$fit)
  pred$se.fit <- as.numeric(p$se.fit)
  pred$lowerCI <- mod$family$linkinv(pred$fit - 1.96 * p$se.fit)
  pred$upperCI <- mod$family$linkinv(pred$fit + 1.96 * p$se.fit)
  pred$fit     <- mod$family$linkinv(pred$fit)
  pred$col <- cols[pred$sex]
  pred
}


#' @title Add error envelopes by sex
#' 
add_error_envelopes_by_sex <- function(pred, predictor) {
  lapply(split(pred, pred$sex), function(d) {
    add_error_envelope(d[, predictor], 
                       ci = list(fit = d$fit, lowerCI = d$lowerCI, upperCI = d$upperCI), 
                       add_fit = list(col = scales::alpha(cols[d$sex[1]], alpha_fit)), 
                       add_ci = list(col = scales::alpha(cols[d$sex[1]], alpha_ci), border = FALSE))
  })
}

#' @title Add observations
#' 
add_obs_by_sex <- function(mframe, predictor, response, ...) {
  points(mframe[, predictor], mframe[, response], 
         pch = 21,
         col = scales::alpha(cols[mframe$sex], alpha_pt), 
         bg = scales::alpha(cols[mframe$sex], alpha_pt), ...)
}


###########################
###########################
#### Statistics helpers

#' @title Compare gam
#' @description This function compares the predictions of a GAM for a hypothetical individual. 
#' @details This function is designed for the model(s) in analyse_h1.R
#' 
compare_gam <- function(model, newdata, 
                        exclude = NULL,
                        newdata.guaranteed = FALSE) {
  newdata$sex <- factor(newdata$sex, levels = c("F", "M"))
  p <- 
    predict(model, 
            newdata = newdata,
            exclude = c("s(stream)", "s(stream,section)"),
            newdata.guaranteed = newdata.guaranteed, 
            se.fit = TRUE) |> 
    list_CIs(inv_link = mod$family$linkinv, plot_suggestions = FALSE)
  do.call(cbind, p)
}


#### End of code.
###########################
###########################