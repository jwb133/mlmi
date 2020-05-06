
#' Simulated example data with continuous outcome measured repeatedly over time
#'
#' A dataset in the wide form containing simulated data with a repeatedly measured
#' outcome. Some outcome values are missing. The missing data pattern is
#' monotone. There are two baseline covariates.
#'
#' @format A data frame with 500 rows and 7 variables:
#' \describe{
#'   \item{id}{ID for individual}
#'   \item{trt}{A numeric 0/1 variable indicating control or active treatment group}
#'   \item{v}{A baseline covariate}
#'   \item{y0}{Baseline measurement of the outcome variable}
#'   \item{y1}{Outcome measurement at visit 1}
#'   \item{y2}{Outcome measurement at visit 2}
#'   \item{y3}{Outcome measurement at visit 3}
#' }
#'
"ctsTrialWide"
