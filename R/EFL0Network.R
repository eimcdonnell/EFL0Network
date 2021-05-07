#' EFL0Network: A package for estimating an Elastic Fuse L0 (EFL0) dynamic
#' Gaussian graphical model
#'
#' The EFL0Network package provides three important functions:
#' EFL0_time_varying_graphical, EFL0_graphical and EFL0.
#' 
#' @docType package
#' @name EFL0Network
#' @importFrom dplyr select filter mutate arrange distinct full_join left_join
#' pull rename bind_rows n
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer pivot_wider nest unnest tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom stats lm loess resid weighted.mean coef density dist time
#' @importFrom utils head
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix Matrix
#' @useDynLib EFL0Network
#' @importFrom Rcpp sourceCpp
NULL
#> NULL
