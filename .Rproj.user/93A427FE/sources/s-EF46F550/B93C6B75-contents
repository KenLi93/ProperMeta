
#' @param yi: study effect size
#' @param vi: variance estimate
#' @param ni: study sample size
#'
test_pwa_meta <- function(yi, vi, ni) {
  ## test: if yi, vi and ni are all numetic vectors
  if ((!is.numeric(yi))|(!is.numeric(vi))|(!is.numeric(ni))) {
    stop("yi, vi, and ni should all be numeric vectors.")
  }

  ## test: if all entries of vi is positive
  if (any(vi <= 0)) {
    stop("Elements in vi have to be positive.")
  }

  ## test: if all entries of n is positive
  if (any(ni <= 0)) {
    stop("Elements in ni have to be positive.")
  }

  ## test: if x and n have the same length
  if ((length(yi) != length(vi)) | (length(yi)!= length(ni))) {
    stop("yi, vi, and ni should all have the same length.")
  }
}

#' Inference for precision weighted average estimator in general meta-analysis
#'
#' @param yi: study effect size
#' @param vi: variance estimate
#' @param ni: study sample size
#'
#' @export
#'
pwa_meta <- function(yi, vi, ni, alternative = "two.sided",
                     conf.level = 0.95) {
  test_pwa_meta(yi, vi, ni)
  K <- length(yi)

  nn <- sum(ni)

  if (K <= 5) {
    warning("the number of studies is small -- variance estimates can be inaccurate!")
  }

  ## constants for computing the overall variance
  ui <- vi * ni
  kappa <- cov(sqrt(ni) * yi, sqrt(ni) * ui)
  eta <- var(sqrt(ni) * ui)
  gamma <- ni / nn

  EST <- sum(yi / vi) / sum(1 / vi)

  VAR_NAIVE <- 1 / sum(1 / vi)
  SE_NAIVE <- sqrt(VAR_NAIVE)

  VAR <- sum(gamma * (1 / ui - 2 * yi * kappa / ui ^ 3 + yi ^ 2 * eta / ui ^ 4)) /
    sum(gamma / ui) ^ 2 -
    2 * sum(gamma * (- kappa / ui ^ 3 + yi * eta / ui ^ 4)) * sum(gamma * yi / ui) /
    sum(gamma / ui) ^ 3 +
    sum(gamma * yi / ui) ^ 2 * sum(gamma * eta / ui ^ 4) /
    sum(gamma / ui) ^ 4
  SE <- sqrt(VAR)

  ## construct confidence intervals
  CI_NAIVE <- wald_ci(EST, SE_NAIVE, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)


  return(list(EST = EST, SE_NAIVE = SE_NAIVE, VAR_NAIVE = VAR_NAIVE, CI_NAIVE = CI_NAIVE,
              SE = SE, VAR = VAR, CI = CI))
}
