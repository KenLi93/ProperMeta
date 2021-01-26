test_binomial <- function(xi, ni) {
  ## test: if x and n are both numetic vectors
  if ((!is.numeric(xi))|(!is.numeric(ni))) {
    stop("xi and ni should both be numeric vectors.")
  }

  ## test: if all entries of x is nonnegative
  if (any(xi < 0)) {
    stop("Elements in xi have to be nonnegative.")
  }

  ## test: if all entries of n is positive
  if (any(ni <= 0)) {
    stop("Elements in ni have to be positive.")
  }

  ## test: if x and n have the same length
  if (length(xi) != length(ni)) {
    stop("xi and ni should have the same length.")
  }

  ## test: if each entry of x and n are integers
  if (any(round(xi) != xi)) {
    warning("Elements in xi are not all integers -- take the closest integers.")
  }

  ## test: if each entry of x is no larger than the corresponding entry in n
  if (any(xi > ni)) {
    stop("Elements of xi should be no larger than the corresponding elements of ni.")
  }
}

#' Inference for overall proportion in proportion meta-analysis
#'
#' @param xi: number of responders in each study
#' @param ni: study sample size
#'
#' @export
#'

overall_proportion <- function(xi, ni, alternative = "two.sided", conf.level = 0.95, ...) {

  test_binomial(xi, ni)
  xi <- round(xi)
  ni <- round(ni)


  ## quantities used in calculation
  pp <- xi / ni
  nn <- sum(ni)
  gamma <- ni / nn

  ## estimator: overall proportion
  EST <- sum(xi) / sum(ni)

  VAR_HOM <- EST * (1 - EST) / (nn - 1)
  SE_HOM <- sqrt(VAR_HOM)

  VAR <- sum(gamma ^ 2 * pp * (1 - pp) / (ni - 1))
  SE <- sqrt(VAR)

  ## construct confidence intervals
  CI_HOM <- wald_ci(EST, SE_HOM, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)

  return(list(EST = EST, SE_HOM = SE_HOM, VAR_HOM = VAR_HOM, CI_HOM = CI_HOM,
              SE = SE, VAR = VAR, CI = CI))

}


#' Inference for precision weighted average proportion in proportion meta-analysis
#'
#' @param xi: number of responders in each study
#' @param ni: study sample size
#'
#' @export
#'

pwa_proportion <- function(xi, ni, alternative = "two.sided", conf.level = 0.95,
                           correction_factor = 0.01, ...) {
  test_binomial(xi, ni)
  xi <- round(xi)
  ni <- round(ni)

  ## zero cell correction: adding a small number to the zero counts,
  ## subtracting a small number to the full counts
  xi[xi == 0] <- xi[xi == 0] + correction_factor
  xi[xi == n] <- xi[xi == n] - correction_factor

  ## quantities used in calculation
  pp <- xi / ni
  nn <- sum(ni)
  gamma <- ni / nn

  ## estimator: overall proportion
  EST <- sum(xi) / sum(ni)

  VAR_NAIVE <- EST * (1 - EST) / (nn - 1)
  SE_NAIVE <- sqrt(VAR_NAIVE)

  VAR <- sum(gamma ^ 2 * pp * (1 - pp) / (n - 1))
  SE <- sqrt(VAR)


  ## construct confidence intervals
  CI_NAIVE <- wald_ci(EST, SE_NAIVE, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)



  return(list(EST = EST, SE_NAIVE = SE_NAIVE, VAR_NAIVE = VAR_NAIVE, CI_NAIVE = CI_NAIVE,
              SE = SE, VAR = VAR, CI = CI))
}
