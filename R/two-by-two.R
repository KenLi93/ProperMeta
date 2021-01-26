test_2by2_tables <- function(ai, n1i, ci, n2i) {
  ## test: if x and n are both numetic vectors
  if ((!is.numeric(ai))|(!is.numeric(n1i))|(!is.numeric(ci))|(!is.numeric(n2i))) {
    stop("The arguments ai, n1i, ci, n2i should be numeric vectors.")
  }

  ## test: if all entries of x is nonnegative
  if (any(c(ai, n1i, ci, n2i) < 0)) {
    stop("Elements in X have to be nonnegative.")
  }


  ## test: if ai, n1i, ci, n2i have the same length
  if ((length(ai) != length(n1i)) | length(ai) != length(ci) |
      length(ai) != length(n2i)) {
    stop("The vectors ai, n1i, ci, n2i should have the same length.")
  }

}

#' Inference for MLE of common log odds ratios in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'

## MLE of the common log odds ratio
mle_logor <- function(ai, n1i, ci, n2i, correction_factor = 0.01,
                      alternative = "two.sided",
                      conf.level = 0.95) {
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  ni <- sum(ni)

  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

  ## zero-cell correction for tables with zeroes
  for (jj in 1:K) {
    if (any0[jj]) {
      ai[jj] <- ai[jj] + correction_factor * n1i[jj] / ni[jj]
      bi[jj] <- bi[jj] + correction_factor * n1i[jj] / ni[jj]
      ci[jj] <- ci[jj] + correction_factor * n2i[jj] / ni[jj]
      di[jj] <- di[jj] + correction_factor * n2i[jj] / ni[jj]
    }
  }




  ## formulate the data for running logistic regression
  dat_reshape <- rbind(data.frame(rr = ai, nr = bi, xx = 1, t = 1:K),
                       data.frame(rr = ci, nr = di, xx = 0, t = 1:K))
  suppressWarnings(
    logit_mod <- glm(cbind(rr, nr) ~ 0 + as.factor(t) + xx, data = dat_reshape, family = "binomial")
  )
  ## MLE <- logit_mod$coefficients["xx"]

  MLE <- logit_mod$coefficients

  ## compute the asymptotic variance-covariance matrix of beta_hat
  aa <- MLE[1:K]
  bb <- MLE[K + 1]

  ## make sandwich estimator of variance
  ## make matrix J.hat as in the document (negtive Fisher Information)
  uu <- - n1i * exp(aa + bb) / ((1 + exp(aa + bb)) ^ 2 * nn)
  vv <- - n2i * exp(aa) / ((1 + exp(aa)) ^ 2 * nn)

  J <- matrix(0, nrow = K + 1, ncol = K + 1)
  diag(J) <- c(uu + vv, sum(uu))
  J[K + 1, 1:K] <- J[1:K, K + 1] <- uu

  ## make matrix U.hat as in the document
  ss <- ai * bi / (n1i * nn)
  tt <- ci * di / (n2i * nn)

  U <- matrix(0, nrow = K + 1, ncol = K + 1)
  diag(U) <- c(ss + tt, sum(ss))
  U[K + 1, 1:K] <- U[1:K, K + 1] <- ss

  ## estimated variance of MLE under heterogeneity: sandwich estimator
  Omega_het <- solve(-J) %*% U %*% solve(-J)

  ## estimated variance of MLE under homogeneity: inverse Fisher Information
  Omega_hom <- solve(-J)

  ## approximate variance of MLE of "common odds ratio" under homogeneity and heterogeneity
  VAR_HOM <- Omega_hom[K + 1, K + 1] / nn
  SE_HOM <- sqrt(VAR_HOM)

  VAR <- Omega_het[K + 1, K + 1] / nn
  SE <- sqrt(VAR)

  ## approximate CI of MLE
  ## construct confidence intervals

  CI_HOM <- wald_ci(EST, SE_HOM, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)

  return(list(EST = bb, SE_HOM = SE_HOM,VAR_HOM = VAR_HOM, CI_HOM = CI_HOM,
              SE = SE, VAR = VAR, CI = CI))
}

#' Inference for Firth mean-bias corrected MLE of common log odds ratios in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'

## Firth's mean-bias corrected logistic regression
firth_logor <- function(ai, n1i, ci, n2i, correction_factor = 0.01,
                        alternative = "two.sided",
                        conf_level = 0.95) {
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  nn <- sum(ni)


  ## objective function of firth logistic regression -- penalized negative log likelihood
  objfunc <- function(param) {
    alpha <- param[1:K]
    psi <- param[K+1]
    negloglikhd <- -sum(ai * (alpha + psi) - n1i * log(1 + exp(alpha + psi)) +
                          ci * alpha - n2i * log(1 + exp(alpha)))

    ## Fisher information
    FI <- matrix(0, nrow = K + 1, ncol = K + 1)

    ## quantities for calculating the Fisher information
    tt1 <- n1i * exp(alpha + psi) / (1 + exp(alpha + psi)) ^ 2
    tt0 <- n2i * exp(alpha) / (1 + exp(alpha + psi)) ^ 2

    FI[1:K, K + 1] <- FI[K + 1, 1:K] <- tt1
    diag(FI)[1:K] <- tt1 + tt0
    FI[K + 1, K + 1] <- sum(tt1)

    negloglikhd - 0.5 * log(det(FI))
  }

  ## use MLE as the initial value for firth's regression



  dat_reshape <- rbind(data.frame(rr = ai, nr = bi, xx = 1, t = 1:K),
                       data.frame(rr = ci, nr = di, xx = 0, t = 1:K))
  suppressWarnings(
    logit_mod <- glm(cbind(rr, nr) ~ 0 + as.factor(t) + xx, data = dat_reshape, family = "binomial")
  )
  mle <- logit_mod$coefficients

  ## point estimates of Firth regression
  firth_est <- as.numeric(optim(mle, objfunc)$par)

  EST <- firth_est[K + 1]

  ## Compute the sandwich covariance matrix for the estimator -- the same as MLE
  ## zero-cell correction for tables with zeroes -- necessary for variance calculation
  aa <- ai; bb <- bi; cc <- ci; dd <- di
  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))
  if (any(any0)) {
    for (jj in 1:K) {
      if (any0[jj]) {
        aa[jj] <- ai[jj] + correction_factor * n1i[jj] / N[jj]
        bb[jj] <- bi[jj] + correction_factor * n1i[jj] / N[jj]
        cc[jj] <- ci[jj] + correction_factor * n2i[jj] / N[jj]
        dd[jj] <- di[jj] + correction_factor * n2i[jj] / N[jj]
      }
    }
    n1i <- aa + bb
    n2i <- cc + dd
    ni <- n1i + n2i
    nn <- sum(nn)
  }

  ## constants that help with the calculation
  gamma <- ni / nn
  delta <- n1i / ni

  p1 <- ai / n1i
  p0 <- ci / n2i

  alpha <- firth_est[1:K]
  psi <- firth_est[K + 1]

  uu <- - gamma * delta * exp(alpha + psi) / (1 + exp(alpha + psi)) ^ 2
  vv <- - gamma * (1 - delta) * exp(alpha) / (1 + exp(alpha)) ^ 2
  ss <- gamma * delta * p1 * (1 - p1)
  tt <- gamma * (1 - delta) * p0 * (1 - p0)

  ## sandwich estimator
  JJ <- UU <- matrix(0, nrow = K + 1, ncol = K + 1)

  JJ[K + 1, 1:K] <- JJ[1:K, K + 1] <- uu
  diag(JJ)[1:K] <- uu + vv
  JJ[K + 1, K + 1] <- sum(uu)

  UU[K + 1, 1:K] <- UU[1:K, K + 1] <- ss
  diag(UU)[1:K] <- ss + tt
  UU[K + 1, K + 1] <- sum(ss)

  ## Fisher Information matrix
  vcov_hom <- solve(-JJ)

  ## sandwich covariance matrix
  vcov_hw <- solve(JJ) %*% UU %*% solve(JJ)

  ## variance and se of the common log or under homo- or heterogeneity
  VAR_HOM <- vcov_hom[K + 1, K + 1] / nn
  VAR <- vcov_hw[K + 1, K + 1] / nn

  SE_HOM <- sqrt(VAR_HOM)
  SE <- sqrt(VAR)



  ## compute the approximate ci
  ## construct confidence intervals

  CI_HOM<- wald_ci(EST, SE_HOM, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)

  list(EST = EST, SE_HOM = SE_HOM, VAR_HOM = VAR_HOM, CI_HOM = CI_HOM,
       SE = SE, VAR = VAR, CI = CI)
}

#' Inference for logarithm of Cochran-Mantel-Haenszel estimator in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'


## Confidence interval of logarithm of Cochran-Mantel-Haenszel statistics with or without assuming homogeneity
cmh_logor <- function(ai, n1i, ci, n2i, alternative = "two.sided", conf.level = 0.95,
                      correction_factor = 0.01, ...) {
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  nn <- sum(ni)


  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

  ## zero-cell correction for tables with zeroes
  if (any(any0)) {
    for (jj in 1:K) {
      if (any0[jj]) {
        ai[jj] <- ai[jj] + correction_factor * n1i[jj] / ni[jj]
        bi[jj] <- bi[jj] + correction_factor * n1i[jj] / ni[jj]
        ci[jj] <- ci[jj] + correction_factor * n2i[jj] / ni[jj]
        di[jj] <- di[jj] + correction_factor * n2i[jj] / ni[jj]
      }
    }
    ni <- n1i + n2i
    nn <- sum(ni)
  }
  ## mantel-haenszel statistics
  CMH <- sum(ai * di / ni) / sum(bi * ci / ni)
  EST <- log(CMH)


  VAR_HOM <- sum((bi * ci / ni) ^ 2 * (1 / ai + 1 / bi + 1 / ci + 1 / di)) /
    sum(bi * ci / ni) ^ 2

  SE_HOM <- sqrt(VAR_HOM)

  VAR <- sum((ai * bi / n1i * (di + ci * CMH) ^ 2 + ci * di / n2i * (ai + bi * CMH) ^ 2) / ni ^ 2) /
    (sum(bi * ci / ni) * CMH) ^ 2

  SE <- sqrt(VAR)

  ## construct confidence intervals
  ## construct confidence intervals
  CI_HOM<- wald_ci(EST, SE_HOM, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)


  return(list(EST = EST, SE_HOM = SE_HOM, VAR_HOM = VAR_HOM, CI_HOM = CI_HOM,
              SE = SE, VAR = VAR, CI = CI))
}


#' Inference for precision weighted average log odds ratios in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'

pwa_logor <- function(ai, n1i, ci, n2i, alternative = "two.sided", conf.level = 0.95,
                correction_factor = 0.01, ...){
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  nn <- sum(ni)


  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

  ## zero-cell correction for tables with zeroes
  if (any(any0)){
    for (jj in 1:K) {
      if (any0[jj]) {
        ai[jj] <- ai[jj] + correction_factor * n1i[jj] / ni[jj]
        bi[jj] <- bi[jj] + correction_factor * n1i[jj] / ni[jj]
        ci[jj] <- ci[jj] + correction_factor * n2i[jj] / ni[jj]
        di[jj] <- di[jj] + correction_factor * n2i[jj] / ni[jj]
      }
    }
    ## column totals
    n1i <- ai + bi
    n2i <- ci + di

    ## table totals
    ni <- n1i + n2i

    nn <- sum(ni)
  }

  gamma <- ni / nn
  delta <- n1i / ni

  p1 <- ai / n1i
  p0 <- ci / n2i

  ## study effect estimate: log odds ratio
  theta <- log((ai * di) / (bi * ci))

  ## scaled variance estimate
  vv <- 1 / (gamma * delta * p1 * (1 - p1)) + 1 / (gamma * (1 - delta) * p0 * (1 - p0))

  ## point estimate: inverse-variance weighted average log odds ratio
  EST <- sum(theta / vv) / sum(1 / vv)

  ## constants for calculating the sandwich covariance matrix
  ww <- - (1 - 2 * p1) / (gamma * delta * p1 * (1 - p1)) ^ 2 +
    (1 - 2 * p0) / (gamma * (1 - delta) * p0 * (1 - p0)) ^ 2

  zz <- (1 - 2 * p1) ^ 2 / (gamma * delta * p1 * (1 - p1)) ^ 3 +
    (1 - 2 * p0) ^ 2 / (gamma * (1 - delta) * p0 * (1 - p0)) ^ 3

  ## common formula for estimated variance of inverse-variance weighted estimators
  ## allowing heterogeneity
  VAR <- (sum(1 / vv - 2 * theta * ww / vv ^ 3 + theta ^ 2 * zz / vv ^ 4) / sum(1 / vv) ^ 2 -
    2 * sum(- ww / vv ^ 3 + theta * zz / vv ^ 4) * sum(theta / vv) / sum(1 / vv) ^ 3 +
    sum(theta / vv) ^ 2 * sum(zz / vv ^ 4) / sum(1 / vv) ^ 4) / nn

  SE <- sqrt(VAR)

  ## common formula for estimated variance of inverse-variance weighted estimators
  ## assuming homogeneity
  VAR_NAIVE <- 1 / (NN * sum(1 / vv))

  SE_NAIVE <- sqrt(VAR_HOM)


  ## construct confidence intervals
  CI_NAIVE <- wald_ci(EST, SE_NAIVE, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)


  return(list(EST = EST, SE_NAIVE = SE_NAIVE, VAR_NAIVE = VAR_NAIVE, CI_NAIVE = CI_NAIVE,
              SE = SE, VAR = VAR, CI = CI))

}

#' Inference for precision weighted average risk difference in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'

## inverse-variance weighted average risk difference, or Woolf's estimator
pwa_rd <- function(X, alternative = "two.sided", conf.level = 0.95,
                      correction_factor = 0.01, ...){
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  nn <- sum(ni)


  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

  ## zero-cell correction for tables with zeroes
  if (any(any0)){
    for (jj in 1:K) {
      if (any0[jj]) {
        ai[jj] <- ai[jj] + correction_factor * n1i[jj] / ni[jj]
        bi[jj] <- bi[jj] + correction_factor * n1i[jj] / ni[jj]
        ci[jj] <- ci[jj] + correction_factor * n2i[jj] / ni[jj]
        di[jj] <- di[jj] + correction_factor * n2i[jj] / ni[jj]
      }
    }
    ## column totals
    n1i <- ai + bi
    n2i <- ci + di

    ## table totals
    ni <- n1i + n2i

    nn <- sum(ni)
  }

  gamma <- ni / nn
  delta <- n1i / ni

  p1 <- ai / n1i
  p0 <- ci / n2i

  ## study effect estimate: risk difference
  theta <- p1 - p0

  ## unbiased variance estimate
  vi <- p1 * (1 - p1) / (m1 - 1) + p0 * (1 - p0) / (m0 - 1)

  ## point estimate: inverse-variance weighted average risk difference
  EST <- sum(theta / vi) / sum(1 / vi)

  ## variance and se estimate assuming homogeneity
  VAR_NAIVE <- 1 / sum(1 / vi)
  SE_NAIVE <- sqrt(VAR_NAIVE)

  ## constants for calculating the sandwich covariance matrix

  vv <- p1 * (1 - p1) / (gamma * delta) + p0 * (1 - p0) / (gamma * (1 - delta))

  ww <- p1 * (1 - p1) * (1 - 2 * p1) / (gamma * delta) ^ 2 -
    p0 * (1 - p0) * (1 - 2 * p0) / (gamma * (1 - delta)) ^ 2

  zz <- p1 * (1 - p1) * (1 - 2 * p1) ^ 2 / (gamma * delta) ^ 3 +
    p0 * (1 - p0) * (1 - 2 * p0) ^ 2 / (gamma * (1 - delta)) ^ 3

  ## common formula for estimated variance of inverse-variance weighted estimators
  ## allowing heterogeneity
  VAR <- (sum(1 / vv - 2 * theta * ww / vv ^ 3 + theta ^ 2 * zz / vv ^ 4) / sum(1 / vv) ^ 2 -
                2 * sum(- ww / vv ^ 3 + theta * zz / vv ^ 4) * sum(theta / vv) / sum(1 / vv) ^ 3 +
                sum(theta / vv) ^ 2 * sum(zz / vv ^ 4) / sum(1 / vv) ^ 4) / nn

  SE <- sqrt(VAR)


  ## construct confidence intervals
  CI_NAIVE <- wald_ci(EST, SE_NAIVE, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)


  return(list(EST = EST, SE_NAIVE = SE_NAIVE, VAR_NAIVE = VAR_NAIVE, CI_NAIVE = CI_NAIVE,
              SE = SE, VAR = VAR, CI = CI))

}

#' Inference for precision weighted average log risk ratio in meta-analysis of two-by-two tables
#'
#' @param ai: number of responders in group 1
#' @param n1i: sample size of group 1
#' @param ci: number of responders in group 2
#' @param n2i: sample size of group 2
#'
#' @export
#'
#'

pwa_lrr <- function(ai, n1i, ci, n2i, alternative = "two.sided", conf.level = 0.95,
                   correction_factor = 0.01, ...){
  test_2by2_tables(ai, n1i, ci, n2i)
  K <- length(ai) ## number of tables

  ## read numbers from the table
  bi <- n1i - ai
  di <- n2i - ci

  ## table totals
  ni <- n1i + n2i

  nn <- sum(ni)


  any0 <- apply(data.frame(ai, bi, ci, di), 1, function(rr) any(rr == 0))

  ## zero-cell correction for tables with zeroes
  if (any(any0)){
    for (jj in 1:K) {
      if (any0[jj]) {
        ai[jj] <- ai[jj] + correction_factor * n1i[jj] / ni[jj]
        bi[jj] <- bi[jj] + correction_factor * n1i[jj] / ni[jj]
        ci[jj] <- ci[jj] + correction_factor * n2i[jj] / ni[jj]
        di[jj] <- di[jj] + correction_factor * n2i[jj] / ni[jj]
      }
    }
    ## column totals
    n1i <- ai + bi
    n2i <- ci + di

    ## table totals
    ni <- n1i + n2i

    nn <- sum(ni)
  }

  gamma <- ni / nn
  delta <- n1i / ni

  p1 <- ai / n1i
  p0 <- ci / n2i


  ## study effect estimate: risk difference
  theta <- log(p1 / p0)

  ## scaled variance estimate
  vv <- (1 - p1) / (gamma * delta * p1) + (1 - p0) / (gamma * (1 - delta) * p0)

  ## point estimate: inverse-variance weighted average risk difference
  EST <- sum(theta / vv) / sum(1 / vv)

  ## variance and se estimate assuming homogeneity
  VAR_NAIVE <- 1 / sum(1 / vv)
  SE_NAIVE <- sqrt(VAR_NAIVE)

  ## constants for calculating the sandwich covariance matrix
  ww <- - (1 - p1) / (gamma * delta * p1) ^ 2 + (1 - p0) / (gamma * (1 - delta) * p0) ^ 2

  zz <- (1 - p1) / (gamma * delta * p1) ^ 3 + (1 - p0) / (gamma * (1 - delta) * p0) ^ 3

  ## common formula for estimated variance of inverse-variance weighted estimators
  ## allowing heterogeneity
  VAR <- (sum(1 / vv - 2 * theta * ww / vv ^ 3 + theta ^ 2 * zz / vv ^ 4) / sum(1 / vv) ^ 2 -
                2 * sum(- ww / vv ^ 3 + theta * zz / vv ^ 4) * sum(theta / vv) / sum(1 / vv) ^ 3 +
                sum(theta / vv) ^ 2 * sum(zz / vv ^ 4) / sum(1 / vv) ^ 4) / nn

  SE <- sqrt(VAR)


  ## construct confidence intervals
  CI_NAIVE <- wald_ci(EST, SE_NAIVE, alternative = alternative, conf.level = conf.level)
  CI <- wald_ci(EST, SE, alternative = alternative, conf.level = conf.level)


  return(list(EST = EST, SE_NAIVE = SE_NAIVE, VAR_NAIVE = VAR_NAIVE, CI_NAIVE = CI_NAIVE,
              SE = SE, VAR = VAR, CI = CI))

}
