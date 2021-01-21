## construct confidence intervals
wald_ci <- function(EST, SE, alternative = "two.sided",
                     conf.level = 0.95) {
  alternative <- char.expand(alternative, c("two.sided", "less", "greater"))

  alpha <- 1 - conf.level
  if (alternative == "less") {
    CI <- c(EST + qnorm(alpha) * SE, Inf)
  } else if (alternative == "greater") {
    CI <- c(-Inf, EST + qnorm(conf.level) * SE)
  } else if (alternative == "two.sided") {
    CI <- EST + qnorm(c(alpha / 2, 1 - alpha / 2)) * SE
  } else {
    stop("'alternative' must be one of 'two.sided', 'less' or 'greater'")
  }

  return(CI)
}
