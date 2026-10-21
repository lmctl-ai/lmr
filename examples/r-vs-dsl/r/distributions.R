# The distribution-function table of a typical introductory R course: d/p/q for the normal, t, chi-square, F,
# binomial, Poisson, exponential, gamma and beta laws.
source("r/_emit.R")
emit("distributions",
     dnorm_0 = dnorm(0), pnorm_196 = pnorm(1.96), qnorm_tails = qnorm(c(0.025, 0.975)), pnorm_upper_5 = pnorm(5, lower.tail = FALSE),
     pt_21_15 = pt(2.1, 15), qt_975 = qt(0.975, c(5, 10, 30)), pchisq_384 = pchisq(3.84, 1),
     qchisq_95 = qchisq(0.95, 1:5), qf_95_3_20 = qf(0.95, 3, 20), pf_3_2_20 = pf(3, 2, 20, lower.tail = FALSE),
     dbinom_10 = dbinom(0:10, 10, 0.5), pbinom_3 = pbinom(3, 10, 0.5), qbinom = qbinom(c(0.1, 0.5, 0.9), 20, 0.3),
     dpois_2 = dpois(0:5, 2), ppois_4 = ppois(4, 2), qpois = qpois(c(0.05, 0.95), 4),
     pexp_1 = pexp(1, 2), qexp_half = qexp(0.5, 0.1), qgamma_half = qgamma(0.5, 2, 1), pbeta_03 = pbeta(0.3, 2, 5),
     dlnorm_1 = dlnorm(1.5, 0, 1), punif = punif(0.3, -1, 1))
