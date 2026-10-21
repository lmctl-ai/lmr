# Power and sample size: power.t.test (two-sample, one-sample, paired) and power.prop.test, plus the noncentral
# t, chi-square and F distributions behind them.
source("r/_emit.R")
a <- power.t.test(n = 20, delta = 1, sd = 1)
b <- power.t.test(power = 0.90, delta = 1)
c1 <- power.t.test(n = 15, delta = 0.8, sd = 1.2, type = "one.sample", alternative = "one.sided")
pp <- power.prop.test(p1 = 0.5, p2 = 0.75, power = 0.9)
print(a); print(b)
emit("power-noncentral",
     power_n20 = a$power, n_for_power90 = b$n, power_one_sample = c1$power, prop_n_for_power90 = pp$n,
     pt_ncp = pt(c(-1, 0, 1, 2, 3), df = 10, ncp = 1.5), pt_ncp_upper = pt(5, 20, ncp = 1, lower.tail = FALSE),
     qt_ncp = qt(c(0.05, 0.5, 0.95), df = 12, ncp = 2), dt_ncp = dt(c(0, 1, 2), df = 8, ncp = 1),
     pchisq_ncp = pchisq(c(1, 5, 10, 20), df = 4, ncp = 3), qchisq_ncp = qchisq(0.95, df = 3, ncp = 5),
     dchisq_ncp = dchisq(4, df = 4, ncp = 2), pf_ncp = pf(c(1, 2, 4), 3, 20, ncp = 4), qf_ncp = qf(0.9, 2, 30, ncp = 3))
