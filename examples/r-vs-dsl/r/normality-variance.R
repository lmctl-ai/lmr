# Normality and variance checks: Shapiro-Wilk on Old Faithful eruptions and on mpg, the F test of equal variances
# for mpg by transmission, and a one-sample Kolmogorov-Smirnov test against N(20, 6).
source("r/_emit.R")
s1 <- shapiro.test(faithful$eruptions)
s2 <- shapiro.test(mtcars$mpg)
v <- var.test(mpg ~ am, data = mtcars)
k <- suppressWarnings(ks.test(mtcars$mpg, "pnorm", 20, 6))
print(s1); print(v); print(k)
emit("normality-variance",
     shapiro_faithful_w = s1$statistic, shapiro_faithful_p = s1$p.value, shapiro_mpg_w = s2$statistic,
     shapiro_mpg_p = s2$p.value, var_f = v$statistic, var_p = v$p.value, var_ci = v$conf.int,
     ks_d = k$statistic, ks_p = k$p.value)
