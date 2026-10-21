# Wilcoxon tests: the exact paired test from ?wilcox.test (depression scores) and the normal approximation with ties
# and continuity correction for mpg by transmission (mtcars).
source("r/_emit.R")
x <- c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
a <- wilcox.test(x, y, paired = TRUE, alternative = "greater")
b <- suppressWarnings(wilcox.test(mpg ~ am, data = mtcars))
c2 <- suppressWarnings(wilcox.test(x, y, conf.int = TRUE))  # ties: normal approximation
print(a); print(b); print(c2)
emit("wilcox-tests",
     paired_v = a$statistic, paired_p = a$p.value, mtcars_w = b$statistic, mtcars_p = b$p.value,
     two_sample_w = c2$statistic, two_sample_p = c2$p.value, two_sample_ci = c2$conf.int,
     two_sample_location_shift = c2$estimate)
