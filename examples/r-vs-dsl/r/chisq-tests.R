# Chi-squared tests from ?chisq.test: independence in a 2 x 3 table of gender by party, and goodness of fit.
source("r/_emit.R")
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
a <- chisq.test(M)
g <- chisq.test(c(89, 37, 30, 28, 2), p = c(40, 20, 20, 15, 5) / 100)
k <- chisq.test(matrix(c(12, 5, 7, 7), 2))
print(a); print(g); print(k)
emit("chisq-tests",
     independence_x2 = a$statistic, independence_df = a$parameter, independence_p = a$p.value, expected = a$expected,
     gof_x2 = g$statistic, gof_p = g$p.value, yates_x2 = k$statistic, yates_p = k$p.value)
