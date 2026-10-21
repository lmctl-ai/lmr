# Exact rank-correlation tests: the Hollander & Wolfe example of ?cor.test (no ties), Kendall and Spearman, one- and
# two-sided.
source("r/_emit.R")
x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c(2.6, 3.1, 2.5, 5.0, 3.6, 4.0, 5.2, 2.8, 3.8)
k <- cor.test(x, y, method = "kendall", alternative = "greater")
s <- cor.test(x, y, method = "spearman", alternative = "greater")
s2 <- cor.test(x, y, method = "spearman")
p <- cor.test(x, y, alternative = "greater")
print(k); print(s)
emit("cor-test-exact",
     kendall_tau = k$estimate, kendall_t = k$statistic, kendall_p = k$p.value, spearman_rho = s$estimate,
     spearman_s = s$statistic, spearman_p = s$p.value, spearman_two_sided_p = s2$p.value, pearson_p = p$p.value,
     pearson_ci = p$conf.int)
