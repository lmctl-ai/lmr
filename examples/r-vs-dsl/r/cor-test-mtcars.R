# Correlation tests of fuel economy against weight (mtcars): Pearson, Spearman and Kendall.
source("r/_emit.R")
p <- cor.test(mtcars$mpg, mtcars$wt)
s <- suppressWarnings(cor.test(mtcars$mpg, mtcars$wt, method = "spearman", exact = FALSE))
k <- suppressWarnings(cor.test(mtcars$mpg, mtcars$wt, method = "kendall", exact = FALSE))
print(p)
emit("cor-test-mtcars",
     pearson_r = p$estimate, pearson_t = p$statistic, pearson_p = p$p.value, pearson_ci = p$conf.int,
     spearman_rho = s$estimate, spearman_s = s$statistic, spearman_p = s$p.value,
     kendall_tau = k$estimate, kendall_z = k$statistic, kendall_p = k$p.value)
