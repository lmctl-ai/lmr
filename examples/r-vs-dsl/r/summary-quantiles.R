# Descriptive statistics of mpg (mtcars): summary(), the nine quantile types, spread measures and a correlation
# matrix with Pearson and Spearman coefficients.
source("r/_emit.R")
x <- mtcars$mpg
print(summary(x))
probs <- c(0.1, 0.25, 0.5, 0.9)
q <- unlist(lapply(1:9, function(t) quantile(x, probs, type = t)))
m <- cor(mtcars[, c("mpg", "wt", "hp")])
sp <- cor(mtcars[, c("mpg", "wt", "hp")], method = "spearman")
emit("summary-quantiles",
     summary = summary(x), quantiles_types_1_to_9 = q, sd = sd(x), var = var(x), iqr = IQR(x), mad = mad(x),
     pearson = m[upper.tri(m)], spearman = sp[upper.tri(sp)], kendall_mpg_wt = cor(x, mtcars$wt, method = "kendall"),
     cov_mpg_wt = cov(x, mtcars$wt))
