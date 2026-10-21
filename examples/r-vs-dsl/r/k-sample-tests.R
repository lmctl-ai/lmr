# k-sample tests: Kruskal-Wallis on ozone by month (airquality, ?kruskal.test), Bartlett and Fligner-Killeen
# homogeneity of variances and the one-way analysis of means (equal and unequal variances) on InsectSprays.
source("r/_emit.R")
k <- kruskal.test(Ozone ~ Month, data = airquality)
b <- bartlett.test(count ~ spray, data = InsectSprays)
f <- fligner.test(count ~ spray, data = InsectSprays)
o1 <- oneway.test(count ~ spray, data = InsectSprays, var.equal = TRUE)
o2 <- oneway.test(count ~ spray, data = InsectSprays)
print(k); print(o2)
emit("k-sample-tests",
     kruskal_chi2 = k$statistic, kruskal_df = k$parameter, kruskal_p = k$p.value,
     bartlett_k2 = b$statistic, bartlett_p = b$p.value, fligner_chi2 = f$statistic, fligner_p = f$p.value,
     oneway_f = o1$statistic, oneway_p = o1$p.value, welch_f = o2$statistic, welch_df = o2$parameter, welch_p = o2$p.value)
