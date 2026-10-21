# Tests for counts from the R help pages: Fisher's Tea Tasting, binom.test(682, 925) (Mendel's peas) and
# prop.test on smokers among four patient groups.
source("r/_emit.R")
TeaTasting <- matrix(c(3, 1, 1, 3), nrow = 2)
f <- fisher.test(TeaTasting, alternative = "greater")
f2 <- fisher.test(TeaTasting)
b <- binom.test(682, 925, p = 0.75)
smokers <- c(83, 90, 129, 70); patients <- c(86, 93, 136, 82)
p <- prop.test(smokers, patients)
p2 <- prop.test(15, 50, p = 0.4)
print(f2); print(b); print(p)
emit("fisher-binom-prop",
     fisher_greater_p = f$p.value, fisher_p = f2$p.value, fisher_odds_ratio = f2$estimate, fisher_ci = f2$conf.int,
     binom_p = b$p.value, binom_ci = b$conf.int, binom_estimate = b$estimate,
     prop_x2 = p$statistic, prop_df = p$parameter, prop_p = p$p.value, prop_estimates = p$estimate,
     one_prop_x2 = p2$statistic, one_prop_p = p2$p.value, one_prop_ci = p2$conf.int)
