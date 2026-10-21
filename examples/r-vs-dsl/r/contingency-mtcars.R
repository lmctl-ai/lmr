# Contingency tables from a data frame: cylinders by gears (chi-squared, small expected counts), transmission by
# engine shape (Fisher's exact test and Yates' chi-squared) on mtcars.
source("r/_emit.R")
t1 <- table(mtcars$cyl, mtcars$gear)
t2 <- table(mtcars$am, mtcars$vs)
c1 <- suppressWarnings(chisq.test(t1))
f2 <- fisher.test(t2)
c2 <- chisq.test(t2)
print(t1); print(c1); print(f2)
emit("contingency-mtcars",
     cyl_gear_counts = as.numeric(t1), cyl_gear_x2 = c1$statistic, cyl_gear_df = c1$parameter, cyl_gear_p = c1$p.value,
     am_vs_counts = as.numeric(t2), fisher_p = f2$p.value, fisher_odds_ratio = f2$estimate, fisher_ci = f2$conf.int,
     yates_x2 = c2$statistic, yates_p = c2$p.value)
