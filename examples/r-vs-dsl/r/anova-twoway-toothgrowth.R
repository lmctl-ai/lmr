# Two-way ANOVA with interaction: tooth length by supplement type and dose (ToothGrowth), the standard example,
# with sequential (type I) sums of squares as summary(aov()) prints them.
source("r/_emit.R")
tg <- ToothGrowth
tg$dose <- factor(tg$dose)
fit <- aov(len ~ supp * dose, data = tg)
tab <- summary(fit)[[1]]
print(tab)
emit("anova-twoway-toothgrowth",
     df = tab$Df, sum_sq = tab$`Sum Sq`, mean_sq = tab$`Mean Sq`, f_value = tab$`F value`[1:3], p_value = tab$`Pr(>F)`[1:3],
     cell_means = as.numeric(tapply(tg$len, list(tg$supp, tg$dose), mean)))
