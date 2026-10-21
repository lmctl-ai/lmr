# One-way ANOVA: plant weight by treatment group (PlantGrowth), as in ?PlantGrowth and most ANOVA tutorials.
source("r/_emit.R")
fit <- aov(weight ~ group, data = PlantGrowth)
tab <- summary(fit)[[1]]
print(tab)
lmfit <- lm(weight ~ group, data = PlantGrowth)
emit("anova-plantgrowth",
     df = tab$Df, sum_sq = tab$`Sum Sq`, mean_sq = tab$`Mean Sq`, f_value = tab$`F value`[1], p_value = tab$`Pr(>F)`[1],
     group_means = tapply(PlantGrowth$weight, PlantGrowth$group, mean), coefficients = coef(lmfit),
     coef_p_values = summary(lmfit)$coefficients[, 4])
