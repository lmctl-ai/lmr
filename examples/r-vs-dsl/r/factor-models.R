# Models with factor predictors as R users write them: lm(len ~ supp * factor(dose)) with its ANOVA table,
# lm(mpg ~ wt + factor(cyl)), a nested-model F test, and a Poisson glm with factors and an analysis of deviance.
source("r/_emit.R")
tg <- lm(len ~ supp * factor(dose), data = ToothGrowth)
at <- anova(tg)
cy <- lm(mpg ~ wt + factor(cyl), data = mtcars)
small <- lm(mpg ~ wt, data = mtcars)
nest <- anova(small, cy)
# From the CSV, as lmr reads it: read.csv gives character columns, so the levels sort alphabetically (tension H, L, M;
# the warpbreaks dataset itself orders them L, M, H, which changes the baseline but not the fit).
wbd <- read.csv("data/warpbreaks.csv")
wb <- glm(breaks ~ wool * tension, family = poisson, data = wbd)
ad <- anova(wb, test = "Chisq")
pg <- lm(weight ~ group, data = PlantGrowth)
print(at); print(ad)
emit("factor-models",
     tooth_coef = coef(tg), tooth_anova_ss = at$`Sum Sq`, tooth_anova_f = at$`F value`[1:3], tooth_anova_p = at$`Pr(>F)`[1:3],
     cyl_coef = coef(cy), cyl_se = summary(cy)$coefficients[, 2], nested_f = nest$F[2], nested_p = nest$`Pr(>F)`[2],
     warp_coef = coef(wb), warp_deviance = ad$Deviance[2:4], warp_resid_dev = ad$`Resid. Dev`, warp_p = ad$`Pr(>Chi)`[2:4],
     plant_coef = coef(pg), plant_anova_f = anova(pg)$`F value`[1])
