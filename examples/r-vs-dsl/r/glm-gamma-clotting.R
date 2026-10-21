# The Gamma example of ?glm (McCullagh & Nelder 1989, p. 300): blood clotting time against log plasma concentration
# for two lots of thromboplastin, with the canonical inverse link.
source("r/_emit.R")
clotting <- data.frame(u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
                       lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
                       lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12))
f1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
f2 <- glm(lot2 ~ log(u), data = clotting, family = Gamma)
s1 <- summary(f1)
print(s1)
emit("glm-gamma-clotting",
     lot1_estimate = coef(f1), lot1_std_error = s1$coefficients[, 2], lot1_p_value = s1$coefficients[, 4],
     lot1_dispersion = s1$dispersion, lot1_deviance = deviance(f1), lot1_aic = AIC(f1),
     lot2_estimate = coef(f2), lot2_dispersion = summary(f2)$dispersion, lot2_aic = AIC(f2))
