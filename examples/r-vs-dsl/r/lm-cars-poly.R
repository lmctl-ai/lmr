# Quadratic regression of stopping distance on speed (cars), the classic polynomial-term example.
source("r/_emit.R")
fit <- lm(dist ~ speed + I(speed^2), data = cars)
s <- summary(fit)
print(s)
pr <- predict(fit, data.frame(speed = c(10, 20, 25)), se.fit = TRUE)
emit("lm-cars-poly",
     estimate = coef(fit), std_error = s$coefficients[, 2], p_value = s$coefficients[, 4], sigma = s$sigma,
     r_squared = s$r.squared, f_statistic = s$fstatistic[1], predict = pr$fit, predict_se = pr$se.fit,
     fitted = fitted(fit))
