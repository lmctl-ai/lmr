# Poisson regression of yarn breaks on wool type and tension (warpbreaks), with factor (dummy) predictors.
source("r/_emit.R")
fit <- glm(breaks ~ wool + tension, family = poisson, data = warpbreaks)
s <- summary(fit)
print(s)
emit("glm-poisson-warpbreaks",
     estimate = coef(fit), std_error = s$coefficients[, 2], p_value = s$coefficients[, 4], deviance = deviance(fit),
     null_deviance = fit$null.deviance, aic = AIC(fit), dispersion_pearson = sum(residuals(fit, "pearson")^2) / fit$df.residual)
