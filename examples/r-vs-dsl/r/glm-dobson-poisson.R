# The Poisson example of ?glm (Dobson 1990, Table 2.2): counts by outcome and treatment, three levels each.
source("r/_emit.R")
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
fit <- glm(counts ~ outcome + treatment, family = poisson())
s <- summary(fit)
print(s)
emit("glm-dobson-poisson",
     estimate = coef(fit), std_error = s$coefficients[, 2], z_value = s$coefficients[, 3], p_value = s$coefficients[, 4],
     deviance = deviance(fit), null_deviance = fit$null.deviance, aic = AIC(fit), fitted = fitted(fit),
     deviance_residuals = residuals(fit, "deviance"))
