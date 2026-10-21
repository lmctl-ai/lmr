# Logistic regression: transmission type from weight and horsepower (mtcars).
source("r/_emit.R")
fit <- glm(am ~ wt + hp, family = binomial, data = mtcars)
s <- summary(fit)
print(s)
emit("glm-logistic-mtcars",
     estimate = coef(fit), std_error = s$coefficients[, 2], z_value = s$coefficients[, 3], p_value = s$coefficients[, 4],
     deviance = deviance(fit), null_deviance = fit$null.deviance, aic = AIC(fit), iterations = fit$iter,
     predict_response = predict(fit, data.frame(wt = c(2.5, 3.5), hp = c(120, 120)), type = "response"),
     fitted = fitted(fit))
