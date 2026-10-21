# Multiple regression on mtcars: the first model of most R regression tutorials.
source("r/_emit.R")
fit <- lm(mpg ~ wt + hp, data = mtcars)
s <- summary(fit)
print(s)
ci <- confint(fit)
new <- data.frame(wt = c(2.5, 3.5), hp = c(100, 200))
pr <- predict(fit, new, se.fit = TRUE)
emit("lm-mtcars",
     estimate = coef(fit), std_error = s$coefficients[, 2], t_value = s$coefficients[, 3],
     p_value = s$coefficients[, 4], sigma = s$sigma, r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
     f_statistic = s$fstatistic[1], f_p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
     conf_low = ci[, 1], conf_high = ci[, 2], loglik = logLik(fit), aic = AIC(fit), bic = BIC(fit),
     predict = pr$fit, predict_se = pr$se.fit, residuals = residuals(fit))
