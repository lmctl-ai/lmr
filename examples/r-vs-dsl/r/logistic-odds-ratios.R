# Reporting a logistic regression as odds ratios: exp(coef) with Wald intervals exp(confint.default()), and the
# likelihood-ratio test against the null model (mtcars, vs ~ mpg + am).
source("r/_emit.R")
fit <- glm(vs ~ mpg + am, family = binomial, data = mtcars)
or <- exp(cbind(OR = coef(fit), confint.default(fit)))
print(or)
lr <- fit$null.deviance - fit$deviance
emit("logistic-odds-ratios",
     odds_ratio = or[, 1], or_lower = or[, 2], or_upper = or[, 3], lr_statistic = lr,
     lr_p = pchisq(lr, fit$df.null - fit$df.residual, lower.tail = FALSE))
