# Confidence and prediction intervals for new observations: predict(..., interval = "confidence"/"prediction"),
# for stopping distance at five speeds (cars).
source("r/_emit.R")
fit <- lm(dist ~ speed, data = cars)
new <- data.frame(speed = c(5, 10, 15, 20, 25))
ci <- predict(fit, new, interval = "confidence")
pi <- predict(fit, new, interval = "prediction", level = 0.9)
print(ci); print(pi)
emit("lm-prediction-intervals",
     fit = ci[, "fit"], conf_lower = ci[, "lwr"], conf_upper = ci[, "upr"], pred90_lower = pi[, "lwr"],
     pred90_upper = pi[, "upr"])
