# Three more lm idioms on R's datasets: an interaction (mpg ~ wt * hp), a log-log model (log(dist) ~ log(speed)) and
# weighted least squares (mpg ~ wt with weights 1 / disp).
source("r/_emit.R")
ia <- lm(mpg ~ wt * hp, data = mtcars)
ll <- lm(log(dist) ~ log(speed), data = cars)
wl <- lm(mpg ~ wt, data = mtcars, weights = 1 / disp)
print(summary(ia))
emit("lm-interaction-weights",
     interaction_coef = coef(ia), interaction_p = summary(ia)$coefficients[, 4], interaction_r2 = summary(ia)$r.squared,
     loglog_coef = coef(ll), loglog_sigma = summary(ll)$sigma, loglog_r2 = summary(ll)$r.squared,
     weighted_coef = coef(wl), weighted_se = summary(wl)$coefficients[, 2], weighted_sigma = summary(wl)$sigma,
     weighted_r2 = summary(wl)$r.squared)
