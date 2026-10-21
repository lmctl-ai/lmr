# The airline model: seasonal ARIMA(0,1,1)(0,1,1)[12] on log(AirPassengers), with a 12-month forecast.
source("r/_emit.R")
fit <- arima(log(AirPassengers), order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
fc <- predict(fit, n.ahead = 12)
print(fit)
emit("arima-airline",
     coef = coef(fit), se = sqrt(diag(fit$var.coef)), sigma2 = fit$sigma2, loglik = fit$loglik, aic = fit$aic,
     forecast_log = fc$pred, forecast_se = fc$se, forecast_passengers = exp(fc$pred))
