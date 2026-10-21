# ARIMA models of the luteinizing hormone series lh, from ?arima and ?predict.Arima.
source("r/_emit.R")
a1 <- arima(lh, order = c(1, 0, 0))
a3 <- arima(lh, order = c(3, 0, 0))
a11 <- arima(lh, order = c(1, 0, 1))
fc <- predict(a3, n.ahead = 12)
print(a1); print(a3)
emit("arima-lh",
     ar1_coef = coef(a1), ar1_se = sqrt(diag(a1$var.coef)), ar1_sigma2 = a1$sigma2, ar1_loglik = a1$loglik, ar1_aic = a1$aic,
     ar3_coef = coef(a3), ar3_se = sqrt(diag(a3$var.coef)), ar3_aic = a3$aic,
     arma11_coef = coef(a11), arma11_aic = a11$aic, forecast = fc$pred, forecast_se = fc$se)
