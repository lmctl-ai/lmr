# Model and check: ARIMA(1,1,1) for the Nile flow, then the Ljung-Box test on the residuals (Box.test with
# fitdf = 2), and the same model by ML and CSS.
source("r/_emit.R")
fit <- arima(Nile, order = c(1, 1, 1))
bt <- Box.test(residuals(fit), lag = 10, type = "Ljung-Box", fitdf = 2)
ml <- arima(Nile, order = c(1, 1, 1), method = "ML")
css <- arima(Nile, order = c(1, 1, 1), method = "CSS")
print(fit); print(bt)
emit("arima-nile-diagnostics",
     coef = coef(fit), se = sqrt(diag(fit$var.coef)), sigma2 = fit$sigma2, aic = fit$aic,
     ljung_box_q = bt$statistic, ljung_box_p = bt$p.value, ml_coef = coef(ml), css_coef = coef(css),
     css_sigma2 = css$sigma2, forecast = predict(fit, n.ahead = 5)$pred)
