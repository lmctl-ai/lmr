# Autocorrelation and partial autocorrelation functions of lh and of the Nile flow series.
source("r/_emit.R")
emit("acf-pacf",
     lh_acf = acf(lh, plot = FALSE)$acf, lh_pacf = pacf(lh, plot = FALSE)$acf,
     nile_acf_20 = acf(Nile, lag.max = 20, plot = FALSE)$acf, nile_pacf_20 = pacf(Nile, lag.max = 20, plot = FALSE)$acf)
