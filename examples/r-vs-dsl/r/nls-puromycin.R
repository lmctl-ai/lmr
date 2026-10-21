# Michaelis-Menten fit of enzyme velocity on substrate concentration for the treated cells (Puromycin), as in ?nls.
source("r/_emit.R")
treated <- subset(Puromycin, state == "treated")
fit <- nls(rate ~ Vm * conc / (K + conc), data = treated, start = c(Vm = 200, K = 0.05))
s <- summary(fit)
print(s)
emit("nls-puromycin",
     estimate = coef(fit), std_error = s$coefficients[, 2], p_value = s$coefficients[, 4], sigma = s$sigma,
     rss = deviance(fit), fitted = fitted(fit))
