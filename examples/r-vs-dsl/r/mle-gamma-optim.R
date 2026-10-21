# Maximum likelihood by hand: fit a gamma distribution to the Old Faithful waiting times by minimising the negative
# log-likelihood with optim (BFGS on log shape and log rate), with standard errors from the Hessian.
source("r/_emit.R")
x <- faithful$waiting
nll <- function(th) -sum(dgamma(x, shape = exp(th[1]), rate = exp(th[2]), log = TRUE))
o <- optim(c(log(10), log(0.1)), nll, method = "BFGS", hessian = TRUE)
se <- sqrt(diag(solve(o$hessian)))
print(exp(o$par))
emit("mle-gamma-optim",
     shape_rate = exp(o$par), log_par = o$par, log_par_se = se, min_nll = o$value, convergence = o$convergence)
