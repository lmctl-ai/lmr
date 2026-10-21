# The examples of ?integrate, ?uniroot and ?optimize: normal probabilities, an improper integral equal to pi,
# the root of x - 1/3 and the minimum of (x - 1/3)^2.
source("r/_emit.R")
i1 <- integrate(dnorm, -1.96, 1.96)
i2 <- integrate(dnorm, -Inf, Inf)
i3 <- integrate(function(x) 1 / ((x + 1) * sqrt(x)), lower = 0, upper = Inf)
u <- uniroot(function(x, a) x - a, c(0, 1), tol = 0.0001, a = 1/3)
o <- optimize(function(x, a) (x - a)^2, c(0, 1), tol = 0.0001, a = 1/3)
print(i3); print(u[1:3])
emit("integrate-roots",
     normal_95 = i1$value, normal_all = i2$value, improper_pi = i3$value, improper_abs_error = i3$abs.error,
     uniroot_root = u$root, uniroot_iter = u$iter, optimize_minimum = o$minimum, optimize_objective = o$objective)
