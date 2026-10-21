# Rosenbrock's banana function from ?optim, minimised from (-1.2, 1) with Nelder-Mead, BFGS and L-BFGS-B
# (numerical gradients).
source("r/_emit.R")
fr <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
nm <- optim(c(-1.2, 1), fr)
bf <- optim(c(-1.2, 1), fr, method = "BFGS")
lb <- optim(c(-1.2, 1), fr, method = "L-BFGS-B")
print(nm[1:3]); print(bf[1:3])
emit("optim-rosenbrock",
     nm_par = nm$par, nm_value = nm$value, nm_fn_count = nm$counts[1],
     bfgs_par = bf$par, bfgs_value = bf$value, bfgs_fn_count = bf$counts[1], bfgs_gr_count = bf$counts[2],
     lbfgsb_par = lb$par)
