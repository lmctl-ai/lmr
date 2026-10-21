# Kernel density estimation of Old Faithful eruption times (density(), ?density): bandwidth rules, the estimate at
# grid points for several kernels, plus the empirical CDF and histogram counts.
source("r/_emit.R")
x <- faithful$eruptions
d <- density(x)
de <- density(x, kernel = "epanechnikov", bw = 0.3)
F <- ecdf(x)
h <- hist(x, plot = FALSE)
idx <- c(1, 100, 200, 256, 300, 400, 512)
emit("density-ecdf",
     bw_nrd0 = bw.nrd0(x), bw_nrd = bw.nrd(x), bw_ucv = bw.ucv(x), bw_bcv = bw.bcv(x), bw_sj = bw.SJ(x),
     grid_from_to = c(d$x[1], d$x[512]), gaussian_y = d$y[idx], epanechnikov_y = de$y[idx],
     ecdf_at = F(c(2, 3, 4, 4.5)), hist_breaks = h$breaks, hist_counts = h$counts)
