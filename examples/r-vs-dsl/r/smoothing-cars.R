# Smoothing stopping distance on speed (cars): a smoothing spline with 5 degrees of freedom, loess with the default
# span (R's default interpolated surface and the exact "direct" surface), and interpolating splines from ?splinefun.
source("r/_emit.R")
ss <- smooth.spline(cars$speed, cars$dist, df = 5)
lo_direct <- loess(dist ~ speed, cars, control = loess.control(surface = "direct"))
lo_default <- loess(dist ~ speed, cars)
x <- 1:10; y <- sin(x)
xx <- c(1.5, 2.25, 5.5, 9.9)
emit("smoothing-cars",
     ss_fit_at = predict(ss, c(5, 10, 15, 20, 25))$y, ss_df = ss$df,
     loess_direct_fitted = fitted(lo_direct), loess_default_fitted = fitted(lo_default),
     spline_fmm = splinefun(x, y, "fmm")(xx), spline_natural = splinefun(x, y, "natural")(xx),
     approx_linear = approx(x, y, xout = xx)$y)
