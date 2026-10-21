# Seasonal decomposition by loess, the two examples of ?stl: nottem with a periodic seasonal component and
# log(co2) with s.window = 21.
source("r/_emit.R")
a <- stl(nottem, "periodic")
b <- stl(log(co2), s.window = 21)
emit("stl-decomposition",
     nottem_seasonal = a$time.series[, "seasonal"], nottem_trend = a$time.series[, "trend"],
     nottem_remainder = a$time.series[, "remainder"], co2_seasonal = b$time.series[, "seasonal"],
     co2_trend = b$time.series[, "trend"])
