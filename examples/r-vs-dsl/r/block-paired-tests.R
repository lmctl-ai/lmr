# Tests for blocked and paired designs from the R help pages: Friedman (RoundingTimes, ?friedman.test), Quade
# (dataFreq, ?quade.test) and McNemar (Performance, ?mcnemar.test).
source("r/_emit.R")
RoundingTimes <- matrix(c(5.40, 5.50, 5.55, 5.85, 5.70, 5.75, 5.20, 5.60, 5.50, 5.55, 5.50, 5.40, 5.90, 5.85, 5.70,
                          5.45, 5.55, 5.60, 5.40, 5.40, 5.35, 5.45, 5.50, 5.35, 5.25, 5.15, 5.00, 5.85, 5.80, 5.70,
                          5.25, 5.20, 5.10, 5.65, 5.55, 5.45, 5.60, 5.35, 5.45, 5.05, 5.00, 4.95, 5.50, 5.50, 5.40,
                          5.45, 5.55, 5.50, 5.55, 5.55, 5.35, 5.45, 5.50, 5.55, 5.50, 5.45, 5.25, 5.65, 5.60, 5.40,
                          5.70, 5.65, 5.55, 6.30, 6.30, 6.25), nrow = 22, byrow = TRUE)
fr <- friedman.test(RoundingTimes)
dataFreq <- matrix(c( 5, 4, 7, 10, 12, 1, 3, 1, 0, 2, 16, 12, 22, 22, 35, 5, 4, 3, 5, 4, 10, 9, 7, 13, 10,
                     19, 18, 28, 37, 58, 10, 7, 6, 8, 7), nrow = 7, byrow = TRUE)
qu <- quade.test(dataFreq)
Performance <- matrix(c(794, 86, 150, 570), nrow = 2)
mc <- mcnemar.test(Performance)
mc0 <- mcnemar.test(Performance, correct = FALSE)
print(fr); print(qu); print(mc)
emit("block-paired-tests",
     friedman_chi2 = fr$statistic, friedman_df = fr$parameter, friedman_p = fr$p.value,
     quade_f = qu$statistic, quade_df = qu$parameter, quade_p = qu$p.value,
     mcnemar_chi2 = mc$statistic, mcnemar_p = mc$p.value, mcnemar_uncorrected = mc0$statistic,
     mcnemar_uncorrected_p = mc0$p.value)
