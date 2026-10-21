# Student's sleep data: Welch, pooled-variance and paired t tests (as in ?sleep and ?t.test).
source("r/_emit.R")
x <- sleep$extra[sleep$group == 1]; y <- sleep$extra[sleep$group == 2]
w <- t.test(extra ~ group, data = sleep)
v <- t.test(extra ~ group, data = sleep, var.equal = TRUE)
p <- t.test(x, y, paired = TRUE)
print(w); print(p)
emit("t-test-sleep",
     welch_t = w$statistic, welch_df = w$parameter, welch_p = w$p.value, welch_ci = w$conf.int,
     pooled_t = v$statistic, pooled_df = v$parameter, pooled_p = v$p.value, pooled_ci = v$conf.int,
     paired_t = p$statistic, paired_df = p$parameter, paired_p = p$p.value, paired_ci = p$conf.int,
     paired_mean_diff = p$estimate)
