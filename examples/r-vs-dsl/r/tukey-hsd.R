# Tukey's honest significant differences after a one-way ANOVA (warpbreaks: breaks by tension), as in ?TukeyHSD,
# and the studentized range distribution ptukey/qtukey.
source("r/_emit.R")
fit <- aov(breaks ~ tension, data = warpbreaks)
th <- TukeyHSD(fit)$tension
print(th)
emit("tukey-hsd",
     diff = th[, "diff"], lwr = th[, "lwr"], upr = th[, "upr"], p_adj = th[, "p adj"],
     ptukey = ptukey(c(1, 2, 3.5, 5), nmeans = 3, df = 51), ptukey_inf_df = ptukey(3, 4, Inf),
     qtukey = qtukey(c(0.9, 0.95, 0.99), nmeans = 3, df = 51), qtukey_k5 = qtukey(0.95, 5, 20))
