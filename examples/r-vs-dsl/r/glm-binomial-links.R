# One binary outcome, three links: transmission type on weight (mtcars) with logit, probit and complementary log-log,
# compared by AIC, plus a dose-response with grouped binomial data and trials as weights (?glm style).
source("r/_emit.R")
lg <- glm(am ~ wt, family = binomial(link = "logit"), data = mtcars)
pr <- glm(am ~ wt, family = binomial(link = "probit"), data = mtcars)
cl <- glm(am ~ wt, family = binomial(link = "cloglog"), data = mtcars)
ldose <- 0:5; dead <- c(1, 4, 9, 13, 18, 20); total <- rep(20, 6)
gr <- glm(dead / total ~ ldose, family = binomial, weights = total)
print(summary(pr))
emit("glm-binomial-links",
     logit_coef = coef(lg), logit_aic = AIC(lg), probit_coef = coef(pr), probit_se = sqrt(diag(vcov(pr))),
     probit_aic = AIC(pr), cloglog_coef = coef(cl), cloglog_aic = AIC(cl), cloglog_iter = cl$iter,
     grouped_coef = coef(gr), grouped_se = sqrt(diag(vcov(gr))), grouped_deviance = deviance(gr),
     grouped_ld50 = -coef(gr)[1] / coef(gr)[2])
