# Pairwise comparisons after a one-way ANOVA: pairwise.t.test(weight, group) with the pooled SD and Holm's
# adjustment (the R default), plus Bonferroni, and the unadjusted p values (PlantGrowth).
source("r/_emit.R")
h <- pairwise.t.test(PlantGrowth$weight, PlantGrowth$group)
b <- pairwise.t.test(PlantGrowth$weight, PlantGrowth$group, p.adjust.method = "bonferroni")
n <- pairwise.t.test(PlantGrowth$weight, PlantGrowth$group, p.adjust.method = "none")
print(h)
# Lower triangle in R's order: trt1-ctrl, trt2-ctrl, trt2-trt1.
lt <- function(m) c(m[1, 1], m[2, 1], m[2, 2])
emit("pairwise-t-plantgrowth", holm = lt(h$p.value), bonferroni = lt(b$p.value), unadjusted = lt(n$p.value))
