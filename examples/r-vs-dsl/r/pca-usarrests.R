# Principal component analysis: prcomp(USArrests, scale. = TRUE), the classic PCA tutorial, prcomp of iris
# measurements without scaling, and the eigenvalues of cor(mtcars). Component signs are arbitrary, so loadings and
# scores are compared in absolute value.
source("r/_emit.R")
p <- prcomp(USArrests, scale. = TRUE)
s <- summary(p)$importance
q <- prcomp(iris[, 1:4])
e <- eigen(cor(mtcars), symmetric = TRUE)
print(summary(p))
emit("pca-usarrests",
     sdev = p$sdev, abs_rotation = abs(p$rotation), abs_scores_first5 = abs(p$x[1:5, ]),
     proportion = s["Proportion of Variance", ], cumulative = s["Cumulative Proportion", ],
     iris_sdev = q$sdev, iris_abs_pc1 = abs(q$rotation[, 1]), mtcars_cor_eigenvalues = e$values)
