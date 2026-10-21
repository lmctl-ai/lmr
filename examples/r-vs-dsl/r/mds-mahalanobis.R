# Classical multidimensional scaling of road distances between 21 European cities (cmdscale(eurodist), ?cmdscale)
# and Mahalanobis distances of the iris setosa measurements from their mean.
source("r/_emit.R")
m <- cmdscale(eurodist, k = 2, eig = TRUE)
se <- iris[iris$Species == "setosa", 1:4]
d2 <- mahalanobis(se, colMeans(se), cov(se))
print(head(m$points))
emit("mds-mahalanobis",
     abs_points = abs(m$points), eig_top5 = m$eig[1:5], gof = m$GOF,
     mahalanobis = d2, mahalanobis_max = max(d2), cor_matrix = cor(se))
