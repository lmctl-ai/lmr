# k-means on the four iris measurements with k = 3 and 25 random starts, the standard clustering tutorial.
# Random starts differ between R and lmr, so only label-free results are compared (sorted sizes and
# within-cluster sums of squares, total and between sums of squares).
source("r/_emit.R")
set.seed(42)
km <- kmeans(iris[, 1:4], centers = 3, nstart = 25)
print(km$size); print(km$tot.withinss)
o <- order(km$size)
emit("kmeans-iris",
     sizes_sorted = sort(km$size), withinss_by_size = km$withinss[o], tot_withinss = km$tot.withinss,
     betweenss = km$betweenss, totss = km$totss,
     sepal_length_centre_by_size = km$centers[o, "Sepal.Length"])
