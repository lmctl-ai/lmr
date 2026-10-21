# Hierarchical clustering of the 50 US states on scaled USArrests (complete and average linkage), cut into 4 groups.
source("r/_emit.R")
d <- dist(scale(USArrests))
hc <- hclust(d, "complete")
ha <- hclust(d, "average")
print(table(cutree(hc, 4)))
emit("hclust-usarrests",
     complete_heights = hc$height, complete_order = hc$order, complete_cut4 = cutree(hc, 4),
     complete_cut4_sizes = as.numeric(table(cutree(hc, 4))), average_heights = ha$height,
     average_cut4 = cutree(ha, 4))
