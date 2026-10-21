# Exports the R datasets the recipes use to data/*.csv, so the lmr side runs without R.
# Run from examples/r-vs-dsl: Rscript r/export-data.R
dir.create("data", showWarnings = FALSE)
put <- function(df, name) write.csv(df, file.path("data", paste0(name, ".csv")), row.names = FALSE, quote = TRUE)
ts_frame <- function(x) data.frame(time = as.numeric(time(x)), value = as.numeric(x))

put(cbind(model = rownames(mtcars), mtcars), "mtcars")
put(data.frame(extra = sleep$extra, group = as.character(sleep$group), ID = as.integer(sleep$ID)), "sleep")
put(cars, "cars")
put(data.frame(weight = PlantGrowth$weight, group = as.character(PlantGrowth$group)), "PlantGrowth")
put(data.frame(breaks = warpbreaks$breaks, wool = as.character(warpbreaks$wool),
               tension = as.character(warpbreaks$tension)), "warpbreaks")
put(data.frame(iris[, 1:4], Species = as.character(iris$Species)), "iris")
put(cbind(state = rownames(USArrests), USArrests), "USArrests")
put(faithful, "faithful")
put(data.frame(conc = Puromycin$conc, rate = Puromycin$rate, state = as.character(Puromycin$state)), "Puromycin")
put(ts_frame(lh), "lh")
put(ts_frame(AirPassengers), "AirPassengers")
put(ts_frame(nottem), "nottem")
put(ts_frame(co2), "co2")
put(ts_frame(Nile), "Nile")
