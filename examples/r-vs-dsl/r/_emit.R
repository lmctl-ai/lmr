# Writes a recipe's results to expected/NAME.json: every value as a number or an array of numbers, printed with 17
# significant digits so the comparison sees R's exact doubles. Base R only (no jsonlite).
emit <- function(name, ...) {
  v <- list(...)
  num <- function(x) {
    s <- sprintf("%.17g", x)
    s[is.na(x)] <- "null"
    s[is.infinite(x) & x > 0] <- "1e999"
    s[is.infinite(x) & x < 0] <- "-1e999"
    s
  }
  enc <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 1) num(x) else paste0("[", paste(num(x), collapse = ", "), "]")
  }
  body <- paste0('    "', names(v), '": ', vapply(v, enc, ""), collapse = ",\n")
  out <- c("{", sprintf('  "recipe": "%s",', name), sprintf('  "r_version": "%s",', R.version.string),
           '  "values": {', body, "  }", "}")
  writeLines(out, file.path("expected", paste0(name, ".json")))
}
