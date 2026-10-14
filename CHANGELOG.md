# Changelog

## Unreleased

- Stream large quantile probabilities and inline model columns through temporary numeric files, load them with vectorized `scan()`, and encode numeric result arrays in one vectorized pass. The 200,000-probability smoke test enforces a five-second ceiling.
