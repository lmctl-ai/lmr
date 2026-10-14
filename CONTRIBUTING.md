# Contributing

Keep the R boundary process-based. Do not include `Rmath.h`, link `-lR` or `-lRmath`, or link libR by another name. Do not copy R implementation source into this project. Calls must use the fixed operation catalogue and must not evaluate caller-provided R code.

Run `scripts/check-no-libR.sh`, build with the `rcli-dev` image, and run CTest before committing. Commit locally; this repository is private and the maintainer pushes changes.
