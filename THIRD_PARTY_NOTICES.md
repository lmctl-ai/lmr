# Third-party notices

## R

This project invokes the separately installed `Rscript` executable as a child process. It does not link to `libR`, `libRmath`, or embed R in the `rcli` binary. R is distributed under the GNU General Public License, version 2 or version 3, at the recipient's choice. Consult the R installation's `COPYING` file and the R Project licensing page for the terms applicable to that installation.

The worker uses base R and the recommended `stats` and `utils` packages shipped with the R distribution. It does not install or require third-party R packages.
