# Third-party notices

## R

This project invokes the separately installed `Rscript` executable as a child process. It does not link to `libR`, `libRmath`, or embed R in the `rcli` binary. R is distributed under the GNU General Public License, version 2 or version 3, at the recipient's choice. Consult the R installation's `COPYING` file and the R Project licensing page for the terms applicable to that installation.

The worker uses base R and the recommended `stats` and `utils` packages shipped with the R distribution. It does not install or require third-party R packages.

## Native `rcli-dsl` binaries

The `rcli-dsl` release binaries contain no R code and do not use R at run time. They contain only this project's
code, which is Copyright (c) 2026 Lattice Systems LLC under the MIT license (see LICENSE). That includes the
vendored `native/cpp/dslvm.hpp` runtime from the cpp-dsl-sdk project (MIT).

- **Ported public-domain algorithms**: the numerical integration (`native/src/num/quad.cpp`) follows QUADPACK
  (Piessens, de Doncker-Kapenga, Überhuber and Kahaner, 1983), and the STL decomposition follows the netlib STL
  of Cleveland et al. (1990). Both are in the public domain.
- **Random numbers**: xoshiro256** and splitmix64 (Blackman and Vigna; Steele, Lea and Flood) are in the public
  domain.
- **Other algorithms**: implemented from the published papers cited in the source comments.
- **C and C++ runtime**: on Linux, libstdc++ and libgcc are linked statically under the GCC Runtime Library
  Exception, and glibc (LGPL-2.1) is linked dynamically. The Windows build links the MinGW-w64 runtime statically
  (permissive licences). macOS uses its system libraries.
