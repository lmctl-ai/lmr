# Runtime and contribution parity

The production runtime is Ubuntu 24.04 with `r-base-core=4.3.3-2build2`. `docker/Dockerfile` pins the Ubuntu 24.04 multi-architecture manifest digest and the public `lmcpp` source commit. Its `test` target runs the R-oracle CTest suite; its `runtime` target contains the rcli executable and R runtime without a compiler. Build on the host architecture (linux/arm64 or linux/amd64) with Docker's native platform selection.

The `bin/rcli` host shim starts one isolated runtime container per call. It requires the current directory to be below `$HOME`, mounts only that directory read-write, uses `--network=none`, passes through uid/gid and only the `RCLI_TIMEOUT_S` control, and places child temporary files in a private `.rcli-tmp` directory inside the mounted working directory. Docker startup failures return `SERVICE_UNAVAILABLE`; a container crash without a complete result returns `INTERNAL_ERROR` with `TOOL_CRASH`. Exit status and valid tool envelopes pass through.

## Adding an operation

For each new operation:

1. Add its schema entry with descriptions for every parameter, defaults, units, and the result shape.
2. Validate types, domains, input sizes, computational complexity, and output headroom in C++ before launching R.
3. Add a fixed worker implementation; caller text must never be parsed or evaluated as R code.
4. Add an R-computed oracle test, using same-environment expected values for version-sensitive behavior, plus validation and boundary cases.
5. Update the README operation table and limits, and add or update an lmscript recipe that imports the published schema and checks numeric output.
6. Regenerate `cli.schema.dev.json`, parse it with lmts, then run CTest, the host smoke, and all recipes.
