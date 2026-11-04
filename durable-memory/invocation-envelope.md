# Invocation envelope for the native host (2026-11-04 port)

The native `lmr` host now answers `--lmcode` with the SDK one-envelope frame
(`LMCODE_RETURN_BEGIN` / `END` around one JSON result), replacing one result line per file.
The vendored SDK moved to lmcode-cpp origin/main (invocation.hpp and result_contract.hpp joined
the vendored set; native/sdk.lock pin, sync via native/scripts/sync-sdk.sh).

Shape (same family as lmweather/lmphysics/lmfinance):

- `--lmcode FILE [FILE ...]` (plain runs) and `--lmcode FILE [--lmcode FILE ...]` (builder slot
  fragments) are both legal argv; slot order comes from the invocation metadata.
- Every file is loaded and verified before anything runs; failures are `results.fail(index, phase, e)`,
  phases map to the exit codes: cli 2, load/run 3, app 4 (`APP_FAIL` from out.fail/out.check).
  Load/run errors with several files keep the located `file: message` text (stderr and the frame).
- The success payload keeps the pre-envelope shape exactly: `{values: {...in first-set order...},
  checks: N, failures: N}` per slot (builder frames keep one slot result per part plus `parts: N`).
  The host summary is text-built by `result_json(host.out)`, parsed with the SDK's own Json parser
  and handed to the slot's ResultSink as a typed `Json` value, so no textual frame assembly.
- Dry mode is invocation-wide from LMCODE_INVOCATION_FILE metadata (`lmscript --dry-run`); a dry
  invocation loads and validates every chunk but runs none; sinks get the empty summary (there is
  no effect-free stats run).
- lmscript side: examples/dsl/lib/lmr.lms builds one-off runs with `lmcode_build([host argv0])` +
  `b.add(program)` + `run_lmcode(b)`, so a fixture compiler is only needed for native/tests
  (native/scripts/compile-fixture.mjs on the compiler library, same as lmok's).
- Version: lmr 0.5.2 (VERSION file; CMake and host print follow it).
