# rcli status and hand-off (2026-10-19)

## State
- `main`: C++17 lmcpp front end, spawns `Rscript` per call in a pinned offline docker image, hardened host shim `bin/rcli`, 28 ops, lmapi schema via `--schema` / `--schema --json`.
- `examples/`: 14 `tut-*.lms` tutorial conversions (`tut-vehicles.csv`, `scripts/generate-tutorial-data.py`, each recipe asserts a known answer). Three `oracle-*` recipes need `statsrust` on PATH and fail without it (expected). Known limitation: the lm "power" parameter is unidentifiable. Run as `lmscript examples/x.lms` with the shim on PATH (never `lmscript run`).
- Docs: `docs/parity.md`, `docs/oracle-findings.md`. R is GPL but is only spawned as a process; never read or translate R/CRAN source.
- Repo is PRIVATE. Per the CI rule, private repos are verified locally; keep paid macOS/Windows legs gated with `github.event.repository.private == false`.

## Not done (deliberately)
- The re-engineering requested by the owner: one static binary with no Python runtime or library dependency. A separate team owns this. Start from `durable-memory/native-plan.md` (multicall Rust binary `lmnum`, invoked as `rcli` it serves this tool's byte-identical schema; self re-exec worker with rlimits + seccomp; docker image frozen as test oracle with tolerance tiers; ~38-45 engineer-weeks). Blockers before linking upstream code publicly: lmda Tukey GPL-derived code, lmstats Numerical Recipes routines, WMA/Erlang-A provenance.
- RNG decision: R-identical streams are NOT required; PCG64-style with a statistical guarantee is fine. Do not copy GPL R/Rmath source (clean-room; scipy/statsmodels/sklearn are BSD and may be consulted with notices kept).

## Rules that apply
- Commit identity Mike Ma <mike.ma@lmctl.com>, one Co-Authored-By per contributing model, explicit `git add`, push after every commit. Staging/deployment belongs to another team.
