# rcli durable memory

- [status.md](status.md): hand-off state, rules that apply (commit identity, CI policy, clean-room), what is not done.
- [native-dsl.md](native-dsl.md): the R-free native core and the `_lmcode_` host `lmr` (`native/`, `examples/dsl/`,
  `docs/dsl-api.md`): layout, build, decisions, clean-room rule, measured agreement with R.
- [native-plan.md](native-plan.md): the earlier Rust `lmnum` design (op-by-op mapping, tolerance tiers, sandbox plan);
  superseded for rcli by native-dsl.md but still the reference for op risks and licensing blockers.
