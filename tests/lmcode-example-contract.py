#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[1]
examples=(
    'examples/linear-fit.lms',
    'examples/normal-quantile.lms',
    'examples/rosenbrock-minimum.lms',
    'examples/dsl/install-and-run.lms',
)
for relative in examples:
    source=(root/relative).read_text()
    assert 'run_lmcode(' in source, f'{relative} must invoke through the lmcode manifest helper'
    for retired in ('savedsl(', 'dsl_cli('):
        assert retired not in source, f'{relative} retains unsupported helper {retired}'
print('PASS lmr examples use the lmcode invocation manifest')
