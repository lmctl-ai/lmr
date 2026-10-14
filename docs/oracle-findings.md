# lmstats cross-tool oracle findings

Runs used lmstats `statsrust` built at lmstats commit
`3e17ca02fac4a28723d71f9e6f47d9f0c626b2a2` and host R 4.3.3 through rcli.

No comparison exceeded its documented tolerance. The measured differences were:

| Operation / quantity | Input | lmstats | R | Absolute difference | Tolerance / diagnosis |
|---|---|---:|---:|---:|---|
| OLS intercept | x=1..8, y=[2.1,4.2,5.8,8.3,9.9,12.2,13.8,16.1] | 0.114285714286 | 0.114285714286 | 0 | 1e-9; agreement |
| OLS slope | same | 1.985714285714 | 1.985714285714 | 0 | 1e-9; agreement |
| ACF lags 0..8 | [3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4] | see recipe output | see recipe output | < 1e-12 each | 1e-9; same biased autocovariance normalization |
| AR(1) phi | 120 values from rcli `rng_arima_sim`, seed 20261014, AR=.6 | 0.558542944600 | 0.558542833630 | 0.000000110970 | 1e-3; small optimizer/numerical difference between R ML and statsmodels state-space MLE |
| AR(1) sigma2 | same | 0.903489269531 | 0.903489228244 | 0.000000041287 | 1e-3; same variance quantity, numerical fit difference |

The ARIMA differences are below tolerance and are not classified as a discrepancy. The ACF
values match because both calls use the biased (divisor n) autocorrelation estimate. This file
should be extended with input, both values, difference, and a diagnosis if a future run exceeds a
recipe tolerance; tolerances must not be raised to hide a mismatch.
