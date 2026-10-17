#!/usr/bin/env python3
"""Regenerate the tutorial CSV fixtures under examples/data (seeded, original synthetic data).

tut-cars.csv   : 60 rows, mpg = 42 - 5.0*weight - 0.015*power + N(0, 1.2^2); weight ~ U(1.5, 4.5) tonnes,
                 power = 60 + 40*weight + N(0, 15^2) hp (correlated predictors). Seed 2718.
tut-counts.csv : 80 rows, dose on a 0..4 grid (20 doses, 4 replicates), count ~ Poisson(exp(0.3 + 0.5*dose)).
                 Seed 2719 (inverse-CDF sampling, so only random.random is used).

Usage: python3 scripts/generate-tutorial-data.py   (from the repository root; rewrites both files)
"""
import math
import pathlib
import random

ROOT = pathlib.Path(__file__).resolve().parents[1] / "examples" / "data"


def cars(seed=2718, n=60):
    rng = random.Random(seed)
    rows = ["mpg,weight,power"]
    for _ in range(n):
        weight = rng.uniform(1.5, 4.5)
        power = 60 + 40 * weight + rng.gauss(0, 15)
        mpg = 42 - 5.0 * weight - 0.015 * power + rng.gauss(0, 1.2)
        rows.append(f"{mpg:.6f},{weight:.6f},{power:.6f}")
    return "\n".join(rows) + "\n"


def poisson(rng, lam):
    # Inverse-CDF sampling; exact for the small means used here.
    u = rng.random()
    k, p, cdf = 0, math.exp(-lam), math.exp(-lam)
    while u > cdf:
        k += 1
        p *= lam / k
        cdf += p
    return k


def counts(seed=2719):
    rng = random.Random(seed)
    rows = ["dose,count"]
    for i in range(20):
        dose = i * 4 / 19
        for _ in range(4):
            rows.append(f"{dose:.6f},{poisson(rng, math.exp(0.3 + 0.5 * dose))}")
    return "\n".join(rows) + "\n"


if __name__ == "__main__":
    (ROOT / "tut-cars.csv").write_text(cars())
    (ROOT / "tut-counts.csv").write_text(counts())
    print("wrote", ROOT / "tut-cars.csv", "and", ROOT / "tut-counts.csv")
