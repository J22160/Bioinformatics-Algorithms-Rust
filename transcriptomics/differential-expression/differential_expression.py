#!/usr/bin/env python3
"""Differential Expression — Negative Binomial Model (Python)"""

import math
from statistics import median


def estimate_size_factors(counts):
    n_genes, n_samples = len(counts), len(counts[0])
    geo_means = []
    for i in range(n_genes):
        if all(c > 0 for c in counts[i]):
            lm = sum(math.log(c) for c in counts[i]) / n_samples
            geo_means.append((i, math.exp(lm)))

    size_factors = []
    for j in range(n_samples):
        ratios = [counts[i][j] / gm for i, gm in geo_means if gm > 0]
        size_factors.append(median(ratios) if ratios else 1.0)
    return size_factors


def wald_test(counts, size_factors, group_a, group_b):
    results = []
    for gene in counts:
        na, nb = len(group_a), len(group_b)
        norm_a = [gene[j] / size_factors[j] for j in group_a]
        norm_b = [gene[j] / size_factors[j] for j in group_b]
        ma = sum(norm_a) / na
        mb = sum(norm_b) / nb
        lfc = math.log2(mb / ma) if ma > 0 and mb > 0 else 0
        va = sum((x - ma)**2 for x in norm_a) / max(na - 1, 1)
        vb = sum((x - mb)**2 for x in norm_b) / max(nb - 1, 1)
        se = math.sqrt(va / na + vb / nb) if (va + vb) > 0 else 1
        stat = (mb - ma) / se if se > 0 else 0
        p = 2 * (1 - 0.5 * (1 + math.erf(abs(stat) / math.sqrt(2))))
        results.append((lfc, p))
    return results


if __name__ == "__main__":
    counts = [
        [150, 130, 160, 450, 420, 480],
        [500, 520, 490, 510, 530, 505],
        [20, 25, 18, 200, 180, 220],
    ]
    names = ["BRCA1", "TP53", "IL6"]
    sf = estimate_size_factors(counts)
    results = wald_test(counts, sf, [0, 1, 2], [3, 4, 5])
    for name, (lfc, p) in zip(names, results):
        print(f"{name}: log2FC={lfc:.2f}, p={p:.2e}")