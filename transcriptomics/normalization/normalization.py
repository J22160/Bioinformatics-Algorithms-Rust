#!/usr/bin/env python3
"""RNA-seq Normalization (Python)"""
import math
from statistics import median

def median_of_ratios(counts):
    n_genes, n_samples = len(counts), len(counts[0])
    geo_means = []
    for i in range(n_genes):
        if all(c > 0 for c in counts[i]):
            lm = sum(math.log(c) for c in counts[i]) / n_samples
            geo_means.append((i, math.exp(lm)))
    sf = []
    for j in range(n_samples):
        ratios = [counts[i][j] / gm for i, gm in geo_means]
        sf.append(median(ratios))
    return sf

if __name__ == "__main__":
    counts = [[100, 200, 150], [500, 1000, 750], [50, 100, 75]]
    print("Size factors:", median_of_ratios(counts))
