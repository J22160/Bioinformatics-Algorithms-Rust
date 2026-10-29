#!/usr/bin/env python3
"""GSEA (Python)"""

import random

def enrichment_score(mask, scores):
    n = len(scores)
    n_hit = sum(mask)
    n_miss = n - n_hit
    nr = sum(abs(s) for m, s in zip(mask, scores) if m)
    if nr == 0 or n_miss == 0:
        return 0.0
    current = 0.0
    max_dev = 0.0
    for i in range(n):
        if mask[i]:
            current += abs(scores[i]) / nr
        else:
            current -= 1.0 / n_miss
        if abs(current) > abs(max_dev):
            max_dev = current
    return max_dev

if __name__ == "__main__":
    scores = [3.2, 2.8, 2.5, 2.1, 1.9, 0.1, -0.05, 0.02, -0.1, 0.0]
    mask = [True]*5 + [False]*5
    es = enrichment_score(mask, scores)
    print(f"Enrichment Score: {es:.4f}")
