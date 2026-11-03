#!/usr/bin/env python3
"""Bayesian Genotyping with Substitution Matrix (Python)"""
import math

BASES = "ACGT"

def bayesian_genotype(counts, ref_idx=0, prior_het=0.001):
    sub = [[0.995,0.002,0.001,0.002],[0.001,0.995,0.002,0.002],
           [0.002,0.002,0.995,0.001],[0.002,0.001,0.002,0.995]]

    results = []
    for a1 in range(4):
        for a2 in range(a1, 4):
            f = [0.0]*4
            if a1 == a2: f[a1] = 1.0
            else: f[a1], f[a2] = 0.5, 0.5

            ll = sum(counts[o] * math.log(max(sum(f[b]*sub[b][o] for b in range(4)), 1e-300))
                     for o in range(4) if counts[o] > 0)

            prior = (1 - 3*prior_het) if a1 == a2 and a1 == ref_idx else \
                    prior_het*0.05 if a1 == a2 else prior_het
            results.append((a1, a2, ll + math.log(max(prior, 1e-15))))

    max_lp = max(r[2] for r in results)
    log_sum = math.log(sum(math.exp(r[2]-max_lp) for r in results)) + max_lp
    for a1, a2, lp in results:
        post = math.exp(lp - log_sum)
        if post > 0.01:
            print(f"  {BASES[a1]}/{BASES[a2]}: {post:.4f}")

if __name__ == "__main__":
    bayesian_genotype([15, 0, 12, 1])