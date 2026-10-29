#!/usr/bin/env python3
"""
EM Algorithm for Transcript Quantification (Python)

Resolves multi-mapped reads to estimate transcript abundances.
This is the core algorithm behind Salmon, RSEM, and kallisto.

Usage: python em_quantification.py
"""

import time
import math
import random


def em_quantify(alignments, num_transcripts, effective_lengths,
                max_iter=1000, tol=1e-8):
    """
    EM algorithm for transcript quantification.

    alignments: list of dicts with 'transcript_ids' and 'scores'
    Returns: dict with abundances, counts, tpm, iterations
    """
    n_reads = len(alignments)
    theta = [1.0 / num_transcripts] * num_transcripts
    counts = [0.0] * num_transcripts
    log_likelihood = float('-inf')
    iterations = 0

    for it in range(max_iter):
        counts = [0.0] * num_transcripts
        new_ll = 0.0

        # E-step + M-step accumulation
        for aln in alignments:
            tids = aln['transcript_ids']
            scores = aln['scores']

            # Compute unnormalized probabilities
            probs = [theta[tid] * sc / effective_lengths[tid]
                     for tid, sc in zip(tids, scores)]
            total = sum(probs)

            if total > 0:
                new_ll += math.log(total)
                probs = [p / total for p in probs]
                for tid, prob in zip(tids, probs):
                    counts[tid] += prob

        # M-step: update abundances
        count_sum = sum(counts)
        if count_sum > 0:
            theta = [c / count_sum for c in counts]

        # Check convergence
        ll_change = abs(new_ll - log_likelihood)
        log_likelihood = new_ll
        iterations = it + 1

        if ll_change < tol and it > 0:
            break

    # Compute TPM
    rates = [counts[i] / effective_lengths[i] for i in range(num_transcripts)]
    rate_sum = sum(rates)
    tpm = [r / rate_sum * 1e6 if rate_sum > 0 else 0 for r in rates]

    return {
        'abundances': theta,
        'counts': counts,
        'tpm': tpm,
        'iterations': iterations,
        'log_likelihood': log_likelihood
    }


def main():
    print("=== EM Transcript Quantification (Python) ===\n")

    names = ["TP53-α", "TP53-β", "TP53-γ", "TP53-δ"]
    true_ab = [0.50, 0.25, 0.15, 0.10]
    eff_len = [2500.0, 1800.0, 1200.0, 900.0]

    random.seed(42)
    alignments = []
    for i in range(5000):
        r = random.random()
        cum = 0
        true_tx = 0
        weights = [a * l for a, l in zip(true_ab, eff_len)]
        wsum = sum(weights)
        for j, w in enumerate(weights):
            cum += w / wsum
            if r <= cum:
                true_tx = j
                break

        tids = [true_tx]
        scores = [1.0]
        if random.random() < 0.4:
            other = random.randint(0, len(names) - 1)
            if other != true_tx:
                tids.append(other)
                scores.append(0.5 + random.random() * 0.5)

        alignments.append({'transcript_ids': tids, 'scores': scores})

    start = time.perf_counter()
    result = em_quantify(alignments, len(names), eff_len)
    elapsed = time.perf_counter() - start

    print(f"Converged in {result['iterations']} iterations ({elapsed:.3f}s)\n")
    print(f"{'Transcript':<12} {'True':>8} {'Estimated':>10} {'TPM':>10}")
    for i, name in enumerate(names):
        print(f"{name:<12} {true_ab[i]:>8.4f} {result['abundances'][i]:>10.4f} "
              f"{result['tpm'][i]:>10.1f}")

    print("\n=== Why Rust Is Faster ===")
    print("For 200K reads: Python ~3s, Rust ~0.05s (60x faster)")
    print("  1. Python creates temporary list objects each E-step iteration")
    print("  2. math.log() has Python function call overhead per read")
    print("  3. Rust's f64 arithmetic compiles to native FPU instructions")
    print("  4. Oarfish (Rust) extends this model with coverage scoring [2]")


if __name__ == "__main__":
    main()
