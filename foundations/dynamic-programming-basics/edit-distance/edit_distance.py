#!/usr/bin/env python3
"""
Edit Distance (Levenshtein Distance) — CRISPR Off-Target Analysis (Python)

Usage: python edit_distance.py
"""

import time


def edit_distance(a: str, b: str) -> int:
    """Full O(n*m) dynamic programming edit distance."""
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,      # deletion
                dp[i][j - 1] + 1,      # insertion
                dp[i - 1][j - 1] + cost # substitution
            )

    return dp[n][m]


def edit_distance_linear_space(a: str, b: str) -> int:
    """O(min(n,m)) space version."""
    if len(a) < len(b):
        a, b = b, a
    n, m = len(a), len(b)

    prev = list(range(m + 1))
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        curr[0] = i
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev, curr = curr, prev

    return prev[m]


def main():
    print("=== Edit Distance: CRISPR Off-Target (Python) ===\n")

    guide = "GAGTCCGAGCAGAAGAAGAA"
    offtarget = "GAGTCCAAGCAGAAGATGAA"

    d = edit_distance(guide, offtarget)
    print(f"Guide:      {guide}")
    print(f"Off-target: {offtarget}")
    print(f"Edit dist:  {d}")

    if d <= 5:
        print(f"\n⚠  Potential off-target site (within 5 edits)")

    # Benchmark
    print("\n=== Benchmark ===")
    seq_a = "ACGT" * 250  # 1000 bp
    seq_b = "TGCA" * 250

    start = time.perf_counter()
    d = edit_distance_linear_space(seq_a, seq_b)
    elapsed = time.perf_counter() - start
    print(f"  1000x1000: distance={d}, time={elapsed:.3f}s")

    print()
    print("=== Why Rust Is Faster ===")
    print("For 1000x1000 matrix: Python ~1-2s, Rust ~2ms (500-1000x faster)")
    print("Reasons:")
    print("  1. Python list-of-lists: each cell is a 28-byte int object on the heap")
    print("     Rust flat Vec<u32>: each cell is exactly 4 bytes, contiguous in memory")
    print("  2. Python's min() creates temporary tuple objects each iteration")
    print("     Rust's .min() chain compiles to two cmov instructions")
    print("  3. Python has per-iteration bytecode dispatch overhead (~100ns)")
    print("     Rust inner loop runs at ~1-2ns per cell")
    print("  4. Sassy (a CRISPR off-target tool) uses SIMD-accelerated edit distance")
    print("     in Rust, achieving 2 Gbp/s throughput for pattern searching")


if __name__ == "__main__":
    main()