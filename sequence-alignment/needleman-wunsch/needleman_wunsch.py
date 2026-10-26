#!/usr/bin/env python3
"""
Needleman-Wunsch Global Alignment (Python)

Aligns two sequences end-to-end using dynamic programming.

Usage: python needleman_wunsch.py
"""

import time


def needleman_wunsch(a: str, b: str, match=1, mismatch=-1, gap=-2):
    """Full Needleman-Wunsch with traceback. O(nm) time and space."""
    n, m = len(a), len(b)

    # Build DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i * gap
    for j in range(1, m + 1):
        dp[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if a[i-1] == b[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

    # Traceback
    aligned_a, aligned_b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match if a[i-1] == b[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + score:
                aligned_a.append(a[i-1])
                aligned_b.append(b[j-1])
                i -= 1; j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            aligned_a.append(a[i-1])
            aligned_b.append('-')
            i -= 1
        else:
            aligned_a.append('-')
            aligned_b.append(b[j-1])
            j -= 1

    return dp[n][m], ''.join(reversed(aligned_a)), ''.join(reversed(aligned_b))


def main():
    print("=== Needleman-Wunsch Global Alignment (Python) ===\n")

    human = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDL"
    gorilla = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFDSFGDL"

    score, aln_a, aln_b = needleman_wunsch(human, gorilla, match=2, mismatch=-1, gap=-3)
    match_line = ''.join('|' if a == b else '.' if a != '-' and b != '-' else ' '
                         for a, b in zip(aln_a, aln_b))
    print(f"  {aln_a}")
    print(f"  {match_line}")
    print(f"  {aln_b}")
    print(f"Score: {score}\n")

    # Benchmark
    print("=== Benchmark ===")
    for n in [100, 500, 1000]:
        a = "ACGT" * (n // 4)
        b = "CGTA" * (n // 4)
        start = time.perf_counter()
        s, _, _ = needleman_wunsch(a, b)
        elapsed = time.perf_counter() - start
        print(f"  {n}x{n}: score={s}, time={elapsed:.3f}s")

    print("\n=== Why Rust Is Faster ===")
    print("For 1000x1000: Python ~2-5s, Rust ~5ms (400-1000x faster)")
    print("  1. Flat Vec<i32> (4 bytes/cell) vs list-of-lists (28+ bytes/cell)")
    print("  2. Compiled to native code with LLVM optimizations")
    print("  3. No GC pauses; deterministic memory management")
    print("  4. max() on i32 compiles to cmov (branchless)")


if __name__ == "__main__":
    main()
