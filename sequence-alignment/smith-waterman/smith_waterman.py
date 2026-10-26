#!/usr/bin/env python3
"""Smith-Waterman Local Alignment (Python)"""

import time


def smith_waterman(a: str, b: str, match=2, mismatch=-1, gap=-2):
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    max_score, max_i, max_j = 0, 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if a[i-1] == b[j-1] else mismatch
            dp[i][j] = max(0,
                           dp[i-1][j-1] + score,
                           dp[i-1][j] + gap,
                           dp[i][j-1] + gap)
            if dp[i][j] > max_score:
                max_score, max_i, max_j = dp[i][j], i, j

    # Traceback
    aligned_a, aligned_b = [], []
    i, j = max_i, max_j
    while i > 0 and j > 0 and dp[i][j] > 0:
        score = match if a[i-1] == b[j-1] else mismatch
        if dp[i][j] == dp[i-1][j-1] + score:
            aligned_a.append(a[i-1]); aligned_b.append(b[j-1])
            i -= 1; j -= 1
        elif dp[i][j] == dp[i-1][j] + gap:
            aligned_a.append(a[i-1]); aligned_b.append('-')
            i -= 1
        else:
            aligned_a.append('-'); aligned_b.append(b[j-1])
            j -= 1

    return max_score, ''.join(reversed(aligned_a)), ''.join(reversed(aligned_b))


if __name__ == "__main__":
    contig = "TTGCAATGCAATCGATCGAATGGCATGCTAGCTAGCATCGATATCGCAATGCAATGCAAGC"
    gene = "GCATGCTAGCTAGCATCGA"
    score, aln_a, aln_b = smith_waterman(contig, gene)
    print(f"Score: {score}")
    print(f"  {aln_a}")
    print(f"  {aln_b}")
