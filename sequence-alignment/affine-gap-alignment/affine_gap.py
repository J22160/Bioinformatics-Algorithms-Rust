#!/usr/bin/env python3
"""Affine Gap Alignment (Python) — Gotoh's three-matrix approach."""

NEG_INF = float('-inf')

def affine_gap_align(a, b, match=2, mismatch=-3, gap_open=-5, gap_extend=-1):
    n, m = len(a), len(b)
    M = [[NEG_INF]*(m+1) for _ in range(n+1)]
    X = [[NEG_INF]*(m+1) for _ in range(n+1)]
    Y = [[NEG_INF]*(m+1) for _ in range(n+1)]

    M[0][0] = 0
    for i in range(1, n+1):
        Y[i][0] = gap_open + i * gap_extend
    for j in range(1, m+1):
        X[0][j] = gap_open + j * gap_extend

    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if a[i-1] == b[j-1] else mismatch
            M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + s
            X[i][j] = max(M[i][j-1] + gap_open + gap_extend, X[i][j-1] + gap_extend)
            Y[i][j] = max(M[i-1][j] + gap_open + gap_extend, Y[i-1][j] + gap_extend)

    return max(M[n][m], X[n][m], Y[n][m])

if __name__ == "__main__":
    score = affine_gap_align("ATGCTAGCTAGCATCGATCGATCGATCGATCGATCG", "ATGCTAGCTAGCATCGATCG")
    print(f"Score: {score}")
