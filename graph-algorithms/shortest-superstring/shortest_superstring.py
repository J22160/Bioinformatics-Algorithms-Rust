#!/usr/bin/env python3
"""Shortest Superstring — Greedy (Python)"""


def overlap_len(a, b):
    for k in range(min(len(a), len(b)), 0, -1):
        if a.endswith(b[:k]):
            return k
    return 0


def greedy_superstring(strings):
    strings = list(strings)
    while len(strings) > 1:
        best_i, best_j, best_ol = 0, 1, -1
        for i in range(len(strings)):
            for j in range(len(strings)):
                if i != j:
                    ol = overlap_len(strings[i], strings[j])
                    if ol > best_ol:
                        best_i, best_j, best_ol = i, j, ol
        merged = strings[best_i] + strings[best_j][best_ol:]
        strings = [s for idx, s in enumerate(strings)
                   if idx not in (best_i, best_j)] + [merged]
    return strings[0]


if __name__ == "__main__":
    reads = ["ATGCCA", "GCCATG", "CATGCC", "TGCCAT"]
    result = greedy_superstring(reads)
    print(f"Superstring: {result} (len={len(result)})")
    for r in reads:
        assert r in result
    print("✓ All reads present")
