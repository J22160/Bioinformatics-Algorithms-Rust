#!/usr/bin/env python3
"""Multiple Sequence Alignment — Star Method (Python)"""

def edit_distance(a, b):
    n, m = len(a), len(b)
    prev = list(range(m + 1))
    for i in range(1, n + 1):
        curr = [i] + [0] * m
        for j in range(1, m + 1):
            cost = 0 if a[i-1] == b[j-1] else 1
            curr[j] = min(prev[j-1] + cost, prev[j] + 1, curr[j-1] + 1)
        prev = curr
    return prev[m]

def star_msa(sequences):
    k = len(sequences)
    # Find center
    best_center, best_total = 0, float('inf')
    for i in range(k):
        total = sum(edit_distance(sequences[i], sequences[j])
                    for j in range(k) if j != i)
        if total < best_total:
            best_total, best_center = total, i

    print(f"Center sequence: S{best_center + 1}")
    return best_center

if __name__ == "__main__":
    seqs = ["ATGCTAGCTA", "ATGCAAGCTA", "ATGCTAGTA", "ATCTAGCTA", "ATGCTAGCTAG"]
    center = star_msa(seqs)
    print(f"Selected S{center + 1} as center (minimizes total edit distance)")
