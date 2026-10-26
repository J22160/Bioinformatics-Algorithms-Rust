#!/usr/bin/env python3
"""Hirschberg's Linear-Space Alignment (Python)"""

def nw_last_row(a, b, match=1, mismatch=-1, gap=-2):
    m = len(b)
    prev = [j * gap for j in range(m + 1)]
    for i in range(1, len(a) + 1):
        curr = [i * gap] + [0] * m
        for j in range(1, m + 1):
            s = match if a[i-1] == b[j-1] else mismatch
            curr[j] = max(prev[j-1] + s, prev[j] + gap, curr[j-1] + gap)
        prev = curr
    return prev

def hirschberg(a, b, match=1, mismatch=-1, gap=-2):
    if len(a) == 0:
        return '-' * len(b), b
    if len(b) == 0:
        return a, '-' * len(a)
    if len(a) == 1:
        if a[0] in b:
            idx = b.index(a[0])
            return '-' * idx + a[0] + '-' * (len(b) - idx - 1), b
        return a + '-' * (len(b) - 1), '-' + b[:len(b)-1] if len(b) > 1 else b

    mid = len(a) // 2
    row_l = nw_last_row(a[:mid], b, match, mismatch, gap)
    row_r = nw_last_row(a[mid:][::-1], b[::-1], match, mismatch, gap)
    best_j = max(range(len(b) + 1), key=lambda j: row_l[j] + row_r[len(b) - j])

    a1, b1 = hirschberg(a[:mid], b[:best_j], match, mismatch, gap)
    a2, b2 = hirschberg(a[mid:], b[best_j:], match, mismatch, gap)
    return a1 + a2, b1 + b2

if __name__ == "__main__":
    a, b = hirschberg("AGTACGCATAG", "TATGCCATG")
    print(f"  {a}\n  {b}")
