#!/usr/bin/env python3
"""Four Russians Speedup — Conceptual Demo (Python)"""

def lcs_standard(a, b):
    n, m = len(a), len(b)
    prev = [0] * (m + 1)
    for i in range(1, n + 1):
        curr = [0] * (m + 1)
        for j in range(1, m + 1):
            if a[i-1] == b[j-1]:
                curr[j] = prev[j-1] + 1
            else:
                curr[j] = max(prev[j], curr[j-1])
        prev = curr
    return prev[m]

if __name__ == "__main__":
    a = "ACGTACGTACGTACGT"
    b = "AGCATCGATCGACGT"
    print(f"LCS length: {lcs_standard(a, b)}")
    print("Note: Myers 1999 bitpacking is the practical O(nm/w) approach.")
