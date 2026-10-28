#!/usr/bin/env python3
"""Suffix Arrays (Python)"""

def suffix_array(text):
    return sorted(range(len(text)), key=lambda i: text[i:])

def sa_search(text, sa, pattern):
    import bisect
    n, m = len(sa), len(pattern)
    lo = bisect.bisect_left(range(n), True,
        key=lambda i: text[sa[i]:sa[i]+m] >= pattern)
    hi = bisect.bisect_right(range(n), True,
        key=lambda i: text[sa[i]:sa[i]+m] <= pattern)
    return [sa[i] for i in range(lo, hi)]

if __name__ == "__main__":
    text = b"ATCGATCGATCG$"
    sa = suffix_array(text)
    print(f"SA: {sa}")
    # Simple search
    pattern = b"ATCG"
    positions = [sa[i] for i in range(len(sa))
                 if text[sa[i]:sa[i]+len(pattern)] == pattern]
    print(f"'{pattern.decode()}' at positions: {positions}")