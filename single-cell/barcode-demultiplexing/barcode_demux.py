#!/usr/bin/env python3
"""Barcode Demultiplexing for scRNA-seq (Python)"""

from collections import Counter
import time

def match_barcode(observed, whitelist_set, qualities=None):
    if observed in whitelist_set:
        return ("exact", observed)
    bases = "ACGT"
    candidates = []
    for pos in range(len(observed)):
        for alt in bases:
            if alt == observed[pos]:
                continue
            variant = observed[:pos] + alt + observed[pos+1:]
            if variant in whitelist_set:
                candidates.append((variant, pos))
    if len(candidates) == 1:
        return ("corrected", candidates[0][0])
    elif len(candidates) > 1:
        return ("ambiguous", observed)
    return ("nomatch", observed)

if __name__ == "__main__":
    whitelist = {"ACGTACGTACGTACGT", "TGCATGCATGCATGCA", "GGCCAATTGGCCAATT"}
    reads = ["ACGTACGTACGTACGT", "ACGTACGTTCGTACGT", "AAAAAAAAAAAAAAAA"]
    for r in reads:
        status, bc = match_barcode(r, whitelist)
        print(f"  {r} -> {status}: {bc}")
