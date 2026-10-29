#!/usr/bin/env python3
"""Pseudoalignment (Python)"""

from collections import defaultdict

class KmerIndex:
    def __init__(self, transcripts, k):
        self.k = k
        self.index = defaultdict(set)
        for tid, seq in enumerate(transcripts):
            for i in range(len(seq) - k + 1):
                self.index[seq[i:i+k]].add(tid)

    def pseudoalign(self, read):
        compat = None
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i+self.k]
            if kmer in self.index:
                tids = self.index[kmer]
                compat = tids if compat is None else compat & tids
                if not compat:
                    return None
        return sorted(compat) if compat else None

if __name__ == "__main__":
    transcripts = ["ATGCTAGCTAGCATCG", "ATGCTAGCTAGC", "CCGGAATTCCGG"]
    names = ["GAPDH-001", "GAPDH-002", "ACTB-001"]
    idx = KmerIndex(transcripts, 5)

    reads = ["CTAGCTAGCATCG", "GAATTCCGG"]
    for r in reads:
        ec = idx.pseudoalign(r)
        matched = [names[i] for i in ec] if ec else []
        print(f"  '{r}' -> {matched}")
