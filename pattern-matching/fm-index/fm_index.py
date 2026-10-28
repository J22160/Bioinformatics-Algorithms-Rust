#!/usr/bin/env python3
"""FM-Index (Python)"""

class FMIndex:
    def __init__(self, text: bytes):
        self.n = len(text)
        self.sa = sorted(range(self.n), key=lambda i: text[i:])
        self.bwt = bytes(text[i-1] if i > 0 else text[-1] for i in self.sa)

        # C array
        counts = [0] * 256
        for b in self.bwt:
            counts[b] += 1
        self.c = [0] * 256
        total = 0
        for i in range(256):
            self.c[i] = total
            total += counts[i]

        # Occurrence table
        self.occ = [[0]*256 for _ in range(self.n + 1)]
        for i in range(self.n):
            self.occ[i+1] = self.occ[i][:]
            self.occ[i+1][self.bwt[i]] += 1

    def count(self, pattern: bytes) -> int:
        lo, hi = 0, self.n
        for c in reversed(pattern):
            lo = self.c[c] + self.occ[lo][c]
            hi = self.c[c] + self.occ[hi][c]
            if lo >= hi:
                return 0
        return hi - lo

    def locate(self, pattern: bytes) -> list:
        lo, hi = 0, self.n
        for c in reversed(pattern):
            lo = self.c[c] + self.occ[lo][c]
            hi = self.c[c] + self.occ[hi][c]
            if lo >= hi:
                return []
        return sorted(self.sa[lo:hi])


if __name__ == "__main__":
    fm = FMIndex(b"ATGCTAGCTAGCATCGATCGATCG$")
    for pat in [b"ATCG", b"TAG", b"CGA"]:
        print(f"'{pat.decode()}': {fm.count(pat)} occurrences at {fm.locate(pat)}")
