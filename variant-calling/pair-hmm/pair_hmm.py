#!/usr/bin/env python3
"""Pair HMM (Python) — Simplified read-haplotype likelihood."""
import math

def phred_to_prob(q): return 10**(-q/10)

def pair_hmm(read, quals, haplotype):
    r, h = len(read), len(haplotype)
    M = [[0.0]*(h+1) for _ in range(r+1)]
    I = [[0.0]*(h+1) for _ in range(r+1)]
    D = [[0.0]*(h+1) for _ in range(r+1)]
    for j in range(1, h+1): D[0][j] = 1.0/h
    gop, gep = 0.001, 0.1
    for i in range(1, r+1):
        e = phred_to_prob(quals[i-1])
        for j in range(1, h+1):
            p_emit = (1-e) if read[i-1]==haplotype[j-1] else e/3
            M[i][j] = p_emit*(M[i-1][j-1]*(1-2*gop)+I[i-1][j-1]*(1-gep)+D[i-1][j-1]*(1-gep))
            I[i][j] = M[i-1][j]*gop + I[i-1][j]*gep
            D[i][j] = M[i][j-1]*gop + D[i][j-1]*gep
    return sum(M[r][j]+I[r][j] for j in range(1,h+1))

if __name__ == "__main__":
    hap = "ACGTACGTACGTACGTACGT"
    read = "ACGTACGTACGT"
    lik = pair_hmm(read, [30]*len(read), hap)
    print(f"P(read|hap) = {lik:.6e}")