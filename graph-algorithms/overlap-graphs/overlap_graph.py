#!/usr/bin/env python3
"""Overlap Graph (Python)"""


def overlap_length(a, b, min_ol=3):
    for k in range(min(len(a), len(b)), min_ol - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return 0


def build_overlap_graph(reads, min_ol=3):
    n = len(reads)
    return [[overlap_length(reads[i], reads[j], min_ol)
             if i != j else 0 for j in range(n)] for i in range(n)]


if __name__ == "__main__":
    reads = ["ATGCCA", "GCCATG", "CATGCC", "TGCCAT", "CCATGC"]
    graph = build_overlap_graph(reads)
    print("Overlap matrix:")
    for row in graph:
        print(row)
