#!/usr/bin/env python3
"""de Bruijn Graph Genome Assembly (Python)"""

from collections import defaultdict


class DeBruijnGraph:
    def __init__(self):
        self.adj = defaultdict(list)
        self.in_deg = defaultdict(int)
        self.out_deg = defaultdict(int)

    @classmethod
    def from_reads(cls, reads, k):
        graph = cls()
        for read in reads:
            for i in range(len(read) - k + 1):
                kmer = read[i:i+k]
                prefix, suffix = kmer[:-1], kmer[1:]
                graph.adj[prefix].append(suffix)
                graph.out_deg[prefix] += 1
                graph.in_deg[suffix] += 1
                if suffix not in graph.adj:
                    graph.adj[suffix] = []
        return graph

    def find_eulerian_path(self):
        # Find start
        start = None
        for n in self.adj:
            if self.out_deg.get(n, 0) - self.in_deg.get(n, 0) == 1:
                start = n
                break
        if start is None:
            start = next(iter(self.adj))

        adj_copy = {n: list(v) for n, v in self.adj.items()}
        stack, path = [start], []

        while stack:
            curr = stack[-1]
            if adj_copy.get(curr):
                stack.append(adj_copy[curr].pop())
            else:
                path.append(stack.pop())

        path.reverse()
        return path

    def reconstruct(self, path):
        return path[0] + "".join(n[-1] for n in path[1:])


if __name__ == "__main__":
    reads = ["AAGCTGCA", "GCTGCAAG", "TGCAAGCT"]
    print(f"Reads: {reads}")

    graph = DeBruijnGraph.from_reads(reads, 4)
    path = graph.find_eulerian_path()
    seq = graph.reconstruct(path)
    print(f"Assembled: {seq}")

    for read in reads:
        circular = seq + seq
        assert read in circular, f"Missing {read}"
    print("✓ All reads verified")
