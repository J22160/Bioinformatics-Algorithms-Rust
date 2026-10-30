#!/usr/bin/env python3
"""UMI Deduplication — Directional Adjacency Method (Python)"""

from collections import deque

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def directional_dedup(umis):
    """umis: list of (sequence, count) tuples. Returns list of clusters."""
    n = len(umis)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if hamming(umis[i][0], umis[j][0]) == 1:
                if umis[i][1] >= umis[j][1]: adj[i].append(j)
                if umis[j][1] >= umis[i][1]: adj[j].append(i)

    visited = [False]*n
    order = sorted(range(n), key=lambda i: -umis[i][1])
    clusters = []
    for seed in order:
        if visited[seed]: continue
        cluster, queue = [], deque([seed])
        visited[seed] = True
        while queue:
            cur = queue.popleft()
            cluster.append(cur)
            for nb in adj[cur]:
                if not visited[nb] and umis[nb][1] <= umis[cur][1]:
                    visited[nb] = True; queue.append(nb)
        clusters.append(cluster)
    return clusters

if __name__ == "__main__":
    umis = [("ACGTACGTAC",150),("ACGTACGTAT",3),("TGCATGCATG",200),("AACCGGTTAA",120)]
    clusters = directional_dedup(umis)
    print(f"Molecules: {len(clusters)}")
    for i, c in enumerate(clusters):
        print(f"  Mol {i+1}: {umis[c[0]][0]} ({len(c)} UMIs merged)")
