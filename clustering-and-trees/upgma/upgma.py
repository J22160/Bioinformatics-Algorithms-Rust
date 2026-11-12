#!/usr/bin/env python3
"""UPGMA — Ultrametric Phylogenetic Tree Construction (Python)"""

def upgma(dist, labels):
    n = len(labels)
    clusters = {i: [i] for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    d = {(i,j): dist[i][j] for i in range(n) for j in range(n)}
    active = set(range(n))
    next_id = n
    tree = {}

    for _ in range(n - 1):
        # Find closest pair
        min_d, mi, mj = float('inf'), -1, -1
        for i in active:
            for j in active:
                if i < j and d.get((i,j), float('inf')) < min_d:
                    min_d, mi, mj = d[(i,j)], i, j

        height = min_d / 2.0
        tree[next_id] = (mi, mj, height)

        # Recompute distances
        si, sj = len(clusters[mi]), len(clusters[mj])
        for k in active:
            if k == mi or k == mj: continue
            nd = (si * d.get((min(mi,k),max(mi,k)), 0) +
                  sj * d.get((min(mj,k),max(mj,k)), 0)) / (si + sj)
            d[(min(next_id,k), max(next_id,k))] = nd

        clusters[next_id] = clusters[mi] + clusters[mj]
        heights[next_id] = height
        active.discard(mi)
        active.discard(mj)
        active.add(next_id)
        next_id += 1

    return tree, heights, labels

def print_newick(tree, heights, labels, node):
    if node < len(labels):
        return labels[node]
    left, right, h = tree[node]
    lh = heights.get(left, 0)
    rh = heights.get(right, 0)
    return f"({print_newick(tree,heights,labels,left)}:{h-lh:.4f},{print_newick(tree,heights,labels,right)}:{h-rh:.4f})"

if __name__ == "__main__":
    labels = ["E_coli", "S_enterica", "B_subtilis", "S_aureus", "M_tuberculosis"]
    dist = [
        [0.0, 0.8, 3.0, 3.2, 4.5],
        [0.8, 0.0, 2.8, 3.0, 4.3],
        [3.0, 2.8, 0.0, 1.5, 4.0],
        [3.2, 3.0, 1.5, 0.0, 3.8],
        [4.5, 4.3, 4.0, 3.8, 0.0],
    ]
    tree, heights, labels = upgma(dist, labels)
    root = max(tree.keys())
    print(f"Newick: {print_newick(tree, heights, labels, root)}")
    print(f"Root height: {heights[root]:.4f}")
