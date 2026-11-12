#!/usr/bin/env python3
"""Neighbor Joining (Python)"""

def neighbor_joining(dist, labels):
    import copy
    n = len(labels)
    d = [row[:] for row in dist]
    active = list(range(n))
    all_labels = list(labels)
    edges = []
    next_id = n

    while len(active) > 2:
        m = len(active)
        u = {}
        for i in active:
            u[i] = sum(d[i][j] for j in active if j != i) / max(m - 2, 1)

        best, bi, bj = float('inf'), 0, 0
        for a in range(m):
            for b in range(a+1, m):
                i, j = active[a], active[b]
                val = d[i][j] - u[i] - u[j]
                if val < best:
                    best, bi, bj = val, i, j

        dij = d[bi][bj]
        li = 0.5*dij + 0.5*(u[bi] - u[bj])
        lj = 0.5*dij + 0.5*(u[bj] - u[bi])

        nid = next_id; next_id += 1
        all_labels.append(f"Int_{nid}")

        # Expand matrix
        for row in d: row.extend([0.0] * (nid + 1 - len(row)))
        while len(d) <= nid: d.append([0.0] * (nid + 1))

        for k in active:
            if k in (bi, bj): continue
            nd = (d[bi][k] + d[bj][k] - dij) / 2.0
            d[nid][k] = nd; d[k][nid] = nd

        edges.append((nid, bi, li))
        edges.append((nid, bj, lj))
        active = [x for x in active if x not in (bi, bj)] + [nid]

    if len(active) == 2:
        edges.append((active[0], active[1], d[active[0]][active[1]]))

    return edges, all_labels

if __name__ == "__main__":
    labels = ["Human", "Chimp", "Gorilla", "Mouse", "Rat", "Dog"]
    dist = [[0,1.2,1.8,6,6.2,5],[1.2,0,1.6,5.8,6,4.8],[1.8,1.6,0,6.2,6.4,5.2],
            [6,5.8,6.2,0,1,5.5],[6.2,6,6.4,1,0,5.7],[5,4.8,5.2,5.5,5.7,0]]
    edges, all_labels = neighbor_joining(dist, labels)
    for f, t, l in edges:
        print(f"  {all_labels[f]} -- {all_labels[t]}: {l:.4f}")
