#!/usr/bin/env python3
"""K-Means Clustering (Python)"""
import math, random

def euclidean(a, b):
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def kmeans(data, k, max_iter=100, seed=42):
    random.seed(seed)
    n, d = len(data), len(data[0])
    centers = random.sample(data, k)
    assignments = [0]*n
    for it in range(max_iter):
        changed = False
        for i in range(n):
            best = min(range(k), key=lambda c: euclidean(data[i], centers[c]))
            if assignments[i] != best: changed = True; assignments[i] = best
        if not changed: break
        centers = [[0]*d for _ in range(k)]
        counts = [0]*k
        for i in range(n):
            c = assignments[i]; counts[c] += 1
            for j in range(d): centers[c][j] += data[i][j]
        for c in range(k):
            if counts[c] > 0:
                centers[c] = [v/counts[c] for v in centers[c]]
    wcss = sum(euclidean(data[i], centers[assignments[i]])**2 for i in range(n))
    return assignments, centers, it+1, wcss

if __name__ == "__main__":
    data = [[10,1,1],[9.5,1.5,0.5],[1,10,1],[0.5,9.5,1.5],[1,1,10],[1.5,0.5,9.5]]
    asgn, centers, iters, wcss = kmeans(data, 3)
    print(f"Converged in {iters} iterations, WCSS={wcss:.2f}")
    for i, a in enumerate(asgn): print(f"  Point {i}: cluster {a}")
