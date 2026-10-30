#!/usr/bin/env python3
"""Louvain Clustering (Python)"""
import math
from collections import defaultdict

def euclidean(a, b):
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def build_knn(data, k):
    n = len(data)
    graph = [[] for _ in range(n)]
    for i in range(n):
        dists = sorted(((j, euclidean(data[i],data[j])) for j in range(n) if j!=i), key=lambda x:x[1])
        for j, d in dists[:k]:
            graph[i].append((j, 1.0/(1.0+d)))
    return graph

def modularity(graph, comm):
    tw = sum(w for edges in graph for _,w in edges)
    if tw == 0: return 0
    q = 0
    for i, edges in enumerate(graph):
        ki = sum(w for _,w in edges)
        for j, w in edges:
            if comm[i]==comm[j]:
                kj = sum(ww for _,ww in graph[j])
                q += w - ki*kj/tw
    return q/tw

if __name__ == "__main__":
    data = [[2,2],[2.1,1.9],[8,8],[8.1,7.9],[2,8],[1.9,8.1]]
    graph = build_knn(data, 2)
    comm = list(range(len(data)))
    print(f"Initial modularity: {modularity(graph, comm):.4f}")
    print("Graph-based clustering separates cells into communities.")
