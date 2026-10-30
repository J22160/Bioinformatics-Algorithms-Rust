#!/usr/bin/env python3
"""PCA for scRNA-seq Dimensionality Reduction (Python)"""
import math

def center(data):
    n, d = len(data), len(data[0])
    means = [sum(row[j] for row in data)/n for j in range(d)]
    return [[row[j]-means[j] for j in range(d)] for row in data]

def dot(a, b):
    return sum(x*y for x,y in zip(a,b))

def normalize(v):
    norm = math.sqrt(sum(x**2 for x in v))
    return [x/norm for x in v] if norm > 0 else v

def pca(data, k=2, max_iter=100):
    centered = center(data)
    n, d = len(centered), len(centered[0])
    components, variances = [], []
    residual = [row[:] for row in centered]

    for _ in range(k):
        pc = normalize([1.0]*d)
        for _ in range(max_iter):
            proj = [dot(row, pc) for row in residual]
            new_pc = [sum(residual[i][j]*proj[i] for i in range(n)) for j in range(d)]
            pc = normalize(new_pc)
        proj = [dot(row, pc) for row in residual]
        var = sum(p**2 for p in proj)/n
        variances.append(var)
        for i in range(n):
            for j in range(d):
                residual[i][j] -= proj[i]*pc[j]
        components.append(pc)
    return components, variances

if __name__ == "__main__":
    data = [[10,0.5,0.3],[9.8,0.3,0.4],[0.4,9.0,10.0],[0.5,9.5,9.8],[0.3,0.4,10.0],[0.4,0.5,9.8]]
    comps, vars = pca(data, 2)
    total = sum(vars)
    for i, v in enumerate(vars):
        print(f"PC{i+1}: {v/total*100:.1f}% variance")
