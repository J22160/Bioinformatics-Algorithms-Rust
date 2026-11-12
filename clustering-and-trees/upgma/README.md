# UPGMA — Ultrametric Tree Construction

## What Is This?

UPGMA builds an evolutionary tree from a distance matrix by
repeatedly merging the two closest clusters. It assumes a
**molecular clock** — that all species evolve at the same rate —
which means every leaf is equidistant from the root.

## How It Works

1. Start with each species as its own cluster
2. Find the two clusters with the smallest average distance
3. Merge them into a new cluster
4. Set the height of the new node to half the merge distance
5. Recompute distances from the new cluster to all others
   using the average of all pairwise distances
6. Repeat until one cluster remains

## Distance Recomputation

When clusters C1 and C2 merge into C, the distance to any
other cluster C* is computed as the weighted average [1]:

    D(C, C*) = (|C1|×D(C1,C*) + |C2|×D(C2,C*)) / (|C1|+|C2|)

## Complexity
- **Time**: O(n³) for n taxa
- **Space**: O(n²) for the distance matrix

## Bioinformatics Problem

**Scenario**: You have 16S rRNA sequences from 5 bacterial
isolates and have computed pairwise evolutionary distances.
Build a UPGMA tree to estimate when each pair of species
diverged, assuming a constant molecular clock.

## Limitations

UPGMA produces incorrect topologies when the molecular clock
assumption is violated — when different lineages evolve at
different rates. The neighbor joining algorithm handles
this case better [1].

UPGMA is arguably the most popular hierarchical clustering
algorithm, especially for gene expression data and protein
sequence clustering [12]. However, standard UPGMA requires
the entire dissimilarity matrix in memory, which limits
scalability to very large datasets [12]. Memory-constrained
variants like MC-UPGMA have been developed to address this [12].
