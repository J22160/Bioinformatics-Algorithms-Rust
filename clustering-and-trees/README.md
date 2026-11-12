# Clustering and Trees

This section covers algorithms for two deeply connected tasks:
grouping biological data into clusters and reconstructing
evolutionary trees from distance or character data.

## Why Clustering and Trees Matter

- **Gene expression analysis**: Clustering groups genes with
  similar expression patterns to reveal functional relationships
- **Phylogenetics**: Distance-based trees reveal evolutionary
  relationships among species from DNA or protein sequences
- **Single-cell analysis**: Clustering assigns cells to types;
  trees reveal differentiation trajectories
- **Metagenomics**: Clustering groups organisms by sequence similarity

## Algorithms

| Algorithm | Type | Purpose | Used By |
|-----------|------|---------|---------|
| K-Means | Partitional | Group data into k clusters | Gene expression analysis |
| Hierarchical (Agglomerative) | Hierarchical | Build cluster dendrograms | Heatmap generation |
| UPGMA | Distance-based tree | Ultrametric phylogenetic trees | Molecular clock studies |
| Neighbor Joining | Distance-based tree | Non-ultrametric trees | Phylogenetic inference |
| Small Parsimony (Fitch/Sankoff) | Character-based tree | Label internal nodes | Ancestral sequence reconstruction |

