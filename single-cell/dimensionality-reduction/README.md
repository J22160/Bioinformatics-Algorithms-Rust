# Dimensionality Reduction: PCA and UMAP

## What Is This?

A single-cell count matrix might have 30,000 genes (dimensions)
measured across 10,000 cells. Most of these dimensions are noise.
Dimensionality reduction projects cells into a much smaller
space (50 dimensions for PCA, 2 for UMAP) while preserving
the meaningful biological variation.

## PCA (Principal Component Analysis)

PCA finds the directions of maximum variance in the data.
It is a linear method — fast and interpretable but cannot
capture nonlinear relationships.

In scRNA-seq, PCA typically reduces from ~2000 variable genes
to 50 principal components. These 50 PCs are then used as
input to UMAP and clustering [10].

## UMAP (Uniform Manifold Approximation and Projection)

UMAP is a nonlinear method that preserves both local and
global structure [10][13]:

1. Build a weighted k-nearest-neighbor graph in high-dimensional space
2. Find a low-dimensional (2D) embedding that best preserves
   the graph structure
3. Optimize the embedding using stochastic gradient descent

UMAP produces the 2D plots you see in every scRNA-seq paper.

## Complexity
- **PCA**: O(n × d × k) for n cells, d genes, k components
- **UMAP**: O(n × k × iterations) for graph construction + optimization

## Connection to Textbook

PCA uses the same eigenvalue decomposition concepts that appear
in spectral clustering (Chapter 10) [1]. UMAP's k-nearest-neighbor
graph construction is related to the neighborhood graphs used
in the clustering algorithms of Chapter 10 [1]. The optimization
of the UMAP embedding uses techniques related to the gradient
descent methods that underpin the Baum-Welch algorithm for
HMM parameter estimation in Chapter 11 [1].

t-SNE, another popular visualization method, achieved the best
overall accuracy in benchmark comparisons of dimensionality
reduction methods, while UMAP exhibited the highest stability
and better global structure preservation [13].
