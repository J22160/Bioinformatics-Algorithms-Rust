# Single-Cell RNA-seq Algorithms

This section covers the core computational algorithms that power
modern single-cell RNA sequencing analysis pipelines such as
Cell Ranger, Scanpy, and Seurat.

## Why Single-Cell?

Bulk RNA-seq measures the average gene expression across millions
of cells. Single-cell RNA-seq (scRNA-seq) measures expression
in individual cells, revealing cellular heterogeneity invisible
to bulk methods.

## The scRNA-seq Analysis Pipeline

```
Raw Reads (FASTQ)
|
Barcode Demultiplexing <-- barcode-demultiplexing/
|
UMI Deduplication <-- umi-deduplication/
|
Count Matrix (genes x cells)
|
Dimensionality Reduction <-- dimensionality-reduction/
|
Clustering <-- clustering/
|
Cell Type Annotation
```


## Algorithms

| Algorithm | Purpose | Used By |
|-----------|---------|---------|
| Barcode Demultiplexing | Assign reads to cells | Cell Ranger, STARsolo |
| UMI Deduplication | Remove PCR duplicates | UMI-tools, Cell Ranger |
| Dimensionality Reduction (PCA/UMAP) | Reduce noise, visualize | Scanpy, Seurat |
| Graph Clustering (Leiden/Louvain) | Identify cell types | Scanpy, Seurat |


- **Barcode matching** uses approximate pattern matching 
- **UMI deduplication** uses graph algorithms and clustering 
- **Leiden clustering** uses community detection on graphs 
- **PCA** is a linear algebra foundation; **UMAP** uses graph
  construction similar to k-nearest-neighbor graphs