# Graph-Based Clustering (Louvain and Leiden)

## What Is Graph-Based Clustering?

After dimensionality reduction, such as Principal Component Analysis (PCA) or Uniform Manifold Approximation and Projection (UMAP), each cell in a single-cell RNA-seq (scRNA-seq) dataset is represented by a numerical vector describing its position in a lower-dimensional space.

The next step is to group cells with similar molecular profiles into clusters. Ideally, these clusters correspond to biologically meaningful cell populations, such as:

- **T cells**
- **B cells**
- **Monocytes**
- **Macrophages**
- **Endothelial cells**
- **Fibroblasts**

Traditional clustering algorithms, such as $k$-means, directly partition points in a feature space. However, modern single-cell analysis commonly uses **graph-based community detection**, particularly the Louvain and Leiden algorithms.

Instead of clustering cells directly based on their coordinates, graph-based methods first represent relationships between cells as a network. They then identify groups of cells that are more strongly connected to one another than to the rest of the network.

This approach is particularly useful for single-cell data because cell populations can have different sizes, complex shapes, and gradual transitions between biological states.

---

## 1. Building a Cell Similarity Graph

The first step is to construct a graph representing similarities between cells.

A graph is defined as:

$$G = (V,E)$$

where:
- $V$ represents the set of vertices (nodes), corresponding to individual cells.
- $E$ represents the edges connecting similar cells.
- Each edge can have a weight representing the strength of similarity between two cells.

For a dataset containing $N$ cells:

$$V = \{v_1, v_2, \ldots, v_N\}$$

Each cell becomes one node in the graph.

### Step 1: Calculate Cell-to-Cell Similarity

Before constructing the graph, cells are typically represented using a reduced-dimensional representation, often obtained through PCA.

For example, each cell might be represented by its coordinates across 30 principal components:

$$X_i = (PC_{i1}, PC_{i2}, \ldots, PC_{i30})$$

Similarity between cells can then be estimated using a distance metric. A commonly used metric is Euclidean distance:

$$d(i,j) = \sqrt{\sum_{p=1}^{D} (X_{ip} - X_{jp})^2}$$

where:
- $i, j$ represent two cells.
- $D$ is the number of dimensions.
- $X_{ip}$ is the coordinate of cell $i$ in dimension $p$.

Cells with smaller distances have more similar representations. Other distance or similarity measures can also be used, depending on the data and analysis workflow.

### Step 2: Construct a $k$-Nearest-Neighbor Graph

For each cell, identify its $k$ nearest neighbors based on the calculated distances. For example, if $k=5$, each cell is connected to its five most similar cells.

Consider a simplified dataset containing six cells:

| Cell | Five Nearest Neighbors (Illustrative) |
| :--- | :--- |
| **A** | B, C, D, E, F |
| **B** | A, C, D, E, F |
| **C** | A, B, D, E, F |
| **D** | A, B, C, E, F |
| **E** | A, B, C, D, F |
| **F** | A, B, C, D, E |

In a real dataset, nearest-neighbor relationships are not necessarily symmetric. Cell A may identify B as a neighbor without B identifying A as one of its own nearest neighbors.

Two common graph construction strategies are:
- **Directed KNN Graph:** An edge from cell A to cell B exists if B is among A's nearest neighbors.
- **Undirected KNN Graph:** Directed relationships are combined to create an undirected graph.

For undirected graphs, two common approaches are:
- **Union KNN:** Connect two cells if either cell identifies the other as a neighbor.
- **Mutual KNN:** Connect two cells only if they identify each other as neighbors.

The choice affects graph connectivity and can influence downstream clustering.

### Step 3: Assign Edge Weights

Not all neighboring cells are equally similar. Therefore, edges are often assigned weights representing the strength of the relationship between two cells.

One possible approach is to transform distances into similarities using a Gaussian kernel:

$$w_{ij} = \exp\left(-\frac{d(i,j)^2}{2\sigma^2}\right)$$

where:
- $d(i,j)$ is the distance between cells $i$ and $j$.
- $\sigma$ controls the distance scale.
- $w_{ij}$ is the similarity weight.

Smaller distances produce larger weights, indicating stronger relationships.

In practice, single-cell analysis tools often use locally scaled distances, shared-nearest-neighbor (SNN) information, or other graph-weighting strategies rather than this simple Gaussian transformation.

The resulting graph contains:
- **Nodes** representing cells.
- **Edges** representing cell-to-cell relationships.
- **Edge weights** representing the strength of those relationships.

> **Important:** The graph is typically constructed using PCA or another selected feature representation, not necessarily UMAP coordinates. UMAP is commonly used for visualization after clustering. Clustering directly on UMAP coordinates is generally not recommended because UMAP is designed primarily to preserve selected neighborhood relationships for visualization rather than exact global distances.

---

## 2. What Is Community Detection?

Once the cell similarity graph has been constructed, the next step is to identify communities.

A **community** is a group of nodes that have relatively strong connections to one another and fewer connections to nodes outside the group. In single-cell analysis, these communities are candidate cell populations.

For example, a graph might contain three densely connected regions:
- **Community 1:** Mostly T cells
- **Community 2:** Mostly B cells
- **Community 3:** Mostly Monocytes

The goal of community detection is to partition the graph into groups that capture this connectivity structure. One widely used objective for community detection is **modularity optimization**.

### Understanding Modularity

Modularity measures how much more strongly connected nodes are within the same community compared with what would be expected under a degree-preserving random graph model.

The modularity $Q$ of a partition is defined as:

$$Q = \frac{1}{2m} \sum_{i,j} \left( A_{ij} - \frac{k_i k_j}{2m} \right) \delta(c_i, c_j)$$

where:
- $A_{ij}$ is the edge weight between nodes $i$ and $j$.
- $k_i$ is the weighted degree (strength) of node $i$.
- $k_j$ is the weighted degree of node $j$.
- $m$ is the total edge weight, where $2m = \sum_{i,j} A_{ij}$ for an undirected graph.
- $c_i$ and $c_j$ are the community assignments of nodes $i$ and $j$.
- $\delta(c_i, c_j) = 1$ if both nodes belong to the same community, and $0$ otherwise.

The expression:

$$\frac{k_i k_j}{2m}$$

represents the expected connection strength between nodes $i$ and $j$ under the configuration-model null expectation.

The modularity equation therefore compares the observed connections with the connections expected under this null model.

**Interpretation:**
- **Positive modularity ($Q > 0$):** Indicates that the partition has more within-community connectivity than expected under the null model.
- **Zero modularity ($Q \approx 0$):** Indicates connectivity similar to the null expectation.
- **Negative modularity ($Q < 0$):** Indicates weaker within-community connectivity than expected.

The objective is to find a partition that maximizes modularity:

$$\max Q$$

Finding the globally optimal modularity partition is computationally NP-hard for large graphs. Louvain and Leiden use efficient heuristic optimization strategies to find high-quality partitions.

---

## 3. The Louvain Algorithm

The **Louvain algorithm** is a greedy, multilevel community detection algorithm introduced by Blondel et al. in 2008. It identifies communities by repeatedly moving nodes between communities to improve modularity, followed by graph aggregation.

The algorithm consists of two main alternating phases:

```
+-------------------------------------------------------+
|                 Initialize Communities                |
|           (Each node in its own community)            |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|            Step 1: Local Optimization                 |
|  - Pass through nodes sequentially                    |
|  - Evaluate modularity gain delta Q for neighbor moves|
|  - Move node to neighbor community if delta Q > 0     |
|  - Repeat until no further local gain                 |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|               Step 2: Graph Aggregation              |
|  - Collapse detected communities into super-nodes     |
|  - Sum edge weights between communities               |
|  - Convert internal community edges to self-loops     |
+-------------------------------------------------------+
                           |
         +-----------------+-----------------+
         |                                   |
         v                                   v
[ Modularity Improved ]            [ No Further Improvement ]
         |                                   |
         +--> Repeat Step 1                 +--> Terminate
```

### Step 1: Local Modularity Optimization

Initially, every node is assigned to its own community. For a graph containing $N$ cells:

$$C_i = \{v_i\}$$

The algorithm then visits nodes sequentially and evaluates whether moving a node into the community of one of its neighboring nodes improves modularity.

For example, consider a cell $A$ connected to cells $B$, $C$, and $D$. If $B$ and $C$ already belong to the same community, the algorithm evaluates whether moving $A$ into that community increases the modularity score. If the move improves modularity, it is accepted.

This process continues until no further local node movement produces an improvement. The algorithm primarily considers communities represented among a node's immediate neighbors, making the optimization scalable for sparse graphs.

### Step 2: Graph Aggregation

After local optimization, each detected community is collapsed into a single super-node.

**Before Aggregation:**
```text
A --- B
| \   |
C --- D       E --- F
```

Suppose the algorithm identifies:
- **Community 1:** A, B, C, D
- **Community 2:** E, F

**After Aggregation:**
```text
[ Super-Node 1 ] -------- [ Super-Node 2 ]
```

- Edges between original nodes in different communities become weighted edges between the corresponding super-nodes.
- Edges connecting nodes within the same community contribute to the self-loop weight of the super-node.

The aggregated graph is smaller than the original graph, allowing the algorithm to repeat local optimization at a higher level.

### Limitations of Louvain

Although Louvain is computationally efficient, it has an important limitation: **it can produce communities that are internally poorly connected or disconnected.**

For example, a community may contain two dense groups connected through only a narrow bridge. Louvain can assign both groups to the same community even when their internal connectivity is weak.

Louvain can also be sensitive to:
- Random node visitation order.
- Graph construction parameters.
- Resolution settings.
- Structural density variations in the input graph.

---

## 4. The Leiden Algorithm

The **Leiden algorithm**, introduced by Traag et al. in 2019, addresses the limitations of Louvain by incorporating a dedicated refinement phase that guarantees communities are well-connected.

```
+-------------------------------------------------------+
|                 Initialize Communities                |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|              Step 1: Local Node Movement             |
|  - Fast greedy local modularity optimization          |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|                 Step 2: Refinement                   |
|  - Identify well-connected subcommunities             |
|  - Allow splitting of poorly connected partitions     |
|  - Randomize node movement within subcommunities      |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|               Step 3: Graph Aggregation              |
|  - Aggregate based on REFINED subcommunities          |
|  - Ensures non-disconnected super-nodes               |
+-------------------------------------------------------+
                           |
                           v
         +-----------------+-----------------+
         |                                   |
[ Graph Changes ]                   [ Convergence Reached ]
         |                                   |
         +--> Repeat Process                 +--> Final Partition
```

### How Leiden Works

1. **Step 1: Local Moving:** Similar to Louvain, Leiden begins by moving nodes locally to optimize the quality function.
2. **Step 2: Refinement:** Rather than directly aggregating communities from Phase 1, Leiden refines them into smaller, well-connected subcommunities. It evaluates sub-partitions and can split communities if internal connectivity is insufficient.
3. **Step 3: Aggregation:** The refined communities—rather than the raw local communities—are aggregated into super-nodes. This guarantees that super-nodes remain well-connected at higher levels.

### Why Leiden Is Preferred

- **Connectivity Guarantees:** Ensures that all communities are internally connected.
- **Higher Partition Quality:** Consistently finds partitions with equal or higher modularity/quality values.
- **Fast Convergence:** Refinement avoids bad local minima, often leading to faster runtimes in practice.

---

## 5. Louvain vs. Leiden Comparison

| Feature | Louvain | Leiden |
| :--- | :--- | :--- |
| **Optimization Strategy** | Greedy, multilevel | Greedy, multilevel with refinement |
| **Local Node Movement** | Yes | Yes |
| **Refinement Phase** | No dedicated refinement | Yes (guarantees connectivity) |
| **Graph Aggregation** | Aggregates raw communities | Aggregates refined subcommunities |
| **Connectivity Guarantees** | Can yield disconnected/bridge-like clusters | Strong internal connectivity guarantees |
| **Computational Efficiency** | Efficient $\mathcal{O}(N \log N)$ average | Fast, often converges quicker than Louvain |
| **Biological Validation Needed** | Yes | Yes |

---

## 6. Understanding the Resolution Parameter

The **resolution parameter** ($\gamma$) controls the scale or granularity of the detected communities.

- **Lower Resolution ($\gamma < 1.0$):** Favors larger, broader clusters.
- **Higher Resolution ($\gamma > 1.0$):** Favors smaller, fine-grained subclusters.

A modularity objective incorporating resolution $\gamma$ is:

$$Q(\gamma) = \frac{1}{2m} \sum_{i,j} \left( A_{ij} - \gamma \frac{k_i k_j}{2m} \right) \delta(c_i, c_j)$$

When $\gamma = 1$, this equals standard modularity.

### Granularity Example (10,000 Immune Cells)

```text
Low Resolution (e.g., gamma = 0.2)
├── T Cells
├── B Cells
└── Monocytes

High Resolution (e.g., gamma = 1.2)
├── CD4+ T Cells
├── CD8+ T Cells
├── Regulatory T Cells (Tregs)
├── Naive B Cells
├── Memory B Cells
├── Classical Monocytes
└── Non-classical Monocytes
```

> **Note:** A higher resolution parameter does not automatically mean better biological fidelity. Oversplitting can partition homogenous populations due to random noise. Multiple resolutions should always be evaluated alongside marker gene profiles.

---

## 7. Worked Example: Single-Cell RNA-seq Workflow

Consider a dataset of 5,000 peripheral blood mononuclear cells (PBMCs) processed through standard preprocessing (Normalization, Highly Variable Gene Selection, PCA).

```
[ Raw Counts ] -> [ QC & Norm ] -> [ HVG Selection ] -> [ PCA (30 PCs) ]
                                                               |
                                                               v
                                                    [ Build KNN Graph (k=20) ]
                                                               |
                                                               v
                                                    [ Leiden Clustering ]
                                                               |
                                                               v
                                                    [ Marker Gene Annotation ]
```

### Step 1: Construct the KNN Graph
- $N = 5,000$ cells
- $k = 20$ nearest neighbors calculated on 30 Principal Components.

### Step 2: Apply Leiden Clustering
Running Leiden community detection yields 4 clusters:

| Cluster | Cell Count | Marker Genes | Biological Annotation |
| :---: | :---: | :--- | :--- |
| **0** | 1,800 | *CD3D*, *CD3E*, *IL7R* | T Cells |
| **1** | 1,200 | *MS4A1*, *CD79A*, *CD74* | B Cells |
| **2** | 1,400 | *LST1*, *S100A8*, *LYZ* | Monocytes |
| **3** | 600 | *NKG7*, *GNLY*, *PRF1* | Natural Killer (NK) Cells |

### Step 3: Validate Biological Identities
Annotations are confirmed by visualizing gene expression using dot plots, feature plots, or differential expression tests against canonical cell markers.

---

## 8. Why Graph-Based Clustering Over $k$-Means?

$k$-means minimizes the within-cluster sum of squared distances:

$$J = \sum_{i=1}^{N} \| X_i - \mu_{c_i} \|^2$$

where $X_i$ is the cell coordinate, $\mu_{c_i}$ is the cluster centroid, and $c_i$ is the cluster assignment.

| Feature / Problem | $k$-Means Clustering | Graph-Based Clustering |
| :--- | :--- | :--- |
| **Number of Clusters** | Requires pre-specifying $k$ | Determined automatically via graph topology & $\gamma$ |
| **Cluster Geometry** | Assumes spherical, equal-sized clusters | Handles arbitrary shapes, manifold geometries |
| **Density Variations** | Sensitive to density differences | Robust to complex density distributions |
| **Neighborhood Modeling** | Uses global centroid distances | Uses local graph connectivity |

---

## 9. Computational Complexity

Let:
- $N$ = Number of cells (nodes)
- $E$ = Number of edges
- $k$ = Number of nearest neighbors per cell
- $D$ = Dimensionality (e.g., number of PCs)
- $I$ = Number of optimization iterations

### Algorithmic Complexity Table

| Step | Naive Complexity | Optimized / Sparse Complexity | Memory Complexity |
| :--- | :--- | :--- | :--- |
| **Pairwise Distance** | $\mathcal{O}(N^2 D)$ | $\mathcal{O}(N \log N \cdot D)$ (approx. search) | $\mathcal{O}(N^2)$ dense $\rightarrow \mathcal{O}(Nk)$ sparse |
| **KNN Graph Build** | $\mathcal{O}(N^2)$ | $\mathcal{O}(Nk)$ | $\mathcal{O}(Nk)$ |
| **Louvain Clustering**| Unbounded | $\mathcal{O}(I \cdot (N + E))$ | $\mathcal{O}(N + E)$ |
| **Leiden Clustering** | Unbounded | $\mathcal{O}(I \cdot (N + E))$ | $\mathcal{O}(N + E)$ |

> **Takeaway:** For modern datasets ($N > 10^5$), constructing the KNN graph via approximate nearest neighbor search (e.g., Annoy, HNSW) dominates computation time rather than the clustering algorithms.

---

## 10. Practical Applications in Single-Cell Analysis

### Cell-Type Identification
Primary markers commonly used across major lineages:

- **T Cells:** *CD3D*, *CD3E*, *TRAC*
- **B Cells:** *MS4A1*, *CD79A*, *CD74*
- **Monocytes:** *LST1*, *LYZ*, *S100A8*
- **NK Cells:** *NKG7*, *GNLY*, *PRF1*
- **Endothelial Cells:** *PECAM1*, *VWF*, *EMCN*
- **Fibroblasts:** *COL1A1*, *COL1A2*, *DCN*

### Subclustering
Subclustering isolates specific broad groups (e.g., T cells) to re-run PCA, KNN construction, and Leiden clustering to resolve finer subtypes (CD4+ Naive, CD8+ Cytotoxic, Tregs).

### Spatial Transcriptomics
In spatial platforms (e.g., 10x Visium, Stereo-seq), graph construction can combine expression vectors with physical spatial coordinates $(x, y)$ to identify spatial domain niches.

---

## 11. Common Tools and Code Implementations

### Python (Scanpy)

```python
import scanpy as sc

# Assume adata has been preprocessed (normalized, PCA computed)

# 1. Construct KNN graph
sc.pp.neighbors(
    adata,
    n_neighbors=20,
    n_pcs=30
)

# 2. Perform Leiden clustering
sc.tl.leiden(
    adata,
    resolution=0.8,
    key_added="leiden"
)

# 3. Compute UMAP and visualize
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data"
)
```

### R (Seurat)

```r
library(Seurat)

# Assume seurat_obj has completed RunPCA

# 1. Construct Shared Nearest Neighbor (SNN) graph
seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = 1:30,
    k.param = 20
)

# 2. Perform Leiden clustering (algorithm = 4)
seurat_obj <- FindClusters(
    seurat_obj,
    resolution = 0.8,
    algorithm = 4
)

# 3. Visualize on UMAP
DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE
)
```

---

## 12. Important Considerations & Best Practices

1. **Choice of $k$ (Neighbors):** Small $k$ leads to fragmented graphs and over-clustering; large $k$ smooths over rare biological populations.
2. **Resolution Tuning:** Test multiple resolutions (e.g., $0.2 \to 1.5$) and evaluate stability via trees (e.g., `clustree` in R).
3. **Batch Effects:** Graph construction on uncorrected data groups cells by technical batch rather than biology. Perform integration (e.g., Harmony, BBKNN, Scanorama) prior to graph building.
4. **Continuous Trajectories:** Differentiation pathways or cell cycles represent continuous states. Partitioning continuous graphs into discrete clusters can be artificial; supplement with trajectory inference (e.g., PAGA, Slingshot, Monocle).

---

## 13. Related Computer Science Concepts

- **Graph Theory:** Adjacency matrix representations, sparse graphs, graph partitioning.
- **Greedy Optimization:** Heuristic decision-making at local nodes without full search-space evaluation.
- **Multilevel Graph Partitioning:** Coarsening networks iteratively to solve optimization problems across spatial scales.
- **Combinatorial Optimization:** Solving NP-hard partition problems on graphs using heuristic approximations.

## 14. Key Takeaways

- Graph-based clustering represents cells as nodes in a similarity network.
- KNN graph construction connects each cell to its most similar neighbors.
- Edge weights encode the strength of cell-to-cell relationships.
- Community detection identifies groups of densely interconnected cells.
- Modularity measures within-community connectivity relative to a null-model expectation.
- Louvain uses local modularity optimization followed by graph aggregation.
- Leiden adds a refinement phase that improves community connectivity.
- The resolution parameter controls the granularity of the detected communities.
- Graph-based clustering is widely used for cell-type identification, subclustering, and cell-state exploration.
- Clusters require biological validation; community detection alone does not establish cell identity.

## References

1. Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*.  
   https://doi.org/10.1088/1742-5468/2008/10/P10008

2. Traag VA, Waltman L, van Eck NJ. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*.  
   https://doi.org/10.1038/s41598-019-41695-z

3. Waltman L, van Eck NJ. (2013). A smart local moving algorithm for large-scale modularity-based community detection. *The European Physical Journal B*.  
   https://doi.org/10.1140/epjb/e2013-40829-0

4. Wolf FA, Angerer P, Theis FJ. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*.  
   https://doi.org/10.1186/s13059-017-1382-0

5. Stuart T, Butler A, Hoffman P, et al. (2019). Comprehensive integration of single-cell data. *Cell*.  
   https://doi.org/10.1016/j.cell.2019.05.031

6. McInnes L, Healy J, Melville J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.  
   https://arxiv.org/abs/1802.03426