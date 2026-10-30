# Dimensionality Reduction: PCA and UMAP

## What Is Dimensionality Reduction?

Single-cell RNA sequencing (scRNA-seq) generates high-dimensional gene expression data. A typical dataset may contain thousands of cells and tens of thousands of genes.

For example, consider a dataset containing:

* **10,000 cells**
* **30,000 measured genes**
* **One expression value for every gene in every cell**

The resulting expression matrix $X$ has dimensions:

$$X \in \mathbb{R}^{10000 \times 30000}$$

Each row represents a cell, and each column represents a gene.

Although this matrix contains a large amount of information, not every gene contributes equally to meaningful biological variation. Some genes:

* Have very low expression across most cells.
* Show little variation between cells.
* Reflect technical noise or sequencing artifacts.
* Are strongly correlated with other genes.
* Capture biological differences such as cell identity, activation, or differentiation.

Working directly with all measured genes can be computationally expensive and may allow noisy or redundant features to influence downstream analyses.

**Dimensionality reduction** transforms the original high-dimensional data into a smaller number of informative features while attempting to preserve important patterns of variation or similarity.

Two widely used methods in single-cell analysis are:

1. **Principal Component Analysis (PCA):** A linear dimensionality reduction method that identifies directions of maximum variance.
2. **Uniform Manifold Approximation and Projection (UMAP):** A nonlinear manifold-learning method that constructs a low-dimensional representation using relationships between neighboring observations.

These methods are complementary rather than interchangeable. PCA is commonly used to generate a compact representation for downstream analysis, while UMAP is frequently used to visualize cell populations in two dimensions.

---

## 1. Why Is Dimensionality Reduction Necessary?

Consider a dataset containing 2,000 highly variable genes across 10,000 cells. Each cell can be represented as a point in a 2,000-dimensional space:

$$X_i = (x_{i1}, x_{i2}, \ldots, x_{i2000})$$

where $x_{ij}$ represents the expression of gene $j$ in cell $i$.

Comparing cells directly in this high-dimensional space presents several key challenges:


```

+-----------------------------------------------------------------------+
|                    Challenges in High-Dimensional Space               |
+-----------------------------------------------------------------------+
|  1. High Computational Cost: Distance metrics scale as O(N^2 * D)     |
|  2. Redundant Features: Pathway-correlated genes inflate variance     |
|  3. Curse of Dimensionality: Distance metrics lose discriminating power|
|  4. Noise Overfitting: Non-informative genes hide biological signal   |
+-----------------------------------------------------------------------+

```

### 1.1 High Computational Cost

Calculating pairwise distances between every pair of cells requires comparing their expression profiles across thousands of genes.

For $N$ cells and $D$ features, a straightforward all-pairs distance calculation has complexity:

$$\mathcal{O}(N^2 D)$$

For large datasets, this can become computationally expensive. Reducing the number of dimensions makes downstream neighbor search and clustering significantly more efficient.

### 1.2 Redundant Information

Genes are not statistically independent. For example, genes involved in the same biological pathway may be strongly correlated and show similar expression patterns. Including many correlated genes effectively gives certain biological signals disproportionate influence in distance calculations. Dimensionality reduction combines correlated features into a smaller set of composite components.

### 1.3 The Curse of Dimensionality

In high-dimensional spaces, observations become increasingly sparse. As dimensionality increases:

* Distances between observations become less informative.
* The difference between nearest and farthest neighbors becomes less pronounced.
* Estimating local relationships becomes more difficult.
* Statistical and computational requirements increase exponentially.

Reducing the feature space restores the effectiveness of distance metrics for identifying cell-to-cell relationships.

### 1.4 Biological Interpretation

Gene expression variation is often driven by a smaller number of underlying biological processes. For example, variation across a PBMC dataset may reflect:

* T-cell versus B-cell identity
* Myeloid versus lymphoid lineage
* Immune activation status
* Cell-cycle activity
* Interferon signaling
* Differentiation states

Dimensionality reduction attempts to represent important patterns in the original expression matrix using fewer numerical features. However, reduced dimensions do not automatically correspond directly to individual biological pathways or single cell types.

---

## 2. Principal Component Analysis (PCA)

### What Is PCA?

**Principal Component Analysis (PCA)** is a linear dimensionality reduction method that transforms correlated variables into a smaller set of uncorrelated variables called **principal components (PCs)**.

Each principal component is a linear combination of the original features. PCA identifies orthogonal directions in the data that explain the greatest amount of variance.

For single-cell RNA-seq:

$$2000 \text{ genes} \longrightarrow 50 \text{ principal components}$$

Each cell is represented by 50 numerical values (PC scores) instead of 2,000 gene expression measurements.

### 2.1 Understanding Variance

PCA is based on the principle that directions with greater variance capture primary structural patterns in the data.

Consider a simplified dataset with two genes: **Gene A** and **Gene B**. Suppose their expression levels co-vary across cells:

```text
Gene B
  |
  |                 * *   (High Variance Axis -> PC1)
  |              * *
  |           * *
  |        * *
  |     * *
  |  * *
  |________________________ Gene A

```

Most of the variation occurs along the diagonal. PCA identifies this diagonal as the **first principal component (PC1)**. The perpendicular direction captures the remaining variation and becomes the **second principal component (PC2)**:

```text
             PC2 (Orthogonal / Low Variance)
              |
              |
        * *   |
      * *     |
    * *       |
  * *         |
 /________________________ PC1

```

* **PC1:** Captures the maximum possible variance along a single direction.
* **PC2:** Captures the maximum remaining variance subject to being orthogonal (uncorrelated) to PC1.
* **Subsequent PCs:** Follow the same orthogonal maximization principle.

### 2.2 Mathematical Foundation of PCA

Let the expression matrix be $X \in \mathbb{R}^{N \times D}$, where $N$ is the number of cells and $D$ is the number of genes.

1. **Mean Center the Features:**

$$X_c = X - \bar{X}$$

where $\bar{X}$ contains the gene-wise means.

2. **Compute Covariance Matrix:**

$$C = \frac{1}{N-1} X_c^T X_c$$

The covariance matrix $C \in \mathbb{R}^{D \times D}$ describes how features vary together.

3. **Eigendecomposition:**

PCA solves for the eigenvectors $v_i$ and eigenvalues $\lambda_i$ of $C$:

$$C v_i = \lambda_i v_i$$

where:

* $v_i$ is an eigenvector representing a principal component direction.
* $\lambda_i$ is the eigenvalue representing the variance captured along $v_i$.

Eigenvectors are ordered by decreasing eigenvalue:

$$\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_D \ge 0$$

4. **Project Data onto Top $K$ Components:**

Construct the loading matrix $V_K = [v_1, v_2, \ldots, v_K] \in \mathbb{R}^{D \times K}$. The projected scores matrix $Z$ is:

$$Z = X_c V_K, \quad Z \in \mathbb{R}^{N \times K}$$

Each row of $Z$ represents a cell in the reduced $K$-dimensional space.

### 2.3 PCA Loadings and Scores

#### PCA Loadings

Loadings represent the weights assigned to original features in constructing each component:

| Gene | PC1 Loading | Interpretation Example |
| --- | --- | --- |
| **CD3D** | `+0.42` | Strong positive contribution (T-cell program) |
| **CD3E** | `+0.39` | Positive contribution |
| **TRAC** | `+0.36` | Positive contribution |
| **MS4A1** | `-0.28` | Negative contribution (B-cell program) |
| **LYZ** | `-0.31` | Negative contribution (Monocyte program) |

A larger absolute loading value indicates stronger feature contribution. The sign indicates opposing directions of expression variation along the component.

#### PCA Scores

Scores represent the coordinate values of each cell in the transformed space:

| Cell ID | PC1 | PC2 | PC3 |
| --- | --- | --- | --- |
| **Cell 1** | `+2.4` | `+0.8` | `-0.3` |
| **Cell 2** | `+2.1` | `+0.6` | `+0.1` |
| **Cell 3** | `-1.8` | `+1.2` | `+0.5` |
| **Cell 4** | `-2.2` | `-0.4` | `+0.2` |

Cells with similar scores across components share similar global expression profiles and group together in graph construction.

### 2.4 Explained Variance

The **Explained Variance Ratio (EVR)** for component $i$ is defined as:

$$\text{EVR}_i = \frac{\lambda_i}{\sum_{j=1}^{D} \lambda_j}$$

The cumulative variance explained by the first $K$ components is:

$$\text{Cumulative EVR}(K) = \frac{\sum_{i=1}^{K} \lambda_i}{\sum_{j=1}^{D} \lambda_j}$$

**Example Variance Distribution:**

| Component | Explained Variance | Cumulative Variance |
| --- | --- | --- |
| **PC1** | 18% | 18% |
| **PC2** | 12% | 30% |
| **PC3** | 8% | 38% |
| **PC4** | 5% | 43% |
| **PC5** | 4% | 47% |

An **elbow plot** (scree plot) displays explained variance per component. In single-cell workflows, retaining 20–50 PCs is standard, as low-variance components can still carry critical signal for small or subtle cell populations.

### 2.5 PCA in Single-Cell Preprocessing Pipelines

PCA is not applied directly to raw gene count matrices. A standard single-cell preprocessing workflow follows this structure:

```text
  [ Raw Count Matrix ]
           |
           v
  [ Quality Control & Filtering ]
           |
           v
  [ Library-Size Normalization ]
           |
           v
  [ Log Transformation (log1p) ]
           |
           v
  [ Highly Variable Gene (HVG) Selection ] (e.g., top 2,000 genes)
           |
           v
  [ Feature Scaling / Z-score Transformation ]
           |
           v
  [ Principal Component Analysis ]
           |
           v
  [ PC Scores Matrix (N x K) ]

```

#### Why Select Highly Variable Genes (HVGs)?

Selecting 2,000–3,000 HVGs eliminates non-informative genes dominated by technical noise or baseline house-keeping functions, optimizing both signal-to-noise ratio and computational performance.

#### Why Scale the Data?

Scaling centers feature expression to a mean of 0 and variance of 1:

$$X'_{ij} = \frac{X_{ij} - \mu_j}{\sigma_j}$$

This prevents highly expressed genes with massive raw variances from single-handedly dominating all top principal components.

### 2.6 Choosing the Number of PCs

To select an appropriate number of principal components ($K$), evaluate:

* Elbow plots (scree plots) showing inflection points in variance explained.
* PC loading heatmaps to check if components encode clear gene programs vs. random noise.
* Downstream cluster stability across different PC thresholds.

---

## 3. Uniform Manifold Approximation and Projection (UMAP)

### What Is UMAP?

**Uniform Manifold Approximation and Projection (UMAP)** is a nonlinear dimensionality reduction technique based on manifold learning and Riemannian geometry.

In single-cell RNA-seq, UMAP transforms high-dimensional representations (e.g., 30–50 PC scores) into 2 (or 3) dimensions for visualization:

$$Z \in \mathbb{R}^{N \times 50} \longrightarrow Y \in \mathbb{R}^{N \times 2}$$

Where each cell $i$ gets a 2D embedding coordinate:

$$Y_i = (\text{UMAP1}_i, \text{UMAP2}_i)$$

### 3.1 Why Use UMAP?

While linear methods like PCA preserve overall global structure, they cannot effectively project complex, nonlinear manifold geometries into 2D spaces without overlapping distinct groups. UMAP generates visually clear 2D embeddings that preserve local neighborhood relationships.

```
High-Dimensional PC Space                 2D UMAP Embedding
  (Complex Manifold)                       (Visually Resolved)

   .  :  *  .  :                            +------------------+
 :  * T-Cells *  :                          | [T-Cells]        |
  *  .  :  .  *                             |         [B-Cells]|
    --- Bridge ---                          | [NK]             |
 :  * Monocytes * :                         |      [Monocytes] |
   .  :  *  .  :                            +------------------+

```

### 3.2 How UMAP Works: Stage-by-Stage

UMAP executes in three core stages:

```
+-----------------------------------------------------------------------+
|  Stage 1: High-Dimensional Graph Construction                          |
|  - Find k-nearest neighbors in PC space                               |
|  - Compute directional local edge weights p_(j|i) using exponential   |
|    distance decay relative to local distance scale sigma_i           |
+-----------------------------------------------------------------------+
                                   |
                                   v
+-----------------------------------------------------------------------+
|  Stage 2: Fuzzy Set Union                                             |
|  - Combine asymmetric edge weights: p_ij = p_(j|i) + p_(i|j) - p_(j|i)p_(i|j) |
|  - Ensures global graph connectivity across varying densities        |
+-----------------------------------------------------------------------+
                                   |
                                   v
+-----------------------------------------------------------------------+
|  Stage 3: Low-Dimensional Optimization                                |
|  - Initialize 2D coordinates (e.g., Spectral Initialization)          |
|  - Compute low-dimensional similarity weights q_ij                   |
|  - Minimize Fuzzy Set Cross-Entropy via Stochastic Gradient Descent   |
+-----------------------------------------------------------------------+

```

#### Step 1: High-Dimensional Neighborhood Graph

For each observation $i$, find its $k$-nearest neighbors in PCA space using Euclidean or Cosine metrics.

#### Step 2: Local Distance Scaling

UMAP adapts to variable local data densities by calibrating per-node scales $\sigma_i$ and local connectivity offsets $\rho_i$.

The asymmetric strength of connection from node $i$ to node $j$ is calculated as:

$$p_{j\vert{}i} = \begin{cases} 1, & \text{if } d(i,j) \le \rho_i \\ \exp\left(-\frac{d(i,j) - \rho_i}{\sigma_i}\right), & \text{if } d(i,j) > \rho_i \end{cases}$$

where $\rho_i$ is the distance from node $i$ to its immediate nearest neighbor.

Directional weights are combined into a symmetric fuzzy union:

$$p_{ij} = p_{j\vert{}i} + p_{i\vert{}j} - p_{j\vert{}i} p_{i\vert{}j}$$

#### Step 3: Low-Dimensional Similarity Representation

In the 2D target space, the similarity weight $q_{ij}$ between low-dimensional coordinates $y_i$ and $y_j$ is modeled using a Student-$t$-like distribution curve:

$$q_{ij} = \frac{1}{1 + a \Vert{}y_i - y_j\Vert{}^{2b}}$$

where $a$ and $b$ are user/default curve-fitting parameters tied to `min_dist` and `spread`.

#### Step 4: Optimization via Fuzzy Set Cross-Entropy

UMAP matches $q_{ij}$ to $p_{ij}$ by minimizing the fuzzy set cross-entropy loss function $C$:

$$C = \sum_{i,j} \left[ p_{ij} \log\left(\frac{p_{ij}}{q_{ij}}\right) + (1 - p_{ij}) \log\left(\frac{1 - p_{ij}}{1 - q_{ij}}\right) \right]$$

* The left term **$p_{ij} \log(p_{ij}/q_{ij})$** pulls similar points together (attractive force).
* The right term **$(1-p_{ij}) \log((1-p_{ij})/(1-q_{ij}))$** pushes dissimilar points apart (repulsive force).

Minimization is executed efficiently via Stochastic Gradient Descent (SGD) with negative sampling.

---

## 4. Key UMAP Hyperparameters

Adjusting UMAP hyperparameters changes the structure and appearance of the final visualization:

| Hyperparameter | Typical Range | Effect of Smaller Values | Effect of Larger Values |
| --- | --- | --- | --- |
| **`n_neighbors`** | `10` – `50` | Preserves fine local structure; can fragment continuous clusters. | Captures broader global topology; may merge subtle sub-populations. |
| **`min_dist`** | `0.001` – `0.5` | Packs points tightly into dense clusters; sharpens boundaries. | Spreads points out; preserves smoother visual transitions between states. |
| **`n_components`** | `2` or `3` | Standard 2D output for plotting/visualization. | 3D interactive projections or embedding inputs for downstream tasks. |
| **`metric`** | `euclidean`, `cosine` | Measures absolute geometric distances (standard for PC scores). | Measures directional angular separation. |
| **`spread`** | `1.0` (default) | Decreases spatial range of low-dimensional clusters. | Increases overall structural scale of points. |
| **`random_state`** | Integer seed | Ensures reproducible stochastic optimization outputs. | Non-deterministic layouts across independent runs. |

---

## 5. Method Comparisons

### 5.1 PCA vs. UMAP

| Feature | Principal Component Analysis (PCA) | Uniform Manifold Approximation (UMAP) |
| --- | --- | --- |
| **Linearly Nature** | Strictly Linear | Nonlinear Manifold Learning |
| **Primary Objective** | Maximize total explained variance | Preserve fuzzy nearest-neighbor relationships |
| **Mathematical Engine** | Eigendecomposition / SVD | Fuzzy set theory & SGD cross-entropy minimization |
| **Input Data** | Scaled/Normalized expression matrix | High-dimensional PC scores matrix |
| **Output Data** | PC score matrix ($N \times K$) | Embedding coordinates ($N \times 2$) |
| **Typical Dimension Count** | 20 to 50 dimensions | 2 (or 3) dimensions |
| **Feature Interpretability** | High (loadings quantify feature weights) | None (embedding coordinate axes have no biological unit) |
| **Distance Preservation** | Preserves global Euclidean distance | Preserves local neighborhood graph distances |
| **Primary Single-Cell Role** | Noise reduction & input for KNN graphs | 2D Data visualization |

### 5.2 UMAP vs. t-SNE

| Feature | UMAP | t-SNE (t-Distributed Stochastic Neighbor Embedding) |
| --- | --- | --- |
| **Optimization Target** | Fuzzy set cross-entropy | Kullback-Leibler (KL) divergence |
| **KL Loss Formula** | $-$ | $KL(P \Vert{} Q) = \sum_{i,j} p_{ij} \log \left( \frac{p_{ij}}{q_{ij}} \right)$ |
| **Global Structure** | Better retention of inter-cluster relationships | Primarily preserves local neighborhoods; distorts inter-cluster distances |
| **Computational Speed** | Fast (scales linearly $\mathcal{O}(N)$ with approx. neighbors) | Can be slow on large datasets without FFT-accelerated approximations |
| **Cluster Density Representation** | Relatively stable | Clusters inflate to fill space regardless of underlying variance |
| **New Point Projection** | Can transform new query points into existing embeddings | Requires re-running entire optimization |

---

## 6. Worked Example: PBMC Single-Cell Workflow

Consider analyzing a dataset of 10,000 Peripheral Blood Mononuclear Cells (PBMCs) containing T cells, B cells, Monocytes, and NK cells.

```
Raw Expression Data: 10,000 cells x 30,000 genes
       |
       v
Normalize & Log-Transform -> Select 2,000 HVGs -> Scale Matrix
       |
       v
Compute PCA (50 components)
       |
       v
Construct KNN Graph (k=20, using Top 30 PCs)
       |
       +------------------------------------+
       |                                    |
       v                                    v
Perform Leiden Clustering           Compute UMAP (2D)
(Identify discrete biological groups)  (Visualize global layout)
       |                                    |
       +-----------------+------------------+
                         |
                         v
              Annotated 2D Visualization

```

### Execution Steps:

1. **Preprocessing:** Filter low-quality cells, normalize total UMI counts to $10,000$, perform $\log(1+x)$ transformation, and select top 2,000 HVGs.
2. **PCA Computation:** Compute top 50 principal components on scaled HVG data.
3. **Graph Construction:** Construct a $k$-nearest neighbor graph ($k=20$) based on the first 30 PC scores.
4. **Graph Partitioning:** Run Leiden community detection on the KNN graph to assign discrete cluster IDs.
5. **UMAP Visualization:** Generate 2D coordinates via UMAP from the same 30 PCs to visualize cell clusters.

```
       UMAP Visualization of Annotated PBMC Clusters
  +---------------------------------------------------+
  |                                                   |
  |    (Cluster 0: T-Cells)     (Cluster 1: B-Cells)  |
  |       CD3D+/TRAC+               MS4A1+/CD79A+     |
  |         * * *                      # # #          |
  |        * * * *                    # # # #         |
  |                                                   |
  |    (Cluster 3: NK Cells)   (Cluster 2: Monocytes) |
  |        NKG7+/GNLY+               LYZ+/S100A8+     |
  |         @ @ @                      % % %          |
  |                                                   |
  +---------------------------------------------------+

```

---

## 7. Practical Implementation (Scanpy Code Example)

```python
import scanpy as sc

# 1. Load preprocessed single-cell AnnData object
adata = sc.datasets.pbmc3k_processed()

# 2. Normalize and identify highly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

# Filter matrix for highly variable genes
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# 3. Scale expression data to unit variance and zero mean
sc.pp.scale(adata_hvg, max_value=10)

# 4. Compute PCA
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

# Plot explained variance ratio per principal component
sc.pl.pca_variance_ratio(adata_hvg, log=True)

# 5. Build Neighborhood Graph using top 30 PCs
sc.pp.neighbors(adata_hvg, n_neighbors=20, n_pcs=30)

# 6. Run Leiden Graph Clustering
sc.tl.leiden(adata_hvg, resolution=0.8, key_added="leiden")

# 7. Compute UMAP Embedding
sc.tl.umap(adata_hvg, min_dist=0.3, spread=1.0, random_state=42)

# 8. Visualize UMAP embedding colored by Leiden Clusters
sc.pl.umap(adata_hvg, color=["leiden"], legend_loc="on data")

```

---

## 8. Computational Complexity & Scalability

Let:

* $N$ = Number of cells (observations)
* $D$ = Number of original features (genes)
* $K$ = Number of principal components retained
* $k$ = Number of nearest neighbors per cell
* $T$ = Number of optimization iterations/epochs

| Algorithm Phase | Method / Variant | Time Complexity | Memory Complexity | Notes |
| --- | --- | --- | --- | --- |
| **Exact PCA** | Full SVD | $\mathcal{O}(N D^2 + D^3)$ | $\mathcal{O}(ND)$ | Impractical for huge scRNA-seq matrices. |
| **Truncated PCA** | Randomized / ARPACK SVD | $\mathcal{O}(N D K)$ | $\mathcal{O}(NK)$ | Standard single-cell implementation. |
| **KNN Graph** | Exact Brute-Force | $\mathcal{O}(N^2 K)$ | $\mathcal{O}(N^2)$ | Poor scaling for big data ($N > 100,000$). |
| **KNN Graph** | Approx. Nearest Neighbor (HNSW/Annoy) | $\mathcal{O}(N \log N \cdot K)$ | $\mathcal{O}(Nk)$ | Highly scalable sparse graph construction. |
| **UMAP Optimization** | SGD + Negative Sampling | $\mathcal{O}(T \cdot N \cdot k)$ | $\mathcal{O}(Nk)$ | Scales linearly with sample size $N$. |

---

## 9. Limitations & Best Practices

```
+-----------------------------------------------------------------------+
|                       CRITICAL PITFALL WARNINGS                       |
+-----------------------------------------------------------------------+
|  1. DO NOT cluster on 2D UMAP coordinates!                            |
|     Clustering MUST be performed on high-dimensional PC scores / KNN  |
|     graphs to avoid embedding compression artifacts.                  |
|                                                                       |
|  2. DO NOT over-interpret inter-cluster distances on UMAP!            |
|     Distances between far-apart islands on UMAP plots do NOT          |
|     linearly reflect real biological distance or evolutionary scale.  |
|                                                                       |
|  3. Inspect batch effects on PC1/2 prior to downstream tasks!         |
|     If technical batch dominates top PCs, apply integration tools     |
|     (e.g., Harmony, scVI, BBKNN) before constructing KNN graphs.      |
+-----------------------------------------------------------------------+

```

### Key Analytical Rules:

* **PCA Identifies Statistical Variance, Not Direct Biological Meaning:** High-variance components can reflect technical noise, cell cycle state, or batch differences rather than cell type identity. Inspect PC loadings to confirm biological drivers.
* **Batch Correction:** Uncorrected batch effects dominate top PCs, causing cells to cluster by sequencing run instead of biological identity. Perform batch integration prior to KNN graph construction and UMAP projection.
* **UMAP Axes Have No Absolute Scale:** Axis coordinates ($\text{UMAP1}, \text{UMAP2}$) are arbitrary embedding parameters without physical units.

---

## 10. References

1. **Jolliffe IT, Cadima J.** (2016). Principal component analysis: a review and recent developments. *Philosophical Transactions of the Royal Society A*. [https://doi.org/10.1098/rsta.2015.0202](https://doi.org/10.1098/rsta.2015.0202)
2. **McInnes L, Healy J, Melville J.** (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv preprint*. [https://arxiv.org/abs/1802.03426](https://arxiv.org/abs/1802.03426)
3. **van der Maaten L, Hinton G.** (2008). Visualizing Data using t-SNE. *Journal of Machine Learning Research*, 9:2579-2605. [https://www.jmlr.org/papers/v9/vandermaaten08a.html](https://www.jmlr.org/papers/v9/vandermaaten08a.html)
4. **Wolf FA, Angerer P, Theis FJ.** (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19:15. [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
5. **Stuart T, Butler A, Hoffman P, et al.** (2019). Comprehensive Integration of Single-Cell Data. *Cell*, 177(7):1888-1902. [https://doi.org/10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
6. **Kobak D, Linderman GC.** (2021). Initialization is critical for preserving global data structure in both t-SNE and UMAP. *Nature Biotechnology*, 39:156-157. [https://doi.org/10.1038/s41587-020-00809-z](https://doi.org/10.1038/s41587-020-00809-z)
7. **Becht E, McInnes L, Healy J, et al.** (2019). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, 37:38-44. [https://doi.org/10.1038/nbt.4314](https://doi.org/10.1038/nbt.4314)
8. **Luecken MD, Theis FJ.** (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6):e8746. [https://doi.org/10.15252/msb.20188746](https://doi.org/10.15252/msb.20188746)