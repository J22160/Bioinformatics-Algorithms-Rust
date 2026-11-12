# K-Means Clustering

## Table of Contents

- [What Is K-Means Clustering?](#what-is-k-means-clustering)
- [1. How K-Means Works](#1-how-k-means-works)
- [2. Mathematical Objective Function](#2-mathematical-objective-function)
- [3. Worked Example](#3-worked-example)
- [4. K-Means for Gene Expression Analysis](#4-k-means-for-gene-expression-analysis)
- [5. Choosing the Optimal Number of Clusters](#5-choosing-the-optimal-number-of-clusters)
- [6. Initialization Strategies](#6-initialization-strategies)
- [7. Computational Complexity](#7-computational-complexity)
- [8. Advantages of K-Means](#8-advantages-of-k-means)
- [9. Limitations of K-Means](#9-limitations-of-k-means)
- [10. K-Means vs. Other Clustering Algorithms](#10-k-means-vs-other-clustering-algorithms)
- [11. K-Means and Deep Learning](#11-k-means-and-deep-learning)
- [12. Practical Applications in Bioinformatics](#12-practical-applications-in-bioinformatics)
- [13. Best Practices for Gene Expression Clustering](#13-best-practices-for-gene-expression-clustering)
- [14. Connection to Statistical and Computational Concepts](#14-connection-to-statistical-and-computational-concepts)
- [References and Further Reading](#references-and-further-reading)

---

## What Is K-Means Clustering?

**K-means** is an unsupervised machine-learning algorithm that partitions a dataset into a predefined number of clusters ($K$), where each data point is assigned to the cluster with the nearest centroid. The objective is to group data points so that observations within the same cluster are more similar to one another than to observations in different clusters.

### Common Applications in Bioinformatics

- Clustering genes based on expression patterns.
- Identifying groups of samples with similar molecular profiles.
- Exploring cell populations using gene expression measurements.
- Grouping proteins based on sequence-derived or structural features.
- Identifying patterns in metabolomics and proteomics datasets.
- Discovering potential biological subgroups in high-dimensional datasets.

For gene expression analysis, K-means can group genes that exhibit similar expression patterns across multiple experimental conditions. For example, genes that are highly expressed in the same tissue or respond similarly to a treatment may be assigned to the same cluster. These genes may participate in related biological processes or regulatory pathways, although shared expression patterns do not prove that they have the same function or are directly regulated by the same mechanism.

---

## 1. How K-Means Works

K-means uses an iterative optimization procedure called **Lloyd's algorithm**. The algorithm repeatedly performs two main operations until the assignments stabilize or a stopping criterion is reached:

1. **Assignment:** Assign each data point to its nearest cluster centroid.
2. **Update:** Update each centroid based on the points assigned to that cluster.

### Step 1: Initialize Cluster Centroids

Choose the number of clusters, $K$. The algorithm initializes $K$ centroids, which represent the centers of the clusters.

For example, if the dataset contains $200$ genes and the goal is to identify five expression groups:

$$K = 5$$

The initial centroids may be selected randomly from the observed data points or generated using another initialization strategy.

### Step 2: Assign Data Points to the Nearest Centroid

For each data point, calculate its distance to every centroid. The most commonly used distance measure is **Euclidean distance**:

$$d(x_i, \mu_j) = \sqrt{\sum_{l=1}^{D}(x_{il} - \mu_{jl})^2}$$

Where:
- $x_i$ = data point $i$
- $\mu_j$ = centroid of cluster $j$
- $D$ = number of features or dimensions
- $x_{il}$ = value of feature $l$ for data point $i$

Each data point is assigned to the cluster with the smallest distance using the following assignment rule:

$$c_i = \underset{j \in \{1, \ldots, K\}}{\operatorname{argmin}} \ d(x_i, \mu_j)$$

Here, $c_i$ represents the cluster assigned to data point $i$.

### Step 3: Update Cluster Centroids

After assigning all data points, calculate the new centroid for each cluster. The centroid is the arithmetic mean of all points assigned to that cluster:

$$\mu_j = \frac{1}{|C_j|} \sum_{x_i \in C_j} x_i$$

Where:
- $C_j$ = set of data points assigned to cluster $j$
- $|C_j|$ = number of points in cluster $j$
- $\mu_j$ = updated centroid

Each feature is averaged independently.

### Step 4: Repeat Until Convergence

Repeat the assignment and centroid-update steps until one of the following conditions is met:

- Cluster assignments no longer change.
- Centroid movement becomes smaller than a predefined tolerance.
- The maximum number of iterations is reached.
- The objective function stops improving meaningfully.

The final result is a set of $K$ clusters and their corresponding centroids.

> **Important:** K-means converges to a **local optimum** of its objective function, not necessarily the globally optimal clustering. Different initial centroid positions can produce different results.

---

## 2. Mathematical Objective Function

K-means minimizes the within-cluster sum of squared distances, also known as the **within-cluster sum of squares (WCSS)** or **inertia**.

The objective function is defined as:

$$J = \sum_{j=1}^{K} \sum_{x_i \in C_j} \|x_i - \mu_j\|^2$$

Where:
- $K$ = number of clusters
- $C_j$ = data points assigned to cluster $j$
- $x_i$ = individual data point
- $\mu_j$ = centroid of cluster $j$
- $\|x_i - \mu_j\|^2$ = squared Euclidean distance between the point and its centroid

The algorithm attempts to minimize $J$.

### Interpreting the Objective Function

A smaller WCSS indicates that data points are, on average, closer to their assigned centroids:
- **Low WCSS:** Indicates relatively compact clusters.
- **High WCSS:** Indicates greater variation within clusters.

However, WCSS always decreases or stays the same as $K$ increases because adding more clusters gives the algorithm additional flexibility. Therefore, a low WCSS alone does not establish that a clustering solution is biologically meaningful or that the selected number of clusters is appropriate.

---

## 3. Worked Example

Consider five genes measured across three experimental conditions. Suppose their normalized expression values are:

| Gene | Control | Treatment A | Treatment B |
| :--- | :---: | :---: | :---: |
| **Gene A** | 2 | 8 | 9 |
| **Gene B** | 3 | 7 | 8 |
| **Gene C** | 8 | 2 | 3 |
| **Gene D** | 9 | 3 | 2 |
| **Gene E** | 5 | 5 | 6 |

Assume we want to divide these genes into two clusters ($K = 2$). For illustration, choose the initial centroids:

$$\mu_1 = (2, 8, 9)$$

$$\mu_2 = (8, 2, 3)$$

### Assignment Step

Calculate the Euclidean distance between each gene and the two centroids. For **Gene B**:

$$d(B, \mu_1) = \sqrt{(3 - 2)^2 + (7 - 8)^2 + (8 - 9)^2} = \sqrt{3} \approx 1.73$$

$$d(B, \mu_2) = \sqrt{(3 - 8)^2 + (7 - 2)^2 + (8 - 3)^2} = \sqrt{75} \approx 8.66$$

Gene B is assigned to **Cluster 1** because its distance to $\mu_1$ is smaller ($1.73 < 8.66$). Repeating this calculation for every gene gives the following grouping:

- **Cluster 1:** Gene A, Gene B, Gene E
- **Cluster 2:** Gene C, Gene D

### Update Step

Recalculate the centroid of Cluster 1:

$$\mu_1 = \frac{(2, 8, 9) + (3, 7, 8) + (5, 5, 6)}{3} = (3.33, 6.67, 7.67)$$

Recalculate the centroid of Cluster 2:

$$\mu_2 = \frac{(8, 2, 3) + (9, 3, 2)}{2} = (8.5, 2.5, 2.5)$$

The algorithm repeats the assignment and update steps using these new centroids until convergence. In this example, the two groups represent genes with broadly similar expression patterns:

- **Cluster 1:** Higher expression in Treatment A and Treatment B relative to Control.
- **Cluster 2:** Higher expression in Control relative to the treatment conditions.

This provides an initial hypothesis about treatment-responsive gene groups.

> **Biological Interpretation:** Genes grouped together may share expression behavior, but additional functional enrichment, pathway analysis, and experimental validation are needed to determine whether they have related biological roles.

---

## 4. K-Means for Gene Expression Analysis

K-means is frequently used to identify groups of genes with similar expression patterns across experimental conditions.

### Bioinformatics Scenario

Suppose you have RNA-seq expression profiles for $200$ genes across $10$ experimental conditions, and you want to find $5$ desired clusters. The expression matrix is:

$$X \in \mathbb{R}^{200 \times 10}$$

Where:
- **Rows** represent genes.
- **Columns** represent experimental conditions.
- **Each entry** represents the expression level of a gene under a particular condition.

### Example Workflow

```text
RNA-seq raw read data
        |
        v
Read quality control
        |
        v
Read alignment / quantification
        |
        v
Gene-level count matrix
        |
        v
Normalization
        |
        v
Log transformation / variance stabilization
        |
        v
Feature selection
        |
        v
Gene expression matrix
        |
        v
K-means clustering
        |
        v
Cluster assignment
        |
        v
Expression heatmap and cluster profiles
        |
        v
Functional enrichment analysis
        |
        v
Biological interpretation
```

### Step 1: Normalize Gene Expression Data

Raw RNA-seq read counts should not generally be used directly for Euclidean-distance-based K-means clustering, as they are influenced by sequencing depth and differences in library composition.

A normalization method such as **DESeq2 size-factor normalization** or **edgeR TMM** can help account for these differences. A variance-stabilizing transformation (**VST**) or regularized log transformation (**rlog**) can also be useful for exploratory expression analysis.

```r
library(DESeq2)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = ~ condition
)

# Normalize and transform
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)

# Extract transformed expression matrix
expression_matrix <- assay(vsd)
```

### Step 2: Select Genes for Clustering

Clustering all genes may introduce substantial noise from genes with little or no variation across conditions. Highly variable genes or genes selected based on a specific biological question provide more informative clusters.

```r
# Calculate variance for each gene
gene_variance <- apply(expression_matrix, 1, var)

# Select the 200 most variable genes
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:200]
selected_matrix <- expression_matrix[top_genes, ]
```

The selection strategy should match the research question. For example, selecting genes based on differential-expression results can be appropriate when clustering treatment-responsive genes, though it may introduce selection bias if the same data are later used for independent statistical claims.

### Step 3: Scale Expression Profiles

Genes can have substantially different average expression levels and variances. If the goal is to group genes by the **shape** of their expression profiles rather than absolute expression magnitude, row-wise standardization (Z-score transformation) is often useful:

$$Z_{ij} = \frac{X_{ij} - \bar{X}_i}{s_i}$$

Where:
- $X_{ij}$ = expression of gene $i$ in condition $j$
- $\bar{X}_i$ = mean expression of gene $i$ across conditions
- $s_i$ = standard deviation of gene $i$ across conditions

This centers each gene around zero and scales it by its variability:

```r
# Standardize each gene across conditions
scaled_matrix <- t(scale(t(selected_matrix)))
```

After scaling:
- **Positive values** indicate expression above that gene's average.
- **Negative values** indicate expression below that gene's average.
- **Zero** indicates expression at that gene's average.

> **Important:** Row-wise scaling changes the clustering objective. It emphasizes similarity in expression patterns rather than absolute expression levels, and it can amplify noise in genes with very low variability. Consider filtering low-information genes before scaling.

### Step 4: Apply K-Means Clustering

Use the `kmeans()` function in R:

```r
set.seed(42)

kmeans_result <- kmeans(
    scaled_matrix,
    centers = 5,
    nstart = 25,
    iter.max = 100
)

# Cluster assignment for each gene
gene_clusters <- kmeans_result$cluster

# Cluster centroids
cluster_centers <- kmeans_result$centers

# Within-cluster sum of squares
wcss <- kmeans_result$tot.withinss
```

#### Key Parameters
- `centers = 5`: Specifies five clusters.
- `nstart = 25`: Runs the algorithm using $25$ different initializations and retains the best solution.
- `iter.max = 100`: Maximum number of iterations per initialization.
- `set.seed(42)`: Makes the random initialization reproducible.

The result object contains cluster membership for each gene, final cluster centroids, within-cluster sum of squares, between-cluster sum of squares, and cluster sizes.

### Step 5: Visualize the Clusters

A heatmap can reveal the expression patterns associated with each cluster:

```r
library(pheatmap)

annotation_row <- data.frame(
    Cluster = factor(gene_clusters)
)
rownames(annotation_row) <- rownames(scaled_matrix)

pheatmap(
    scaled_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    show_rownames = FALSE,
    main = "K-means Gene Expression Clusters"
)
```

You can also calculate average expression profiles for each cluster and plot them to compare cluster-level trends.

### Step 6: Perform Functional Enrichment Analysis

After identifying clusters, perform functional enrichment analysis to investigate their potential biological roles.

Common approaches include:
- Gene Ontology (GO) enrichment
- KEGG pathway enrichment
- Reactome pathway enrichment
- Gene set enrichment analysis (GSEA)

For example, individual clusters may be enriched for immune response, cell-cycle regulation, mitochondrial function, synaptic signaling, or inflammatory pathways, helping generate hypotheses about the biological processes represented.

---

## 5. Choosing the Optimal Number of Clusters

K-means requires the number of clusters, $K$, to be specified beforehand. When the optimal $K$ is unknown, several quantitative methods can guide selection.

### A. Elbow Method

The elbow method evaluates the within-cluster sum of squares across different values of $K$:

$$\text{WCSS}(K) = \sum_{j=1}^{K} \sum_{x_i \in C_j} \|x_i - \mu_j\|^2$$

As $K$ increases, WCSS decreases. The objective is to identify the "elbow"—the point where increasing $K$ yields diminishing returns in WCSS reduction.

```r
set.seed(42)
wcss <- numeric(10)

for (k in 1:10) {
    result <- kmeans(
        scaled_matrix,
        centers = k,
        nstart = 25
    )
    wcss[k] <- result$tot.withinss
}

plot(
    1:10,
    wcss,
    type = "b",
    xlab = "Number of clusters (K)",
    ylab = "Within-cluster sum of squares",
    main = "Elbow Method"
)
```

*Note: The elbow method is useful but subjective; some datasets lack a sharp inflection point.*

### B. Silhouette Score

The silhouette score measures how well each data point fits within its assigned cluster compared with neighboring clusters. For a data point $i$:

$$s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}$$

Where:
- $a(i)$ = average distance between point $i$ and other points in its own cluster
- $b(i)$ = lowest average distance between point $i$ and points in another cluster

| Silhouette Score | Interpretation |
| :--- | :--- |
| **Close to +1** | Point is well matched to its cluster |
| **Around 0** | Point lies near a cluster boundary |
| **Negative** | Point may be assigned to an inappropriate cluster |

A higher average silhouette score across all points generally indicates better cluster separation, though it should be interpreted alongside biological context.

### C. Gap Statistic

The gap statistic compares the observed within-cluster dispersion with the expected dispersion under a reference null distribution. A larger gap suggests that the observed clustering structure is stronger than expected by chance, making it helpful when the elbow method is ambiguous.

### Important Considerations

No single metric guarantees that the selected $K$ is biologically correct. The final choice should account for:
- Cluster separation
- Cluster stability across initializations
- Expression-pattern interpretability
- Biological relevance
- Reproducibility across datasets or experimental batches

---

## 6. Initialization Strategies

Because K-means is sensitive to initial centroid placement, poor initialization can result in suboptimal clustering.

### Random Initialization

Selects $K$ observations randomly as initial centroids.
- **Advantages:** Simple, fast, and easy to implement.
- **Disadvantages:** Can produce inconsistent results and converge to poor local optima.

### K-Means++ Initialization

K-means++ selects initial centroids using a probability distribution that favors points far away from already selected centroids:

1. Choose the first centroid randomly from the data points.
2. Calculate each point's squared distance to its nearest selected centroid.
3. Select the next centroid with probability proportional to that squared distance.
4. Repeat steps 2–3 until $K$ centroids have been selected.
5. Run the standard K-means algorithm.

K-means++ generally yields better, more consistent initializations than naive random selection and is the default in most modern implementations.

### Multiple Initializations

Running K-means multiple times with different initial centroids reduces dependence on a single random seed. The algorithm retains the solution with the lowest total inertia:

```python
from sklearn.cluster import KMeans

kmeans = KMeans(
    n_clusters=5,
    init="k-means++",
    n_init=25,
    max_iter=300,
    random_state=42
)

labels = kmeans.fit_predict(expression_matrix)
```

---

## 7. Computational Complexity

Let:
- $n$ = number of data points
- $K$ = number of clusters
- $d$ = number of features or dimensions
- $I$ = number of iterations

### Time Complexity

During each iteration, every data point is compared with every cluster centroid. Calculating one Euclidean distance requires $O(d)$ operations:
- **Assignment step:** $O(nKd)$
- **Centroid-update step:** $O(nd)$

Combining both steps across $I$ iterations gives the standard practical complexity for Lloyd's algorithm on dense data:

$$O(nKdI)$$

Actual runtime depends on convergence speed, initialization strategy, data representation, and low-level implementation.

### Space Complexity

- **Input dataset:** $O(nd)$
- **Centroids:** $O(Kd)$
- **Cluster assignments:** $O(n)$

The overall space complexity is:

$$O(nd + Kd + n)$$

For datasets where $d > 1$, this simplifies to $O(nd + Kd)$ (excluding temporary buffers used by specific implementations).

### Practical Example

Suppose a dataset contains $10,000$ genes, $20$ experimental conditions, $5$ clusters, and runs for $100$ iterations. The approximate number of distance calculations is:

$$n \times K \times I = 10,000 \times 5 \times 100 = 5,000,000$$

Each distance calculation evaluates $20$ dimensions, illustrating why K-means scales efficiently for moderate-sized gene expression datasets.

---

## 8. Advantages of K-Means

1. **Simple implementation:** Straightforward to understand, implement, and debug.
2. **Computationally efficient:** Scales linearly with data size and iterations for dense datasets.
3. **Interpretable results:** Cluster centroids provide direct, representative expression profiles.
4. **Flexible:** Applicable across transcriptomics, proteomics, and metabolomics datasets.
5. **Scalable:** Handles large datasets effectively with optimized or mini-batch implementations.
6. **Easy visualization:** Cluster assignments integrate cleanly into heatmaps and dimensionality-reduction plots.

---

## 9. Limitations of K-Means

Despite its utility, K-means has several well-known constraints:

- **A. Requires $K$ in advance:** The number of clusters must be specified prior to running the algorithm, and heuristic selection methods can remain ambiguous.
- **B. Sensitive to initialization:** Different starting centroids can yield different cluster partitions (mitigated via K-means++ and multiple starts).
- **C. Sensitive to outliers:** Because centroids are arithmetic means, extreme values or technical artifacts can pull centroids away from true cluster centers.
- **D. Assumes spherical clusters:** Minimizing Euclidean distance assumes clusters are convex, compact, and roughly spherical with similar variance and density. It struggles with irregular shapes or heavily overlapping populations.
- **E. Sensitive to feature scaling:** Features with larger numerical ranges dominate Euclidean distances; preprocessing must align with biological goals.
- **F. Does not model biological relationships directly:** Grouping relies strictly on numerical similarity and ignores gene regulatory networks, protein-protein interactions, and pathway topology.
- **G. Arbitrary cluster labels:** Cluster numbers (e.g., Cluster 1 vs. Cluster 2) are nominal identifiers with no inherent ordering or rank.

---

## 10. K-Means vs. Other Clustering Algorithms

| Feature | K-Means | Hierarchical Clustering | DBSCAN | Gaussian Mixture Model |
| :--- | :--- | :--- | :--- | :--- |
| **Number of clusters** | Specified in advance | Selected by cutting the dendrogram | Not specified directly | Usually specified in advance |
| **Cluster shape** | Compact, approximately spherical | Depends on linkage and distance | Arbitrary / irregular shapes | Elliptical / Gaussian shapes |
| **Assignment** | Hard assignment | Hard assignment at tree cut | Density-based; identifies noise | Probabilistic (soft) membership |
| **Outlier sensitivity** | Relatively high | Depends on linkage | High robustness (labels noise) | Moderate (depends on covariance) |
| **Scalability** | High (linear in $n$) | Lower ($O(n^2)$ or $O(n^3)$) | Moderate (depends on indexing) | Moderate |
| **Main strength** | Simplicity and speed | Nested hierarchical structure | Noise handling & irregular shapes | Probabilistic cluster modeling |

- **K-Means vs. Hierarchical Clustering:** Hierarchical clustering constructs a nested tree (dendrogram) without requiring $K$ upfront, making it ideal for examining multi-level relationships.
- **K-Means vs. DBSCAN:** DBSCAN groups points based on local density, detecting irregular shapes and isolating outliers as noise without needing $K$, though it requires tuning neighborhood radius ($\epsilon$) and minimum points.
- **K-Means vs. Gaussian Mixture Models (GMMs):** GMMs model data as a mixture of Gaussian distributions, assigning soft probabilities of cluster membership—useful when biological populations overlap.

---

## 11. K-Means and Deep Learning

While traditional K-means operates directly on raw or standardized input features, deep learning can learn compressed, non-linear representations of complex, high-dimensional biological data prior to clustering.

```text
High-dimensional gene expression matrix
                  |
                  v
               Encoder
                  |
                  v
      Low-dimensional latent space
                  |
                  v
               K-means
                  |
                  v
         Cluster assignments
                  |
                  v
       Biological interpretation
```

### Autoencoder + K-Means

An autoencoder consists of an **encoder** (which compresses input data into a low-dimensional latent space) and a **decoder** (which reconstructs the original input). Applying K-means to the latent representations allows the algorithm to capture non-linear biological patterns that standard Euclidean distances in high-dimensional space might miss.

### Deep Embedded Clustering (DEC)

**Deep Embedded Clustering (DEC)** jointly optimizes the low-dimensional representation and cluster assignments using a unified clustering objective. While powerful for complex datasets, deep clustering introduces additional challenges:
- Architecture and hyperparameter selection
- Training stability and sensitivity to initialization
- Risk of overfitting to technical artifacts rather than true biological signals

Deep learning does not automatically guarantee biologically superior clusters; results still require rigorous quality metrics and biological validation.

---

## 12. Practical Applications in Bioinformatics

- **Gene Expression Clustering:** Identifying co-expressed genes, treatment-responsive modules, developmental time-course patterns, and candidate genes for pathway analysis.
- **Single-Cell RNA Sequencing (scRNA-seq):** Grouping cells with similar transcriptional profiles. While graph-based methods (such as Louvain and Leiden) are standard for handling scRNA-seq sparsity and complex manifolds, K-means remains useful for exploratory comparisons and sub-clustering.
- **Proteomics:** Clustering protein abundance profiles to identify co-regulated proteins across experimental states.
- **Metabolomics:** Grouping metabolite abundance patterns associated with disease progression, drug response, or environmental exposure.
- **Multi-Omics Analysis:** Clustering integrated representations across multiple molecular modalities (requires careful cross-modality scaling so one data type does not dominate).

---

## 13. Best Practices for Gene Expression Clustering

1. **Normalize expression data:** Never cluster raw RNA-seq counts directly using Euclidean distance.
2. **Apply an appropriate transformation:** Use VST, `rlog`, or log-transformed normalized counts.
3. **Select informative genes:** Filter out invariant or low-expression genes to reduce noise.
4. **Choose scaling based on the biological question:** Use row-wise Z-scores to group by expression pattern shape rather than absolute abundance.
5. **Use multiple initializations:** Rely on K-means++ and multiple random starts (`nstart` / `n_init`) to avoid poor local optima.
6. **Evaluate multiple values of $K$:** Combine the elbow method, silhouette score, and biological interpretability.
7. **Check cluster stability:** Verify reproducibility across different random seeds or bootstrapped subsamples.
8. **Visualize cluster profiles:** Inspect results using clustered heatmaps and mean profile line plots.
9. **Perform functional enrichment analysis:** Test clusters for enrichment in GO terms, KEGG, or Reactome pathways.
10. **Account for batch effects:** Ensure technical batch variation is corrected so it does not drive cluster separation.
11. **Avoid overinterpreting clusters:** Co-expression suggests correlation, not direct transcriptional regulation or identical function.
12. **Document preprocessing decisions:** Record normalization methods, transformations, gene filters, scaling choices, $K$, initialization strategies, and random seeds.

---

## 14. Connection to Statistical and Computational Concepts

| Concept | Connection to K-Means |
| :--- | :--- |
| **Euclidean geometry** | Measures straight-line distance between observations in feature space |
| **Optimization** | Minimizes the within-cluster sum of squares (WCSS) objective function |
| **Iterative algorithms** | Alternates between expectation (assignment) and maximization (centroid update) steps |
| **Vector arithmetic** | Computes cluster centroids as component-wise arithmetic means |
| **Unsupervised learning** | Discovers latent group structure without labeled training data |
| **Dimensionality reduction** | Frequently paired with PCA or autoencoders to mitigate the curse of dimensionality |
| **Functional enrichment** | Bridges numerical clustering with biological pathway interpretation |

> **Key Conceptual Distinction:** K-means identifies **numerical similarity**, not biological causality. Two genes may share nearly identical expression profiles without directly regulating each other or participating in the same molecular pathway.

---

## References and Further Reading

1. **MacQueen, J. (1967).** Some Methods for Classification and Analysis of Multivariate Observations. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 281–297.
2. **Lloyd, S. (1982).** Least Squares Quantization in PCM. *IEEE Transactions on Information Theory*, 28(2), 129–137.
3. **Arthur, D. and Vassilvitskii, S. (2007).** k-means++: The Advantages of Careful Seeding. *Proceedings of the Eighteenth Annual ACM-SIAM Symposium on Discrete Algorithms*, 1027–1035.
4. **Tibshirani, R., Walther, G., and Hastie, T. (2001).** Estimating the Number of Clusters in a Data Set Via the Gap Statistic. *Journal of the Royal Statistical Society: Series B*, 63(2), 411–423.
5. **Rousseeuw, P. J. (1987).** Silhouettes: A Graphical Aid to the Interpretation and Validation of Cluster Analysis. *Journal of Computational and Applied Mathematics*, 20, 53–65.
6. **Eisen, M. B. et al. (1998).** Cluster analysis and display of genome-wide expression patterns. *Proceedings of the National Academy of Sciences*, 95(25), 14863–14868.
7. **Xie, J., Girshick, R., and Farhadi, A. (2016).** Unsupervised Deep Embedding for Clustering Analysis. *Proceedings of the 33rd International Conference on Machine Learning*, 478–487.
8. **DESeq2 Documentation:** [https://bioconductor.org/packages/DESeq2/](https://bioconductor.org/packages/DESeq2/)
9. **Scikit-learn K-Means Documentation:** [https://scikit-learn.org/stable/modules/clustering.html#k-means](https://scikit-learn.org/stable/modules/clustering.html#k-means)
10. **R `stats` Documentation — `kmeans()`:** [https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html)