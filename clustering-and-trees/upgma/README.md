# UPGMA — Unweighted Pair Group Method with Arithmetic Mean

## Table of Contents

* [What Is UPGMA?](#what-is-upgma)
* [Understanding the Molecular Clock Assumption](#understanding-the-molecular-clock-assumption)
* [How UPGMA Works](#how-upgma-works)
* [Distance Recomputation: Why Is It Weighted?](#distance-recomputation-why-is-it-weighted)
* [Worked Example: Five Bacterial Isolates](#worked-example-five-bacterial-isolates)
* [Python Implementation](#python-implementation)
* [Generating Evolutionary Distances from 16S rRNA Sequences](#generating-evolutionary-distances-from-16s-rrna-sequences)
* [Computational Complexity](#computational-complexity)
* [Advantages of UPGMA](#advantages-of-upgma)
* [Limitations](#limitations)
* [UPGMA vs. Neighbor Joining](#upgma-vs-neighbor-joining)
* [UPGMA Beyond Phylogenetics](#upgma-beyond-phylogenetics)
* [Scalability and Memory-Constrained Variants](#scalability-and-memory-constrained-variants)
* [Bootstrap Support and Tree Reliability](#bootstrap-support-and-tree-reliability)
* [Key Takeaways](#key-takeaways)
* [References](#references)

## What Is UPGMA?

**UPGMA (Unweighted Pair Group Method with Arithmetic Mean)** is a hierarchical clustering algorithm commonly used to construct rooted phylogenetic trees from pairwise evolutionary distance matrices.

The algorithm progressively merges the two closest clusters and recalculates distances between the newly formed cluster and all remaining clusters. UPGMA assumes a **molecular clock**, meaning that all evolutionary lineages accumulate changes at the same average rate. Under this assumption, the resulting phylogenetic tree is **ultrametric**: every leaf (taxon) has the exact same distance from the root.

For example, suppose we have 16S rRNA sequences from five bacterial isolates:

* Isolate A
* Isolate B
* Isolate C
* Isolate D
* Isolate E

After multiple sequence alignment, we calculate pairwise evolutionary distances. UPGMA uses these distances to construct a hierarchical tree representing the relationships among the isolates.

### Key Characteristics

* **Input:** A symmetric pairwise distance matrix.
* **Output:** A rooted, ultrametric tree.
* **Clustering strategy:** Average linkage.
* **Evolutionary assumption:** Constant evolutionary rate across lineages (molecular clock).
* **Branch lengths:** Derived from cluster merge heights.
* **Objective:** Repeatedly merge the closest clusters based on average pairwise distances.

---

## Understanding the Molecular Clock Assumption

The **molecular clock hypothesis** assumes that molecular substitutions accumulate at approximately constant rates over evolutionary time.

Suppose two bacterial lineages diverged from a common ancestor $10$ million years ago:

* **Lineage A** accumulates $10$ substitutions per million years.
* **Lineage B** accumulates $10$ substitutions per million years.

After $10$ million years, both lineages have accumulated approximately $100$ substitutions. Their evolutionary distances from the common ancestor are therefore equal under this simplified model. UPGMA relies on this assumption to estimate the relative timing of evolutionary divergence.

### What Is an Ultrametric Tree?

An **ultrametric tree** is a rooted tree in which every leaf has the same total branch length from the root. For any three taxa $A$, $B$, and $C$, their pairwise distances satisfy the three-point condition:

$$
D(A, B) \leq \max(D(A, C), D(B, C))
$$

Equivalently, among the three pairwise distances, the two largest values are equal. For example:

$$
D(A, B) = 4, \quad D(A, C) = 8, \quad D(B, C) = 8
$$

This distance matrix is consistent with an ultrametric tree. The equality of root-to-leaf distances allows UPGMA to interpret cluster heights as relative divergence times, provided the molecular-clock assumptions and distance calibration are appropriate.

> **Important:** An ultrametric tree is not simply a tree with equal branch lengths. Individual branches can have different lengths; what must be equal is the total distance from the root to every leaf.

---

## How UPGMA Works

UPGMA begins by treating every taxon as an independent cluster. At each iteration, it identifies the two clusters with the smallest average pairwise distance and merges them into a new cluster until all taxa belong to a single cluster.

### Step 1: Initialize Clusters

Suppose there are five bacterial isolates:

$$
C = \{A, B, C, D, E\}
$$

Initially, each isolate is its own singleton cluster:

$$
\{A\}, \{B\}, \{C\}, \{D\}, \{E\}
$$

### Step 2: Find the Closest Pair

Examine the pairwise distance matrix and identify the two clusters with the smallest distance. For example, if $D(A, B) = 4$ is the smallest non-zero distance in the matrix, merge $A$ and $B$.

### Step 3: Create a New Cluster

Create a new cluster $C_{AB} = \{A, B\}$. The height of the new internal node is half the distance between the two merged clusters:

$$
h(C_{AB}) = \frac{D(A, B)}{2}
$$

For a distance of $4$:

$$
h(C_{AB}) = \frac{4}{2} = 2
$$

This height determines the branch lengths from the new node to its children.

### Step 4: Recalculate Distances

Calculate the distance between the new cluster and every remaining cluster using **average linkage**. If clusters $C_1$ and $C_2$ merge into cluster $C$, then the distance to another cluster $C^*$ is:

$$
D(C, C^*) = \frac{\vert{}C_1\vert{} D(C_1, C^*) + \vert{}C_2\vert{} D(C_2, C^*)}{\vert{}C_1\vert{} + \vert{}C_2\vert{}}
$$

Where:

* $\vert{}C_1\vert{}$ = number of original taxa in cluster $C_1$
* $\vert{}C_2\vert{}$ = number of original taxa in cluster $C_2$
* $D(C_1, C^*)$ = average distance between clusters $C_1$ and $C^*$
* $D(C_2, C^*)$ = average distance between clusters $C_2$ and $C^*$

This is a **weighted average of the cluster distances**, where the weights are the numbers of original taxa in the clusters.

### Step 5: Repeat Clustering

After updating the distance matrix:

1. Find the closest pair of remaining clusters.
2. Merge them.
3. Calculate the new node height.
4. Update distances using average linkage.

Repeat until only one cluster remains.

### Step 6: Construct the Rooted Tree

The final cluster becomes the root of the tree. Each internal node represents a cluster-merging event, and its height represents the evolutionary distance from the leaves to that node under the molecular-clock interpretation.

---

## Distance Recomputation: Why Is It Weighted?

UPGMA uses the arithmetic mean of all pairwise distances between two clusters. Suppose cluster $C_1 = \{A, B\}$ contains two taxa and cluster $C_2 = \{C\}$ contains one taxon. The distance between the clusters is:

$$
D(C_1, C_2) = \frac{D(A, C) + D(B, C)}{2}
$$

For example, if $D(A, C) = 8$ and $D(B, C) = 10$:

$$
D(C_1, C_2) = \frac{8 + 10}{2} = 9
$$

If another cluster contains three taxa, its contribution must be weighted by three when calculating the average distance to the new cluster. This ensures that UPGMA calculates the average distance across all pairs of original taxa rather than simply averaging cluster distances without considering cluster sizes.

> **Important:** UPGMA is called *unweighted* because each original taxon contributes equally to the average, not because all clusters receive equal weight during distance updates.

---

## Worked Example: Five Bacterial Isolates

### Scenario

Suppose you have 16S rRNA sequences from five bacterial isolates. After multiple sequence alignment and evolutionary distance estimation, you obtain the following illustrative distance matrix:

| Isolate | A | B | C | D | E |
| :--- | ---: | ---: | ---: | ---: | ---: |
| **A** | 0 | 4 | 8 | 10 | 10 |
| **B** | 4 | 0 | 8 | 10 | 10 |
| **C** | 8 | 8 | 0 | 6 | 10 |
| **D** | 10 | 10 | 6 | 0 | 10 |
| **E** | 10 | 10 | 10 | 10 | 0 |

The objective is to construct an ultrametric tree using UPGMA.

### Step 1: Identify the Closest Pair

The smallest non-zero distance is $D(A, B) = 4$. Merge $A$ and $B$ into $C_{AB} = \{A, B\}$. Its node height is:

$$
h(C_{AB}) = \frac{4}{2} = 2
$$

The branch lengths from the new node to $A$ and $B$ are both $2$.

### Step 2: Update Distances

Calculate the distance between cluster $C_{AB}$ and the remaining isolates:

$$
D(C_{AB}, C) = \frac{D(A, C) + D(B, C)}{2} = \frac{8 + 8}{2} = 8
$$

$$
D(C_{AB}, D) = \frac{10 + 10}{2} = 10
$$

$$
D(C_{AB}, E) = \frac{10 + 10}{2} = 10
$$

The updated distance matrix is:

| Cluster | AB | C | D | E |
| :--- | ---: | ---: | ---: | ---: |
| **AB** | 0 | 8 | 10 | 10 |
| **C** | 8 | 0 | 6 | 10 |
| **D** | 10 | 6 | 0 | 10 |
| **E** | 10 | 10 | 10 | 0 |

### Step 3: Merge $C$ and $D$

The smallest remaining non-zero distance is $D(C, D) = 6$. Merge $C$ and $D$ into $C_{CD} = \{C, D\}$. The new node height is:

$$
h(C_{CD}) = \frac{6}{2} = 3
$$

The branch lengths from the new node to $C$ and $D$ are both $3$.

### Step 4: Update Distances Again

Calculate the distance between clusters $AB$ (size $2$) and $CD$ (size $2$):

$$
D(C_{AB}, C_{CD}) = \frac{2 D(C_{AB}, C) + 2 D(C_{AB}, D)}{4} = \frac{2(8) + 2(10)}{4} = \frac{36}{4} = 9
$$

The distances to $E$ are:

$$
D(C_{AB}, E) = 10
$$

$$
D(C_{CD}, E) = \frac{D(C, E) + D(D, E)}{2} = \frac{10 + 10}{2} = 10
$$

The updated matrix becomes:

| Cluster | AB | CD | E |
| :--- | ---: | ---: | ---: |
| **AB** | 0 | 9 | 10 |
| **CD** | 9 | 0 | 10 |
| **E** | 10 | 10 | 0 |

### Step 5: Merge $AB$ and $CD$

The smallest remaining non-zero distance is $D(AB, CD) = 9$. Merge the two clusters into $C_{ABCD} = \{A, B, C, D\}$. The new node height is:

$$
h(C_{ABCD}) = \frac{9}{2} = 4.5
$$

Because the child nodes have different heights, their branch lengths to $C_{ABCD}$ differ:

$$
L(AB, C_{ABCD}) = 4.5 - 2 = 2.5
$$

$$
L(CD, C_{ABCD}) = 4.5 - 3 = 1.5
$$

### Step 6: Connect the Final Isolate

The distance between cluster $ABCD$ and $E$ is:

$$
D(C_{ABCD}, E) = \frac{10 + 10 + 10 + 10}{4} = 10
$$

The final root height is:

$$
h(\text{Root}) = \frac{10}{2} = 5
$$

The branch lengths to the root are:

$$
L(ABCD, \text{Root}) = 5 - 4.5 = 0.5
$$

$$
L(E, \text{Root}) = 5 - 0 = 5
$$

### Final UPGMA Tree

```text
                    Root (5)
                   /        \
             ABCD (4.5)      E
               /      \
         AB (2)       CD (3)
          / \          / \
         A   B        C   D
```

| Branch | Length |
| :--- | ---: |
| **Root to E** | 5.0 |
| **Root to ABCD** | 0.5 |
| **ABCD to AB** | 2.5 |
| **ABCD to CD** | 1.5 |
| **AB to A** | 2.0 |
| **AB to B** | 2.0 |
| **CD to C** | 3.0 |
| **CD to D** | 3.0 |

Every leaf is now at a total distance of $5$ from the root:

* **A:** $0.5 + 2.5 + 2 = 5$
* **B:** $0.5 + 2.5 + 2 = 5$
* **C:** $0.5 + 1.5 + 3 = 5$
* **D:** $0.5 + 1.5 + 3 = 5$
* **E:** $5$

This demonstrates the ultrametric property of the UPGMA tree.

---

## Python Implementation

UPGMA can be implemented from scratch using a distance matrix and cluster-size tracking.

### Educational Implementation from Scratch

```python
def upgma(distance_matrix, labels):
    """
    Construct a rooted UPGMA tree from a distance matrix.

    Parameters
    ----------
    distance_matrix : list of lists
        Symmetric pairwise distance matrix.
    labels : list
        Taxon names.

    Returns
    -------
    dict
        Tree representation containing cluster relationships,
        node heights, and branch lengths.
    """
    n = len(labels)

    # Initialize singleton clusters
    clusters = {
        i: {
            "name": labels[i],
            "size": 1,
            "height": 0.0,
            "children": None
        }
        for i in range(n)
    }

    # Store pairwise distances for upper-triangular pairs
    distances = {}
    for i in range(n):
        for j in range(i + 1, n):
            distances[(i, j)] = float(distance_matrix[i][j])

    next_id = n

    while len(clusters) > 1:
        active = list(clusters.keys())

        # Find closest pair of active clusters
        i, j = min(distances, key=lambda pair: distances[pair])
        merge_distance = distances[(i, j)]

        # Calculate new node height
        new_height = merge_distance / 2.0

        cluster_i = clusters[i]
        cluster_j = clusters[j]

        # Calculate branch lengths from new internal node to children
        branch_i = new_height - cluster_i["height"]
        branch_j = new_height - cluster_j["height"]

        # Create merged cluster
        new_cluster = {
            "name": f"Cluster_{next_id}",
            "size": cluster_i["size"] + cluster_j["size"],
            "height": new_height,
            "children": [
                (cluster_i, branch_i),
                (cluster_j, branch_j)
            ]
        }

        # Update distances to remaining clusters using weighted average
        new_distances = {}
        for k in active:
            if k == i or k == j:
                continue

            d_ik = distances[tuple(sorted((i, k)))]
            d_jk = distances[tuple(sorted((j, k)))]

            size_i = cluster_i["size"]
            size_j = cluster_j["size"]

            new_d = (size_i * d_ik + size_j * d_jk) / (size_i + size_j)
            new_distances[(min(next_id, k), max(next_id, k))] = new_d

        # Remove merged clusters and their old distances
        del clusters[i]
        del clusters[j]

        distances = {
            pair: value
            for pair, value in distances.items()
            if i not in pair and j not in pair
        }

        # Add new cluster and updated distances
        clusters[next_id] = new_cluster
        distances.update(new_distances)
        next_id += 1

    return next(iter(clusters.values()))


def print_tree(node, prefix="", branch_length=None):
    """Print a simple indented text representation of the tree."""
    if branch_length is not None:
        print(f"{prefix}{node['name']} ({branch_length:.2f})")
    else:
        print(f"{prefix}{node['name']}")

    if node["children"]:
        for child, length in node["children"]:
            print_tree(child, prefix + "  ", length)


# Example distance matrix
distance_matrix = [
    [0, 4, 8, 10, 10],
    [4, 0, 8, 10, 10],
    [8, 8, 0, 6, 10],
    [10, 10, 6, 0, 10],
    [10, 10, 10, 10, 0]
]

labels = ["A", "B", "C", "D", "E"]
tree = upgma(distance_matrix, labels)
print_tree(tree)
```

*Note: This is an educational implementation. For production use, validate the input matrix, handle ties explicitly, and use a well-tested phylogenetic or clustering library.*

### Using Biopython

Biopython provides a UPGMA tree constructor through its `DistanceTreeConstructor` class:

```python
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix,
    DistanceTreeConstructor
)
from Bio import Phylo

names = ["A", "B", "C", "D", "E"]

# Lower-triangular distance matrix
matrix = [
    [0],
    [4, 0],
    [8, 8, 0],
    [10, 10, 6, 0],
    [10, 10, 10, 10, 0]
]

dm = DistanceMatrix(names, matrix)
constructor = DistanceTreeConstructor()

# Construct UPGMA tree
tree = constructor.upgma(dm)

# Write tree in Newick format
Phylo.write(tree, "upgma_tree.nwk", "newick")

# Display tree
Phylo.draw(tree)
```

The resulting Newick file stores the tree topology and branch lengths.

---

## Generating Evolutionary Distances from 16S rRNA Sequences

UPGMA requires a pairwise distance matrix rather than raw sequences. A typical phylogenetic workflow for bacterial isolates is:

1. Collect homologous 16S rRNA sequences.
2. Perform multiple sequence alignment.
3. Inspect alignment quality.
4. Trim poorly aligned or unreliable regions when appropriate.
5. Estimate pairwise evolutionary distances.
6. Construct a distance matrix.
7. Apply UPGMA.
8. Evaluate the tree and molecular-clock assumption.

### Step 1: Multiple Sequence Alignment

Common tools include **MAFFT**, **MUSCLE**, and **Clustal Omega**. The sequences should represent homologous regions so that each alignment column corresponds to the same evolutionary position.

### Step 2: Estimate Evolutionary Distances

A simple distance measure is the proportion of differing aligned sites ($p$-distance):

$$
p = \frac{\text{Number of differing sites}}{\text{Number of comparable sites}}
$$

For example, if two sequences differ at $20$ out of $1,000$ aligned positions, $p = \frac{20}{1000} = 0.02$ ($2\%$ divergence). However, raw $p$-distance does not account for multiple substitutions at the same site. Model-based distances, such as the **Jukes-Cantor** correction for nucleotide data, account for unobserved substitutions:

$$
d = -\frac{3}{4} \ln\left(1 - \frac{4p}{3}\right)
$$

The Jukes-Cantor model assumes equal nucleotide frequencies and equal substitution rates; more complex models (e.g., Kimura 2-parameter, GTR) may be appropriate when these assumptions are not reasonable.

> **Important:** The choice of distance model affects the resulting UPGMA tree. The molecular-clock assumption is separate from the substitution model used to estimate pairwise distances.

---

## Computational Complexity

Let:

* $n$ = number of taxa
* $D$ = pairwise distance matrix
* $k$ = number of active clusters at an iteration

### Time Complexity: $O(n^3)$

For a straightforward implementation, each iteration performs:

1. Search for the smallest distance between active clusters: $O(k^2)$
2. Merge the selected clusters: $O(1)$
3. Recalculate distances to the remaining clusters: $O(k)$

Across $n - 1$ merges, the minimum-distance search dominates:

$$
\sum_{k=2}^{n} O(k^2) = O(n^3)
$$

Optimized implementations using priority queues or nearest-neighbor chains can reduce practical runtime, particularly for repeated minimum-pair searches.

### Space Complexity: $O(n^2)$

The algorithm stores a pairwise distance matrix containing up to $\frac{n(n - 1)}{2}$ unique distances, alongside cluster membership and node heights ($O(n)$). The distance matrix dominates memory consumption.

---

## Advantages of UPGMA

* **Simple implementation:** The algorithm is straightforward and easy to understand.
* **Efficient for moderate datasets:** Conventional $O(n^3)$ complexity runs quickly on small-to-moderate taxon counts.
* **Produces rooted trees:** The final merge naturally defines a root under the molecular-clock assumption.
* **Branch-length estimation:** Merge heights directly provide branch lengths and relative divergence estimates.
* **Useful for hierarchical clustering:** The average-linkage principle is widely used across general data science and bioinformatics.
* **Interpretable results:** Produces clear dendrograms showing nested clusters and relative distances.

## Limitations

1. **Requires a molecular clock for phylogenetic interpretation:** Assumes approximately equal evolutionary rates across all lineages. If some species evolve faster than others, UPGMA can infer an incorrect topology.
2. **Forces an ultrametric tree:** All leaves are constrained to have the exact same distance from the root, which is unrealistic when evolutionary rates vary.
3. **Sensitive to distance estimates:** Errors in sequence alignment, substitution models, and pairwise distance estimation directly affect clustering.
4. **Greedy clustering:** Makes locally optimal merging decisions at each step without searching all possible tree topologies for a globally optimal explanation.
5. **Requires a complete distance matrix:** Quadratic memory usage ($O(n^2)$) becomes a bottleneck for very large datasets.
6. **Divergence times require calibration:** Cluster heights are measured in units of sequence divergence, not years; estimating absolute divergence times requires an external molecular clock rate or fossil calibration.

---

## UPGMA vs. Neighbor Joining

UPGMA and Neighbor Joining are both distance-based phylogenetic reconstruction methods, but they make fundamentally different assumptions:

| Feature | UPGMA | Neighbor Joining |
| :--- | :--- | :--- |
| **Input** | Distance matrix | Distance matrix |
| **Tree type** | Rooted | Unrooted |
| **Molecular clock** | Assumes a molecular clock | Does not require a molecular clock |
| **Clustering strategy** | Average linkage | Q-matrix criterion |
| **Branch lengths** | Derived from merge heights | Estimated during joining |
| **Rate variation** | Assumes equal rates | Allows unequal rates |
| **Conventional time complexity** | $O(n^3)$ | $O(n^3)$ |
| **Conventional space complexity** | $O(n^2)$ | $O(n^2)$ |

### When Should You Use UPGMA?

* The molecular-clock assumption is biologically reasonable.
* The input distances are approximately ultrametric.
* You need a simple rooted hierarchical tree.
* You want to explore clustering relationships using average linkage.

### When Should You Use Neighbor Joining?

* Evolutionary rates may differ among lineages.
* You do not want to impose an ultrametric constraint.
* You need an unrooted distance-based phylogeny.
* The distance matrix is not expected to satisfy the molecular-clock assumption.

---

## UPGMA Beyond Phylogenetics

Outside evolutionary biology, UPGMA is known as **average-linkage hierarchical clustering** and is widely applied to non-phylogenetic biological datasets:

1. **Gene expression clustering:** Grouping genes with similar expression profiles across tissue samples or experimental conditions (via normalization $\to$ feature selection $\to$ pairwise distances $\to$ average-linkage clustering $\to$ heatmap visualization).
2. **Protein sequence clustering:** Grouping homologous proteins based on pairwise sequence similarity to explore protein families.
3. **Microbial comparative genomics:** Exploring genomic similarity among bacterial isolates using ANI (Average Nucleotide Identity) or k-mer distances.
4. **Hierarchical sample clustering:** Clustering patient or experimental samples based on transcriptomic, methylation, proteomic, or metabolomic profiles.

> **Important:** In gene expression and sample clustering, UPGMA acts purely as a hierarchical clustering method rather than a molecular-clock model; the ultrametric evolutionary interpretation does not apply.

---

## Scalability and Memory-Constrained Variants

Because a complete pairwise distance matrix stores $\frac{n(n - 1)}{2}$ entries, memory grows quadratically with the number of taxa ($n$):

* **$1,000$ taxa:** $\approx 500,000$ unique pairwise distances
* **$10,000$ taxa:** $\approx 50,000,000$ unique pairwise distances
* **$100,000$ taxa:** $\approx 5,000,000,000$ unique pairwise distances

Memory-constrained approaches, such as **MC-UPGMA**, have been developed to reduce memory demands in large-scale clustering settings by avoiding full in-memory distance matrices. However, exact memory savings and computational trade-offs depend on implementation details and input sparsity. For massive phylogenetic datasets, alternatives include representative taxon selection, approximate clustering, or fast sequence-based phylogenetic tools.

---

## Bootstrap Support and Tree Reliability

A UPGMA tree can be evaluated using **bootstrap resampling** of the sequence alignment:

1. Start with a multiple sequence alignment containing $L$ columns.
2. Randomly sample $L$ columns with replacement.
3. Recalculate pairwise distances.
4. Reconstruct a UPGMA tree.
5. Repeat the process many times (e.g., $1,000$ replicates).
6. Calculate how frequently each cluster appears across replicate trees.

For example, if a cluster is recovered in $900$ out of $1,000$ bootstrap trees:

$$
\text{Bootstrap Support} = \frac{900}{1000} \times 100 = 90\%
$$

A bootstrap value of $90\%$ means that the cluster was recovered in $90\%$ of the resampled datasets under the chosen analysis procedure—it does **not** mean that the cluster has a $90\%$ probability of being historically true, nor does it validate whether the molecular-clock assumption holds.

---

## Key Takeaways

* **UPGMA** is an agglomerative hierarchical clustering algorithm that constructs rooted phylogenetic trees from pairwise distances.
* It iteratively merges the two closest clusters using **average linkage**, recalculating distances via a weighted average based on the number of original taxa in each cluster.
* The height of each internal node is half the distance at which its two child clusters merge ($h = D / 2$).
* UPGMA assumes a **molecular clock** and produces an **ultrametric tree** where every leaf has the exact same total distance from the root.
* Conventional implementations have $O(n^3)$ **time complexity** and $O(n^2)$ **space complexity**.
* UPGMA is suitable when evolutionary rates are approximately constant, whereas **Neighbor Joining** is preferred when evolutionary rates vary across lineages.
* Beyond phylogenetics, UPGMA is widely used as **average-linkage hierarchical clustering** without a molecular-clock interpretation.

## References

1. **Sokal, R. R., & Michener, C. D. (1958).** A statistical method for evaluating systematic relationships. *University of Kansas Science Bulletin*, 38, 1409–1438.
2. **Sneath, P. H. A., & Sokal, R. R. (1973).** *Numerical Taxonomy: The Principles and Practice of Numerical Classification*. W. H. Freeman.
3. **Felsenstein, J. (2004).** *Inferring Phylogenies*. Sinauer Associates.
4. **Biopython Documentation:** <https://biopython.org/>