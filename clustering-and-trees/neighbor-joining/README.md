# Neighbor Joining

## Table of Contents

* [What Is Neighbor Joining?](#what-is-neighbor-joining)
* [How Neighbor Joining Differs from UPGMA](#how-neighbor-joining-differs-from-upgma)
* [The Key Insight: Corrected Evolutionary Distances](#the-key-insight-corrected-evolutionary-distances)
* [Step-by-Step Neighbor Joining Algorithm](#step-by-step-neighbor-joining-algorithm)
* [Worked Example: Six Mammalian Species](#worked-example-six-mammalian-species)
* [Python Implementation](#python-implementation)
* [Generating Evolutionary Distances from Protein Sequences](#generating-evolutionary-distances-from-protein-sequences)
* [Computational Complexity](#computational-complexity)
* [Advantages of Neighbor Joining](#advantages-of-neighbor-joining)
* [Limitations](#limitations)
* [Neighbor Joining vs. Other Phylogenetic Methods](#neighbor-joining-vs-other-phylogenetic-methods)
* [Bootstrap Support for Neighbor Joining Trees](#bootstrap-support-for-neighbor-joining-trees)
* [Applications in Bioinformatics](#applications-in-bioinformatics)
* [Key Takeaways](#key-takeaways)
* [References](#references)

---

## What Is Neighbor Joining?

**Neighbor Joining (NJ)** is a distance-based phylogenetic reconstruction algorithm that builds an **unrooted phylogenetic tree** from a matrix of pairwise evolutionary distances between taxa.

Introduced by Naruya Saitou and Masatoshi Nei in 1987, the algorithm progressively joins pairs of taxa or existing clusters that are likely to be neighboring branches in the underlying evolutionary tree.

Unlike UPGMA, Neighbor Joining does not assume that all lineages evolve at the same rate. This makes it particularly useful for reconstructing phylogenies when evolutionary rates differ among species.

For example, consider protein sequences from six mammalian species:

* Human
* Chimpanzee
* Gorilla
* Orangutan
* Mouse
* Rat

After aligning their homologous protein sequences, we calculate pairwise evolutionary distances. NJ uses this distance matrix to infer a tree representing their evolutionary relationships.

### Key Characteristics

* **Input:** A symmetric pairwise distance matrix.
* **Output:** An unrooted phylogenetic tree.
* **Evolutionary assumption:** Does not require a molecular clock.
* **Approach:** Iteratively joins neighboring taxa or clusters.
* **Objective:** Minimize the total branch length of the tree under the algorithm's distance-based criterion.

---

## How Neighbor Joining Differs from UPGMA

Both Neighbor Joining and UPGMA reconstruct phylogenetic trees using pairwise evolutionary distances. However, their assumptions and joining strategies differ significantly:

| Feature | UPGMA | Neighbor Joining |
| :--- | :--- | :--- |
| **Molecular clock** | Assumes a constant evolutionary rate | Does not assume a constant rate |
| **Tree type** | Rooted, ultrametric tree | Unrooted tree |
| **Joining strategy** | Merges the closest clusters | Uses corrected distances to identify neighboring clusters |
| **Branch lengths** | Based on cluster merging heights | Estimated from pairwise distances and corrected distances |
| **Evolutionary rates** | Assumes equal rates across lineages | Allows unequal rates |
| **Typical application** | Data approximately satisfying a molecular clock | Data with potentially different evolutionary rates |

### Why Does the Molecular Clock Matter?

A **molecular clock** assumes that evolutionary changes accumulate at approximately the same rate across different lineages.

For example, suppose two species diverged from a common ancestor $10$ million years ago:

* **Lineage A** accumulates substitutions at $2$ substitutions per million years.
* **Lineage B** accumulates substitutions at $5$ substitutions per million years.

Their evolutionary distances from the common ancestor will differ even though they diverged at the exact same time. UPGMA assumes equal rates and can therefore produce an incorrect topology when this assumption is violated. Neighbor Joining accommodates unequal rates by considering how each taxon is related to all other taxa when deciding which pair to join.

> **Important:** NJ does not estimate divergence times directly. Because its output is unrooted, an external rooting method such as an outgroup or a suitable evolutionary model is needed to infer the direction of evolution.

---

## The Key Insight: Corrected Evolutionary Distances

A naive approach to building a phylogenetic tree would be to repeatedly join the pair of taxa with the smallest pairwise distance. However, this can be misleading:

* Two taxa may have a relatively large distance between them yet still be immediate evolutionary neighbors if their lineages have accumulated substitutions at different rates.
* Conversely, two taxa may have a small pairwise distance simply because both are relatively distant from the other taxa in the dataset.

Neighbor Joining addresses this by correcting pairwise distances using each taxon's total distance to all other taxa.

### 1. Pairwise Distance Matrix

Let:

* $n$ = number of taxa or current clusters
* $D(i, j)$ = evolutionary distance between taxa or clusters $i$ and $j$

The distance matrix is symmetric ($D(i, j) = D(j, i)$) with zero along the diagonal ($D(i, i) = 0$). For example, a distance matrix for four taxa might look like this:

| Taxon | A | B | C | D |
| :--- | ---: | ---: | ---: | ---: |
| **A** | 0 | 5 | 9 | 9 |
| **B** | 5 | 0 | 10 | 10 |
| **C** | 9 | 10 | 0 | 8 |
| **D** | 9 | 10 | 8 | 0 |

*Note: These distances represent estimated evolutionary divergence, not necessarily raw sequence mismatch percentages.*

### 2. Calculate the Total Distance for Each Taxon

For each taxon $i$, calculate its total distance ($r_i$) to all other taxa:

$$
r_i = \sum_{k=1}^{n} D(i, k)
$$

This provides information about how distant that taxon is from the rest of the dataset. For the example matrix:

$$
r_A = 5 + 9 + 9 = 23
$$

$$
r_B = 5 + 10 + 10 = 25
$$

$$
r_C = 9 + 10 + 8 = 27
$$

$$
r_D = 9 + 10 + 8 = 27
$$

### 3. Calculate the Q-Matrix

Neighbor Joining uses a corrected distance criterion called the **Q-matrix**. For every pair of taxa $i$ and $j$:

$$
Q(i, j) = (n - 2)D(i, j) - r_i - r_j
$$

Where:

* $n$ is the current number of taxa or clusters.
* $D(i, j)$ is their pairwise distance.
* $r_i$ and $r_j$ are their total distances to all other taxa.

The pair with the **smallest (most negative) Q-value** is selected for joining. For four taxa ($n = 4$), the Q-values are:

$$
Q(A, B) = 2(5) - 23 - 25 = -38
$$

$$
Q(A, C) = 2(9) - 23 - 27 = -32
$$

$$
Q(A, D) = 2(9) - 23 - 27 = -32
$$

$$
Q(B, C) = 2(10) - 25 - 27 = -32
$$

$$
Q(B, D) = 2(10) - 25 - 27 = -32
$$

$$
Q(C, D) = 2(8) - 27 - 27 = -38
$$

The smallest values are $Q(A, B) = Q(C, D) = -38$. Because there is a tie, either pair can be selected first, illustrating that multiple joining orders can be equally preferred by the criterion at a given step.

> **Why is this correction useful?** The Q-matrix considers both the distance between two taxa and their distances to the rest of the dataset, avoiding reliance on the smallest raw pairwise distance alone. Note that the Q-matrix is strictly a selection criterion it is not an evolutionary distance or a branch length.

---

## Step-by-Step Neighbor Joining Algorithm

Neighbor Joining repeatedly joins pairs of neighboring taxa until only two clusters remain.

### Step 1: Initialize the Distance Matrix

Start with a symmetric matrix of pairwise evolutionary distances between all $n$ taxa (e.g., six species produce a $6 \times 6$ distance matrix).

### Step 2: Calculate Row Sums

For each active taxon or cluster $i$, calculate its total distance to all other active taxa:

$$
r_i = \sum_{k=1}^{n} D(i, k)
$$

### Step 3: Construct the Q-Matrix

Calculate the Q-value for every pair $(i, j)$:

$$
Q(i, j) = (n - 2)D(i, j) - r_i - r_j
$$

Select the pair $(i, j)$ with the smallest Q-value.

### Step 4: Join the Selected Pair

Create a new internal node $u$ connecting the selected taxa $i$ and $j$. Estimate the branch lengths from the new node $u$ to each selected taxon using:

$$
L(i, u) = \frac{1}{2}D(i, j) + \frac{r_i - r_j}{2(n - 2)}
$$

$$
L(j, u) = D(i, j) - L(i, u)
$$

These equations preserve the pairwise distance between the selected taxa ($L(i, u) + L(j, u) = D(i, j)$). The branches need not have equal lengths because NJ does not assume a molecular clock.

### Step 5: Update the Distance Matrix

Replace taxa $i$ and $j$ with their newly formed cluster $u$. Calculate the distance between the new cluster $u$ and every remaining taxon $k$:

$$
D(u, k) = \frac{D(i, k) + D(j, k) - D(i, j)}{2}
$$

The number of active taxa or clusters decreases by one ($n \leftarrow n - 1$).

### Step 6: Repeat

Recalculate row sums and the Q-matrix using the updated distance matrix. Continue selecting and joining pairs until only two clusters remain.

### Step 7: Connect the Final Two Clusters

When only two clusters remain, connect them with a final branch whose length equals their remaining pairwise distance. The resulting structure is an unrooted phylogenetic tree.

---

## Worked Example: Six Mammalian Species

### Scenario

Suppose you have aligned homologous protein sequences from six mammalian species and calculated their pairwise evolutionary distances. The objective is to reconstruct their phylogenetic relationships without assuming a molecular clock.

Consider the following illustrative distance matrix:

| Species | Human | Chimp | Gorilla | Orangutan | Mouse | Rat |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **Human** | 0 | 5 | 7 | 9 | 20 | 21 |
| **Chimp** | 5 | 0 | 6 | 9 | 20 | 21 |
| **Gorilla** | 7 | 6 | 0 | 10 | 20 | 21 |
| **Orangutan** | 9 | 9 | 10 | 0 | 21 | 22 |
| **Mouse** | 20 | 20 | 20 | 21 | 0 | 4 |
| **Rat** | 21 | 21 | 21 | 22 | 4 | 0 |

*Note: These values are illustrative evolutionary distances, not measured distances from a specific protein dataset.*

### Step 1: Calculate Row Sums

Sum the distances in each row:

$$
r_{\text{Human}} = 5 + 7 + 9 + 20 + 21 = 62
$$

$$
r_{\text{Chimp}} = 5 + 6 + 9 + 20 + 21 = 61
$$

$$
r_{\text{Gorilla}} = 7 + 6 + 10 + 20 + 21 = 64
$$

$$
r_{\text{Orangutan}} = 9 + 9 + 10 + 21 + 22 = 71
$$

$$
r_{\text{Mouse}} = 20 + 20 + 20 + 21 + 4 = 85
$$

$$
r_{\text{Rat}} = 21 + 21 + 21 + 22 + 4 = 89
$$

### Step 2: Calculate Q-Values

Initially, $n = 6$, so $Q(i, j) = 4D(i, j) - r_i - r_j$:

$$
Q(\text{Human}, \text{Chimp}) = 4(5) - 62 - 61 = -103
$$

$$
Q(\text{Mouse}, \text{Rat}) = 4(4) - 85 - 89 = -158
$$

The smallest Q-value in the matrix is for **Mouse** and **Rat** ($-158$). Therefore, Mouse and Rat are selected as the first neighboring pair.

### Step 3: Calculate Branch Lengths

Given $D(\text{Mouse}, \text{Rat}) = 4$, calculate the branch lengths to the new internal node $u$:

$$
L(\text{Mouse}, u) = \frac{4}{2} + \frac{85 - 89}{2(6 - 2)} = 2 - 0.5 = 1.5
$$

$$
L(\text{Rat}, u) = 4 - 1.5 = 2.5
$$

The new internal node $u$ connects Mouse and Rat with branch lengths of $1.5$ and $2.5$. Notice that the branch lengths are unequal, reflecting different evolutionary rates permitted under NJ.

### Step 4: Update the Distance Matrix

Calculate the distance between the new cluster $u$ (Mouse-Rat) and each remaining species:

$$
D(u, \text{Human}) = \frac{D(\text{Mouse}, \text{Human}) + D(\text{Rat}, \text{Human}) - D(\text{Mouse}, \text{Rat})}{2} = \frac{20 + 21 - 4}{2} = 18.5
$$

Similarly:

$$
D(u, \text{Chimp}) = \frac{20 + 21 - 4}{2} = 18.5
$$

$$
D(u, \text{Gorilla}) = \frac{20 + 21 - 4}{2} = 18.5
$$

$$
D(u, \text{Orangutan}) = \frac{21 + 22 - 4}{2} = 19.5
$$

The updated matrix now contains five active clusters:

* Human
* Chimp
* Gorilla
* Orangutan
* Mouse-Rat ($u$)

### Step 5: Continue Until the Tree Is Complete

At each subsequent iteration:

1. Calculate row sums for the current matrix.
2. Construct the Q-matrix.
3. Select the pair with the smallest Q-value.
4. Estimate branch lengths.
5. Merge the selected pair.
6. Update distances.

This procedure continues until only two clusters remain, yielding an unrooted tree containing all six species.

> **Interpretation:** The algorithm evaluates the entire distance matrix to identify likely neighboring branches rather than greedily joining the pair with the smallest raw distance at every iteration.

---

## Python Implementation

The `scikit-bio` and Biopython ecosystems provide tools for working with phylogenetic distance matrices and trees. For a straightforward NJ implementation, Biopython provides `DistanceTreeConstructor`.

### Example Using Biopython

```python
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix,
    DistanceTreeConstructor
)
from Bio import Phylo

# Lower-triangular distance matrix.
# Each row contains distances to preceding taxa,
# including the diagonal zero.
names = [
    "Human",
    "Chimp",
    "Gorilla",
    "Orangutan",
    "Mouse",
    "Rat"
]

matrix = [
    [0],
    [5, 0],
    [7, 6, 0],
    [9, 9, 10, 0],
    [20, 20, 20, 21, 0],
    [21, 21, 21, 22, 4, 0]
]

dm = DistanceMatrix(names, matrix)

# Construct Neighbor Joining tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Export tree in Newick format
Phylo.write(tree, "mammalian_nj_tree.nwk", "newick")

# Display tree
Phylo.draw(tree)
```

The Newick output represents the tree topology and branch lengths. For example, a simplified Newick representation looks like:

```text
((Mouse:1.5,Rat:2.5):..., ...);
```

---

## Generating Evolutionary Distances from Protein Sequences

Because NJ requires a precomputed distance matrix, the first phase of a phylogenetic analysis involves estimating evolutionary distances from homologous sequences.

### Typical Workflow

1. Obtain homologous protein sequences.
2. Perform multiple sequence alignment (MSA).
3. Inspect and trim poorly aligned or gapped regions when appropriate.
4. Estimate pairwise evolutionary distances using a suitable amino acid substitution model.
5. Construct the symmetric distance matrix.
6. Apply Neighbor Joining.
7. Evaluate branch support and root the resulting tree.

### Common Tools and Models

* **Multiple Sequence Alignment Tools:** MAFFT, MUSCLE, Clustal Omega.
* **Protein Substitution Models:** JTT, WAG, LG, Dayhoff.

Raw percentage differences (p-distances) can underestimate true evolutionary divergence because multiple substitutions may occur at the same site over long timescales (saturation). Model-based distance estimates account for these unobserved substitutions.

> **Important:** NJ does not select or apply the substitution model itself it relies entirely on the distance matrix supplied to it. Errors or biases in sequence alignment or distance estimation directly affect the resulting tree.

---

## Computational Complexity

Let:

* $n$ = total number of taxa
* $D$ = pairwise distance matrix
* $k$ = number of active taxa or clusters during a given iteration

### Time Complexity: $O(n^3)$

A standard NJ implementation has cubic time complexity. At each iteration:

1. Calculate row sums: $O(k^2)$
2. Construct or evaluate the Q-matrix: $O(k^2)$
3. Find the minimum Q-value: $O(k^2)$
4. Update distances after merging: $O(k)$

Across approximately $n - 2$ joining iterations:

$$
\sum_{k=3}^{n} O(k^2) = O(n^3)
$$

While optimized variants (such as RapidNJ or FastTree-style heuristics) improve practical runtime using specialized data structures, the canonical algorithm runs in $O(n^3)$ time.

### Space Complexity: $O(n^2)$

A symmetric distance matrix stores $\frac{n(n - 1)}{2}$ unique pairwise distances, so memory usage grows quadratically with the number of taxa. The resulting tree structure requires only $O(n)$ nodes and branches, making the distance matrix the dominant memory cost.

---

## Advantages of Neighbor Joining

* **No molecular clock assumption:** Accommodates varying evolutionary rates across different lineages.
* **Computationally efficient:** $O(n^3)$ time complexity is fast and practical for moderate-sized datasets.
* **Unrooted tree reconstruction:** Infers evolutionary relationships without requiring prior knowledge of the root.
* **Branch-length estimation:** Estimates evolutionary divergence along individual branches.
* **Broad applicability:** Works with nucleotide, protein, morphological, or any valid metric/evolutionary distance matrix.
* **Well established:** Widely implemented and validated across standard bioinformatics packages.

---

## Limitations

1. **Depends on distance estimates:** NJ is only as accurate as the input distance matrix; poor alignments or misspecified substitution models can yield misleading topologies.
2. **Information loss:** Summarizing aligned sequences into pairwise distances discards site-specific character information utilized directly by Maximum Likelihood and Bayesian methods.
3. **Sensitivity at deep divergences:** Long evolutionary distances suffer from site saturation and high variance in distance estimation.
4. **Produces an unrooted tree:** Requires an outgroup, midpoint rooting, or molecular clock model to determine evolutionary direction.
5. **No inherent statistical support:** Produces a single point-estimate tree without confidence intervals; requires resampling methods like bootstrapping to assess branch support.
6. **Negative branch lengths:** Can produce negative branch lengths when input distances deviate from strict tree additivity due to estimation noise. Software implementations may retain, truncate to zero, or adjust adjacent branches to handle these values.

---

## Neighbor Joining vs. Other Phylogenetic Methods

| Method | Input | Main Approach | Molecular Clock |
| :--- | :--- | :--- | :--- |
| **UPGMA** | Distance matrix | Agglomerative hierarchical clustering | Required |
| **Neighbor Joining** | Distance matrix | Corrected-distance (Q-matrix) joining | Not required |
| **Maximum Likelihood** | Sequence alignment | Optimizes tree topology and branch lengths to maximize sequence likelihood | Not required by default |
| **Bayesian Phylogenetics** | Sequence alignment | Estimates posterior probability distributions of trees and parameters | Depends on model |
| **Minimum Evolution** | Distance matrix | Searches tree space to minimize total branch length | Not required |

### When Should You Use Neighbor Joining?

NJ is well suited when:

* You have a reliable pairwise distance matrix.
* You need a fast initial phylogenetic reconstruction or a starting tree for heuristic searches.
* Evolutionary rates vary across lineages (violating the molecular clock).
* You are exploring relationships in large datasets before running computationally intensive Maximum Likelihood or Bayesian analyses.

---

## Bootstrap Support for Neighbor Joining Trees

A phylogenetic tree should not be interpreted solely from its branching topology. **Bootstrap analysis** estimates the stability of inferred branches by repeatedly resampling alignment columns.

### Standard Bootstrap Procedure

1. Start with a multiple sequence alignment containing $L$ columns.
2. Randomly sample $L$ columns **with replacement** to create a pseudoreplicate alignment.
3. Reconstruct an NJ tree from the resampled alignment.
4. Repeat steps 2–3 many times (typically $100$ to $1,000+$ replicates).
5. Calculate the percentage of replicate trees in which each branch (clade/bipartition) appears.

For example, if a clade appears in $950$ out of $1,000$ replicate trees:

$$
\text{Bootstrap Support} = \frac{950}{1000} \times 100 = 95\%
$$

> **Interpreting Bootstrap Values:** A $95\%$ bootstrap value indicates that the branch was recovered in $95\%$ of resampled datasets under the chosen distance and clustering procedure. It is **not** a $95\%$ probability that the clade is biologically true; it measures branch stability against sampling variation.

---

## Applications in Bioinformatics

1. **Comparative genomics:** Reconstructing evolutionary relationships across species using orthologous genes or proteins.
2. **Protein evolution:** Analyzing evolutionary divergence among protein families, superfamilies, and functional domains.
3. **Microbial phylogenetics:** Typing and tracing relationships among bacterial or viral isolates using sequence-derived distances.
4. **Molecular evolution:** Studying patterns of sequence divergence and rate variation across gene families.
5. **Initial phylogenetic exploration:** Generating rapid baseline trees or starting topologies prior to Maximum Likelihood or Bayesian inference.

---

## Key Takeaways

* Neighbor Joining reconstructs an **unrooted phylogenetic tree** from a symmetric pairwise distance matrix.
* It does **not** assume a constant evolutionary rate (no molecular clock required).
* It identifies neighboring taxa using a corrected-distance criterion called the **Q-matrix**, which balances pairwise distance against each taxon's total distance to all other taxa.
* Individual branch lengths are estimated at each joining step and can differ between sister lineages.
* The canonical algorithm has **$O(n^3)$ time complexity** and **$O(n^2)$ space complexity**.
* **Bootstrap resampling** is used to evaluate the reproducibility and stability of inferred branches.
* While fast and effective for exploratory analysis, NJ discards site-level sequence information and is often complemented by Maximum Likelihood or Bayesian methods for high-confidence inference.

---

## References

1. **Saitou, N., & Nei, M. (1987).** The neighbor-joining method: A new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406–425. <https://doi.org/10.1093/oxfordjournals.molbev.a040454>
2. **Studier, J. A., & Keppler, K. J. (1988).** A note on the neighbor-joining algorithm of Saitou and Nei. *Molecular Biology and Evolution*, 5(6), 729–731.
3. **Felsenstein, J. (2004).** *Inferring Phylogenies*. Sinauer Associates.
4. **Biopython Documentation:** <https://biopython.org/>