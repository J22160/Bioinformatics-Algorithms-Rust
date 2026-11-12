# Small Parsimony — Fitch's and Sankoff's Algorithms

## Table of Contents

* [What Is the Small Parsimony Problem?](#what-is-the-small-parsimony-problem)
* [What Is an Evolutionary Character?](#what-is-an-evolutionary-character)
* [Two Algorithms for Small Parsimony](#two-algorithms-for-small-parsimony)
* [1. Fitch's Algorithm — Unweighted Parsimony](#1-fitchs-algorithm--unweighted-parsimony)
  * [What Is Fitch's Algorithm?](#what-is-fitchs-algorithm)
  * [Core Idea](#core-idea)
  * [Fitch's Algorithm: Step-by-Step](#fitchs-algorithm-step-by-step)
  * [Worked Example: Fitch's Algorithm](#worked-example-fitchs-algorithm)
* [2. Sankoff's Algorithm — Weighted Parsimony](#2-sankoffs-algorithm--weighted-parsimony)
  * [What Is Sankoff's Algorithm?](#what-is-sankoffs-algorithm)
  * [Core Idea](#core-idea-1)
  * [Weighted Substitution Cost Matrix](#weighted-substitution-cost-matrix)
  * [Worked Example: Sankoff's Algorithm](#worked-example-sankoffs-algorithm)
* [3. Fitch's vs. Sankoff's Algorithm: Key Differences](#3-fitchs-vs-sankoffs-algorithm-key-differences)
* [4. Backtracking: Recovering Ancestral Character States](#4-backtracking-recovering-ancestral-character-states)
* [5. Extending Small Parsimony to Multiple Sequence Alignments](#5-extending-small-parsimony-to-multiple-sequence-alignments)
* [6. Python Implementation](#6-python-implementation)
* [7. Computational Complexity](#7-computational-complexity)
* [8. Advantages and Limitations](#8-advantages-and-limitations)
* [9. Small Parsimony vs. Large Parsimony](#9-small-parsimony-vs-large-parsimony)
* [10. Applications in Bioinformatics](#10-applications-in-bioinformatics)
* [11. When Should You Use Fitch's or Sankoff's Algorithm?](#11-when-should-you-use-fitchs-or-sankoffs-algorithm)
* [Key Takeaways](#key-takeaways)
* [References](#references)

---

## What Is the Small Parsimony Problem?

The **Small Parsimony problem** is a computational problem in phylogenetics that aims to reconstruct ancestral character states at the internal nodes of a known phylogenetic tree.

Given:
* A phylogenetic tree with known topology.
* Character states assigned to the leaves (observed taxa).
* A cost associated with changes between character states.

The objective is to assign character states to the internal nodes such that the **total evolutionary cost across all tree edges is minimized**.

In molecular evolution, a character can represent a nucleotide at a particular alignment position:

$$
\{A, C, G, T\}
$$

For example, consider four species with the following observed nucleotides at one aligned position:

| Species | Observed Nucleotide |
| :--- | :---: |
| **Human** | A |
| **Chimpanzee** | A |
| **Gorilla** | G |
| **Orangutan** | G |

The tree topology is assumed to be known:

```text
             Root
            /    \
           X      Y
          / \    / \
      Human Chimp Gorilla Orangutan
```

The objective is to determine the most parsimonious nucleotide assignment at internal nodes $X$, $Y$, and the root. A parsimonious reconstruction minimizes the number of evolutionary changes required to explain the observed data.

> **Important:** Parsimony does not necessarily identify the true historical ancestral sequence. It identifies the reconstruction with the minimum cost under the chosen tree and character-change model.

---

## What Is an Evolutionary Character?

A **character** is an inherited feature that can take different states across taxa. In molecular phylogenetics, examples include:

* A nucleotide at a particular DNA alignment position.
* An amino acid at a protein alignment position.
* An insertion or deletion state.
* A binary presence/absence character.

For nucleotide reconstruction, the possible character states are:

$$
\Sigma = \{A, C, G, T\}
$$

An evolutionary change occurs when the state at one end of a tree edge differs from the state at the other end:

```text
A ---- A   (Cost = 0 changes)
A ---- G   (Cost = 1 change under an unweighted model)
```

The **total parsimony score** is the sum of the costs across all edges in the tree.

---

## Two Algorithms for Small Parsimony

Two commonly used algorithms solve the Small Parsimony problem:

1. **Fitch's Algorithm:** Assumes all character substitutions have equal cost.
2. **Sankoff's Algorithm:** Supports different costs for different substitutions.

Both use dynamic programming and solve the Small Parsimony problem optimally for a fixed tree topology and a specified cost model.

| Feature | Fitch's Algorithm | Sankoff's Algorithm |
| :--- | :--- | :--- |
| **Substitution costs** | Equal | Can be different |
| **Character states** | Usually unordered | Arbitrary finite state set |
| **Main operation** | Set intersection and union | Minimum-cost dynamic programming |
| **Data structure** | Sets of possible states | Cost vectors |
| **Complexity per node** | $O(1)$ for a fixed four-state alphabet | $O(k^2)$ for $k$ states |
| **Typical application** | Standard nucleotide parsimony | Weighted or custom evolutionary models |

---

## 1. Fitch's Algorithm — Unweighted Parsimony

### What Is Fitch's Algorithm?

**Fitch's algorithm** is a bottom-up dynamic programming algorithm for finding the minimum number of character-state changes on a fixed phylogenetic tree, introduced by Walter M. Fitch in 1971.

The algorithm assumes that:
* All substitutions between different character states have equal cost.
* Keeping the same state across an edge has zero cost.
* The tree topology is already known.
* Each leaf has an observed character state.

For nucleotide data ($\Sigma = \{A, C, G, T\}$), the cost of a transition is:

$$
\delta(i, j) =
\begin{cases}
0 & \text{if } i = j \\
1 & \text{if } i \neq j
\end{cases}
$$

This means that an $A \leftrightarrow G$ substitution costs the same as an $A \leftrightarrow T$ substitution.

### Core Idea

Fitch's algorithm processes the tree from the leaves toward the root. At each internal node, it determines which character states can be assigned to that node while minimizing the number of changes in its descendant subtree.

For a binary tree, let:
* $u$ = left child of node $v$
* $w$ = right child of node $v$
* $S(u)$ = set of possible states at child $u$
* $S(w)$ = set of possible states at child $w$

There are two cases:

#### Case 1: The Child State Sets Overlap

If the intersection is non-empty ($S(u) \cap S(w) \neq \varnothing$):

$$
S(v) = S(u) \cap S(w)
$$

No additional substitution is required at this node (**additional parsimony cost = $0$**).

* **Example:** If $S(u) = \{A, G\}$ and $S(w) = \{A, T\}$, their intersection is $S(v) = \{A\}$. The algorithm assigns the shared possible state set $\{A\}$ to the parent with $0$ additional cost.

#### Case 2: The Child State Sets Do Not Overlap

If the intersection is empty ($S(u) \cap S(w) = \varnothing$):

$$
S(v) = S(u) \cup S(w)
$$

The parsimony score increases by one (**additional parsimony cost = $1$**).

* **Example:** If $S(u) = \{A\}$ and $S(w) = \{G\}$, their intersection is empty, so $S(v) = \{A, G\}$ with an additional parsimony cost of $1$. The parent can be assigned either $A$ or $G$ during subsequent reconstruction with the same minimum score.

### Fitch's Algorithm: Step-by-Step

1. Assign each leaf its observed character state as a singleton set.
2. Traverse the tree from the leaves toward the root.
3. At each internal node, compare the state sets of its children.
4. If the sets overlap, assign their intersection to the parent.
5. If the sets do not overlap, assign their union and increment the parsimony score by $1$.
6. Continue until the root is reached.
7. Perform a top-down traversal to assign specific states to internal nodes.

The total number of union operations equals the minimum parsimony score for a binary tree under the unweighted model.

> **Important:** The bottom-up pass calculates the minimum score and possible states. It does not necessarily determine one unique ancestral reconstruction; a top-down pass is needed to select a consistent assignment.

### Worked Example: Fitch's Algorithm

#### Scenario

Consider four primate species with the following nucleotide states at one aligned position:

| Species | Nucleotide |
| :--- | :---: |
| **Human** | A |
| **Chimpanzee** | A |
| **Gorilla** | G |
| **Orangutan** | G |

Assume the tree topology is:

```text
                 Root
                /    \
               X      Y
              / \    / \
          Human Chimp Gorilla Orangutan
```

We want to reconstruct the ancestral nucleotide states at nodes $X$, $Y$, and $\text{Root}$.

#### Step 1: Initialize the Leaf State Sets

$$
S(\text{Human}) = \{A\}, \quad S(\text{Chimp}) = \{A\}, \quad S(\text{Gorilla}) = \{G\}, \quad S(\text{Orangutan}) = \{G\}
$$

#### Step 2: Process Node $X$

The children of $X$ are Human and Chimpanzee:

$$
S(X) = S(\text{Human}) \cap S(\text{Chimp}) = \{A\} \cap \{A\} = \{A\}
$$

No substitution is required (**score contribution = $0$**).

#### Step 3: Process Node $Y$

The children of $Y$ are Gorilla and Orangutan:

$$
S(Y) = S(\text{Gorilla}) \cap S(\text{Orangutan}) = \{G\} \cap \{G\} = \{G\}
$$

No substitution is required (**score contribution = $0$**).

#### Step 4: Process the Root

The children of the root are $X$ and $Y$. Since their intersection is empty ($\{A\} \cap \{G\} = \varnothing$):

$$
S(\text{Root}) = S(X) \cup S(Y) = \{A, G\}
$$

The algorithm increments the score by one (**score contribution = $1$**).

#### Step 5: Calculate the Minimum Parsimony Score

$$
\text{Score} = 0 + 0 + 1 = 1
$$

**Minimum parsimony score = $1$**

#### Step 6: Reconstruct Ancestral States Using a Top-Down Pass

Since $S(\text{Root}) = \{A, G\}$, we can choose either state as the root assignment:

* **Reconstruction 1 ($\text{Root} = A$):** Node $X$ must be $A$, and node $Y$ is $G$, requiring one substitution on the $\text{Root} \to Y$ edge.

```text
                 A
                / \
               A   G
              / \ / \
             A  A G  G
```

* **Reconstruction 2 ($\text{Root} = G$):** Node $Y$ must be $G$, and node $X$ is $A$, requiring one substitution on the $\text{Root} \to X$ edge.

```text
                 G
                / \
               A   G
              / \ / \
             A  A G  G
```

Both reconstructions yield a parsimony score of $1$. Therefore, the ancestral state at the root is **ambiguous** under this model—demonstrating that Fitch's algorithm can find the minimum score without uniquely determining every ancestral nucleotide.

---

## 2. Sankoff's Algorithm — Weighted Parsimony

### What Is Sankoff's Algorithm?

**Sankoff's algorithm** is a generalized dynamic programming approach to the Small Parsimony problem that allows different evolutionary changes to have different costs.

For example, nucleotide substitutions can be weighted according to biological properties:
* **Transitions:** $A \leftrightarrow G$ and $C \leftrightarrow T$ (often more frequent).
* **Transversions:** $A \leftrightarrow C$, $A \leftrightarrow T$, $G \leftrightarrow C$, and $G \leftrightarrow T$ (often less frequent).

A weighted parsimony model can therefore assign:

$$
\text{Cost}(\text{Transition}) = 1, \quad \text{Cost}(\text{Transversion}) = 2
$$

This allows the algorithm to favor evolutionary reconstructions involving less expensive substitutions.

> **Important:** These costs are parsimony penalties, not automatically empirical substitution probabilities or evolutionary rates.

### Core Idea

Instead of maintaining a set of possible states at each internal node, Sankoff's algorithm maintains a **cost vector**. For every node $v$ and possible character state $t$, define:

$$
s_t(v)
$$

as the minimum cost of explaining all observed characters in the subtree rooted at $v$, conditional on node $v$ having state $t$.

For a leaf $v$ with observed nucleotide $x$:

$$
s_t(v) =
\begin{cases}
0 & \text{if } t = x \\
\infty & \text{if } t \neq x
\end{cases}
$$

This ensures that a leaf retains its observed state. For an internal node $v$ with children $u$ and $w$:

$$
s_t(v) = \min_i \left[ s_i(u) + \delta(i, t) \right] + \min_j \left[ s_j(w) + \delta(j, t) \right]
$$

Where:
* $t$ is the state assigned to the current node $v$.
* $i$ is a possible state at child $u$.
* $j$ is a possible state at child $w$.
* $\delta(i, t)$ is the cost of changing from state $i$ to state $t$.

For a general tree with multiple children:

$$
s_t(v) = \sum_{u \in \text{Children}(v)} \min_i \left[ s_i(u) + \delta(i, t) \right]
$$

The minimum parsimony score is obtained at the root:

$$
\text{Score} = \min_t s_t(\text{Root})
$$

The algorithm then uses backtracking to identify the ancestral states that achieve this minimum.

### Weighted Substitution Cost Matrix

A substitution cost matrix defines the penalty of changing from one nucleotide to another:

| From / To | A | C | G | T |
| :--- | :---: | :---: | :---: | :---: |
| **A** | 0 | 2 | 1 | 2 |
| **C** | 2 | 0 | 2 | 1 |
| **G** | 1 | 2 | 0 | 2 |
| **T** | 2 | 1 | 2 | 0 |

This matrix assigns cost $0$ for no substitution, cost $1$ for transitions ($\delta(A, G) = 1$, $\delta(C, T) = 1$), and cost $2$ for transversions ($\delta(A, C) = 2$, $\delta(G, T) = 2$). While symmetric here ($\delta(i, j) = \delta(j, i)$), Sankoff's algorithm also supports asymmetric cost matrices when required.

### Worked Example: Sankoff's Algorithm

#### Scenario

Consider two observed nucleotides at the leaves of a small tree:

```text
          Root
          /  \
         A    G
```

We want to determine the minimum-cost ancestral nucleotide at the root using the weighted substitution matrix above.

#### Step 1: Initialize the Leaf Cost Vectors

Ordering vector entries as $[A, C, G, T]$:

$$
s(A) = [0, \infty, \infty, \infty]
$$

$$
s(G) = [\infty, \infty, 0, \infty]
$$

#### Step 2: Calculate the Cost of Each Possible Root State

* **Root = A:** $s_A(\text{Root}) = \delta(A, A) + \delta(G, A) = 0 + 1 = 1$
* **Root = C:** $s_C(\text{Root}) = \delta(A, C) + \delta(G, C) = 2 + 2 = 4$
* **Root = G:** $s_G(\text{Root}) = \delta(A, G) + \delta(G, G) = 1 + 0 = 1$
* **Root = T:** $s_T(\text{Root}) = \delta(A, T) + \delta(G, T) = 2 + 2 = 4$

| Root State | Cost to Left Leaf (A) | Cost to Right Leaf (G) | Total Cost |
| :---: | :---: | :---: | :---: |
| **A** | 0 | 1 | **1** |
| **C** | 2 | 2 | **4** |
| **G** | 1 | 0 | **1** |
| **T** | 2 | 2 | **4** |

#### Step 3: Determine the Minimum Score

$$
\text{Score} = \min(1, 4, 1, 4) = 1
$$

The optimal ancestral root states are $\{A, G\}$, both yielding the same minimum cost of $1$.

---

## 3. Fitch's vs. Sankoff's Algorithm: Key Differences

The main difference lies in how each algorithm represents possible ancestral states and computes evolutionary costs:

* **Fitch's Algorithm:**
  * Stores sets of possible states at each node.
  * Uses set intersection and union operations.
  * Adds $1$ to the parsimony score when child sets are disjoint.
  * Assumes all substitutions have equal cost.
  * Computationally fast for small, fixed alphabets.
* **Sankoff's Algorithm:**
  * Stores a minimum cost vector across every possible state at every node.
  * Uses dynamic programming and minimum-cost transitions.
  * Supports weighted and asymmetric substitutions.
  * Supports arbitrary finite character-state alphabets.
  * Requires more computation per node than Fitch's algorithm for general cost matrices.

*Note: Fitch's algorithm is a special case of Sankoff's algorithm where $\delta(i, j) = 1$ for all $i \neq j$ and $\delta(i, i) = 0$.*

---

## 4. Backtracking: Recovering Ancestral Character States

The bottom-up dynamic programming pass calculates the minimum score and possible/optimal states at each node. To produce a complete ancestral reconstruction, a top-down pass is required.

### Fitch Backtracking

1. Select one state from the root's possible state set.
2. Traverse from the root toward the leaves.
3. At each child node, retain the parent's assigned state if it is present in the child's state set; otherwise, select a state from the child's set.
4. Repeat until all internal nodes have been assigned.

### Sankoff Backtracking

For a child node $u$ whose parent has already been assigned state $t$, choose a child state $i$ that minimizes:

$$
s_i(u) + \delta(i, t)
$$

This selects a state consistent with the optimal cost computed during the bottom-up pass. If multiple states achieve the same minimum, the reconstruction is ambiguous and multiple valid optimal assignments exist.

> **Important:** Choosing one state arbitrarily when ties occur is acceptable for producing a single optimal reconstruction, but it should not be interpreted as evidence that the selected state was historically more likely.

---

## 5. Extending Small Parsimony to Multiple Sequence Alignments

For DNA or protein sequence data, ancestral sequences are reconstructed by applying the algorithm independently to each alignment column.

### Example Alignment

| Species | Position 1 | Position 2 | Position 3 | Position 4 |
| :--- | :---: | :---: | :---: | :---: |
| **Human** | A | C | G | T |
| **Chimp** | A | T | G | T |
| **Gorilla** | G | C | G | A |
| **Orangutan** | G | T | A | A |

For each alignment column:
1. Extract the observed character states across taxa.
2. Apply Fitch's or Sankoff's algorithm using the known tree topology.
3. Calculate the minimum parsimony score for that site.
4. Reconstruct possible ancestral states at every internal node.
5. Concatenate the reconstructed states across all columns to form ancestral sequences:

```text
Ancestral node X: A   C   G   T
Ancestral node Y: G   T   A   A
Root:             A/T C/T G/A T/A
```

*(Here, slash notation at the root illustrates ambiguous possible states at individual sites rather than a single unique sequence.)*

### Gaps and Ambiguous Nucleotides

Real sequence alignments frequently contain gaps (`-`), IUPAC ambiguity codes (`N`, `R`, `Y`, `W`), and missing data (`?`), which require explicit handling:
* **Ambiguous nucleotides:** Can be initialized as a set of compatible states (e.g., `R` $= \{A, G\}$).
* **Missing data (`?`):** Often treated as compatible with every state in $\Sigma$.
* **Gaps (`-`):** Can be treated as a fifth character state, as missing data, or modeled under specialized indel penalties—each choice can alter the resulting ancestral reconstruction.

---

## 6. Python Implementation

The following examples demonstrate the core bottom-up dynamic programming logic for a rooted binary tree.

### Fitch's Algorithm

```python
def fitch(left, right):
    """
    Perform one Fitch operation on two child state sets.

    Returns:
        parent_states (set): Possible states at the parent node.
        score (int): Additional parsimony cost incurred at this node.
    """
    intersection = left & right

    if intersection:
        return intersection, 0

    return left | right, 1


# Example: four observed nucleotide states
human = {"A"}
chimp = {"A"}
gorilla = {"G"}
orangutan = {"G"}

# Bottom-up traversal
node_x, score_x = fitch(human, chimp)
node_y, score_y = fitch(gorilla, orangutan)
root, score_root = fitch(node_x, node_y)

total_score = score_x + score_y + score_root

print("Node X:", node_x)
print("Node Y:", node_y)
print("Root:", root)
print("Parsimony score:", total_score)
```

**Expected Output:**

```text
Node X: {'A'}
Node Y: {'G'}
Root: {'A', 'G'}
Parsimony score: 1
```

### Sankoff's Algorithm

```python
import math

states = ["A", "C", "G", "T"]

# Weighted substitution cost matrix
cost = {
    "A": {"A": 0, "C": 2, "G": 1, "T": 2},
    "C": {"A": 2, "C": 0, "G": 2, "T": 1},
    "G": {"A": 1, "C": 2, "G": 0, "T": 2},
    "T": {"A": 2, "C": 1, "G": 2, "T": 0},
}


def leaf_cost(observed):
    """Initialize cost vector for an observed leaf."""
    return {
        state: 0 if state == observed else math.inf
        for state in states
    }


def sankoff(left, right):
    """
    Calculate the minimum cost for each possible
    parent state using two child cost vectors.
    """
    parent = {}

    for parent_state in states:
        left_cost = min(
            left[child_state] + cost[child_state][parent_state]
            for child_state in states
        )
        right_cost = min(
            right[child_state] + cost[child_state][parent_state]
            for child_state in states
        )
        parent[parent_state] = left_cost + right_cost

    return parent


# Example: observed A and G
left = leaf_cost("A")
right = leaf_cost("G")

root_costs = sankoff(left, right)
minimum_score = min(root_costs.values())
optimal_root_states = [
    state
    for state, value in root_costs.items()
    if value == minimum_score
]

print("Root costs:", root_costs)
print("Minimum score:", minimum_score)
print("Optimal root states:", optimal_root_states)
```

**Expected Output:**

```text
Root costs: {'A': 1, 'C': 4, 'G': 1, 'T': 4}
Minimum score: 1
Optimal root states: ['A', 'G']
```

---

## 7. Computational Complexity

Let:
* $n$ = number of leaves (taxa)
* $k$ = number of possible character states
* $V = 2n - 1$ = number of nodes in a binary tree
* $E = 2n - 2$ = number of edges in a binary tree

### Fitch's Algorithm

At each internal node, the algorithm performs set intersection or union operations over an alphabet of size $k$:
* **Time Complexity:** $O(nk)$ per column ($O(n)$ for a fixed four-state nucleotide alphabet).
* **Space Complexity:** $O(Vk) = O(nk)$ to store state sets for backtracking.

### Sankoff's Algorithm

For each internal node and each possible parent state ($k$), the algorithm evaluates all possible states ($k$) at each child:
* **Time Complexity:** $O(Ek^2) = O(nk^2)$ per column.
* **Space Complexity:** $O(Vk) = O(nk)$ to store cost vectors at all nodes.

### Multiple Sequence Alignment ($L$ Columns)

| Algorithm | Time Complexity | Space Complexity |
| :--- | :---: | :---: |
| **Fitch's Algorithm** | $O(Lnk)$ | $O(nk)$ per column |
| **Sankoff's Algorithm** | $O(Lnk^2)$ | $O(nk)$ per column |

*Note: These complexities assume a fixed tree topology. Searching across alternative tree topologies is the separate, NP-hard Large Parsimony problem.*

---

## 8. Advantages and Limitations

### Advantages

* **Computationally efficient:** Uses dynamic programming to find the exact minimum parsimony score on a fixed tree in linear time with respect to the number of taxa.
* **Easy to interpret:** Directly minimizes the total number or weighted cost of evolutionary changes.
* **Ancestral reconstruction:** Infers possible ancestral states at all internal nodes.
* **Flexible cost models:** Sankoff's algorithm accommodates custom transition/transversion or step matrices.
* **Broadly applicable:** Works across DNA, RNA, protein, and discrete morphological characters.

### Limitations

1. **Not a probabilistic evolutionary model:** The minimum-cost reconstruction is not guaranteed to be the most probable historical sequence, as parsimony does not explicitly account for branch lengths or multiple substitutions along long branches.
2. **Long-branch attraction (LBA):** Rapidly evolving lineages can independently accumulate identical character states (homoplasy), causing parsimony to underestimate changes or misinterpret relationships.
3. **Depends on a known input tree:** Small parsimony evaluates a fixed topology; it does not infer the tree itself.
4. **Multiple optimal reconstructions:** Ties frequently occur, leaving ancestral states at certain nodes ambiguous.
5. **Model dependence:** Sankoff's output depends directly on the user-defined substitution cost matrix.
6. **Sensitive to alignment errors:** Misaligned columns introduce false homologies and distort ancestral state inference.

---

## 9. Small Parsimony vs. Large Parsimony

| Feature | Small Parsimony | Large Parsimony |
| :--- | :--- | :--- |
| **Tree topology** | Given (fixed) | Unknown (must be inferred) |
| **Objective** | Find ancestral states minimizing cost | Find both tree topology and ancestral states minimizing cost |
| **Computational difficulty** | Polynomial-time dynamic programming | NP-hard combinatorial tree search |
| **Algorithms** | Fitch's and Sankoff's algorithms | Heuristic search (NNI, SPR, TBR) scored via Small Parsimony |

---

## 10. Applications in Bioinformatics

1. **Ancestral sequence reconstruction:** Inferring ancient nucleotide or amino acid sequences at internal nodes of a phylogeny.
2. **Comparative genomics:** Tracing the evolutionary history of genes, regulatory elements, and conserved regions across species.
3. **Protein evolution:** Studying historical amino acid replacements to guide ancestral protein resurrection experiments.
4. **Molecular evolution:** Quantifying the minimum number of substitutions required to explain observed sequence variation.
5. **Viral evolution:** Mapping mutations across known viral transmission trees (often paired with model-based methods for rapidly evolving pathogens).
6. **Morphological evolution:** Reconstructing ancestral anatomical, behavioral, or developmental traits using discrete character matrices.

---

## 11. When Should You Use Fitch's or Sankoff's Algorithm?

* **Use Fitch's algorithm when:**
  * All character-state changes should carry equal cost.
  * You are analyzing standard unordered nucleotide or binary morphological characters.
  * You need fast ancestral-state scoring and reconstruction.
* **Use Sankoff's algorithm when:**
  * Different substitutions require different penalties (e.g., transitions vs. transversions, or amino acid step matrices).
  * Character transitions follow custom or asymmetric rules.
* **Use Maximum Likelihood or Bayesian methods when:**
  * You need explicit probabilistic models incorporating branch lengths, rate heterogeneity, and statistical uncertainty (posterior probabilities) for ancestral states.

---

## Key Takeaways

* The **Small Parsimony problem** reconstructs ancestral character states on a known phylogenetic tree while minimizing total evolutionary change.
* **Fitch's algorithm** uses bottom-up set intersections and unions under an equal-cost model.
* **Sankoff's algorithm** uses dynamic programming cost vectors to support weighted substitution matrices.
* Both algorithms guarantee the optimal parsimony score for a fixed tree and cost model, followed by a top-down **backtracking pass** to assign internal states.
* Multiple equally optimal ancestral reconstructions may exist at ambiguous nodes.
* For a fixed binary tree, Fitch's algorithm runs in $O(nk)$ time per alignment column ($O(n)$ for DNA), whereas Sankoff's algorithm runs in $O(nk^2)$ time.

---

## References

1. **Fitch, W. M. (1971).** Toward defining the course of evolution: Minimum change for a specific tree topology. *Systematic Zoology*, 20(4), 406–416. <https://doi.org/10.2307/2412116>
2. **Sankoff, D. (1975).** Minimal mutation trees of sequences. *SIAM Journal on Applied Mathematics*, 28(1), 35–42. <https://doi.org/10.1137/0128004>
3. **Felsenstein, J. (2004).** *Inferring Phylogenies*. Sinauer Associates.
4. **Yang, Z. (2014).** *Molecular Evolution: A Statistical Approach*. Oxford University Press.