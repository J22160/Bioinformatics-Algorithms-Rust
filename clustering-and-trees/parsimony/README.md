# Small Parsimony — Fitch's and Sankoff's Algorithms

## What Is the Small Parsimony Problem?

The **Small Parsimony problem** is a computational problem in phylogenetics that aims to reconstruct ancestral character states at the internal nodes of a known phylogenetic tree.

Given:
- A phylogenetic tree with known topology.
- Character states assigned to the leaves (observed taxa).
- A cost associated with changes between character states.

The objective is to assign character states to the internal nodes such that the **total evolutionary cost across all tree edges is minimized**.

In molecular evolution, a character can represent a nucleotide at a particular alignment position:

\[
\{A,C,G,T\}
\]

For example, consider four species with the following observed nucleotides at one aligned position:

| Species | Observed nucleotide |
|---|---|
| Human | A |
| Chimpanzee | A |
| Gorilla | G |
| Orangutan | G |

The tree topology is assumed to be known:

```text
              Root
             /    \
            X      Y
           / \    / \
        Human Chimp Gorilla Orangutan
```

The objective is to determine the most parsimonious nucleotide assignment at internal nodes \(X\), \(Y\), and the root.

A parsimonious reconstruction minimizes the number of evolutionary changes required to explain the observed data.

**Important:** Parsimony does not necessarily identify the true historical ancestral sequence. It identifies the reconstruction with the minimum cost under the chosen tree and character-change model.

---

## What Is an Evolutionary Character?

A character is an inherited feature that can take different states across taxa.

In molecular phylogenetics, examples include:

- A nucleotide at a particular DNA alignment position.
- An amino acid at a protein alignment position.
- An insertion or deletion state.
- A binary presence/absence character.

For nucleotide reconstruction, the possible character states are:

\[
\Sigma=\{A,C,G,T\}
\]

An evolutionary change occurs when the state at one end of a tree edge differs from the state at the other end.

For example:

```text
A ---- A
```

Cost = 0 changes.

```text
A ---- G
```

Cost = 1 change under an unweighted model.

The total parsimony score is the sum of the costs across all edges in the tree.

---

## Two Algorithms for Small Parsimony

Two commonly used algorithms solve the Small Parsimony problem:

1. **Fitch's Algorithm:** Assumes all character substitutions have equal cost.
2. **Sankoff's Algorithm:** Supports different costs for different substitutions.

Both use dynamic programming and solve the small parsimony problem optimally for a fixed tree topology and a specified cost model.

| Feature | Fitch's Algorithm | Sankoff's Algorithm |
|---|---|---|
| Substitution costs | Equal | Can be different |
| Character states | Usually unordered | Arbitrary finite state set |
| Main operation | Set intersection and union | Minimum-cost dynamic programming |
| Data structure | Sets of possible states | Cost vectors |
| Complexity per node | \(O(1)\) for a fixed four-state alphabet | \(O(k^2)\) for \(k\) states |
| Typical application | Standard nucleotide parsimony | Weighted or custom evolutionary models |

---

# 1. Fitch's Algorithm — Unweighted Parsimony

## What Is Fitch's Algorithm?

Fitch's algorithm is a bottom-up dynamic programming algorithm for finding the minimum number of character-state changes on a fixed phylogenetic tree.

It was introduced by Walter M. Fitch in 1971.

The algorithm assumes that:

- All substitutions between different character states have equal cost.
- Keeping the same state across an edge has zero cost.
- The tree topology is already known.
- Each leaf has an observed character state.

For nucleotide data:

\[
\Sigma=\{A,C,G,T\}
\]

The cost of a transition is:

\[
\delta(i,j)=
\begin{cases}
0 & \text{if } i=j\\
1 & \text{if } i\neq j
\end{cases}
\]

This means that an A-to-G substitution costs the same as an A-to-T substitution.

## Core Idea

Fitch's algorithm processes the tree from the leaves toward the root.

At each internal node, it determines which character states can be assigned to that node while minimizing the number of changes in its descendant subtree.

For a binary tree, let:

- \(u\) = left child of node \(v\).
- \(w\) = right child of node \(v\).
- \(S(u)\) = set of possible states at child \(u\).
- \(S(w)\) = set of possible states at child \(w\).

There are two cases.

### Case 1: The child state sets overlap

If the intersection is non-empty:

\[
S(u)\cap S(w)\neq\varnothing
\]

Then:

\[
S(v)=S(u)\cap S(w)
\]

No additional substitution is required at this node.

**Example:**

\[
S(u)=\{A,G\}
\]

\[
S(w)=\{A,T\}
\]

Their intersection is:

\[
S(v)=\{A\}
\]

The algorithm assigns the shared possible state set \(\{A\}\) to the parent.

Additional parsimony cost:

\[
0
\]

### Case 2: The child state sets do not overlap

If:

\[
S(u)\cap S(w)=\varnothing
\]

Then:

\[
S(v)=S(u)\cup S(w)
\]

The parsimony score increases by one.

**Example:**

\[
S(u)=\{A\}
\]

\[
S(w)=\{G\}
\]

Their intersection is empty, so:

\[
S(v)=\{A,G\}
\]

Additional parsimony cost:

\[
1
\]

The parent can be assigned either A or G during the subsequent reconstruction, with the same minimum score.

## Fitch's Algorithm: Step-by-Step

1. Assign each leaf its observed character state as a singleton set.
2. Traverse the tree from the leaves toward the root.
3. At each internal node, compare the state sets of its children.
4. If the sets overlap, assign their intersection to the parent.
5. If the sets do not overlap, assign their union and increment the parsimony score by one.
6. Continue until the root is reached.
7. Perform a top-down traversal to assign specific states to internal nodes.

The total number of union operations is the minimum parsimony score for a binary tree under the unweighted model.

**Important:** The bottom-up pass calculates the minimum score and possible states. It does not necessarily determine one unique ancestral reconstruction. A top-down pass is needed to select a consistent assignment.

---

## Worked Example: Fitch's Algorithm

### Scenario

Consider four primate species with the following nucleotide states at one aligned position:

| Species | Nucleotide |
|---|---|
| Human | A |
| Chimpanzee | A |
| Gorilla | G |
| Orangutan | G |

Assume the tree topology is:

```text
                 Root
                /    \
               X      Y
              / \    / \
          Human Chimp Gorilla Orangutan
```

We want to reconstruct the ancestral nucleotide states at nodes X, Y, and Root.

### Step 1: Initialize the leaf state sets

Each leaf receives its observed nucleotide:

\[
S(Human)=\{A\}
\]

\[
S(Chimp)=\{A\}
\]

\[
S(Gorilla)=\{G\}
\]

\[
S(Orangutan)=\{G\}
\]

### Step 2: Process node X

The children of X are Human and Chimpanzee.

\[
S(Human)\cap S(Chimp)=\{A\}
\]

Therefore:

\[
S(X)=\{A\}
\]

No substitution is required.

Score contribution:

\[
0
\]

### Step 3: Process node Y

The children of Y are Gorilla and Orangutan.

\[
S(Gorilla)\cap S(Orangutan)=\{G\}
\]

Therefore:

\[
S(Y)=\{G\}
\]

No substitution is required.

Score contribution:

\[
0
\]

### Step 4: Process the root

The children of the root are X and Y.

\[
S(X)=\{A\}
\]

\[
S(Y)=\{G\}
\]

Their intersection is empty:

\[
\{A\}\cap\{G\}=\varnothing
\]

Therefore:

\[
S(Root)=\{A,G\}
\]

The algorithm increments the score by one.

Score contribution:

\[
1
\]

### Step 5: Calculate the minimum parsimony score

The total number of union operations is:

\[
Score=0+0+1=1
\]

<badge color="success">Minimum parsimony score = 1</badge>

### Step 6: Reconstruct ancestral states using a top-down pass

The root has two possible states:

\[
S(Root)=\{A,G\}
\]

We can choose either state as the root assignment.

**Reconstruction 1:**

Assign:

\[
Root=A
\]

Then:
- Node X must be A.
- Node Y can be G, requiring one substitution on the root-to-Y edge.

```text
                 A
                / \
               A   G
              / \ / \
             A  A G  G
```

**Reconstruction 2:**

Assign:

\[
Root=G
\]

Then:
- Node Y must be G.
- Node X can be A, requiring one substitution on the root-to-X edge.

```text
                 G
                / \
               A   G
              / \ / \
             A  A G  G
```

Both reconstructions have a parsimony score of one.

Therefore, the ancestral state at the root is **ambiguous** under this model.

This is an important distinction: Fitch's algorithm can identify the minimum score without uniquely identifying every ancestral nucleotide.

---

# 2. Sankoff's Algorithm — Weighted Parsimony

## What Is Sankoff's Algorithm?

Sankoff's algorithm is a generalized dynamic programming approach to the Small Parsimony problem.

Unlike Fitch's algorithm, it allows different evolutionary changes to have different costs.

For example, nucleotide substitutions can be assigned different costs depending on their biological or evolutionary properties.

A common distinction is between:

- **Transitions:** A ↔ G and C ↔ T.
- **Transversions:** A ↔ C, A ↔ T, G ↔ C, and G ↔ T.

Transitions are often more frequent than transversions, although their relative rates depend on the organism and evolutionary context.

A weighted parsimony model can therefore assign:

\[
Cost(Transition)=1
\]

\[
Cost(Transversion)=2
\]

This allows the algorithm to favor evolutionary reconstructions involving less expensive substitutions.

**Important:** These costs are parsimony penalties, not automatically empirical substitution probabilities or evolutionary rates.

## Core Idea

Instead of maintaining a set of possible states at each internal node, Sankoff's algorithm maintains a cost vector.

For every node \(v\) and possible character state \(t\), define:

\[
s_t(v)
\]

as the minimum cost of explaining all observed characters in the subtree rooted at \(v\), conditional on node \(v\) having state \(t\).

For a leaf \(v\) with observed nucleotide \(x\):

\[
s_t(v)=
\begin{cases}
0 & \text{if } t=x\\
\infty & \text{if } t\neq x
\end{cases}
\]

This ensures that a leaf retains its observed state.

For an internal node \(v\) with children \(u\) and \(w\):

\[
s_t(v)=
\min_i\left[s_i(u)+\delta(i,t)\right]
+
\min_j\left[s_j(w)+\delta(j,t)\right]
\]

where:

- \(t\) is the state assigned to the current node.
- \(i\) is a possible state at child \(u\).
- \(j\) is a possible state at child \(w\).
- \(\delta(i,t)\) is the cost of changing from state \(i\) to state \(t\).

For a general tree with multiple children:

\[
s_t(v)=
\sum_{u\in Children(v)}
\min_i\left[s_i(u)+\delta(i,t)\right]
\]

The minimum parsimony score is obtained at the root:

\[
Score=\min_t s_t(Root)
\]

The algorithm then uses backtracking to identify the ancestral states that achieve this minimum.

---

## Weighted Substitution Cost Matrix

For nucleotide data, a substitution cost matrix defines the cost of changing from one nucleotide to another.

Consider the following example:

| From / To | A | C | G | T |
|---|---:|---:|---:|---:|
| A | 0 | 2 | 1 | 2 |
| C | 2 | 0 | 2 | 1 |
| G | 1 | 2 | 0 | 2 |
| T | 2 | 1 | 2 | 0 |

This matrix assigns:

- Cost 0 for no substitution.
- Cost 1 for transitions.
- Cost 2 for transversions.

For example:

\[
\delta(A,G)=1
\]

\[
\delta(A,C)=2
\]

\[
\delta(C,T)=1
\]

\[
\delta(G,T)=2
\]

The matrix is symmetric in this example:

\[
\delta(i,j)=\delta(j,i)
\]

However, Sankoff's algorithm can also support asymmetric costs when the evolutionary model requires them.

---

## Worked Example: Sankoff's Algorithm

### Scenario

Consider two observed nucleotides at the leaves of a small tree:

```text
          Root
          /  \
         A    G
```

We want to determine the minimum-cost ancestral nucleotide at the root.

Using the weighted substitution matrix above:

| Root state | Cost to A | Cost to G | Total cost |
|---|---:|---:|---:|
| A | 0 | 1 | 1 |
| C | 2 | 2 | 4 |
| G | 1 | 0 | 1 |
| T | 2 | 2 | 4 |

### Step 1: Initialize the leaf cost vectors

For the left leaf, observed as A:

\[
s(A)=[0,\infty,\infty,\infty]
\]

For the right leaf, observed as G:

\[
s(G)=[\infty,\infty,0,\infty]
\]

The vector entries correspond to:

\[
[A,C,G,T]
\]

### Step 2: Calculate the cost of each possible root state

For a root assigned A:

\[
s_A(Root)=\delta(A,A)+\delta(G,A)
\]

\[
s_A(Root)=0+1=1
\]

For a root assigned C:

\[
s_C(Root)=\delta(A,C)+\delta(G,C)
\]

\[
s_C(Root)=2+2=4
\]

For a root assigned G:

\[
s_G(Root)=\delta(A,G)+\delta(G,G)
\]

\[
s_G(Root)=1+0=1
\]

For a root assigned T:

\[
s_T(Root)=\delta(A,T)+\delta(G,T)
\]

\[
s_T(Root)=2+2=4
\]

### Step 3: Determine the minimum score

\[
Score=\min(1,4,1,4)=1
\]

The possible ancestral root states are:

\[
\{A,G\}
\]

Both states produce the same minimum cost.

This example demonstrates that weighted parsimony can produce different results from a different substitution cost matrix, even when the tree topology and observed nucleotides remain unchanged.

---

# 3. Fitch's vs. Sankoff's Algorithm: Key Differences

The main difference is how each algorithm represents the possible ancestral states and calculates evolutionary costs.

### Fitch's algorithm

- Stores sets of possible states.
- Uses set intersection and union.
- Adds one to the score when child sets are disjoint.
- Assumes all substitutions have equal cost.
- Is computationally efficient for small fixed alphabets.

### Sankoff's algorithm

- Stores a minimum cost for every possible state at every node.
- Uses dynamic programming and minimum-cost calculations.
- Supports weighted substitutions.
- Supports arbitrary finite character-state alphabets.
- Requires more computation than Fitch's algorithm for a general cost matrix.

Fitch's algorithm can be viewed as a special case of Sankoff's algorithm when all changes between distinct states have equal cost.

---

# 4. Backtracking: Recovering Ancestral Character States

The bottom-up dynamic programming pass calculates the minimum score and the possible or optimal states at each node.

To produce a complete ancestral reconstruction, a top-down pass is required.

### Fitch backtracking

1. Select one state from the root's possible state set.
2. Move from the root toward the leaves.
3. At each child, retain the parent state if it is present in the child's state set.
4. Otherwise, choose a state from the child's set.
5. Repeat until all internal nodes have been assigned.

When multiple states are possible, different valid optimal reconstructions may exist.

### Sankoff backtracking

For a child node \(u\) whose parent has already been assigned state \(t\), choose a child state \(i\) minimizing:

\[
s_i(u)+\delta(i,t)
\]

This identifies a state consistent with the optimal cost calculated during the bottom-up pass.

If multiple states attain the same minimum, the reconstruction is ambiguous and more than one optimal assignment may exist.

**Important:** Choosing one state arbitrarily is acceptable for producing a single optimal reconstruction, but it should not be interpreted as evidence that the selected state was historically more likely.

---

# 5. Extending Small Parsimony to Multiple Sequence Alignments

The algorithms above describe reconstruction at a single character position.

For DNA or protein sequence data, ancestral sequences are reconstructed by applying the algorithm independently to each alignment column.

### Example

Suppose the following multiple sequence alignment is available:

| Species | Position 1 | Position 2 | Position 3 | Position 4 |
|---|---|---|---|---|
| Human | A | C | G | T |
| Chimp | A | T | G | T |
| Gorilla | G | C | G | A |
| Orangutan | G | T | A | A |

For each alignment column:

1. Extract the observed character states.
2. Apply Fitch's or Sankoff's algorithm using the known tree.
3. Calculate the minimum parsimony score.
4. Reconstruct possible ancestral states at every internal node.
5. Store the reconstructed states for that alignment position.

After processing all columns, concatenate the reconstructed states to obtain ancestral sequences.

For example:

```text
Ancestral node X: A C G T
Ancestral node Y: G T A A
Root:             A/T C/T G/A T/A
```

The root notation here illustrates possible states at individual positions, not a single uniquely inferred sequence. The exact states depend on the tree topology and the reconstruction at each column.

### Gaps and ambiguous nucleotides

Real sequence alignments may contain:

- Gaps (`-`).
- Ambiguous nucleotide codes such as N, R, Y, and W.
- Missing data (`?`).

These require an explicit treatment.

For example, an ambiguous nucleotide can be represented as a set of allowed states. A missing character may be treated as compatible with every state, but this depends on the implementation and analysis assumptions.

Gaps can be treated as a fifth state, as missing data, or under a more specialized indel model. These choices can change the inferred ancestral reconstruction.

---

# 6. Python Implementation

The following examples demonstrate the core dynamic programming logic for a rooted binary tree.

## Fitch's Algorithm

```python
def fitch(left, right):
    """
    Perform one Fitch operation on two child state sets.

    Returns:
        parent_states: possible states at the parent
        score: additional parsimony cost
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

Expected output:

```text
Node X: {'A'}
Node Y: {'G'}
Root: {'A', 'G'}
Parsimony score: 1
```

This implementation demonstrates the bottom-up scoring pass. A separate top-down traversal is needed to assign a consistent state to every internal node.

## Sankoff's Algorithm

```python
import math

states = ["A", "C", "G", "T"]

# Weighted substitution cost matrix
cost = {
    "A": {"A": 0, "C": 2, "G": 1, "T": 2},
    "C": {"A": 2, "C": 0, "G": 2, "T": 1},
    "G": {"A": 1, "C": 2, "G": 0, "T": 2},
    "T": {"A": 2, "C": 1, "G": 2, "T": 0}
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

Expected output:

```text
Root costs: {'A': 1, 'C': 4, 'G': 1, 'T': 4}
Minimum score: 1
Optimal root states: ['A', 'G']
```

This example calculates the optimal cost at the root. A complete implementation for a larger tree must also store or recompute the optimal child-state choices to perform top-down backtracking.

---

# 7. Computational Complexity

Let:

- \(n\) = number of leaves (taxa).
- \(k\) = number of possible character states.
- \(V\) = number of nodes in the tree.
- \(E\) = number of edges.

For a binary tree:

\[
V=2n-1
\]

and:

\[
E=2n-2
\]

### Fitch's Algorithm

At each internal node, the algorithm performs set intersection or union operations.

For a fixed alphabet, such as four nucleotide states, these operations take constant time.

Therefore:

**Time Complexity:**

\[
O(nk)
\]

For a fixed nucleotide alphabet:

\[
O(n)
\]

**Space Complexity:**

\[
O(Vk)=O(nk)
\]

This includes storing possible state sets at the nodes. A streaming or carefully implemented traversal can reduce auxiliary storage, although backtracking generally requires retaining or recomputing relevant information.

### Sankoff's Algorithm

For each internal node and each possible parent state, the algorithm considers all possible states at each child.

This requires \(O(k^2)\) operations per edge for a general cost matrix.

**Time Complexity:**

\[
O(Ek^2)=O(nk^2)
\]

**Space Complexity:**

\[
O(Vk)=O(nk)
\]

For a fixed four-state nucleotide alphabet, the computation is also linear in the number of taxa, although it has a larger constant factor than Fitch's algorithm.

### Multiple Sequence Alignment

For an alignment with:

- \(n\) taxa.
- \(L\) alignment columns.
- \(k\) character states.

The overall complexity becomes:

| Algorithm | Time complexity | Space complexity |
|---|---|---|
| Fitch | \(O(Lnk)\) | \(O(nk)\) per column |
| Sankoff | \(O(Lnk^2)\) | \(O(nk)\) per column |

These estimates assume a binary tree, a fixed cost matrix, and standard bottom-up dynamic programming.

The tree topology is assumed to be known. Searching across alternative tree topologies is a separate and generally much harder problem.

---

# 8. Advantages and Limitations

## Advantages

- **Computationally efficient:** Both algorithms use dynamic programming to find the minimum parsimony score on a fixed tree.
- **Easy to interpret:** The objective is to minimize the total number or cost of character changes.
- **Ancestral reconstruction:** Identifies possible ancestral states at internal nodes.
- **Flexible cost models:** Sankoff's algorithm supports different substitution penalties.
- **Useful for evolutionary studies:** Can be applied to DNA, RNA, protein, and morphological characters.

## Limitations

### 1. Parsimony is not a probabilistic evolutionary model

The minimum-cost reconstruction is not necessarily the most probable ancestral history.

A substitution that occurs frequently under a realistic evolutionary process may be treated as equally costly or more costly under a simple parsimony model.

### 2. Long-branch attraction

Parsimony methods can be vulnerable to long-branch attraction.

This occurs when rapidly evolving lineages independently accumulate similar character states, making them appear more closely related than they actually are.

### 3. Depends on the input tree

Small parsimony assumes that the tree topology is already known.

It reconstructs ancestral states for that topology; it does not determine the best tree topology.

### 4. Multiple optimal reconstructions

Different ancestral assignments may produce exactly the same minimum parsimony score.

Therefore, an ancestral nucleotide may not be uniquely identifiable.

### 5. Model dependence

Sankoff's results depend on the chosen substitution cost matrix.

Different cost matrices can produce different ancestral reconstructions.

### 6. Alignment errors

Incorrectly aligned sequences can produce misleading character states and therefore incorrect ancestral reconstructions.

---

# 9. Small Parsimony vs. Large Parsimony

These two problems are closely related but have different objectives.

| Feature | Small Parsimony | Large Parsimony |
|---|---|---|
| Tree topology | Given | Unknown |
| Objective | Find ancestral states minimizing cost | Find the tree and ancestral states minimizing cost |
| Computational difficulty | Efficient dynamic programming | Computationally difficult tree search |
| Algorithms | Fitch and Sankoff | Heuristic tree search using parsimony scores |

### Small Parsimony

Given a fixed tree, find the minimum-cost ancestral character assignments.

### Large Parsimony

Given observed character states, find the phylogenetic tree with the minimum total parsimony score.

Large parsimony must search among possible tree topologies, making it substantially more computationally challenging.

---

# 10. Applications in Bioinformatics

### 1. Ancestral sequence reconstruction

Infer possible ancestral nucleotides or amino acids at internal nodes of a phylogenetic tree.

### 2. Comparative genomics

Study the evolutionary history of genes and conserved genomic regions across species.

### 3. Protein evolution

Reconstruct ancestral amino acid states to investigate how protein sequences may have changed over evolutionary time.

### 4. Molecular evolution

Identify the minimum number of substitutions required to explain observed nucleotide variation.

### 5. Viral evolution

Explore possible ancestral character states across a known viral phylogeny.

For rapidly evolving viruses, however, model-based methods and careful treatment of temporal and sampling information are particularly important.

### 6. Morphological evolution

Apply parsimony to discrete characteristics such as anatomical features or developmental traits.

---

# 11. When Should You Use Fitch's or Sankoff's Algorithm?

Use **Fitch's algorithm** when:

- All character-state changes should have equal cost.
- You are analyzing standard unordered nucleotide or morphological characters.
- You need a simple and efficient ancestral-state reconstruction.

Use **Sankoff's algorithm** when:

- Different substitutions should have different costs.
- You want to penalize transitions and transversions differently.
- You are analyzing characters with more complex state-transition rules.
- You need a general dynamic programming framework for weighted parsimony.

For probabilistic ancestral reconstruction, consider maximum-likelihood or Bayesian phylogenetic methods instead.

These methods use explicit evolutionary models and can estimate probabilities for ancestral states rather than only reporting minimum-cost assignments.

---

## Key Takeaways

- The Small Parsimony problem reconstructs ancestral character states on a known phylogenetic tree while minimizing evolutionary change.
- Fitch's algorithm uses bottom-up set operations and assumes equal substitution costs.
- Sankoff's algorithm uses dynamic programming cost vectors and supports weighted substitution models.
- Both algorithms calculate the optimal parsimony score for a fixed tree and cost model.
- A top-down backtracking pass is needed to reconstruct specific internal-node states.
- Multiple equally optimal ancestral reconstructions may exist.
- For a fixed binary tree and nucleotide alphabet, Fitch's algorithm runs in \(O(n)\) time per alignment column, while Sankoff's algorithm runs in \(O(nk^2)\).
- Small parsimony assumes a known tree; finding the optimal tree is the separate Large Parsimony problem.

## References

1. Fitch, W. M. (1971). Toward defining the course of evolution: Minimum change for a specific tree topology. *Systematic Zoology*, 20(4), 406–416. https://doi.org/10.2307/2412116

2. Sankoff, D. (1975). Minimal mutation trees of sequences. *SIAM Journal on Applied Mathematics*, 28(1), 35–42. https://doi.org/10.1137/0128004

3. Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.

4. Yang, Z. (2014). *Molecular Evolution: A Statistical Approach*. Oxford University Press.
