# Hirschberg's Linear-Space Alignment Algorithm

## What Is Hirschberg's Algorithm?

**Hirschberg's algorithm** is a divide-and-conquer algorithm that reconstructs an optimal global sequence alignment using **linear auxiliary memory** instead of storing the entire dynamic programming matrix.

It is based on the **Needleman–Wunsch algorithm**, which finds an optimal global alignment between two DNA, RNA, or protein sequences.

The standard Needleman–Wunsch algorithm creates a matrix with $(n+1) \times (m+1)$ cells for sequences of lengths $n$ and $m$. Every cell stores the best alignment score for a pair of sequence prefixes. To reconstruct the alignment, a conventional implementation retains the entire matrix or equivalent traceback matrices.

This becomes computationally expensive for long sequences because the memory requirement grows quadratically ($\mathcal{O}(nm)$).

### Key Insight
> **We do not need to retain the entire dynamic programming matrix to recover an optimal alignment.** 

By calculating selected rows of the matrix, we can identify where an optimal alignment crosses a midpoint, split the problem at that intersection, and recursively solve the smaller subproblems. 

This technique reduces auxiliary memory from $\mathcal{O}(nm)$ to $\mathcal{O}(\min(n,m))$, while preserving the exact optimal alignment score under the same scoring model.

---

## Why Is Linear-Space Alignment Important?

Consider two DNA sequences that are each $10\text{,}000$ bases long. The full dynamic programming matrix contains:

$$(10\text{,}000 + 1)^2 = 100\text{,}020\text{,}001 \text{ cells}$$

If each cell stores a $32$-bit ($4$-byte) integer, the score matrix alone requires approximately **$400 \text{ MB}$** of memory. A conventional implementation storing traceback directions or additional affine state matrices requires even more.

For sequences of $500\text{,}000$ bases each, the matrix contains approximately $250 \text{ billion cells}$. Storing a single $32$-bit score per cell would require about **$1 \text{ TB}$** of memory.

Hirschberg's algorithm avoids retaining these billions of cells. It keeps only the 1D score vectors required to find a split point and then reconstructs the alignment recursively.

### Memory & Complexity Comparison

| Feature | Standard Needleman–Wunsch | Hirschberg's Algorithm |
| :--- | :--- | :--- |
| **Alignment Type** | Global | Global |
| **Optimality** | Guaranteed optimal score | Guaranteed optimal score |
| **Time Complexity** | $\mathcal{O}(nm)$ | $\mathcal{O}(nm)$ |
| **Auxiliary Memory** | $\mathcal{O}(nm)$ | $\mathcal{O}(\min(n,m))$ |
| **Reconstruction** | Reads stored matrix traceback | Divide-and-conquer recursion |
| **Main Trade-off** | Simple logic, quadratic memory | Slightly higher CPU overhead, linear memory |

---

## The Core Idea: Divide and Conquer

The algorithm divides sequence $A$ into two halves, identifies the exact index where an optimal alignment crosses that division in sequence $B$, and recursively aligns the resulting subsequences.

Let the input sequences be:
- $A = a_1 a_2 \ldots a_n$
- $B = b_1 b_2 \ldots b_m$

### Step 1: Split Sequence $A$ at its Midpoint
Choose the midpoint index of sequence $A$:

$$k = \left\lfloor \frac{n}{2} \right\rfloor$$

This divides sequence $A$ into two parts:
- **Left half:** $A_L = A[0..k]$
- **Right half:** $A_R = A[k..n]$

*(Note: While sequence $A$ is split evenly, sequence $B$ will be split at an optimal calculated position $j^*$.)*

### Step 2: Calculate Forward Alignment Scores
Perform a Needleman–Wunsch dynamic programming calculation between $A_L$ and the complete sequence $B$. 

Instead of storing the entire $2\text{D}$ matrix, retain only the current and previous rows. This yields a score vector $F$, where:

$$F[j] = \operatorname{score}(A_L, B[0..j])$$

Each entry $F[j]$ represents the optimal global alignment score between $A_L$ and the prefix $B[0..j]$.

### Step 3: Calculate Reverse Alignment Scores
Next, calculate alignment scores between the reversed right half $A_R$ and the reversed sequence $B$. This yields a score vector $R$, where:

$$R[j] = \operatorname{score}(\text{Rev}(A_R), \text{Rev}(B[j..m]))$$

Conceptually, $R[j]$ provides the optimal score for aligning the right half $A_R$ with the suffix $B[j..m]$. 

> **Implementation Note:** A reverse calculation can operate on reversed sequence views or iterators directly rather than copying strings in memory.

### Step 4: Find the Optimal Split in Sequence $B$
For each possible split position $j \in [0, m]$, sum the forward and reverse scores:

$$S[j] = F[j] + R[j]$$

Find the split position $j^*$ that maximizes this total score:

$$j^* = \underset{0 \le j \le m}{\operatorname{argmax}} \left( F[j] + R[j] \right)$$

The coordinate pair $(k, j^*)$ defines a midpoint through which an optimal global alignment must pass. The original problem is now divided into two independent subproblems:
1. Align $A[0..k]$ with $B[0..j^*]$
2. Align $A[k..n]$ with $B[j^*..m]$

### Step 5: Recurse on Subproblems
Recursively apply Hirschberg's algorithm to both subproblems. Continue dividing until a subproblem reaches a base case (e.g., when a sequence length is $0$ or $1$). Base cases are solved directly using standard DP.

### Step 6: Combine Results
Concatenate the alignments returned from the left and right subproblems. The combined result is a full, optimal global alignment.

---

## Mathematical Justification of the Split

An optimal global alignment can be split at any row $k$ into two optimal sub-alignments.

For any arbitrary split point $j$ in sequence $B$, an alignment passing through $(k, j)$ consists of:
1. An alignment of prefix $A[0..k]$ with prefix $B[0..j]$
2. An alignment of suffix $A[k..n]$ with suffix $B[j..m]$

The score of the best alignment passing through $(k, j)$ is strictly additive:

$$S[j] = F[j] + R[j]$$

Maximizing $S[j]$ over all $0 \le j \le m$ guarantees that $(k, j^*)$ lies on at least one optimal alignment path. Because each subproblem is solved optimally, the concatenated alignment is globally optimal.

*If multiple split positions yield the same maximum score, any of them may be chosen. Different tie-breaking logic will yield valid co-optimal alignments.*

---

## Worked Example

### Sequence Data
- **Sequence A:** `ACGT` ($n = 4$)
- **Sequence B:** `AGT` ($m = 3$)

### Parameters
- Match: $+2$
- Mismatch: $-1$
- Gap: $-2$

---

### Step-by-Step Execution

1. **Split Sequence A:**  
   $$k = \left\lfloor \frac{4}{2} \right\rfloor = 2 \implies A_L = \text{`AC`}, \quad A_R = \text{`GT`}$$

2. **Compute Forward Scores ($F$):**  
   Align `AC` against prefixes of `AGT` (`""`, `"A"`, `"AG"`, `"AGT"`):
   $$F = [-4, -1, 0, -2]$$

3. **Compute Reverse Scores ($R$):**  
   Align `GT` against suffixes of `AGT` (`"AGT"`, `"GT"`, `"T"`, `""`):
   $$R = [0, 2, 4, -4]$$

4. **Combine Scores ($S = F + R$):**
   $$S[0] = F[0] + R[0] = -4 + 0 = -4$$
   $$S[1] = F[1] + R[1] = -1 + 2 = 1$$
   $$S[2] = F[2] + R[2] = 0 + 4 = 4 \quad \star \text{ (Maximum)}$$
   $$S[3] = F[3] + R[3] = -2 + (-4) = -6$$

   $$\implies j^* = 2$$

5. **Subproblems Formed:**
   - **Left Subproblem:** Align $A[0..2]$ (`AC`) with $B[0..2]$ (`AG`)
   - **Right Subproblem:** Align $A[2..4]$ (`GT`) with $B[2..3]$ (`T`)

6. **Base Case Reconstruction:**
   - Left sub-alignment: `A C` / `A -` (Score: $2 - 2 = 0$)
   - Right sub-alignment: `G T` / `G T` (Score: $2 + 2 = 4$)

7. **Final Combined Alignment:**
```text
Sequence A:  A  C  G  T
Sequence B:  A  -  G  T

Total Score: 4
```

---

## Complexity Analysis

### Time Complexity: $\mathcal{O}(nm)$
Computing score vectors at the top level takes $c \cdot nm$ work. At the next level, the total area of the subproblems is:

$$k \cdot j^* + (n - k)(m - j^*) \le \frac{nm}{2}$$

Summing over all recursive levels:

$$\text{Total Work} = nm + \frac{nm}{2} + \frac{nm}{4} + \frac{nm}{8} + \ldots \le 2nm$$

Thus, the overall time complexity remains **$\mathcal{O}(nm)$**. Hirschberg is asymptotically as fast as standard Needleman–Wunsch, typically requiring roughly twice the number of low-level matrix cell evaluations.

### Auxiliary Space Complexity: $\mathcal{O}(\min(n,m))$
By orienting the DP calculation so that the shorter sequence forms the columns, the vector length is bounded by $\min(n, m)$.

$$\text{Auxiliary Space} = \boxed{\mathcal{O}(\min(n,m))}$$

The recursion stack depth is $\mathcal{O}(\log n)$, which uses negligible memory.

---

## Bioinformatics Application: Bacterial Contig Alignment

### Scenario
Aligning two bacterial genome contigs, each $500\text{,}000 \text{ bp}$ long.

- **Full DP Matrix:** $500\text{,}000^2 = 2.5 \times 10^{11} \text{ cells} \approx \mathbf{1 \text{ TB RAM}}$.
- **Hirschberg's Algorithm:** Requires a few $1\text{D}$ score vectors of length $500\text{,}001$. A $32$-bit integer vector of this size uses **$\sim 2 \text{ MB}$ of memory**.

### Use Case Fit
Hirschberg is essential when an **exact full-length global alignment** is required and available system RAM cannot hold quadratic matrices.

> **Note:** For very large genomes, modern aligners often utilize heuristics (e.g., seed-and-extend or minimizers) to bypass quadratic time constraints altogether. Hirschberg optimizes **space**, not **time**.

---

## Implementation Considerations (e.g., in Rust)

When implementing Hirschberg's algorithm in systems languages like Rust:

- **1D DP Buffer Reuse:** Allocate two reusable vectors (`prev_row`, `curr_row`) and swap them at each step to avoid repeated dynamic allocations.
- **Sequence Slicing:** Pass sequence slices (`&[u8]`) to recursive functions to prevent copying string data.
- **Base Case Threshold:** Fall back to standard 2D DP traceback when sequence lengths drop below a small threshold (e.g., $n, m \le 4$) for performance.
- **Affine Gaps:** Standard Hirschberg assumes a linear gap model. Extending Hirschberg to Gotoh's **affine gap model** requires tracking 3 state vectors ($M, X, Y$) across the split boundary.

---

## Key Takeaways

1. **Space Savings:** Reduces dynamic programming auxiliary memory from $\mathcal{O}(nm)$ to $\mathcal{O}(\min(n,m))$.
2. **Time Preserved:** Retains $\mathcal{O}(nm)$ time complexity (approx. $2\times$ runtime compared to full-matrix NW).
3. **Exact Optimality:** Guaranteed to yield an optimal global alignment score.
4. **Core Mechanism:** Uses forward and reverse 1D dynamic programming sweeps to locate midpoint intersections recursively.