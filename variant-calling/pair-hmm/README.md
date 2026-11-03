# Pair Hidden Markov Model (Pair HMM) for Read–Haplotype Likelihood

## 1. What Is Pair HMM?

A **Pair Hidden Markov Model (Pair HMM)** is a probabilistic sequence-alignment algorithm used to calculate the likelihood that a sequencing read originated from a particular candidate DNA haplotype.

Its central objective is to calculate the likelihood:

$$P(R \mid H)$$

Where:
* $R$ = Observed sequencing read.
* $H$ = Candidate haplotype sequence.
* $P(R \mid H)$ = Probability of observing the read given the haplotype.

Unlike conventional sequence alignment algorithms that search for a single alignment with the highest score, Pair HMM considers multiple possible alignments and calculates their **combined probability** under a unified probabilistic model.

It accounts for three major sequencing events:
1. **Matches and Substitutions:** A read base corresponds to a haplotype base, accounting for potential sequencing errors.
2. **Insertions:** The read contains an additional base relative to the haplotype.
3. **Deletions:** A haplotype base is absent from the read.

This makes Pair HMM indispensable for evaluating reads in regions containing single-nucleotide variants (SNVs), insertions, deletions (indels), and complex structural variants.

---

### Why is Pair HMM Important?

In tools like **GATK HaplotypeCaller**, reads are not simply compared position-by-position against the reference genome. Instead:
1. **Active regions** with potential variation are identified.
2. **Candidate haplotypes** are constructed locally for those regions.
3. **Pair HMM evaluates** how well the observed sequencing reads support each candidate haplotype.

#### Example Scenario

Consider a genomic region with two candidate haplotypes:

* **Reference Haplotype:** `ACGTACGT`
* **Alternate Haplotype:** `ACGTTACGT` *(contains an extra T)*

Suppose a sequencing read spanning this region is observed:
* **Read:** `ACGTTACGT`

Comparing the read against the reference requires introducing an insertion into the alignment. The alternate haplotype, however, matches the read directly. Pair HMM calculates the probability of generating this read from each candidate haplotype, accounting for both sequencing error rates and gap penalties.

---

## 2. Why Use a Hidden Markov Model?

A conventional sequence alignment algorithm identifies an optimal alignment between two sequences by optimizing a deterministic scoring function (e.g., assigning positive scores for matches, penalties for mismatches, and affine penalties for gap opening/extension).

Pair HMM treats sequence alignment as a **probabilistic process** where the alignment state sequence is **hidden** and the observed read bases are **emitted** by those states.

```
Haplotype:  A  C  G  T  A  C  G  T
Read:       A  C  G  -  A  C  G  T
State:      M  M  M  D  M  M  M  M
```

Here, the read is missing the `T` at the fourth haplotype position. The hidden state sequence describes the alignment geometry, while the read bases constitute the observed data.

### Core Distinction

> **Conventional Alignment:** Selects a single optimal alignment path (e.g., via Needleman-Wunsch or Smith-Waterman).
> 
> **Pair HMM (Forward Algorithm):** Sums the probabilities of **all allowed alignment paths** to compute total read–haplotype likelihood.

---

## 3. The Three Pair HMM States

The model transition system operates across three hidden states:

| State | Consumes Read Base? | Consumes Haplotype Base? | Modeled Event |
| :--- | :---: | :---: | :--- |
| **Match ($M$)** | **Yes** | **Yes** | Aligned bases (exact matches or substitutions) |
| **Insertion ($I$)** | **Yes** | **No** | Extra base present in read |
| **Deletion ($D$)** | **No** | **Yes** | Base missing from read |

---

### 3.1 Match State ($M$)

Represents an alignment between one read base and one haplotype base (which may or may not be identical).

```
Haplotype:  A  C  G  T  A
Read:       A  C  A  T  A
State:      M  M  M  M  M
                  ^ Substitution (G -> A)
```

The probability of observing read base $r_i$ aligned to haplotype base $h_j$ with Phred error probability $e_i$:

$$P(r_i \mid h_j) = \begin{cases} 1 - e_i & \text{if } r_i = h_j \\ \frac{e_i}{3} & \text{if } r_i \neq h_j \end{cases}$$

---

### 3.2 Insertion State ($I$)

Represents a base present in the read but absent from the haplotype.

```
Haplotype:  A  C  G  -  T  A
Read:       A  C  G  T  T  A
State:      M  M  M  I  M  M
```

The insertion state consumes a read base without advancing the haplotype position.

---

### 3.3 Deletion State ($D$)

Represents a base present in the haplotype but absent from the read.

```
Haplotype:  A  C  G  T  A  C
Read:       A  C  G  -  A  C
State:      M  M  M  D  M  M
```

The deletion state advances the haplotype position without consuming a read base.

---

## 4. Transition Probabilities

Transitions between states model the frequency of matching, opening gaps, and extending gaps:

* $g_o$: **Gap-open probability** (transitioning from $M \to I$ or $M \to D$)
* $g_e$: **Gap-extension probability** (remaining in $I \to I$ or $D \to D$)

```
             ┌──────────────┐
             │              │
             ▼              │
          ┌───────┐      ┌───────┐
          │   M   │─────>│   I   │
          └───────┘      └───────┘
             │              │
             │              │
             ▼              ▼
          ┌───────┐      ┌───────┐
          │   D   │<─────│   M   │
          └───────┘      └───────┘
```

### Simplified Transition Matrix

| From State | To State | Transition Probability |
| :--- | :--- | :--- |
| **Match ($M$)** | Match ($M$) | $1 - 2g_o$ |
| **Match ($M$)** | Insertion ($I$) | $g_o$ |
| **Match ($M$)** | Deletion ($D$) | $g_o$ |
| **Insertion ($I$)** | Match ($M$) | $1 - g_e$ |
| **Insertion ($I$)** | Insertion ($I$) | $g_e$ |
| **Deletion ($D$)** | Match ($M$) | $1 - g_e$ |
| **Deletion ($D$)** | Deletion ($D$) | $g_e$ |

#### Total Probability Conservation

$$\sum_{Y \in \{M, I, D\}} P(X \to Y) = 1 \quad \forall X \in \{M, I, D\}$$

*For $M$:* $(1 - 2g_o) + g_o + g_o = 1$  
*For $I$:* $(1 - g_e) + g_e = 1$  
*For $D$:* $(1 - g_e) + g_e = 1$

> **Note:** Production platforms (such as GATK) often utilize base-quality-dependent dynamic transition matrices rather than constant global transition probabilities.

---

## 5. Emission Probabilities

Emission probabilities model the likelihood of generating specific read bases given the state.

### 5.1 Match-State Emissions

Given Phred Quality score $Q_i$, the base error probability $e_i$ is:

$$e_i = 10^{-\frac{Q_i}{10}}$$

$$\text{Emission Probability: } P(r_i \mid h_j) = \begin{cases} 1 - e_i & \text{if } r_i = h_j \\ \frac{e_i}{3} & \text{if } r_i \neq h_j \end{cases}$$

#### Numerical Example
For a $Q30$ base ($e_i = 10^{-3} = 0.001$):
* **Match ($A \mid A$):** $1 - 0.001 = 0.999$
* **Mismatch ($G \mid A$):** $\frac{0.001}{3} \approx 0.000333$

---

### 5.2 Insertion & Deletion Emissions

* **Insertion State ($I$):** Emits a read base independent of the haplotype, typically modeled using uniform background probability:
  $$P(r_i \mid I) = \frac{1}{4} = 0.25$$
* **Deletion State ($D$):** Emits **no** read base (probability factor = 1).

---

## 6. Dynamic Programming Formulation

To evaluate read $R = r_1 r_2 \dots r_N$ against haplotype $H = h_1 h_2 \dots h_M$, we define three dynamic programming matrices:

* $M_{i,j}$: Probability of aligning read prefix $r_{1..i}$ to haplotype prefix $h_{1..j}$, ending in state **$M$**.
* $I_{i,j}$: Probability of aligning read prefix $r_{1..i}$ to haplotype prefix $h_{1..j}$, ending in state **$I$**.
* $D_{i,j}$: Probability of aligning read prefix $r_{1..i}$ to haplotype prefix $h_{1..j}$, ending in state **$D$**.

Let $a_{XY} = P(X \to Y)$ represent state transition probabilities.

---

### 6.1 Recurrence Relations

$$\begin{aligned}
M_{i,j} &= P(r_i \mid h_j) \cdot \left[ M_{i-1,j-1} \, a_{MM} + I_{i-1,j-1} \, a_{IM} + D_{i-1,j-1} \, a_{DM} \right] \\[0.75em]
I_{i,j} &= P(r_i \mid I) \cdot \left[ M_{i-1,j} \, a_{MI} + I_{i-1,j} \, a_{II} \right] \\[0.75em]
D_{i,j} &= M_{i,j-1} \, a_{MD} + D_{i,j-1} \, a_{DD}
\end{aligned}$$

---

### 6.2 Matrix Form Summary

$$\begin{aligned}
M_{i,j} &= P(r_i \mid h_j) \left[ M_{i-1,j-1} (1 - 2g_o) + I_{i-1,j-1} (1 - g_e) + D_{i-1,j-1} (1 - g_e) \right] \\
I_{i,j} &= P(r_i \mid I) \left[ M_{i-1,j} \, g_o + I_{i-1,j} \, g_e \right] \\
D_{i,j} &= M_{i,j-1} \, g_o + D_{i,j-1} \, g_e
\end{aligned}$$

---

## 7. Initialization and Termination

### 7.1 Global Alignment Initialization
For standard full-sequence matching:

$$M_{0,0} = 1, \quad I_{0,0} = 0, \quad D_{0,0} = 0$$

### 7.2 Substring / Local Alignment
Reads usually map to arbitrary substrings of longer haplotypes. Initial probabilities across $j \in [1, M]$ can be uniformly initialized to allow flexible alignment start sites.

### 7.3 Termination
Total read likelihood $P(R \mid H)$ is computed by summing final state probabilities over allowed terminal positions:

$$P(R \mid H) = \sum_{j=1}^{M} \left( M_{N,j} + I_{N,j} \right)$$

---

## 8. Worked Example: Evaluating Support Across Haplotypes

Suppose we have two haplotypes and four sequencing reads:

* **$H_1$ (Reference):** `ACGTACGT`
* **$H_2$ (Alternate):** `ACGTTACGT` *(Insertion of T)*

### Example Likelihood Table

| Read | Sequence | $P(R \mid H_1)$ | $P(R \mid H_2)$ | Supported Haplotype |
| :--- | :--- | :---: | :---: | :---: |
| **Read 1** | `ACGTTACGT` | $0.001$ | $0.950$ | $H_2$ |
| **Read 2** | `ACGTTACGT` | $0.002$ | $0.900$ | $H_2$ |
| **Read 3** | `ACGT-ACGT` | $0.850$ | $0.020$ | $H_1$ |
| **Read 4** | `ACGTTACGT` | $0.001$ | $0.920$ | $H_2$ |

By aggregating likelihood values across all reads, the variant caller establishes strong statistical evidence favoring $H_2$.

---

## 9. Integration in GATK HaplotypeCaller Workflow

```
               [ Aligned Read Data ]
                         │
                         ▼
             [ Active Region Detection ]
                         │
                         ▼
           [ Local Haplotype Assembly ]
                         │
                         ▼
        ┌────────────────────────────────┐
        │   Pair HMM Engine              │
        │   Calculates P(Read | Haplotype)│
        └────────────────────────────────┘
                         │
                         ▼
         [ Read-Haplotype Likelihood Matrix ]
                         │
                         ▼
          [ Bayesian Genotyping Engine ]
                         │
                         ▼
               [ Final VCF / GVCF ]
```

---

## 10. Computational Complexity

Let:
* $N$ = Read length
* $M$ = Haplotype length
* $S$ = Number of states ($S = 3$)
* $R$ = Number of reads
* $H$ = Number of haplotypes

### 10.1 Single Read–Haplotype Pair
* **Time Complexity:** $\mathcal{O}(N \cdot M \cdot S^2) \implies \mathbf{\mathcal{O}(N \cdot M)}$
* **Space Complexity:** $\mathcal{O}(N \cdot M)$ standard, or $\mathbf{\mathcal{O}(M)}$ using two-row rolling memory buffers.

### 10.2 Total Dataset Complexity
Evaluating $R$ reads against $H$ haplotypes:

$$\text{Time Complexity} = \mathbf{\mathcal{O}(R \cdot H \cdot N \cdot M)}$$

Assuming $N \approx M \approx L$:

$$\text{Time Complexity} = \mathbf{\mathcal{O}(R \cdot H \cdot L^2)}$$

---

## 11. Pair HMM vs. Classical Alignment

| Feature | Classical Alignment (e.g., Smith-Waterman) | Pair HMM (Forward Algorithm) |
| :--- | :--- | :--- |
| **Goal** | Single optimal alignment | Integrated overall probability |
| **Scoring** | Arbitrary match/gap scores | Phred-scaled probabilistic emissions |
| **Output** | Alignment path & score | Total likelihood $P(R \mid H)$ |
| **Path Selection** | Selects maximum ($\max$) | Sums overall paths ($\sum$) |
| **Algorithm** | Viterbi / Dynamic Programming | Forward Algorithm |

### Connection to Affine Gap Penalties
Taking the negative logarithm of probability transitions transforms Pair HMM calculations into equivalent additive affine gap penalties:

$$\text{Penalty} = -\log P(\text{path})$$

---

## 12. Forward Algorithm vs. Viterbi Algorithm

$$\begin{aligned}
\text{Forward (Total Likelihood):} \quad P(R \mid H) &= \sum_{\pi} P(R, \pi \mid H) \\[0.5em]
\text{Viterbi (Best Path):} \quad \pi^* &= \arg\max_{\pi} P(R, \pi \mid H)
\end{aligned}$$

---

## 13. High-Performance Optimizations

1. **SIMD Vectorization:** Parallel evaluation of multiple DP matrix cells or read sequences using CPU AVX2/AVX-512 vector operations.
2. **Rolling Memory Buffers:** Reduced memory footprint from $\mathcal{O}(NM)$ to $\mathcal{O}(\min(N,M))$.
3. **Log-Domain Scaling:** Underflow prevention via floating-point log-scale transformations or exponent normalization.
4. **Hardware Acceleration:** Native C++ implementations (e.g., Intel Genomics Kernel Library / OpenMP).

---

## 14. Summary of Key Takeaways

* **Core Function:** Pair HMM computes $P(R \mid H)$, the probabilistic likelihood of observing read $R$ given haplotype $H$.
* **Hidden States:** Match ($M$), Insertion ($I$), and Deletion ($D$).
* **Forward Summing:** Sums probabilities across all valid alignment paths rather than picking a single best path.
* **Primary Complexity:** $\mathcal{O}(N \cdot M)$ time per read–haplotype pair.