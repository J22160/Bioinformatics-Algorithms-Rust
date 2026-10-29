# Expectation-Maximization for Transcript Quantification

A comprehensive guide to understanding and applying the Expectation-Maximization (EM) algorithm for RNA-seq transcript isoform quantification.

---

## Table of Contents

1. [What Is Expectation-Maximization?](#1-what-is-expectation-maximization)
2. [Why Is Transcript Quantification Difficult?](#2-why-is-transcript-quantification-difficult)
3. [Understanding EM with a Simple Analogy](#3-understanding-em-with-a-simple-analogy)
4. [The Mathematics Behind EM](#4-the-mathematics-behind-em)
   - [4.1 The Likelihood Function](#41-the-likelihood-function)
   - [4.2 E-step: Calculate Read Assignment Probabilities](#42-e-step-calculate-read-assignment-probabilities)
   - [4.3 M-step: Update Transcript Abundances](#43-m-step-update-transcript-abundances)
5. [Worked Example: Quantifying Four TP53 Isoforms](#5-worked-example-quantifying-four-tp53-isoforms)
6. [Transcript Abundance Is Not the Same as Read Count](#6-transcript-abundance-is-not-the-same-as-read-count)
7. [Computational Complexity](#7-computational-complexity)
8. [EM in Modern RNA-seq Quantification Tools](#8-em-in-modern-rna-seq-quantification-tools)
9. [Limitations of EM-Based Transcript Quantification](#9-limitations-of-em-based-transcript-quantification)
10. [Why This Matters in Bioinformatics](#10-why-this-matters-in-bioinformatics)
- [Key Takeaways](#key-takeaways)

---

## 1. What Is Expectation-Maximization?

**Expectation-Maximization (EM)** is an iterative statistical algorithm used to estimate unknown parameters when some of the underlying data assignments are uncertain or hidden (latent variables).

In RNA-seq analysis, EM helps estimate the relative abundance of different transcript isoforms when sequencing reads cannot be uniquely assigned to a single transcript.

A single gene may produce multiple transcript isoforms through:
- **Alternative splicing**
- **Alternative transcription start sites (TSS)**
- **Alternative polyadenylation (APA)**

Because these isoforms often share identical exon sequences, short sequencing reads frequently match multiple transcripts.

### Example Isoform Structure

```text
Transcript A: [Exon 1] ── [Exon 2] ── [Exon 3]
Transcript B: [Exon 1] ───────────── [Exon 3]
Transcript C: [Exon 1] ── [Exon 2] ── [Exon 4]
```

- A read originating entirely within **Exon 1** is compatible with all three transcripts.
- A read spanning the **Exon 2–Exon 3** junction distinguishes **Transcript A** from Transcript B and Transcript C.

This introduces two fundamental categories of reads:

| Read Category | Description |
| :--- | :--- |
| **Uniquely Compatible Reads** | Reads that map to only one transcript under the given model. |
| **Ambiguously Compatible Reads** | Reads that map to multiple candidate transcripts. |

The EM algorithm leverages observed reads, transcript compatibility mappings, and current abundance estimates to infer how likely each transcript was to have generated each read.

> **Important:** EM does not permanently assign an ambiguous read to a single transcript, nor does it recover the true physical origin of every read with certainty. Instead, it distributes read contributions probabilistically to estimate overall transcript abundances from available evidence.

---

## 2. Why Is Transcript Quantification Difficult?

Consider a simple gene with two transcript isoforms:

```text
Transcript A: ATGCCAGT
Transcript B: ATGCTAGT
```

Because most of their sequence is identical, short reads from shared regions map equally well to both:

```text
Read 1: ATG  --> Maps to Transcript A, Transcript B
Read 2: GCC  --> Maps to Transcript A
Read 3: CAG  --> Maps to Transcript A, Transcript B
Read 4: CTA  --> Maps to Transcript B
```

### Potential Quantification Pitfalls

1. **Equal Division:** Dividing ambiguous reads equally across all compatible transcripts skews abundance estimates toward shared regions regardless of true expression levels.
2. **Random Assignment:** Assigning ambiguous reads randomly introduces unnecessary variance and non-reproducible noise between runs.
3. **Discarding Ambiguous Reads:** Throwing away ambiguous reads loses significant sequencing data and severely biases estimates against highly conserved regions.

EM offers a mathematically principled alternative: it assigns each read a **probability** of originating from each compatible transcript, using those probabilities to iteratively update abundance estimates.

In real-world workflows, transcript quantifiers also model factors such as:
- **Effective transcript length**
- **Fragment-length distributions**
- **Sequencing biases** (GC content, positional bias)
- **Alignment quality and compatibility scores**

---

## 3. Understanding EM with a Simple Analogy

Imagine three cookie jars representing three transcript isoforms and a collection of cookies representing sequencing reads.

- Some cookies have distinct labels identifying their exact jar.
- Other cookies have ambiguous labels and could have come from multiple jars.

The goal is to estimate the total fraction of cookies produced by each jar.

```
       =================================================
       |                EM Iteration Loop              |
       =================================================
                               |
                               v
               +-------------------------------+
               |    Initialize Abundances      |
               | (e.g., A: 33.3%, B/C: 33.3%)  |
               +-------------------------------+
                               |
                               v
               +-------------------------------+
               |            E-step             |
               | Estimate assignment probabilities |
               |   for all ambiguous reads     |
               +-------------------------------+
                               |
                               v
               +-------------------------------+
               |            M-step             |
               | Recalculate abundance estimates|
               |  from total expected counts   |
               +-------------------------------+
                               |
                               v
                     /-------------------\
                    /   Has the model     \
                   <     converged?        >
                    \  (Threshold met?)   /
                     \-------------------/
                       /               \
                 No   /                 \   Yes
                     /                   \
      +---------------------+     +--------------------+
      | Repeat E-M Cycle    |     | Return Final       |
      |  (Loop back)        |     | Abundance Estimates|
      +---------------------+     +--------------------+
```

### Step 1: Initialize
Start with an initial estimate of the abundance of each transcript. For instance, equal proportions:
$$\theta_A = 33.3\%, \quad \theta_B = 33.3\%, \quad \theta_C = 33.3\%$$

*(Note: Real implementations may initialize using unique read counts or uniform distributions).*

### Step 2: E-step (Expectation)
For each ambiguous read, calculate the probability that it originated from each compatible transcript. If Transcript A is currently estimated to be twice as abundant as Transcript B, an ambiguous read compatible with both will assign a larger expected contribution to Transcript A.

### Step 3: M-step (Maximization)
Re-estimate transcript abundances using the expected assignments from the E-step.
- **Unique reads** contribute full counts to their single transcript.
- **Ambiguous reads** contribute fractional expected counts across compatible transcripts.

### Step 4: Repeat
Iterate the E-step and M-step until the estimates or log-likelihood stabilize (change falls below a specified convergence threshold $\epsilon$) or a maximum iteration limit is reached.

---

## 4. The Mathematics Behind EM

### Parameters and Notation

Let:
- $N$ = total number of observed reads or fragments.
- $M$ = total number of candidate transcripts.
- $r_n$ = the $n$-th observed read ($n = 1, \dots, N$).
- $\theta_i$ = relative abundance of transcript $i$ ($\sum_{i=1}^{M}\theta_i = 1, \theta_i \ge 0$).
- $p(r_n \mid i)$ = likelihood of observing read $r_n$ given it came from transcript $i$.

The mixture model assumes each read originates from exactly one candidate transcript, but the true generating transcript is an unobserved latent variable.

---

### 4.1 The Likelihood Function

The probability of observing a single read $r_n$ is the marginal probability across all candidate transcripts:

$$P(r_n \mid \theta) = \sum_{i=1}^{M} \theta_i \, p(r_n \mid i)$$

Assuming reads are conditionally independent given $\theta$, the total likelihood over all $N$ reads is:

$$L(\theta) = \prod_{n=1}^{N} \left[ \sum_{i=1}^{M} \theta_i \, p(r_n \mid i) \right]$$

In practice, we maximize the **log-likelihood** $\ell(\theta)$ to convert products to sums:

$$\ell(\theta) = \sum_{n=1}^{N} \log \left[ \sum_{i=1}^{M} \theta_i \, p(r_n \mid i) \right]$$

---

### 4.2 E-step: Calculate Read Assignment Probabilities

Let $z_{ni}$ represent the posterior probability that read $r_n$ originated from transcript $i$, given current parameter estimates $\theta$:

$$z_{ni} = \frac{\theta_i \, p(r_n \mid i)}{\sum_{j=1}^{M} \theta_j \, p(r_n \mid j)}$$

> **Property:** If transcript $i$ cannot physically produce read $r_n$, then $p(r_n \mid i) = 0$, leading directly to $z_{ni} = 0$.

#### Simplification for Equal Read Likelihoods
When candidate transcripts explain a read equally well ($p(r_n \mid A) = p(r_n \mid B)$), the formula simplifies to relative transcript abundance:

$$z_{nA} = \frac{\theta_A}{\theta_A + \theta_B}, \qquad z_{nB} = \frac{\theta_B}{\theta_A + \theta_B}$$

**Numerical Example:**
If current abundance parameters are $\theta_A = 0.7$ and $\theta_B = 0.3$, an ambiguous read compatible with both contributes $0.7$ fractional expected counts to Transcript A and $0.3$ to Transcript B.

---

### 4.3 M-step: Update Transcript Abundances

Under the basic mixture model, the updated abundance estimate $\theta_i^{\text{new}}$ is the sum of expected read assignments to transcript $i$ divided by total reads $N$:

$$\theta_i^{\text{new}} = \frac{1}{N} \sum_{n=1}^{N} z_{ni}$$

This formulation guarantees that $\sum_{i=1}^{M} \theta_i^{\text{new}} = 1$.

---

## 5. Worked Example: Quantifying Four TP53 Isoforms

Consider an RNA-seq experiment yielding **500 total reads** compatible with four isoforms of the *TP53* gene.

### Dataset Breakdown
- **Unique reads:** 200
- **Ambiguous reads:** 300 (compatible with all four isoforms)

#### Initial Unique Read Distribution

| Isoform | Unique Reads ($U_i$) |
| :--- | :---: |
| **TP53-A** | 80 |
| **TP53-B** | 50 |
| **TP53-C** | 40 |
| **TP53-D** | 30 |
| **Total** | **200** |

#### Simplifying Assumptions for this Example
1. All four isoforms have equal effective lengths.
2. All 300 ambiguous reads are equally compatible with all four isoforms ($p(r \mid i)$ is constant).
3. No sequencing or positional biases exist.

---

### Step 1: Initialization

Initialize transcript abundance parameters $\theta^{(0)}$ using the proportion of unique reads:

$$\theta_i^{(0)} = \frac{U_i}{\sum_{j} U_j}$$

| Isoform | Unique Reads ($U_i$) | Initial Proportion ($\theta_i^{(0)}$) |
| :--- | :---: | :---: |
| **TP53-A** | 80 | $80 / 200 = 0.40$ |
| **TP53-B** | 50 | $50 / 200 = 0.25$ |
| **TP53-C** | 40 | $40 / 200 = 0.20$ |
| **TP53-D** | 30 | $30 / 200 = 0.15$ |

---

### Step 2: E-step — Allocating Ambiguous Reads

Calculate expected fractional assignments for the 300 ambiguous reads:

$$E[\text{Ambiguous Reads}_i] = 300 \times \theta_i^{(0)}$$

| Isoform | $\theta_i^{(0)}$ | Fractional Assignment Calculation | Expected Ambiguous Reads |
| :--- | :---: | :---: | :---: |
| **TP53-A** | 0.40 | $300 \times 0.40$ | 120 |
| **TP53-B** | 0.25 | $300 \times 0.25$ | 75 |
| **TP53-C** | 0.20 | $300 \times 0.20$ | 60 |
| **TP53-D** | 0.15 | $300 \times 0.15$ | 45 |
| **Total** | **1.00** | — | **300** |

---

### Step 3: M-step — Updating Abundance Estimates

Combine unique read counts with expected ambiguous counts to compute total expected counts and new parameters:

| Isoform | Unique Reads | Expected Ambiguous Reads | Expected Total Count ($C_i$) | Updated Parameter ($\theta_i^{(1)}$) |
| :--- | :---: | :---: | :---: | :---: |
| **TP53-A** | 80 | 120 | 200 | $200 / 500 = 0.40$ |
| **TP53-B** | 50 | 75 | 125 | $125 / 500 = 0.25$ |
| **TP53-C** | 40 | 60 | 100 | $100 / 500 = 0.20$ |
| **TP53-D** | 30 | 45 | 75 | $75 / 500 = 0.15$ |
| **Total** | **200** | **300** | **500** | **1.00** |

*Note: Because all ambiguous reads had identical multi-mapping compatibility across all isoforms, the proportion converged in one step. In real datasets with complex overlap patterns, multiple iterations are required.*

---

### Step 4: Complexity in Real Data

In realistic datasets, reads group into overlapping compatibility subsets:

```text
Read Group 1: Compatible with {TP53-A, TP53-B}
Read Group 2: Compatible with {TP53-B, TP53-C, TP53-D}
Read Group 3: Compatible with {TP53-A, TP53-C}
Read Group 4: Compatible with {TP53-A, TP53-B, TP53-C, TP53-D}
```

Updating $\theta$ dynamically shifts assignment fractions across groups, driving multi-step iterative convergence.

---

## 6. Transcript Abundance Is Not the Same as Read Count

Longer transcript molecules produce more physical RNA-seq fragments than shorter molecules at equal molar concentrations. 

```
Transcript X (1,000 bp): [========]               --> Generates 10 fragments
Transcript Y (3,000 bp): [========================] --> Generates 30 fragments
```

To account for this, quantification models incorporate **Effective Transcript Length** ($\ell_i^{\text{eff}}$), which represents the number of valid starting positions for fragment generation given fragment size distribution constraints.

$$\text{Expected Fragment Count}_i \propto \text{Molar Abundance}_i \times \ell_i^{\text{eff}}$$

### Common Quantification Metrics

1. **Estimated Counts ($C_i$):** Expected total fragment assignments resulting from the EM algorithm.
2. **Transcripts Per Million (TPM):** Relative expression measure scaled to sum to $1,000,000$.
3. **FPKM / RPKM:** Fragments/Reads Per Kilobase of transcript per Million mapped reads (legacy length-normalized metrics).

#### Calculating TPM

First, compute length-normalized rate $q_i$:

$$q_i = \frac{C_i}{\ell_i^{\text{eff}}}$$

Then scale to sum to $10^6$:

$$\text{TPM}_i = 10^6 \times \frac{q_i}{\sum_{j=1}^{M} q_j}$$

---

## 7. Computational Complexity

### Variables
- $N$ = Total number of reads/fragments
- $M$ = Number of target transcripts
- $K$ = Average number of compatible transcripts per read ($K \ll M$)
- $I$ = Number of EM iterations until convergence
- $E$ = Number of distinct equivalence classes ($E \ll N$)

### Time Complexity
- **Naive Per-Read Approach:** $\mathcal{O}(I \cdot N \cdot K)$
- **Equivalence Class Approach:** $\mathcal{O}(I \cdot E \cdot K)$

### Space Complexity
- **Naive Storage:** $\mathcal{O}(N \cdot K + M)$
- **Equivalence Class Storage:** $\mathcal{O}(E \cdot K + M)$

### Equivalence Class Optimization
Modern tools avoid iterating over billions of individual reads by grouping reads with identical compatibility profiles into **Equivalence Classes**. The EM updates are then computed per class rather than per read, reducing runtime by orders of magnitude.

---

## 8. EM in Modern RNA-seq Quantification Tools

| Tool | Mapping Approach | Primary Inference Method | Key Features |
| :--- | :--- | :--- | :--- |
| **RSEM** | Full alignment (Bowtie/STAR) | Standard EM / Variational Bayes | Detailed generative probabilistic modeling |
| **Salmon** | Selective alignment / Lightweight mapping | Online EM / Variational Bayesian EM | High speed, comprehensive sequence bias correction |
| **kallisto** | Pseudoalignment via $k$-mer DAGs | EM over equivalence classes | Extremely fast, lightweight computational footprint |
| **oarfish** | Transcript-aware alignment mapping | Expectation-Maximization in Rust | Optimized for modern sequencing technologies (e.g., long-read/short-read handling) |

---

## 9. Limitations of EM-Based Transcript Quantification

1. **Indistinguishable Isoforms:** If two isoforms share identical sequence across all observed reads, EM cannot resolve their individual abundances uniquely.
2. **Low-Coverage Uncertainty:** Lowly expressed transcripts exhibit high estimator variance; small sample fluctuations drastically alter proportional outcomes.
3. **Reference Annotation Errors:** Reads originating from unannotated transcripts or genomic contamination will be misallocated to annotated isoforms.
4. **Complex Biases:** Non-uniform fragment sampling (e.g., $3'$ bias, GC-content bias) requires detailed modeling to prevent systematic quantification errors.

---

## 10. Why This Matters in Bioinformatics

EM-based quantification provides critical foundation for downstream genomic analysis:

- **Differential Isoform Usage:** Detecting alternative splicing switches between disease states and controls.
- **Gene-Level Summarization:** Aggregating transcript-level TPM estimates yields accurate total gene expression measurements without bias from length variation.
- **Uncertainty Propagation:** Downstream differential expression tools (e.g., sleuth, Swish, DESeq2) utilize quantification uncertainty to control false discovery rates.

---

## Key Takeaways

- **Iterative Nature:** EM iteratively estimates unobserved read origins (E-step) and updates transcript relative abundances (M-step).
- **Probabilistic Allocation:** Ambiguous reads are distributed as fractional expected counts based on current abundance parameter estimates.
- **Length Normalization:** Fragment counts must be adjusted by effective transcript length ($\ell^{\text{eff}}$) to calculate accurate molar abundances (TPM).
- **Efficiency via Equivalence Classes:** Modern quantifiers (Salmon, kallisto, RSEM, oarfish) aggregate reads into equivalence classes to perform EM iterations efficiently in memory.

## References

1. Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). *Maximum likelihood from incomplete data via the EM algorithm*. Journal of the Royal Statistical Society: Series B, 39(1), 1–38.
2. Li, B., Ruotti, V., Stewart, R. M., Thomson, J. A., & Dewey, C. N. (2010). *RNA-Seq gene expression estimation with read mapping uncertainty*. Bioinformatics, 26(4), 493–500.
3. Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). *Salmon provides fast and bias-aware quantification of transcript expression*. Nature Methods, 14, 417–419.
4. Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). *Near-optimal probabilistic RNA-seq quantification*. Nature Biotechnology, 34, 525–527.