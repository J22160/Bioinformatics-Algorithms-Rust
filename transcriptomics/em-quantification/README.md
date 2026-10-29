# Expectation-Maximization for Transcript Quantification

## 1. What Is Expectation-Maximization?

**Expectation-Maximization (EM)** is an iterative statistical algorithm used to estimate unknown parameters when some of the underlying data assignments are uncertain or hidden.

In RNA-seq analysis, EM can help estimate the abundance of different transcripts when sequencing reads cannot be uniquely assigned to a single transcript.

A gene may produce multiple transcript isoforms through alternative splicing, alternative transcription start sites, or alternative polyadenylation. These isoforms can share many exons, meaning that a short sequencing read may be compatible with several transcripts.

For example, consider three transcript isoforms:

```text
Transcript A: [Exon 1] ── [Exon 2] ── [Exon 3]
Transcript B: [Exon 1] ───────────── [Exon 3]
Transcript C: [Exon 1] ── [Exon 2] ── [Exon 4]
```

A read originating from Exon 1 may be compatible with all three transcripts. However, a read spanning the Exon 2–Exon 3 junction may distinguish Transcript A from the other isoforms.

This creates two broad categories of reads:

- **Uniquely compatible reads:** Reads that support one transcript under the quantification model.
- **Ambiguously compatible reads:** Reads that are compatible with multiple transcripts.

The EM algorithm uses the observed reads, their transcript compatibility, and estimated transcript abundances to infer how likely each transcript is to have generated each read.

Instead of permanently assigning every ambiguous read to one transcript, EM distributes its contribution probabilistically and repeatedly updates the abundance estimates.

**Important:** EM does not recover the true origin of every individual read with certainty. It estimates transcript abundances from the available evidence.

---

## 2. Why Is Transcript Quantification Difficult?

Suppose a gene has two transcript isoforms:

```text
Transcript A: ATGCCAGT
Transcript B: ATGCTAGT
```

Most of their sequences are identical. A short read from a shared region may map equally well to both transcripts.

Consider these hypothetical RNA-seq reads:

```text
Read 1: ATG
Read 2: GCC
Read 3: CAG
Read 4: CTA
```

Depending on the read positions and the transcript sequences, some reads may support both isoforms, while others may distinguish them.

If every ambiguous read is assigned equally to every compatible transcript, the estimated abundances can be misleading. If every read is assigned randomly, the result can vary between runs and introduce unnecessary noise.

EM offers a principled alternative: assign each read a probability of originating from each compatible transcript, then use those assignments to update the abundance estimates.

In practice, transcript quantifiers also consider factors such as effective transcript length, fragment-length distributions, sequencing biases, and alignment or compatibility scores.

---

## 3. Understanding EM with a Simple Analogy

Imagine three cookie jars representing three transcript isoforms and a collection of cookies representing sequencing reads.

Some cookies have labels identifying the jar they came from. Others have ambiguous labels and could have come from multiple jars.

The goal is to estimate how many cookies came from each jar.

### Step 1: Initialize

Start with an initial estimate of the abundance of each transcript.

For example:

```text
Transcript A: 33.3%
Transcript B: 33.3%
Transcript C: 33.3%
```

Equal initialization is a simple illustration. Actual implementations may use other initialization strategies.

### Step 2: E-step — Expectation

For each ambiguous read, estimate the probability that it originated from each compatible transcript.

If Transcript A is estimated to be more abundant than Transcript B, and both are equally likely to generate the observed read apart from abundance, the read receives a larger expected contribution from Transcript A.

### Step 3: M-step — Maximization

Use the expected assignments of all reads to update the transcript abundance estimates.

Uniquely compatible reads contribute to their compatible transcripts. Ambiguous reads contribute fractional expected counts to the transcripts that could have generated them.

### Step 4: Repeat

Repeat the E-step and M-step until the abundance estimates or likelihood change by less than a specified convergence threshold, or until a maximum iteration count is reached.

```text
Initialize transcript abundances
             |
             v
    E-step: assign reads
    probabilistically
             |
             v
    M-step: update transcript
    abundance estimates
             |
             v
    Has the model converged?
          /       \
        No         Yes
        |           |
        +-----------+
        Repeat      Return estimates
```

The final result is a set of estimated transcript abundances, not necessarily a definitive assignment of each read to one transcript.

---

## 4. The Mathematics Behind EM

Let:

- \(N\) = number of observed reads or fragments.
- \(M\) = number of candidate transcripts.
- \(r_n\) = the \(n\)-th observed read or fragment.
- \(\theta_i\) = the relative abundance of transcript \(i\).
- \(p(r_n \mid i)\) = the probability of observing read \(r_n\) if it originated from transcript \(i\).

For the simplified mixture model, transcript proportions satisfy:

\[
\theta_i \geq 0
\]

and

\[
\sum_{i=1}^{M}\theta_i=1
\]

The model assumes that each read originates from one of the candidate transcripts, but its exact source is not directly observed.

### 4.1 The likelihood function

The probability of observing a read is the sum of its probabilities under all candidate transcripts, weighted by their abundance estimates:

\[
P(r_n\mid\theta)
=
\sum_{i=1}^{M}
\theta_i p(r_n\mid i)
\]

Assuming reads are conditionally independent under the model, the likelihood is:

\[
L(\theta)
=
\prod_{n=1}^{N}
\left[
\sum_{i=1}^{M}
\theta_i p(r_n\mid i)
\right]
\]

It is generally more convenient to maximize the log-likelihood:

\[
\ell(\theta)
=
\sum_{n=1}^{N}
\log
\left[
\sum_{i=1}^{M}
\theta_i p(r_n\mid i)
\right]
\]

The EM algorithm iteratively improves this likelihood under the specified model.

### 4.2 E-step: Calculate read assignment probabilities

Let \(z_{ni}\) denote the probability that read \(n\) originated from transcript \(i\), given the current abundance estimates.

Using Bayes' rule:

\[
\boxed{
z_{ni}
=
\frac{
\theta_i p(r_n\mid i)
}{
\sum_{j=1}^{M}
\theta_j p(r_n\mid j)
}
}
\]

The probability is zero if transcript \(i\) cannot generate the observed read under the model.

For a read compatible with transcripts A and B, the assignment probabilities depend on both the abundance estimates and the read likelihoods.

If the two transcripts explain the read equally well, then:

\[
z_{nA}
=
\frac{\theta_A}{\theta_A+\theta_B}
\]

\[
z_{nB}
=
\frac{\theta_B}{\theta_A+\theta_B}
\]

For example, if the current abundance estimates are:

\[
\theta_A=0.7,\qquad \theta_B=0.3
\]

and both transcripts have equal read likelihoods, the read contributes an expected 0.7 to Transcript A and 0.3 to Transcript B.

This is a probabilistic allocation, not a claim that 70% of the individual read physically originated from A.

### 4.3 M-step: Update transcript abundances

Under the simplified mixture model, the updated abundance estimate is the expected number of reads assigned to transcript \(i\), divided by the total number of reads:

\[
\boxed{
\theta_i^{\text{new}}
=
\frac{1}{N}
\sum_{n=1}^{N}z_{ni}
}
\]

This update maintains the constraint that transcript proportions sum to one.

The E-step and M-step are repeated until the estimates converge.

**A practical qualification:** This simplified update assumes that the transcript-specific read likelihoods are appropriately normalized and that the model parameter represents the probability that a sampled fragment originates from each transcript. Actual RNA-seq quantifiers may parameterize transcript abundance differently and account for effective length and other biases.

---

## 5. Worked Example: Quantifying Four TP53 Isoforms

Suppose an RNA-seq experiment produces 500 reads compatible with the TP53 gene's four transcript isoforms.

The initial analysis identifies:

- 200 uniquely compatible reads.
- 300 ambiguous reads.

The unique-read evidence is:

| Transcript | Unique reads |
|---|---:|
| TP53-A | 80 |
| TP53-B | 50 |
| TP53-C | 40 |
| TP53-D | 30 |
| **Total** | **200** |

The remaining 300 reads can map to multiple isoforms.

**Important:** Knowing only that 300 reads are ambiguous is not enough to determine their allocation. We also need to know which transcripts each read is compatible with and how likely each transcript is to have generated it.

For a transparent worked example, assume:

1. All four isoforms have equal effective lengths.
2. The 300 ambiguous reads are equally compatible with all four isoforms.
3. All compatible transcripts have equal read likelihood for these ambiguous reads.
4. There are no additional sequence or fragment biases.

These assumptions make this a simplified example rather than a realistic TP53 quantification model.

### Step 1: Initialize using the unique reads

A simple initialization uses the relative unique-read counts:

\[
\theta_i^{(0)}
=
\frac{U_i}{\sum_j U_j}
\]

where \(U_i\) is the number of uniquely compatible reads for transcript \(i\).

The initial estimates are:

| Transcript | Unique reads | Initial proportion |
|---|---:|---:|
| TP53-A | 80 | 0.40 |
| TP53-B | 50 | 0.25 |
| TP53-C | 40 | 0.20 |
| TP53-D | 30 | 0.15 |

These estimates sum to 1.0.

### Step 2: E-step — Allocate ambiguous reads

Since every ambiguous read is assumed equally compatible with all four transcripts, the posterior assignment probabilities are proportional to the current abundance estimates.

For TP53-A:

\[
P(A\mid r)=0.40
\]

Therefore, the expected number of ambiguous reads assigned to TP53-A is:

\[
300\times0.40=120
\]

Similarly:

\[
300\times0.25=75
\]

\[
300\times0.20=60
\]

\[
300\times0.15=45
\]

The expected assignments are:

| Transcript | Probability per ambiguous read | Expected ambiguous reads |
|---|---:|---:|
| TP53-A | 0.40 | 120 |
| TP53-B | 0.25 | 75 |
| TP53-C | 0.20 | 60 |
| TP53-D | 0.15 | 45 |
| **Total** | **1.00** | **300** |

Fractional expected counts are valid in probabilistic quantification. They represent expected contributions rather than literal fractions of physical reads.

### Step 3: M-step — Update the abundance estimates

Add the expected ambiguous-read contributions to the unique-read counts.

| Transcript | Unique reads | Expected ambiguous reads | Expected total |
|---|---:|---:|---:|
| TP53-A | 80 | 120 | 200 |
| TP53-B | 50 | 75 | 125 |
| TP53-C | 40 | 60 | 100 |
| TP53-D | 30 | 45 | 75 |
| **Total** | **200** | **300** | **500** |

The updated proportions are:

\[
\theta_A=\frac{200}{500}=0.40
\]

\[
\theta_B=\frac{125}{500}=0.25
\]

\[
\theta_C=\frac{100}{500}=0.20
\]

\[
\theta_D=\frac{75}{500}=0.15
\]

In this specially constructed example, the estimates do not change after the first update. All ambiguous reads have the same compatibility pattern, so the relative proportions established by the unique reads remain unchanged.

### Step 4: Why real datasets require iteration

In real RNA-seq data, different reads are compatible with different subsets of transcripts.

For example:

```text
Read group 1: TP53-A, TP53-B
Read group 2: TP53-B, TP53-C, TP53-D
Read group 3: TP53-A, TP53-C
Read group 4: TP53-A, TP53-B, TP53-C, TP53-D
```

The assignment probabilities for each group differ. Updating the abundance estimates changes the relative probabilities, which then changes the next round of assignments.

This feedback between read assignments and abundance estimates is why EM is iterative.

---

## 6. Transcript Abundance Is Not the Same as Read Count

A key issue in RNA-seq quantification is that longer transcripts can generate more sequencing fragments than shorter transcripts at the same molecule concentration.

Consequently, raw fragment counts are not directly equivalent to transcript abundance.

Many quantifiers account for this using an **effective transcript length**, which reflects the number of fragment start positions that can produce observable fragments under the model.

A simplified relationship is:

\[
\text{Expected fragment count}_i
\propto
\text{Molecule abundance}_i
\times
\text{Effective length}_i
\]

Effective length can depend on the fragment-length distribution and other modeling assumptions.

This is why transcript quantification commonly distinguishes between:

- **Estimated counts:** Expected fragment assignments to transcripts.
- **TPM (Transcripts Per Million):** A length-normalized relative abundance measure scaled to one million.
- **FPKM/RPKM:** Length- and sequencing-depth-normalized measures, with different conventions for paired-end and single-end data.

For TPM, a simplified calculation is:

\[
q_i=\frac{C_i}{\ell_i^{\mathrm{eff}}}
\]

\[
TPM_i
=
10^6
\frac{q_i}{\sum_j q_j}
\]

where \(C_i\) is the estimated fragment count and \(\ell_i^{\mathrm{eff}}\) is the effective length.

The exact quantification and normalization procedures vary by tool. Transcript-level TPM estimates are relative abundance measures; they are not direct measurements of absolute transcript molecule counts.

---

## 7. Computational Complexity

Let:

- \(N\) = number of reads or fragments.
- \(M\) = number of candidate transcripts.
- \(K\) = average number of compatible transcripts per read.
- \(I\) = number of EM iterations.
- \(E\) = number of distinct read-compatibility classes, when such classes are used.

### Time complexity

For a straightforward implementation that processes every read and each compatible transcript in every iteration, the cost is approximately:

\[
\boxed{O(I N K)}
\]

This assumes the candidate transcript set for each read is already available and that evaluating a read-transcript likelihood is constant-time or has already been computed.

If the likelihood calculations require additional sequence processing, that cost must also be included.

The actual runtime depends on the number of candidate mappings, the number of iterations, the number of transcripts, and the inference strategy.

### Space complexity

A naive implementation that stores every read-to-transcript compatibility can require:

\[
O(NK+M)
\]

space, in addition to the input sequences and any alignment or index data.

However, implementations need not store all per-read assignments. Many quantifiers aggregate reads into equivalence classes: groups of reads compatible with the same set of transcripts and sharing the relevant likelihood information.

When the model permits this aggregation, the iterative computation can operate on the much smaller collection of classes rather than on every read independently. Memory use then depends on the number of classes and their transcript memberships, along with other data structures.

### Convergence

EM is typically iterated until the change in likelihood or abundance estimates falls below a threshold, or a maximum iteration limit is reached.

There is no universal requirement that convergence take 50–200 iterations. The number of iterations depends on the data, initialization, model, and stopping criterion.

EM generally improves or maintains the likelihood at each iteration under its assumptions, but convergence does not guarantee that the individual transcript abundances are uniquely identifiable. Transcripts with identical or nearly indistinguishable sequence evidence may remain difficult or impossible to separate reliably.

---

## 8. EM in Modern RNA-seq Quantification Tools

Several RNA-seq tools use probabilistic models to estimate transcript abundance, but their algorithms are not identical.

| Tool | General approach | Key consideration |
|---|---|---|
| RSEM | Probabilistic transcript quantification using an expectation-maximization framework | Models uncertainty in read assignment and transcript abundance |
| Salmon | Transcript quantification using lightweight mapping and probabilistic inference | Supports bias correction and efficient quantification; inference details depend on configuration |
| kallisto | Pseudoalignment and transcript quantification using a probabilistic model | Groups reads by transcript compatibility and estimates abundance efficiently |
| oarfish | Rust-based transcript quantification | Uses efficient transcript-aware quantification methods designed for modern sequencing data |

These tools share the general goal of estimating transcript abundances despite ambiguous read assignments. They differ in mapping or compatibility detection, statistical models, bias correction, and implementation.

It is therefore more accurate to say that these tools use related probabilistic quantification strategies than to assume every tool implements exactly the same EM equations or update procedure.

---

## 9. Limitations of EM-Based Transcript Quantification

### Shared exons and indistinguishable isoforms

If two isoforms contain the same observable sequence regions, reads from those regions may not distinguish them. EM can allocate expected counts based on the model, but it cannot create information that the sequencing data do not contain.

### Lowly expressed transcripts

Transcripts with few supporting reads have uncertain abundance estimates. Small changes in the assignment of ambiguous reads may substantially change their estimated proportions.

### Sequencing and mapping biases

GC content, fragment-length distributions, positional biases, sequencing errors, and alignment ambiguity can affect the evidence available to the model. Bias correction helps, but cannot eliminate every source of error.

### Model assumptions

Incorrect transcript annotations or missing isoforms can lead to misleading estimates. If a read originates from a transcript absent from the reference, the quantifier may assign it to an imperfectly matching alternative or leave it unassigned.

### Gene-level versus transcript-level quantification

When several isoforms share exons, estimating total gene expression is often more robust than estimating the abundance of each individual isoform. Transcript-level estimates require sufficient distinguishing evidence.

---

## 10. Why This Matters in Bioinformatics

EM-based transcript quantification is useful whenever multiple biological sources can explain the same observation.

In RNA-seq, it enables researchers to estimate the relative expression of transcript isoforms even when many reads map ambiguously. These estimates can support analyses of alternative splicing, isoform switching, transcript expression, and gene regulation.

For downstream analysis, the uncertainty in transcript-level estimates should be considered when comparing conditions or interpreting low-abundance isoforms. Differential transcript expression methods and differential transcript usage methods answer related but distinct questions.

EM is not a substitute for experimental evidence. When isoforms cannot be distinguished by the available reads, longer reads, junction-spanning reads, improved transcript annotations, or additional experimental data may be needed.

---

## Key Takeaways

- EM estimates hidden assignments and model parameters iteratively.
- In transcript quantification, it assigns probabilistic contributions of ambiguous reads to compatible transcripts.
- The E-step calculates read assignment probabilities using transcript abundance and read likelihood.
- The M-step updates abundance estimates from the expected read assignments.
- Effective transcript length and sequencing biases matter when converting fragment evidence into transcript abundance.
- A basic implementation has an approximate per-iteration cost of \(O(NK)\) when read compatibility information is already available.
- Actual tools may aggregate reads into equivalence classes to improve efficiency.
- EM cannot resolve isoforms that are indistinguishable from the available sequencing evidence.
- RSEM, Salmon, kallisto, and oarfish use related probabilistic approaches, but their implementation details differ.
