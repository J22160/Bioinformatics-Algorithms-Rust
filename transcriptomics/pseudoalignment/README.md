# Pseudoalignment for Transcript Quantification

A comprehensive guide to understanding pseudoalignment, k-mer indexing, equivalence classes, and lightweight mapping in RNA-seq transcript quantification.

## Table of Contents

 1. [What Is Pseudoalignment?](#1-what-is-pseudoalignment)
 2. [Why Is Pseudoalignment Useful?](#2-why-is-pseudoalignment-useful)
 3. [How Pseudoalignment Works](#3-how-pseudoalignment-works)
    * [Step 1: Build an index of transcript sequences](#step-1-build-an-index-of-transcript-sequences)
    * [Step 2: Determine transcript compatibility for each read](#step-2-determine-transcript-compatibility-for-each-read)
    * [Step 3: Group reads into equivalence classes](#step-3-group-reads-into-equivalence-classes)
    * [Step 4: Estimate transcript abundance](#step-4-estimate-transcript-abundance)
 4. [Worked Example](#4-worked-example)
 5. [Why Is Pseudoalignment Fast?](#5-why-is-pseudoalignment-fast)
 6. [Computational Complexity](#6-computational-complexity)
 7. [Pseudoalignment vs. Conventional Alignment](#7-pseudoalignment-vs-conventional-alignment)
 8. [Pseudoalignment vs. Quasi-Mapping](#8-pseudoalignment-vs-quasi-mapping)
 9. [Applications in Bioinformatics](#9-applications-in-bioinformatics)
10. [Limitations and Important Considerations](#10-limitations-and-important-considerations)
* [Key Takeaways](#key-takeaways)

---

## 1. What Is Pseudoalignment?

**Pseudoalignment** is a computational technique used in RNA-seq analysis to determine which transcripts are compatible with a sequencing read without performing a conventional base-by-base alignment.

Traditional alignment attempts to place a read onto a reference sequence, often determining the exact position and alignment of individual bases while accounting for mismatches, insertions, and deletions.

Pseudoalignment instead focuses on a more targeted question:

> **Which transcripts could have generated this read?**

For transcript quantification, knowing the set of compatible transcripts is often sufficient to estimate their relative abundances. Exact genomic coordinates and base-level alignments are not required.

Consider three transcript isoforms:

```text
Transcript A: ATGCCAGT
Transcript B: ATGCTAGT
Transcript C: ATGCCATT
```

* A read containing a sequence shared by all three transcripts is compatible with every isoform.
* A read containing a distinguishing sequence is compatible with only one or two isoforms.

Pseudoalignment identifies these compatibility relationships efficiently, passing the resulting evidence to downstream statistical quantification algorithms.

### A Simple Analogy

Imagine sorting incoming mail by postal code rather than reading and verifying every character of each recipient's full street address.

* **Full Alignment:** Read, verify, and map the detailed placement of every letter/base in the address.
* **Pseudoalignment:** Read key routing codes (sequence signatures) to route the letter directly to candidate regional destinations (transcripts).

While not identical, the analogy captures the core principle: **avoid performing computationally expensive work that is unnecessary for the final analytical goal.**

---

## 2. Why Is Pseudoalignment Useful?

Modern RNA-seq experiments routinely generate tens or hundreds of millions of reads per sample. Conventional base-by-base alignment of every read is computationally expensive, especially when the goal is simply quantifying transcript expression.

Pseudoalignment offers significant practical advantages:

* **Speed:** Bypasses dynamic programming matrices and per-base gap/mismatch scoring for compatibility detection.
* **Efficiency:** Rapidly narrows down candidate transcripts for each read using pre-built indices.
* **Natural Ambiguity Handling:** Represents multi-mapping reads directly as compatibility sets.
* **Data Reduction:** Enables aggregation of millions of reads with identical compatibility patterns into compact equivalence classes.
* **Integration with EM:** Output compatibility structures interface directly with Expectation-Maximization (EM) and Variational Bayes abundance estimators.

> **Note:** Pseudoalignment is not a universal replacement for conventional alignment. Tasks requiring precise genomic coordinates, allele-specific variant calls, splice-junction discovery, or indel detection still require traditional full alignment algorithms.

---

## 3. How Pseudoalignment Works

A simplified pseudoalignment workflow consists of four primary stages:

```text
  +-------------------------------------------------------+
  | Step 1: Build Transcriptome Index                     |
  | Parse reference transcripts into k-mer lookup tables  |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  | Step 2: Determine Read Compatibility                  |
  | Intersect k-mer transcript sets for each read         |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  | Step 3: Construct Equivalence Classes                 |
  | Group reads with identical transcript compatibility   |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  | Step 4: Estimate Abundance (e.g., EM / VB)            |
  | Infer transcript proportions from class counts        |
  +-------------------------------------------------------+
```

### Step 1: Build an Index of Transcript Sequences

The input is a reference transcriptome. Sequences are broken down into **k-mers** (subsequences of fixed length $k$).

For example, the sequence `ATGCCAGT` contains the following 3-mers ($k=3$):

```text
ATG, TGC, GCC, CCA, CAG, AGT
```

An index structure (e.g., a hash table or a colored de Bruijn graph) maps each k-mer to the set of transcripts containing it:

| k-mer | Compatible Transcripts |
|---|---|
| **ATG** | Transcript A, Transcript B |
| **TGC** | Transcript A, Transcript B |
| **GCC** | Transcript A |
| **CTA** | Transcript B |
| **CCA** | Transcript A, Transcript C |

### Step 2: Determine Transcript Compatibility for Each Read

For a given read, sequence signatures (k-mers) are queried against the index. The intersection of transcript sets across all k-mers in the read defines its overall compatibility set.

Suppose the index yields:

```text
Read k-mer 1  -->  {Transcript A, Transcript B, Transcript C}
Read k-mer 2  -->  {Transcript A, Transcript B}
Read k-mer 3  -->  {Transcript A}
```

Taking the set intersection narrows the candidate transcripts:

$$
\{A, B, C\} \cap \{A, B\} \cap \{A\} = \{A\}
$$

Under this simplified model, **Transcript A** is the only candidate compatible with all three k-mers.

If a read's sequence matches regions shared across multiple transcripts, its final set will contain multiple candidates (an ambiguous read):

```text
Read 1  -->  {Transcript A}
Read 2  -->  {Transcript A, Transcript B}
Read 3  -->  {Transcript B, Transcript C}
Read 4  -->  {Transcript A, Transcript B, Transcript C}
```

### Step 3: Group Reads into Equivalence Classes

Many reads map to the exact same subset of transcripts. Rather than processing millions of reads individually during abundance estimation, software aggregates them into **Equivalence Classes (ECs)**.

| Equivalence Class (EC) | Compatible Transcripts | Read Count |
|---|---|---:|
| **EC1** | {Transcript A} | 1,200 |
| **EC2** | {Transcript A, Transcript B} | 3,500 |
| **EC3** | {Transcript B, Transcript C} | 800 |
| **EC4** | {Transcript A, Transcript B, Transcript C} | 2,100 |

*Equivalence classes store read count totals for each distinct compatibility pattern, compressing the dataset significantly without losing quantification evidence.*

### Step 4: Estimate Transcript Abundance

Equivalence class counts are passed to a probabilistic quantifier (such as EM). The quantifier calculates fractional expected counts based on current abundance parameters, updating estimates iteratively until convergence.

* **Pseudoalignment:** Determines *which* transcripts are compatible with each read.
* **EM / Inference Engine:** Determines *how much* expression to assign to each transcript.

---

## 4. Worked Example

Consider a sample with four reads and a transcriptome containing three transcript isoforms.

### Reference Transcripts

```text
Transcript A: ATGCCAGT
Transcript B: ATGCTAGT
Transcript C: ATGCCATT
```

Assume the indexing/compatibility step produces the following mapping sets:

```text
Read 1  -->  {A}
Read 2  -->  {A, B}
Read 3  -->  {B, C}
Read 4  -->  {A, B}
```

### Step 1: Construct Equivalence Classes

Reads 2 and 4 share identical transcript compatibility sets `{A, B}` and are grouped together into `EC2`:

| Equivalence Class | Included Reads | Compatible Transcripts | Count |
|---|---|---|---:|
| **EC1** | Read 1 | {A} | 1 |
| **EC2** | Read 2, Read 4 | {A, B} | 2 |
| **EC3** | Read 3 | {B, C} | 1 |

*Four reads are compressed into three equivalence classes.*

### Step 2: Interpret the Compatibility Evidence

* `EC1` provides unique evidence for **Transcript A**.
* `EC2` provides ambiguous evidence for **Transcript A** or **Transcript B**.
* `EC3` provides ambiguous evidence for **Transcript B** or **Transcript C**.

### Step 3: Probabilistic Quantification

A quantification model allocates fractional counts for ambiguous reads based on transcript lengths, relative abundances, and read likelihoods.

Suppose an EM iteration assigns expected contributions as follows:
* The 2 reads in `EC2` contribute $0.6$ each to A and $0.4$ each to B.
* The 1 read in `EC3` contributes $0.7$ to B and $0.3$ to C.

| Transcript | Unique Contribution | Ambiguous Contribution | Expected Total Count |
|---|---:|---:|---:|
| **Transcript A** | $1.0$ (EC1) | $2 \times 0.6 = 1.2$ (EC2) | **2.2** |
| **Transcript B** | $0.0$ | $(2 \times 0.4) + (1 \times 0.7) = 1.5$ (EC2, EC3) | **1.5** |
| **Transcript C** | $0.0$ | $1 \times 0.3 = 0.3$ (EC3) | **0.3** |
| **Total** | **1.0** | **3.0** | **4.0** |

---

## 5. Why Is Pseudoalignment Fast?

1. **Indexed Lookup Tables:** Replaces full-read alignment loops with rapid $k$-mer lookups.
2. **Avoidance of Per-Base Scoring:** Omits Smith-Waterman / Needleman-Wunsch dynamic programming calculations.
3. **Data Compression via Equivalence Classes:** Reduces the parameter estimation space from millions of reads $N$ to thousands of unique equivalence classes $E$.
4. **Optimized Data Structures:** Utilizes compact colored de Bruijn graphs (kallisto) or high-efficiency space-saving indices (Salmon).

---

## 6. Computational Complexity

Let:

* $T$ = Total bases in reference transcriptome
* $N$ = Total number of reads
* $L$ = Average read length
* $k$ = k-mer length
* $M$ = Number of reference transcripts
* $K$ = Average number of compatible transcripts per read
* $E$ = Total number of distinct equivalence classes

### Complexity Summary

| Process | Conceptual Complexity | Description |
|---|---|---|
| **Index Construction** | $\mathcal{O}(T)$ | Linear in transcriptome size for basic structures (graph/suffix construction may vary). |
| **Per-Read Query** | $\mathcal{O}(L - k + 1)$ | Proportional to the number of k-mer windows in a read of length $L$. |
| **All Reads Processing** | $\mathcal{O}(N \cdot (L - k + 1))$ | Total lookup time for $N$ reads. |
| **Quantification Iteration** | $\mathcal{O}(I \cdot E \cdot K)$ | Iterating over $E$ equivalence classes across $I$ EM loops. |

---

## 7. Pseudoalignment vs. Conventional Alignment

| Feature | Conventional Alignment | Pseudoalignment |
|---|---|---|
| **Primary Objective** | Find exact genomic/transcriptomic coordinates for every base | Identify the set of candidate generating transcripts |
| **Base-Level Coordinates** | Provided (CIGAR strings, alignment positions) | Not computed or reported |
| **Indels & Mismatches** | Explicitly modeled via gap penalties and substitution matrices | Managed via k-mer lookup thresholds or graph traversals |
| **Downstream Applications** | Variant calling, splice-junction discovery, allele expression | Rapid transcript abundance estimation |
| **Computational Overhead** | High (memory and CPU intensive) | Low (orders of magnitude faster) |

---

## 8. Pseudoalignment vs. Quasi-Mapping

While *pseudoalignment* and *quasi-mapping* are both lightweight mapping paradigms designed for rapid transcript quantification, they differ in implementation and spatial detail:

```text
   Methodology Comparison:
   
   [ Read Sequence ]
          |
          +---> kallisto (Pseudoalignment)
          |     • Colored de Bruijn Graph index
          |     • Identifies transcript compatibility sets directly
          |
          +---> Salmon (Quasi-Mapping / Selective Alignment)
                • Suffixed-array / Generalized FM-index
                • Tracks candidate mapping positions & strand orientation
```

### Tool Comparison

| Feature | kallisto | Salmon |
|---|---|---|
| **Core Concept** | Pseudoalignment | Quasi-mapping / Selective Alignment |
| **Index Structure** | Colored de Bruijn Graph (cdBG) | Suffix Array / FM-index |
| **Positional Information** | Does not track alignment positions | Retains approximate mapping coordinates & strand orientation |
| **Bias Modeling** | Basic length normalization | Detailed positional, sequence, and GC-content bias modeling |
| **Primary Output** | Transcript abundance (counts, TPM) | Transcript abundance (counts, TPM) |

---

## 9. Applications in Bioinformatics

* **Transcript-Level Expression:** Estimating relative expression across alternative splice isoforms.
* **Large-Scale RNA-Seq Datasets:** Re-analyzing vast public repositories (e.g., TCGA, GTEx, SRA) rapidly.
* **Single-Cell RNA-Seq (scRNA-seq):** Quantifying UMI and cell-barcoded transcriptomes efficiently (e.g., via `kb-python`).
* **Differential Transcript Usage (DTU):** Providing input counts and variance estimates for isoform-level differential analysis tools.

---

## 10. Limitations and Important Considerations

1. **Shared Sequence Ambiguity:** Reads mapping to identical regions shared across multiple isoforms cannot be resolved by pseudoalignment alone; they require statistical estimation during quantification.
2. **Sequencing Errors:** High error rates in reads can interrupt k-mer matches. Production pseudoaligners use skip steps, k-mer voting, or error-tolerant heuristics to mitigate this.
3. **Dependence on Reference Quality:** Pseudoalignment maps strictly to provided transcript sequences. Incomplete or incorrect annotations will misdirect compatibility assignments.
4. **Lack of Coordinates:** Not suitable for analyses requiring genomic coordinates, splice-site discovery, or structural variant detection.

---

## Key Takeaways

* **Goal-Oriented Efficiency:** Pseudoalignment answers *"Which transcripts could have generated this read?"* rather than finding exact base-level coordinates.
* **K-mer Lookups:** Uses indexed sequence signatures ($k$-mers) to rapidly construct transcript compatibility sets.
* **Equivalence Classes:** Groups reads with identical compatibility patterns, vastly reducing data volume prior to abundance estimation.
* **Separation of Concerns:** Pseudoalignment identifies compatibility; probabilistic frameworks (EM/VB) infer abundance.
* **Tool Ecosystem:** `kallisto` uses graph-based pseudoalignment, whereas `Salmon` uses quasi-mapping/selective alignment with mapping-position tracking.