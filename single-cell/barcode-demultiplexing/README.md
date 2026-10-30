# Barcode Demultiplexing and Error Correction in Single-Cell RNA-seq

## Table of Contents
1. [What Is This?](#what-is-this)
2. [The Barcode Demultiplexing Problem](#the-barcode-demultiplexing-problem)
3. [Sources of Barcode Errors](#sources-of-barcode-errors)
4. [The Barcode Whitelist](#the-barcode-whitelist)
5. [How Barcode Correction Works](#how-barcode-correction-works)
6. [Cell Barcode Correction vs. Cell Calling](#cell-barcode-correction-vs-cell-calling)
7. [Computational Complexity](#computational-complexity)
8. [Efficient Data Structures and Algorithms](#efficient-data-structures-and-algorithms)
9. [Bioinformatics Problem](#bioinformatics-problem)
10. [Limitations and Best Practices](#limitations-and-best-practices)
11. [Connection to Textbook Concepts](#connection-to-textbook-concepts)
12. [Key Takeaway](#key-takeaway)
13. [References](#references)

---

## What Is This?

Barcode demultiplexing is the process of identifying the cell or sample of origin associated with each sequencing read using an attached nucleotide barcode.

In droplet-based single-cell RNA sequencing (scRNA-seq), such as the 10x Genomics Chromium platform, thousands of individual cells are encapsulated into tiny oil droplets along with gel beads. Each gel bead contains millions of copies of barcode oligonucleotides.

These oligonucleotides contain several functional components, including:

- **Cell barcode (CB):** Identifies the droplet or cell-associated partition from which a molecule originated.
- **Unique molecular identifier (UMI):** Identifies an individual captured RNA molecule, helping distinguish original molecules from PCR duplicates.
- **Capture sequence:** Facilitates capture and reverse transcription of RNA molecules.
- **Sequencing adapter or primer sequences:** Enable amplification and sequencing.

During library preparation, mRNA molecules are captured and reverse-transcribed into cDNA. The cell barcode and UMI are incorporated into the resulting cDNA molecules.

After sequencing, reads contain information about the cell barcode, UMI, and transcript sequence. The barcode allows reads to be grouped by their originating partition, while the UMI helps estimate the number of original RNA molecules rather than simply counting sequencing reads.

### Simplified sequencing read structure

A simplified 10x Genomics 3' gene expression read structure is:

```text
Read 1:
[ Cell Barcode ][ UMI ]

Read 2:
[ cDNA / Transcript Sequence ]

```

For a typical 10x Chromium 3' gene expression library, the cell barcode is 16 bases long and the UMI is 12 bases long. Exact structures and lengths depend on the assay chemistry and protocol.

For example:

```text
Cell barcode:  ACGTACGTACGTACGT
UMI:           TTGCAAGTCCGA
Transcript:    ATGGCTAGCTAC...

```

The cell barcode identifies the cell-associated partition, while the UMI identifies a captured molecule within that partition.

> **Barcode demultiplexing is an essential early step in single-cell RNA-seq processing.** It assigns reads to valid barcode identities before downstream operations such as UMI processing, transcript alignment, gene quantification, and cell calling.

---

## The Barcode Demultiplexing Problem

Cell barcodes are designed to be selected from a predefined set of valid sequences called a **barcode whitelist**.

For example, a 10x Chromium whitelist contains hundreds of thousands of possible barcode sequences, depending on the chemistry and assay version.

However, sequencing is not error-free. During sequencing, a nucleotide can be incorrectly identified, resulting in an observed barcode that differs from its original valid sequence.

This creates two major challenges:

1. Determine whether an observed barcode corresponds to a valid barcode in the whitelist.
2. If it does not match exactly, determine whether it can be reliably corrected to a valid barcode.

### Example: Single-base substitution

Suppose the original barcode is:

```text
True barcode:      ACGTACGTACGTACGT
Observed barcode:  ACGTACGTTCGTACGT
                             ^
                         Substitution

```

The observed barcode differs from the original barcode at one position.

The Hamming distance between two equal-length sequences is the number of positions at which they differ:

$$d_H(A,B) = \sum_{i=1}^{L} \mathbb{1}(A_i \neq B_i)$$

Where:

* $A$ and $B$ are two barcode sequences.
* $L$ is the barcode length.
* $\mathbb{1}$ is an indicator function that equals 1 when the bases differ and 0 otherwise.

For the example:

$$d_H(A,B) = 1$$

The demultiplexing algorithm can potentially correct the observed barcode to the valid whitelist sequence.

However, not every unmatched sequence should be corrected. Some observed sequences may originate from sequencing errors, barcode synthesis errors, background molecules, or other artifacts.

---

## Sources of Barcode Errors

Barcode mismatches can arise from several sources:

1. **Sequencing errors:** Instruments occasionally report an incorrect nucleotide via substitution.
2. **Insertions:** An extra nucleotide may be introduced into the observed sequence (changes sequence length; not corrected by ordinary Hamming distance).
3. **Deletions:** A nucleotide may be missing from the observed sequence (changes sequence length).
4. **Low-quality base calls:** Mismatches at positions with low base quality scores are more likely to represent sequencing errors.

---

## The Barcode Whitelist

A barcode whitelist is a predefined collection of valid barcode sequences for a particular single-cell assay.

```text
ACGTACGTACGTACGT
TGCATGCATGCATGCA
GATCGATCGATCGATC
CTAGCTAGCTAGCTAG

```

It is important to use the correct whitelist for the sequencing chemistry and assay version. Different assays use different barcode sets.

---

## How Barcode Correction Works

A simplified whitelist-based barcode correction workflow consists of several stages:

### Step 1: Extract the cell barcode

The barcode sequence is extracted from the appropriate read and position according to the assay-specific library structure.

### Step 2: Exact whitelist matching

Check whether the observed barcode exactly matches a valid barcode in the whitelist. A hash set or hash map provides efficient exact lookup:

```python
whitelist = {
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA",
    "GATCGATCGATCGATC"
}

observed = "ACGTACGTACGTACGT"

if observed in whitelist:
    corrected_barcode = observed
else:
    corrected_barcode = None

```

### Step 3: Generate candidate barcodes with one substitution

If an observed barcode does not match the whitelist, candidate sequences are generated by substituting each nucleotide position with the other three possible bases ($3L$ possible variants). For a 16-base barcode, $3 \times 16 = 48$ candidates are generated.

### Step 4: Evaluate base quality scores

Phred quality scores represent the estimated probability of a base-call error:

$$Q = -10 \log_{10}(P_{\text{error}})$$

$$P_{\text{error}} = 10^{-Q/10}$$

### Step 5: Calculate correction probabilities

A simplified Bayesian framework evaluates the posterior probability of a candidate barcode:

$$P(w \mid b, Q) \propto P(b \ w, Q) P(w)$$

### Step 6: Assign, correct, or reject

* **Exact match:** Observed barcode matches a whitelist entry.
* **Corrected match:** A supported sequencing error is corrected to a valid barcode.
* **Ambiguous:** Multiple candidates or insufficient evidence prevent a reliable assignment.
* **Unmatched:** No acceptable whitelist candidate is found.

---

## Cell Barcode Correction vs. Cell Calling

* **Barcode correction:** Determines which valid whitelist barcode most likely produced an observed barcode sequence (corrects sequencing errors).
* **Cell calling:** Determines which barcode-associated partitions represent genuine cells rather than empty droplets or background partitions (e.g., using knee plots and ambient RNA filtering).

---

## Computational Complexity

| Operation | Approximate Cost | Explanation |
| --- | --- | --- |
| Whitelist construction | $O(W L)$ | Store $W$ sequences of length $L$ |
| Exact barcode lookup | $O(1)$ to $O(L)$ | Expected constant time via hash set |
| Single-substitution candidate generation | $O(L^2)$ or $O(L)$ | Testing $3L$ candidates against whitelist |
| Total complexity (worst-case for $N$ reads) | $O(N L^2)$ | Assuming unmatched reads undergo candidate checks |

---

## Efficient Data Structures and Algorithms

1. **Hash Set:** Fast exact lookup $O(1)$.
2. **Trie:** Tree of nucleotide characters supporting bounded Hamming-distance searches.
3. **BK-tree:** Metric tree designed for searching objects within a specified distance.
4. **Packed Nucleotide Encoding:** Representing 2 bits per base (A=00, C=01, G=10, T=11) packs a 16-base barcode into a 32-bit integer for ultra-fast bitwise comparisons.

---

## Bioinformatics Problem

**Scenario:** You have 500 million sequencing reads from a 10x Chromium single-cell RNA-seq experiment with a 16-base cell barcode, 12-base UMI, and a whitelist of 737,280 valid barcodes.

### Recommended processing workflow:

1. Load and encode the barcode whitelist.
2. Extract cell barcode and UMI per read.
3. Perform exact whitelist matching.
4. Fall back to Hamming distance 1 substitution checks for unmatched barcodes.
5. Apply quality score filters and correct supported barcodes.
6. Group reads by cell barcode and UMI for gene quantification and cell calling.

---

## Limitations and Best Practices

* **Avoid correcting ambiguous barcodes:** If multiple whitelist candidates tie, mark as unassigned.
* **Use the correct whitelist:** Mismatched chemistries cause high drop rates.
* **Incorporate base quality:** Guard against high-confidence sequencing errors.
* **Account for multiplets and ambient RNA:** Barcode correction does not eliminate doublets or free-floating RNA.

---

## Key Takeaway

Barcode demultiplexing converts noisy barcode reads into valid cell-associated assignments via efficient exact and approximate dictionary matching. However, barcode correction is distinct from cell calling; filtering empty droplets requires dedicated statistical modeling of UMI distributions.

---

## References

1. Zheng GXY, Terry JM, Belgrader P, et al. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8, 14049.
2. Macosko EZ, Basu A, Satija R, et al. (2015). Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. *Cell*, 161(5), 1202–1214.
3. 10x Genomics documentation: https://www.10xgenomics.com/support
4. Lun ATL, Riesenfeld S, Andrews T, et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology*, 20, 63.
5. Hamming RW. (1950). Error detecting and error correcting codes. *Bell System Technical Journal*, 29(2), 147–160.
6. Burkhard WA, Keller RM. (1973). Some approaches to best-match file searching. *Communications of the ACM*, 16(4), 230–236.

