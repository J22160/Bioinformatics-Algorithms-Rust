# UMI Deduplication and Error Correction

## What Is UMI Deduplication?

Single-cell RNA sequencing (scRNA-seq) measures gene expression by sequencing fragments derived from RNA molecules inside individual cells.

However, the amount of RNA available from a single cell is extremely small. To generate enough material for sequencing, the molecules are typically reverse-transcribed and amplified through polymerase chain reaction (PCR).

PCR amplification produces multiple copies of the same original molecule. Consequently, the number of sequencing reads generated from a gene does not necessarily equal the number of original RNA molecules.

For example, suppose a cell contains 10 original transcripts of a particular gene. After PCR amplification, these molecules might generate 1,000 sequencing reads. If every read is counted independently, the gene would appear to have 1,000 molecules instead of the original 10.

This introduces **PCR amplification bias**, where highly amplified molecules contribute disproportionately to the final read count.

```text
[ Original Single Cell RNA ] (e.g., 10 Molecules)
             |
             v
   [ Reverse Transcription ]
             |
             v
   [ PCR Amplification ]  <-- Introduces non-uniform copy counts per molecule
             |
             v
   [ Sequencing Output ]   (e.g., 1,000 Reads)

```

**Unique Molecular Identifiers (UMIs)** address this problem by assigning a short, random nucleotide sequence to each original molecule before amplification. All amplified copies derived from the same tagged molecule inherit the same UMI.

Instead of counting every sequencing read, the analysis groups reads by their molecular identifiers and estimates the number of original molecules by counting distinct corrected UMIs. This process is called **UMI deduplication**.

UMI deduplication is a fundamental step in many droplet-based single-cell RNA-seq workflows, including 10x Genomics Chromium.

---

## 1. Understanding Unique Molecular Identifiers

A UMI is a short random nucleotide sequence attached to an original RNA-derived molecule before PCR amplification.

For example, consider a transcript from gene *GATA1*. During library preparation, the original molecule receives a UMI:

```text
Original RNA molecule
        |
        v
Gene: GATA1
UMI:  ACGTACGTAC
        |
        v
    PCR Amplification
        |
        +----------------+----------------+----------------+
        |                |                |                |
        v                v                v                v
     Read 1           Read 2           Read 3           Read 4
   ACGTACGTAC       ACGTACGTAC       ACGTACGTAC       ACGTACGTAC

```

All reads derived from this original molecule share the same UMI. Even though four sequencing reads are generated, they represent one original molecule.

Therefore:

$$\text{Read Count} = 4$$

$$\text{Estimated Molecule Count} = 1$$

### 1.1 Why Are UMIs Necessary?

Without UMIs, sequencing reads cannot reliably distinguish between:

* Multiple original RNA molecules from the same gene.
* Multiple PCR copies generated from a single original molecule.

Consider two genes:

| Gene | Original Molecules | Reads After PCR |
| --- | --- | --- |
| **Gene A** | 10 | 1,000 |
| **Gene B** | 50 | 500 |

If raw reads are counted directly:

* **Gene A** appears more highly expressed than **Gene B**.
* **Gene B** appears to have fewer molecules.

After UMI deduplication:

| Gene | Raw Reads | Unique UMIs | Estimated Molecules |
| --- | --- | --- | --- |
| **Gene A** | 1,000 | 10 | 10 |
| **Gene B** | 500 | 50 | 50 |

The corrected counts more closely reflect the number of captured molecules rather than the extent of PCR amplification.

> **Important:** UMI counts estimate *captured* molecules, not necessarily the total number of RNA molecules originally present in the biological cell. Capture efficiency, reverse-transcription efficiency, sequencing depth, and technical factors still influence observed counts.

---

## 2. UMI Structure and Sequencing

In droplet-based single-cell RNA-seq protocols, the sequencing library contains three core pieces of information:

1. **Cell Barcode (CB):** Identifies the cell or droplet associated with the molecule.
2. **UMI:** Identifies the original captured molecule.
3. **cDNA Sequence:** Contains the transcript-derived sequence used to identify the gene.

```text
+----------------+----------------+----------------------+
|  Cell Barcode  |      UMI       | Transcript Sequence  |
+----------------+----------------+----------------------+
|    Cell_001    |   ACGTACGTAC   |  GENE_A-derived cDNA |
+----------------+----------------+----------------------+

```

The cell barcode distinguishes molecules originating from different cells, while the UMI distinguishes molecules originating from the same cell.

For example:

| Read | Cell Barcode | UMI | Assigned Gene |
| --- | --- | --- | --- |
| **R1** | CB001 | ACGTACGTAC | GENE_A |
| **R2** | CB001 | ACGTACGTAC | GENE_A |
| **R3** | CB001 | ACGTACGTAC | GENE_A |
| **R4** | CB001 | TGCATGCATG | GENE_A |
| **R5** | CB001 | TGCATGCATG | GENE_A |
| **R6** | CB002 | ACGTACGTAC | GENE_A |

The interpretation is:

* R1, R2, and R3 represent one molecule in CB001.
* R4 and R5 represent a second molecule in CB001.
* R6 represents a molecule from a different cell (CB002).

The same UMI sequence can legitimately occur in different cells. UMIs must always be evaluated in the context of cell barcode and gene assignment.

---

## 3. The UMI Deduplication Problem

Ideally, every original molecule receives a unique UMI, and all its PCR copies retain the exact same sequence. In practice, sequencing and amplification errors introduce mutations into UMI sequences.

### Example: Sequencing Error in a UMI

Suppose an original molecule receives:

$$\text{Original UMI: } \texttt{ACGTACGTAC}$$

After PCR amplification, several copies are sequenced:

| Read | Observed UMI | Read Count |
| --- | --- | --- |
| **R1** | ACGTACGTAC | 1 |
| **R2** | ACGTACGTAC | 1 |
| **R3** | ACGTACGTAC | 1 |
| **R4** | ACGTACGTAC | 1 |
| **R5** | ACGTACGTAT | 1 |

The two observed UMIs differ by only one base:

```text
ACGTACGTAC
ACGTACGTAT
         ^ (Mutation / Error)

```

If exact UMI matching is used, the algorithm identifies two distinct UMIs, overestimating molecule count. Therefore, robust UMI deduplication requires **error correction**.

---

## 4. Hamming Distance and UMI Error Detection

A common way to measure similarity between UMIs is the **Hamming distance** ($d_H$), which counts the number of position differences between two equal-length sequences.

For two UMIs:

$$U_1 = \texttt{ACGTACGTAC}, \quad U_2 = \texttt{ACGTACGTAT}$$

$$d_H(U_1, U_2) = 1$$

Another example:

$$U_1 = \texttt{ACGTACGTAC}, \quad U_3 = \texttt{ACGTTCGTAT}$$

$$d_H(U_1, U_3) = 2$$

UMIs with small Hamming distances are candidates for representing the same original molecule. However, two independent molecules can randomly receive UMIs that differ by one nucleotide (a **UMI collision** or near-collision). Blindly merging all UMIs within a Hamming distance of $1$ can cause over-deduplication.

---

## 5. Directional UMI Clustering: UMI-tools

The `UMI-tools` package provides a widely used method called the **directional adjacency method**.

It combines two sources of information:

1. **Sequence Similarity:** Hamming distance between UMIs.
2. **Read Abundance:** Number of reads supporting each UMI.

The central heuristic assumes that an abundant UMI represents a true original molecule, whereas a closely related low-abundance UMI arose from a sequencing or PCR error.

```text
+-------------------------------------------------------+
|  Step 1: Group reads by Cell Barcode + Gene Tag       |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|  Step 2: Build UMI Adjacency Graph                    |
|  - Nodes = Observed UMIs                              |
|  - Edges = Connected if d_H(U_i, U_j) <= 1            |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|  Step 3: Apply Directional Criterion                  |
|  - Check if n_A >= 2*n_B - 1                          |
|  - Absorb low-count nodes into parent clusters        |
+-------------------------------------------------------+
                           |
                           v
+-------------------------------------------------------+
|  Step 4: Count Corrected Molecule Clusters            |
+-------------------------------------------------------+

```

### 5.1 Step 1: Group Reads Before UMI Correction

Reads are grouped by molecular context (e.g., Cell Barcode + Gene Tag + Alignment Position):

| UMI | Number of Reads |
| --- | --- |
| **ACGTACGTAC** | 100 |
| **ACGTACGTAT** | 2 |
| **ACGTACGTAG** | 1 |
| **TGCATGCATG** | 45 |
| **TGCATGCATA** | 3 |
| **GATCGATCGA** | 20 |

### 5.2 Step 2: Build the UMI Adjacency Graph

Each observed UMI is a node in graph $G=(V,E)$, where an edge exists if $d_H(U_i, U_j) \le 1$.

```text
             100 reads
            ACGTACGTAC
             /       \
            /         \
           /           \
   2 reads               1 read
ACGTACGTAT             ACGTACGTAG

```

### 5.3 Step 3: Apply Directional Clustering

A directional edge is drawn from a higher-abundance UMI ($A$) to a lower-abundance UMI ($B$) if their read counts satisfy the directional condition:

$$n_A \ge 2n_B - 1$$

where:

* $n_A$ is the read count of the more abundant UMI.
* $n_B$ is the read count of the less abundant UMI.

#### Evaluation Example:

For $n_A = 100$ and $n_B = 2$:

$$100 \ge 2(2) - 1 \implies 100 \ge 3 \quad (\text{Satisfied})$$

For $n_A = 100$ and $n_C = 1$:

$$100 \ge 2(1) - 1 \implies 100 \ge 1 \quad (\text{Satisfied})$$

Both low-abundance UMIs are absorbed into the cluster seeded by `ACGTACGTAC`.

### 5.4 Step 4: Follow Directional Paths

Directional clustering can absorb UMIs transitively through intermediate nodes:

$$\texttt{UMI-A (100 reads)} \longleftrightarrow \texttt{UMI-B (10 reads)} \longleftrightarrow \texttt{UMI-C (2 reads)}$$

Checking relationships:

1. $100 \ge 2(10) - 1 \implies 100 \ge 19$ (Valid)
2. $10 \ge 2(2) - 1 \implies 10 \ge 3$ (Valid)

Even if $d_H(\text{UMI-A}, \text{UMI-C}) = 2$, UMI-C is merged into the cluster via UMI-B.

### 5.5 Step 5: Count Corrected Molecules

Each resolved cluster equals one estimated original transcript:

| Cluster | Observed UMIs | Corrected Molecules |
| --- | --- | --- |
| **1** | ACGTACGTAC, ACGTACGTAT, ACGTACGTAG | 1 |
| **2** | TGCATGCATG, TGCATGCATA | 1 |
| **3** | GATCGATCGA | 1 |

Total estimated molecules:

$$N_{\text{molecules}} = N_{\text{corrected UMI clusters}} = 3$$

---

## 6. Comparison of Deduplication Methods

| Method | Uses Similarity? | Uses Read Abundance? | Key Characteristic |
| --- | --- | --- | --- |
| **Exact** | No | No | Fast; fails to correct UMI sequencing errors. |
| **Adjacent** | Yes ($d_H \le 1$) | No | Merges neighboring UMIs; can over-cluster independent molecules. |
| **Directional** | Yes ($d_H \le 1$) | Yes ($n_A \ge 2n_B - 1$) | Uses abundance to protect against false merging. |
| **Directional-Adjacent** | Yes | Yes | Combines directional rules with adjacency graph connected components. |

---

## 7. UMI Collisions and Limitations

For a UMI of length $L$ composed of 4 nucleotides ($A, C, G, T$), the total sequence space $M$ is:

$$M = 4^L$$

* **10-bp UMI:** $M = 4^{10} = 1,048,576$
* **12-bp UMI:** $M = 4^{12} = 16,777,216$

If $n$ independent molecules are sampled from space $M$, the expected number of colliding pairs is approximated by the Birthday Paradox:

$$E[\text{colliding pairs}] \approx \frac{n(n-1)}{2M}$$

For $n = 1,000$ molecules with $L = 10$:

$$E[\text{colliding pairs}] \approx \frac{1000 \times 999}{2 \times 1,048,576} \approx 0.48$$

---

## 8. Computational Complexity

Let:

* $R$ = Total sequencing reads
* $G$ = Number of grouping contexts (Cell Barcode + Gene)
* $U_g$ = Unique UMIs in group $g$
* $L$ = UMI sequence length
* $E_g$ = Edges in UMI graph $g$

| Step | Algorithm / Task | Time Complexity | Memory Complexity |
| --- | --- | --- | --- |
| **1. Read Hashing** | Count unique UMIs | $\mathcal{O}(R)$ | $\mathcal{O}(U)$ |
| **2. Graph Build (Naive)** | Pairwise Hamming comparisons | $\mathcal{O}(G \cdot U_g^2 L)$ | $\mathcal{O}(U + E)$ |
| **3. Graph Build (Optimized)** | $3L$ candidate neighbor lookup | $\mathcal{O}(G \cdot U_g \cdot 3L)$ | $\mathcal{O}(U + E)$ |
| **4. Directional Traversal** | Graph clustering per group | $\mathcal{O}(U_g + E_g)$ | $\mathcal{O}(U_g + E_g)$ |

---

## 9. Practical Implementation with UMI-tools

### UMI Extraction

Extract cell barcodes and UMIs from raw FASTQ files before alignment:

```bash
umi_tools extract \
    --bc-pattern=CCCCCCCCCCCCNNNNNNNNNNNN \
    --stdin=R1.fastq.gz \
    --stdout=R1_extracted.fastq.gz \
    --read2-in=R2.fastq.gz \
    --read2-out=R2_extracted.fastq.gz

```

### UMI Deduplication

Deduplicate aligned BAM files using the directional method:

```bash
umi_tools dedup \
    --stdin=aligned_reads.bam \
    --stdout=deduplicated_reads.bam \
    --method=directional \
    --per-cell \
    --gene-tag=GX \
    --per-gene

```

---

## 10. From Deduplication to Gene Expression Matrix

After UMI deduplication, corrected molecule counts form the digital expression matrix:

$$X \in \mathbb{N}^{G \times C}$$

| Gene | Cell 1 | Cell 2 | Cell 3 |
| --- | --- | --- | --- |
| **GATA1** | 10 | 0 | 4 |
| **CD3D** | 0 | 25 | 2 |
| **MS4A1** | 1 | 0 | 30 |
| **LYZ** | 5 | 12 | 3 |

This matrix serves as the foundation for library normalization, cell clustering, and differential expression analysis.

---

## 11. Key QC Metrics

* **Reads per Cell:** Measures raw sequencing depth per droplet.
* **UMIs per Cell:** Reflects captured molecular complexity.
* **Genes Detected per Cell:** Number of genes with $\ge 1$ UMI count.
* **Reads per UMI:** Quantifies average PCR duplication level.
* **UMI Correction Rate:** Fraction of observed UMIs merged due to errors.
* **Sequencing Saturation:** Estimates probability of discovering new molecules with additional sequencing reads.

---

## References

1. **Smith T, Heger A, Sudbery I.** (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Research*, 27(3):491–499. [https://doi.org/10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116)
2. **Islam S, Zeisel A, Joost S, et al.** (2014). Quantitative single-cell RNA-seq with unique molecular identifiers. *Nature Methods*, 11(2):163–166. [https://doi.org/10.1038/nmeth.2772](https://doi.org/10.1038/nmeth.2772)
3. **Kivioja T, Vähärautio A, Karlsson K, et al.** (2012). Counting absolute numbers of molecules using unique molecular identifiers. *Nature Methods*, 9(1):72–74. [https://doi.org/10.1038/nmeth.1778](https://doi.org/10.1038/nmeth.1778)
4. **Zheng GXY, Terry JM, Belgrader P, et al.** (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8:14049. [https://doi.org/10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049)
5. **10x Genomics.** Single Cell Gene Expression Documentation. [https://www.10xgenomics.com/support](https://www.10xgenomics.com/support)