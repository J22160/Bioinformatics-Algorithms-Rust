# RNA-seq Normalization

A comprehensive guide to understanding technical variation, scaling methods, composition bias, and implementation of Median-of-Ratios, TMM, and TPM.

## Table of Contents

 1. [What Is This?](#1-what-is-this)
 2. [The Main Sources of Variation](#2-the-main-sources-of-variation)
 3. [Method 1: Median of Ratios (DESeq2)](#3-method-1-median-of-ratios-deseq2)
 4. [Method 2: TMM (Trimmed Mean of M-values)](#4-method-2-tmm-trimmed-mean-of-m-values)
 5. [Method 3: TPM (Transcripts Per Million)](#5-method-3-tpm-transcripts-per-million)
 6. [Comparing the Methods](#6-comparing-the-methods)
 7. [Other Normalization Approaches](#7-other-normalization-approaches)
 8. [Practical Implementation](#8-practical-implementation)
 9. [Common Mistakes to Avoid](#9-common-mistakes-to-avoid)
10. [Bioinformatics Problem](#10-bioinformatics-problem)
11. [References](#11-references)

---

## 1. What Is This?

RNA sequencing (RNA-seq) quantifies gene expression by counting sequencing reads or fragments assigned to genes or transcripts. These counts are influenced not only by biological expression but also by technical factors, including sequencing depth, RNA composition, and transcript length.

**Normalization adjusts for systematic differences between samples so that expression measurements can be compared more appropriately.**

Consider a simple example:

| Gene | Sample A (10M reads) | Sample B (20M reads) | Raw Ratio (B / A) |
|---|---:|---:|---:|
| **Gene X** | 100 | 200 | 2.0x |
| **Gene Y** | 500 | 1,000 | 2.0x |
| **Gene Z** | 50 | 100 | 2.0x |

Sample B has twice as many total reads as Sample A. If the underlying relative expression profiles are identical, the counts for each gene may approximately double simply because more sequencing was performed.

Comparing raw counts directly would suggest that every gene is more highly expressed in Sample B. Normalization accounts for this systematic difference.

```text
Raw Count Data
     │
     ├── Depth Bias? ────────► [Library Size Adjustment]
     │
     ├── Composition Bias? ──► [Median-of-Ratios / TMM]
     │
     └── Length Bias? ───────► [TPM / Effective Length]
     │
     v
Normalized Data (Ready for downstream GLM or Visualization)
```

However, real RNA-seq datasets are more complicated. Some genes may change substantially between conditions, and a small number of highly expressed transcripts can dominate a sample's total reads. Consequently, normalization methods must account for more than total sequencing depth alone.

### Why Normalization Matters

Without appropriate normalization:
- Samples with deeper sequencing may appear to have higher expression across many genes.
- Differences in RNA composition can distort comparisons based on total library size.
- Differential expression tests may report misleading results if systematic sample-level differences are mistaken for biological effects.
- Expression estimates may be difficult to interpret across samples or experimental conditions.

> **Important Distinction:** Normalization does not remove genuine biological variation, correct every batch effect, or turn relative RNA-seq measurements into absolute numbers of RNA molecules per cell. Its purpose depends on the method and its underlying assumptions.

---

## 2. The Main Sources of Variation

### 1. Sequencing Depth
Sequencing depth refers to the total number of reads or fragments obtained from a sample. A sample with more usable reads generally has more opportunities to capture transcripts, leading to higher raw counts for many genes. Library-size normalization adjusts for differences in the amount of sequencing data available per sample.

### 2. RNA Composition
Samples can differ in the relative abundance of their transcripts. 

For example, suppose a treatment causes one transcript to become extremely abundant. That transcript may consume a large fraction of the available sequencing reads, leaving fewer reads for other genes. Even if the absolute abundance of another gene remains unchanged, its fraction of the total library will decrease.

This is why simply dividing counts by total library size can be misleading when RNA composition differs substantially between samples. Methods such as DESeq2's median-of-ratios and edgeR's TMM attempt to estimate appropriate relative scaling factors while reducing the influence of genes with extreme expression patterns.

### 3. Transcript Length
Longer transcripts generally provide more opportunities for sequencing fragments to originate from them than shorter transcripts, assuming similar transcript abundance.

Transcript-length normalization is important when comparing the abundance of different transcripts within a sample. However, for many gene-level differential expression analyses, the same gene is compared across samples. Its effective length is often similar across samples, so transcript length does not need to be explicitly corrected in the same way as library size.

Transcript length can still matter when isoform usage changes between samples, because the effective length of a gene's expressed transcript mixture may change.

### 4. Technical and Biological Variation
Normalization factors do not automatically correct all technical differences. Batch effects, differences in RNA quality, sample handling, and other experimental factors may require quality control, covariates in the statistical model, or other analytical strategies. Biological variability should generally be modeled rather than removed as if it were a technical artifact.

---

## 3. Method 1: Median of Ratios (DESeq2)

The **median-of-ratios method** is the default size-factor estimation approach in DESeq2. It estimates a relative scaling factor for each sample using the expression patterns of genes across the dataset.

Unlike simple total-count normalization, it is designed to reduce the influence of a small number of highly expressed or strongly changing genes.

### How It Works

Let:
- $K_{ij}$ be the raw count for gene $i$ in sample $j$.
- $G$ be the number of genes.
- $S$ be the number of samples.
- $s_j$ be the size factor estimated for sample $j$.

#### Step 1: Calculate a Pseudo-Reference for Each Gene
For each gene, calculate its geometric mean across all samples:

$$
g_i = \left(\prod_{j=1}^{S} K_{ij}\right)^{1/S}
$$

This geometric mean acts as a gene-specific reference value. Genes with a zero count in any sample have a geometric mean of zero and cannot contribute to the standard ratio calculation. DESeq2 handles this through its ratio-based estimation procedure; alternative size-factor estimation methods are available for datasets with many zeros.

#### Step 2: Calculate Gene-Wise Ratios
For each eligible gene and sample, calculate:

$$
r_{ij} = \frac{K_{ij}}{g_i}
$$

These ratios measure how the observed count in each sample compares with the reference for that gene.

#### Step 3: Estimate the Sample Size Factor
For each sample, take the median of the ratios across eligible genes:

$$
s_j = \text{median}_i(r_{ij})
$$

The estimated size factors are scaled so that their geometric mean is one.

#### Step 4: Normalize the Counts
Normalized counts can be calculated as:

$$
K_{ij}^{\text{norm}} = \frac{K_{ij}}{s_j}
$$

These normalized counts are useful for visualization and descriptive comparisons. For differential expression testing, DESeq2 incorporates the size factors directly into its negative binomial model rather than relying on normalized counts as input to a conventional statistical test.

### Worked Example

Suppose the estimated size factors are calculated for a gene with varying raw counts:

| Sample | Raw Count for Gene X | Size Factor ($s_j$) | Normalized Count ($K/s_j$) |
|---|---:|---:|---:|
| **Sample A** | 100 | 1.0 | 100 |
| **Sample B** | 200 | 2.0 | 100 |
| **Sample C** | 150 | 1.5 | 100 |

After normalization, the counts for Gene X are identical. The interpretation is that the observed count differences are consistent with the estimated sample-wide scaling differences rather than a gene-specific expression change.

### Pros & Cons

* **Advantages:**
  - Robust to a subset of genes with extreme counts or large expression changes.
  - Integrates seamlessly into DESeq2's negative binomial GLM.
  - Does not require equal total read counts across samples.
  - Performs well when most genes do not change substantially in the same direction.
* **Limitations:**
  - Depends on assumptions about the distribution of expression changes across genes.
  - Requires alternative estimation approaches when many genes have zero counts.
  - Can be problematic in datasets with global shifts in RNA abundance.
  - Does not correct batch effects or other confounders on its own.

---

## 4. Method 2: TMM (Trimmed Mean of M-values)

**TMM**, or Trimmed Mean of M-values, is the default normalization method used in edgeR.

Like the median-of-ratios method, TMM estimates relative scaling factors from gene expression patterns rather than simply using total library sizes. It is designed to reduce the impact of highly expressed genes and genes with extreme expression differences.

### How It Works

TMM compares each sample with a reference sample and calculates gene-wise measures of relative expression.

For gene $i$, a simplified log-ratio measure ($M_i$) is:

$$
M_i = \log_2\left( \frac{K_{i,\text{sample}} / N_{\text{sample}}}{K_{i,\text{reference}} / N_{\text{reference}}} \right)
$$

Where:
- $K_{i,\text{sample}}$ and $K_{i,\text{reference}}$ are gene counts.
- $N_{\text{sample}}$ and $N_{\text{reference}}$ are the corresponding library sizes.

TMM also calculates an average-expression measure ($A_i$), often called the A-value, for each gene.

### Algorithm Steps

```text
Select Reference Sample
         │
         v
Calculate M-values (log-ratios) & A-values (absolute expression)
         │
         v
Trim Extreme M-values (default top/bottom 30%) 
& Extreme A-values (default top/bottom 5%)
         │
         v
Calculate Precision-Weighted Mean of Remaining M-values
         │
         v
Compute Final Normalization Factor (f_j) -> Effective Library Size (N_j * f_j)
```

1. Select a reference sample (typically using automated selection).
2. Calculate log-ratios ($M$) and average-expression values ($A$) for genes in the sample-reference comparison.
3. Trim genes with extreme log-ratios (default: 30%) and extreme average-expression values (default: 5%).
4. Calculate a precision-weighted mean of the remaining log-ratios.
5. Convert the estimated mean log-ratio into a normalization factor $f_j$.
6. Combine $f_j$ with the original library size $N_j$ to obtain the effective library size:

$$
N_j^{\text{effective}} = N_j \cdot f_j
$$

In edgeR, these effective library sizes are used directly in the statistical model to account for differences in sequencing depth and RNA composition.

### Why Trimming Helps
Suppose a small number of genes become extremely abundant after treatment, dominating raw library sizes. If total-count normalization were used, the increase in those few genes would make the normalized counts of unchanged genes appear lower. TMM eliminates those extreme values prior to calculating the scaling factor.

---

## 5. Method 3: TPM (Transcripts Per Million)

**TPM**, or Transcripts Per Million, is a transcript-abundance measure that accounts for both transcript length and sequencing depth. It is commonly used for reporting and comparing relative transcript abundance within a sample.

Unlike DESeq2 size factors and TMM factors, TPM is not a sample-level scaling factor applied to raw gene counts.

### Why Transcript Length Matters

Suppose two transcripts have the same number of physical RNA molecules in a sample:
- **Transcript A:** 1,000 bases long.
- **Transcript B:** 4,000 bases long.

Transcript B has approximately four times as many possible fragment origin positions. Under comparable fragmentation and sequencing conditions, it will yield roughly four times as many reads despite having identical molecular abundance. TPM corrects for length prior to depth scaling.

### How TPM Is Calculated

Let:
- $C_i$ be the estimated fragment count for transcript $i$.
- $L_i$ be its effective length in kilobases.
- $R_i$ be its length-normalized abundance (Reads Per Kilobase, RPK).

#### Step 1: Calculate Reads Per Kilobase (RPK)

$$
R_i = \frac{C_i}{L_i / 1000}
$$

#### Step 2: Calculate Sum of RPK Values

$$
R_{\text{total}} = \sum_{i=1}^{T} R_i
$$

Where $T$ is the total number of quantified transcripts.

#### Step 3: Scale to One Million

$$
\text{TPM}_i = \left( \frac{R_i}{R_{\text{total}}} \right) \times 10^6
$$

By construction, TPM values sum to exactly one million within every sample.

### Worked Example

Suppose a sample contains two transcripts:

| Transcript | Estimated Fragments | Length (kb) | RPK ($C / L$) | TPM Calculation | Final TPM |
|---|---:|---:|---:|---|---:|
| **Transcript A** | 100 | 1.0 | 100 | $(100 / 200) \times 10^6$ | 500,000 |
| **Transcript B** | 400 | 4.0 | 100 | $(100 / 200) \times 10^6$ | 500,000 |

Both transcripts have equal TPM because Transcript B's fourfold higher read count is entirely explained by its fourfold greater length.

---

## 6. Comparing the Methods

| Feature | DESeq2 Median-of-Ratios | edgeR TMM | TPM |
|---|---|---|---|
| **Main Purpose** | Sample scaling for count-based GLM | Sample scaling for count-based GLM | Relative expression reporting |
| **Accounts for Sequencing Depth?** | Yes | Yes | Yes |
| **Accounts for Transcript Length?** | No (handled at model level if needed) | No (handled at model level if needed) | Yes |
| **Accounts for Composition Bias?** | Yes (robust median assumption) | Yes (robust trimming assumption) | No (forces sum to $10^6$) |
| **Primary Software** | DESeq2 | edgeR | Salmon, kallisto, RSEM |
| **Valid input to DESeq2/edgeR GLM?** | **Yes** (as raw counts + size factors) | **Yes** (as raw counts + norm factors) | **No** (breaks count distribution) |
| **Primary Interpretation** | Relative expression between conditions | Relative expression between conditions | Relative abundance within sample |

---

## 7. Other Normalization Approaches

### Counts Per Million (CPM)
CPM divides each gene's count by total library size and scales to one million:

$$
\text{CPM}_{ij} = \left( \frac{K_{ij}}{N_j} \right) \times 10^6
$$

CPM is useful for simple descriptive filtering (e.g., removing unexpressed genes with $\text{CPM} < 1$). However, CPM does not correct for composition bias or transcript length.

### RPKM and FPKM
- **RPKM:** Reads Per Kilobase of transcript per Million mapped reads (single-end).
- **FPKM:** Fragments Per Kilobase of transcript per Million mapped fragments (paired-end).

RPKM and FPKM divide by library size first and length second. Because the sum of RPKM/FPKM values varies across samples, TPM has largely superseded them for reporting relative abundances. Neither RPKM nor FPKM should be used as input for count-based differential expression tools.

### Single-Cell RNA-seq (scRNA-seq) Normalization
scRNA-seq datasets present unique challenges, such as high dropout rates (zeros), low coverage per cell, and cell-size differences.
- **Pseudobulk Approach:** Aggregate single-cell counts across cells of the same type per biological replicate, then use standard bulk normalization (DESeq2/edgeR).
- **Single-Cell Specific:** Methods like `scran` compute pooled size factors across cell clusters to handle extreme sparsity.

---

## 8. Practical Implementation

### DESeq2 Size Factors (R)

```r
library(DESeq2)

# Input: 'counts' (raw integer matrix) and 'metadata' (sample info dataframe)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Estimate size factors via median-of-ratios
dds <- estimateSizeFactors(dds)

# Inspect sample size factors
sizeFactors(dds)

# Retrieve normalized counts for visualization/clustering
norm_counts <- counts(dds, normalized = TRUE)

# Perform full differential expression workflow
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"))
```

### edgeR TMM Normalization (R)

```r
library(edgeR)

# Create DGEList object from raw counts
y <- DGEList(counts = counts, group = metadata$condition)

# Calculate TMM normalization factors
y <- calcNormFactors(y, method = "TMM")

# Inspect effective library sizes
y$samples

# Proceed with edgeR GLM fitting
design <- model.matrix(~ metadata$condition)
y      <- estimateDisp(y, design)
fit    <- glmQLFit(y, design)
res    <- glmQLFTest(fit, coef = 2)
```

### Manual TPM Calculation (R)

```r
# Calculates TPM from raw counts and feature lengths in base pairs
calculate_tpm <- function(counts_matrix, lengths_bp) {
  # Convert length to kilobases
  lengths_kb <- lengths_bp / 1000
  
  # Step 1: Compute RPK
  rpk <- counts_matrix / lengths_kb
  
  # Step 2 & 3: Scale by total RPK per sample
  tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
  return(tpm)
}
```

---

## 9. Common Mistakes to Avoid

1. **Feeding TPM or FPKM into DESeq2 / edgeR:** Differential expression models expect non-negative raw integer counts to properly parameterize negative binomial noise.
2. **Confusing Normalization with Transformation:** Normalization adjusts for scaling differences (e.g., size factors). Transformations like VST (Variance Stabilizing Transformation) or `rlog` are log-scale transformations intended for PCA, heatmaps, and clustering.
3. **Assuming Normalization Fixes Batch Effects:** Size factors adjust for global sequencing scale, not batch-specific technical noise. Batches must be included as covariates in GLM designs.
4. **Comparing TPM Values Directly Between Different Genes Across Conditions:** Compositional shifts can alter TPM proportions even if absolute expression remains constant.

---

## 10. Bioinformatics Problem

**Scenario:** You receive a raw count matrix from an experiment with six samples (3 Controls vs. 3 Treated). Total library sizes range from 12 million to 28 million reads. Treatment induces massive expression of a specific gene family.

### Recommended Execution Workflow

```text
1. Inspect Raw Counts & QC
   └── Check total library sizes & zero counts per sample.

2. Select Normalization Approach
   ├── For Differential Expression ──► Use DESeq2 (Median-of-Ratios) or edgeR (TMM)
   └── For Sample PCA/Heatmaps     ──► Apply VST or rlog on DESeqDataSet

3. Model Execution
   └── Fit Negative Binomial GLM using raw counts and estimated scaling factors.

4. Validation
   └── Check MA-plots and dispersion fits to confirm composition bias was resolved.
```

---

## 11. References

1. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
2. Robinson MD, McCarthy DJ, Smyth GK. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140. [https://doi.org/10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)
3. Robinson MD, Oshlack A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11, R25. [https://doi.org/10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
4. Wagner GP, Kin K, Lynch VJ. (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory in Biosciences*, 131, 281–285. [https://doi.org/10.1007/s12064-012-0162-3](https://doi.org/10.1007/s12064-012-0162-3)
5. Soneson C, Love MI, Robinson MD. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521. [https://doi.org/10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)