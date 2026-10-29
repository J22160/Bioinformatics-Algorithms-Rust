# Differential Expression: Negative Binomial Generalized Linear Model (GLM)

A comprehensive guide to understanding count data modeling, size-factor normalization, dispersion shrinkage, hypothesis testing, and implementation in DESeq2.

## Table of Contents

 1. [What Is This?](#1-what-is-this)
 2. [Why Not Just Use a t-test?](#2-why-not-just-use-a-t-test)
 3. [The Negative Binomial Distribution](#3-the-negative-binomial-distribution)
 4. [The DESeq2 Model](#4-the-deseq2-model)
 5. [How Size Factor Normalization Works](#5-how-size-factor-normalization-works)
 6. [How Dispersion Estimation Works](#6-how-dispersion-estimation-works)
 7. [Differential Expression Testing](#7-differential-expression-testing)
 8. [Worked Example: Treatment vs. Control](#8-worked-example-treatment-vs-control)
 9. [Practical Implementation with DESeq2](#9-practical-implementation-with-deseq2)
10. [Computational Complexity](#10-computational-complexity)
11. [Bioinformatics Problem](#11-bioinformatics-problem)
12. [Limitations and Best Practices](#12-limitations-and-best-practices)
13. [References](#13-references)

---

## 1. What Is This?

RNA sequencing (RNA-seq) measures gene expression by sequencing RNA-derived fragments and counting how many reads or fragments are assigned to each gene. After quantification, the result is typically represented as a **count matrix**, where rows represent genes, columns represent biological samples, and each entry contains the observed count for a gene in a sample.

For example:

| Gene | Control 1 | Control 2 | Control 3 | Treated 1 | Treated 2 | Treated 3 |
|---|---:|---:|---:|---:|---:|---:|
| **Gene A** | 120 | 135 | 110 | 250 | 275 | 260 |
| **Gene B** | 80 | 75 | 90 | 82 | 78 | 85 |
| **Gene C** | 12 | 18 | 10 | 45 | 50 | 42 |
| **Gene D** | 500 | 620 | 450 | 480 | 700 | 510 |

The goal of differential expression (DE) analysis is to identify genes whose expression differs systematically between experimental conditions, such as treatment and control.

However, observed counts are affected by several factors:

- **Biological condition:** Treatment may increase or decrease gene expression.
- **Sequencing depth:** Samples with more sequenced reads may have higher counts across many genes.
- **Biological variability:** Individuals or biological replicates may express genes differently even under the same condition.
- **Technical variability:** Library preparation, sequencing, and other experimental factors can introduce additional variation.

Differential expression analysis must distinguish changes associated with the condition of interest from these other sources of variation.

A widely used approach is the **negative binomial generalized linear model (NB GLM)**, implemented in tools such as DESeq2 and edgeR.

```
  +-------------------------------------------------------+
  |                   Raw Count Matrix                    |
  |             (Genes x Biological Samples)              |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  |              Size Factor Normalization                |
  |     (Adjusts for sequencing depth & composition)     |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  |          Dispersion Estimation & Shrinkage            |
  |       (Shares information across genes via EB)        |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  |             Negative Binomial GLM Fitting             |
  |               (Log link: ln(u) = Xb + offset)          |
  +-------------------------------------------------------+
                              |
                              v
  +-------------------------------------------------------+
  |          Hypothesis Testing & FDR Correction          |
  |         (Wald Test / LRT -> Benjamini-Hochberg)        |
  +-------------------------------------------------------+
```

---

## 2. Why Not Just Use a t-test?

A conventional t-test compares group means under assumptions about the distribution of the observations and the variability within groups. Applying a t-test directly to raw RNA-seq counts is generally inappropriate because count data have statistical properties that differ from the assumptions of a standard t-test.

### 1. RNA-seq counts are discrete and non-negative

Raw counts are integers greater than or equal to zero. They cannot be negative or fractional. A normal distribution is continuous and can assign probability to negative values, whereas a count distribution naturally represents non-negative integer observations.

### 2. Variance depends on the mean

RNA-seq counts exhibit a mean-variance relationship. Genes with higher expected counts often have greater absolute variability than genes with lower expected counts.

The Poisson distribution assumes that the variance equals the mean:

$$
\text{Var}(K) = \mu
$$

In real RNA-seq data, biological variability often produces more variation than the Poisson model allows. This phenomenon is called **overdispersion**. The negative binomial distribution accommodates this extra variability through a gene-specific dispersion parameter.

### 3. Samples have different sequencing depths

Suppose one sample contains twice as many usable reads as another. Even if the underlying expression proportions are identical, the first sample may have approximately twice as many reads assigned to a gene. A differential expression model must account for differences in sequencing depth and RNA composition before comparing gene expression.

### 4. Biological replicates are essential

A difference between two observed counts does not automatically indicate a biological effect. Replicates help estimate the variability expected within each condition. With few replicates, estimating variability independently for every gene is difficult. Methods such as DESeq2 therefore share information across genes to improve dispersion estimates.

> **Important:** A t-test is not inherently unsuitable for every analysis involving RNA-seq. Tests on appropriately transformed or normalized data can be useful in certain settings. However, applying a standard t-test directly to raw counts does not adequately model their mean-variance relationship, library size differences, and count-generating process.

---

## 3. The Negative Binomial Distribution

The negative binomial (NB) distribution is a probability distribution for non-negative integer counts. In RNA-seq, it models the observed count for a particular gene in a particular sample.

Let:
- $K_{ij}$: observed count for gene $i$ in sample $j$
- $\mu_{ij}$: expected count for gene $i$ in sample $j$
- $\alpha_i$: dispersion parameter for gene $i$

The model is:

$$
K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i)
$$

Its variance is:

$$
\text{Var}(K_{ij}) = \mu_{ij} + \alpha_i \mu_{ij}^{2}
$$

The two terms have distinct roles:
- $\mu_{ij}$ represents the expected count and determines the Poisson-like component of the variance.
- $\alpha_i \mu_{ij}^{2}$ represents additional variation beyond the Poisson assumption.

### Understanding Dispersion

Dispersion determines how much variability exists beyond what is expected from the mean alone.

- When $\alpha_i$ approaches zero, the negative binomial model approaches a Poisson model.
- When $\alpha_i$ is larger, the count distribution permits greater variability around its mean.
- Two genes with the same expected count can have different variances if their dispersion values differ.

For example, consider two genes with an expected count of 100:

| Parameter | Gene A | Gene B |
|---|---:|---:|
| **Mean count ($\mu$)** | 100 | 100 |
| **Dispersion ($\alpha$)** | 0.01 | 0.10 |
| **Variance ($\mu + \alpha\mu^2$)** | 200 | 1,100 |
| **Standard deviation** | 14.1 | 33.2 |

Both genes have the same expected count, but Gene B has substantially greater variability. A differential expression method should account for this difference when determining whether an observed change is statistically convincing.

---

## 4. The DESeq2 Model

DESeq2 models the observed count $K_{ij}$ for gene $i$ in sample $j$ using a negative binomial generalized linear model.

The model has two main components: a normalization offset and a linear predictor describing experimental effects.

### 1. Mean Model

$$
K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i)
$$

$$
\mu_{ij} = s_j q_{ij}
$$

Where:
- $s_j$: sample-specific size factor used to account for differences in sequencing depth and RNA composition
- $q_{ij}$: normalized expected expression for gene $i$ in sample $j$
- $\mu_{ij}$: expected observed count after accounting for the sample-specific size factor
- $\alpha_i$: dispersion parameter for gene $i$

The size factor acts as a multiplicative adjustment. A sample with a larger size factor is expected to have more observed reads for the same underlying normalized expression.

### 2. Linear Predictor

For a design matrix $X$, DESeq2 uses a log link:

$$
\ln(\mu_{ij}) = \ln(s_j) + X_j \beta_i
$$

Equivalently, for normalized expression:

$$
\ln(q_{ij}) = X_j \beta_i
$$

Here:
- $X_j$ is the row of the design matrix corresponding to sample $j$.
- $\beta_i$ is the vector of regression coefficients for gene $i$.
- $\ln(s_j)$ is an offset that accounts for the sample-specific size factor.

**Note on logarithm bases:** DESeq2's GLM is parameterized using the natural logarithm. Its reported $\log_2$ fold changes are obtained by expressing the estimated condition effect on a base-2 scale.

For a simple treatment-versus-control experiment, the design formula is:

```r
design = ~ condition
```

With control as the reference level, the condition coefficient represents the estimated treatment effect relative to control.

If the estimated $\log_2$ fold change coefficient is $\beta_{i,\text{treatment}}$, then:

$$
\text{log}_2\text{FC}_i = \frac{\beta_{i,\text{treatment}}}{\ln(2)}
$$

The corresponding fold change is:

$$
\text{FC}_i = 2^{\text{log}_2\text{FC}_i}
$$

For example:
- $\log_2 \text{fold change} = 0$: no estimated change.
- $\log_2 \text{fold change} = 1$: approximately twofold higher expression.
- $\log_2 \text{fold change} = -1$: approximately twofold lower expression.
- $\log_2 \text{fold change} = 2$: approximately fourfold higher expression.

In a simple two-group design, the condition coefficient estimates the log ratio of normalized expected counts. More complex designs can include additional coefficients for covariates, interactions, batches, or other experimental factors.

---

## 5. How Size Factor Normalization Works

Before comparing expression between samples, DESeq2 estimates a size factor for each sample. The default median-of-ratios method estimates relative sample scaling factors using genes across the dataset.

### Step 1: Calculate a Reference Expression for Each Gene

For each gene $i$, DESeq2 calculates a geometric mean across samples with usable counts:

$$
g_i = \left(\prod_{j=1}^{S} K_{ij}\right)^{1/S}
$$

Genes with a zero count in any sample have a zero geometric mean under this formula and do not contribute to the standard ratio calculation.

### Step 2: Calculate Sample-to-Reference Ratios

For each eligible gene and sample:

$$
r_{ij} = \frac{K_{ij}}{g_i}
$$

### Step 3: Estimate the Sample Size Factor

The size factor is the median of the gene-wise ratios for that sample:

$$
s_j = \text{median}_i(r_{ij})
$$

The factors are relative normalization values, conventionally scaled so that their geometric mean is one.

### Why Use the Median?

The median is robust to a subset of genes that are extremely highly expressed or strongly differentially expressed. It helps estimate the typical sample-wide scaling difference rather than allowing a few genes to dominate normalization.

However, this method relies on assumptions about the overall distribution of expression changes. If a large fraction of genes shifts in the same direction, or the dataset has unusual composition or many zeros, alternative normalization approaches or additional diagnostics may be needed.

> Normalization does not convert raw counts into directly comparable measurements of absolute RNA molecules per cell. It adjusts for relative sample scaling under the method's assumptions.

---

## 6. How Dispersion Estimation Works

Dispersion estimation is one of the key statistical features of DESeq2.

The model needs a dispersion estimate for each gene because biological variability differs among genes. However, estimating dispersion independently from only a few biological replicates can be unreliable.

DESeq2 addresses this problem through **information sharing across genes**, also called empirical Bayes shrinkage.

```
   Gene-wise Dispersions (Noisy)
             \
              \
               +--->  Mean-Dispersion Trend  --->  Shrunk Final Dispersions
              /                                      (Stable Estimates)
             /
   Prior Information (Gene Density)
```

### Step 1: Estimate Gene-wise Dispersion

DESeq2 estimates an initial dispersion for each gene using the observed counts and the fitted mean model. These gene-wise estimates capture gene-specific variability, but they can be noisy, particularly when sample sizes are small.

### Step 2: Fit a Mean-Dispersion Trend

DESeq2 models the general relationship between mean expression and dispersion across genes. The trend captures the typical dispersion expected for genes with similar mean expression.

### Step 3: Shrink Dispersion Estimates Toward the Trend

DESeq2 combines the gene-wise information with the trend to obtain more stable final dispersion estimates. Genes with limited information are generally pulled more strongly toward the trend. Genes with sufficient evidence for dispersion values different from the trend can retain more gene-specific estimates.

This improves stability and statistical power, especially in experiments with few biological replicates.

**Important distinction:** Dispersion shrinkage is not the same as $\log_2$ fold-change shrinkage. DESeq2 can also shrink fold-change estimates using separate methods, such as `apeglm`, to improve their stability when estimates are noisy.

---

## 7. Differential Expression Testing

After estimating normalization factors, dispersions, and regression coefficients, DESeq2 tests whether each gene has evidence of a condition-associated expression change.

### 1. Define the Null Hypothesis

For a treatment-versus-control comparison, the usual null hypothesis is:

$$
H_0: \beta_{i,\text{treatment}} = 0
$$

The alternative hypothesis is:

$$
H_A: \beta_{i,\text{treatment}} \neq 0
$$

The null hypothesis states that the condition effect is zero on the log scale, corresponding to a fold change of one.

### 2. Estimate the Condition Effect

The fitted GLM estimates the treatment coefficient for each gene. The estimated coefficient is converted to a $\log_2$ fold change for interpretation.

### 3. Calculate a Test Statistic

DESeq2 commonly uses a Wald test for the default coefficient-based comparison. Conceptually, the statistic is:

$$
W_i = \frac{\widehat{\beta}_i}{\text{SE}(\widehat{\beta}_i)}
$$

Here, $\widehat{\beta}_i$ is the estimated condition coefficient on the GLM's natural-log scale, and $\text{SE}(\widehat{\beta}_i)$ is its standard error. A large absolute Wald statistic indicates that the estimated effect is large relative to its uncertainty.

### 4. Correct for Multiple Testing

RNA-seq experiments may test thousands or tens of thousands of genes simultaneously. If 20,000 independent null hypotheses were tested at an unadjusted p-value threshold of 0.05, approximately 1,000 false positives would be expected on average under the global null.

To control the false discovery rate (FDR), DESeq2 uses the Benjamini-Hochberg procedure by default. The adjusted p-value is commonly called the **padj** or FDR-adjusted p-value. A frequently used significance threshold is:

```text
padj < 0.05
```

This controls the expected proportion of false discoveries among rejected hypotheses under the procedure's assumptions. It does not mean that every individual significant gene has exactly a 5% probability of being false.

DESeq2 also uses independent filtering by default to improve detection power. Genes with very low mean normalized counts may be filtered from the multiple-testing procedure because they generally provide little power to detect meaningful changes.

### 5. Interpret Effect Size and Significance Together

For complete interpretation, evaluate:
- **$\log_2$ fold change:** Magnitude and direction of the estimated change.
- **Adjusted p-value:** Evidence after multiple-testing correction.
- **Mean normalized count:** Abundance and context for interpreting the estimate.
- **Standard error or confidence interval:** Uncertainty in the effect estimate.
- **Biological relevance:** Whether the magnitude and direction are meaningful for the system being studied.

---

## 8. Worked Example: Treatment vs. Control

Suppose an experiment contains three biological replicates per condition and measures 20,000 genes.

Consider three hypothetical genes:

| Gene | Mean normalized count (control) | Mean normalized count (treated) | $\log_2$ fold change | Adjusted p-value |
|---|---:|---:|---:|---:|
| **Gene A** | 100 | 200 | 1.0 | 0.001 |
| **Gene B** | 100 | 110 | 0.14 | 0.72 |
| **Gene C** | 10 | 40 | 2.0 | 0.08 |

### Interpretation

- **Gene A:** Estimated expression doubles in treatment, with strong evidence after multiple-testing correction. It meets the conventional `padj < 0.05` threshold.
- **Gene B:** The estimated increase is small, and the evidence is insufficient to declare differential expression.
- **Gene C:** The estimated expression increases fourfold, but the adjusted p-value does not meet the conventional significance threshold. Low counts and/or high variability contribute to the uncertainty.

---

## 9. Practical Implementation with DESeq2

```r
library(DESeq2)

# Example setup:
# counts   <- read.csv("gene_counts.csv", row.names = 1, check.names = FALSE)
# metadata <- read.csv("metadata.csv", row.names = 1)

metadata$condition <- factor(
  metadata$condition,
  levels = c("control", "treated")
)

# Confirm sample alignment between count matrix and metadata
stopifnot(
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(rownames(metadata)),
  setequal(colnames(counts), rownames(metadata))
)

metadata <- metadata[colnames(counts), , drop = FALSE]

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Filter low count genes to reduce memory and testing burden
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run the DESeq2 pipeline (Normalization, Dispersion Estimation, GLM Fitting)
dds <- DESeq(dds)

# Extract results for treated vs. control
res <- results(
  dds,
  contrast = c("condition", "treated", "control"),
  alpha    = 0.05
)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj, na.last = TRUE), ]

# Filter for statistically significant genes (FDR < 0.05)
significant_genes <- subset(
  as.data.frame(res_ordered),
  !is.na(padj) & padj < 0.05
)

# Save results
write.csv(
  as.data.frame(res_ordered),
  "deseq2_results.csv",
  row.names = TRUE
)
```

### Important Implementation Notes

1. **Use raw counts:** DESeq2 expects unnormalized integer counts. Do not supply TPM, FPKM, or normalized values.
2. **Biological replicates:** Technical replicates do not substitute for biological replicates and should be collapsed using `collapseReplicates()`.
3. **Complex designs:** Incorporate terms for confounding factors or batches when applicable (e.g., `design = ~ batch + condition`).
4. **Outlier detection:** Review Cooks distance and diagnostic plots (PCA, MA-plot, dispersion plots) to detect sample or gene outliers.

---

## 10. Computational Complexity

Let:
- $G$: number of genes tested
- $S$: number of samples
- $P$: number of coefficients in the design matrix
- $I$: number of iterations used in fitting

| Operation | Approximate Computational Cost | Explanation |
|---|---|---|
| **Size Factor Estimation** | $\mathcal{O}(G \cdot S)$ | Calculates gene geometric means and sample medians. |
| **Initial GLM Fit** | $\mathcal{O}(G \cdot S \cdot P)$ | Fits iterative regression coefficients per gene. |
| **Dispersion Estimation** | $\mathcal{O}(G \cdot S \cdot I)$ | Estimates gene-wise dispersion and fits mean-dispersion trend. |
| **Wald Testing** | $\mathcal{O}(G \cdot P)$ | Evaluates test statistics and raw p-values per gene. |
| **FDR Adjustment** | $\mathcal{O}(G \log G)$ | Sorts p-values for Benjamini-Hochberg adjustment. |
| **Storage Requirement** | $\mathcal{O}(G \cdot S)$ | Stores raw and normalized count matrices alongside model states. |

---

## 11. Bioinformatics Problem

**Scenario:** You have RNA-seq counts for 20,000 genes across six biological samples: three controls and three treated samples. Identify genes that are significantly differentially expressed between conditions while controlling for false discoveries.

### Inputs
- A raw gene-level count matrix ($20,000 \times 6$).
- Sample metadata specifying conditions (`control` vs `treated`).
- Experimental design formulation.

### Execution Workflow

1. Perform sample-level Quality Control (QC) and PCA inspect.
2. Filter unexpressed/lowly expressed genes.
3. Compute size factors via median-of-ratios method.
4. Fit Negative Binomial GLM for each gene.
5. Apply empirical Bayes dispersion shrinkage toward fitted trend.
6. Perform Wald test for condition parameter $\beta_{\text{treatment}}$.
7. Adjust p-values using Benjamini-Hochberg FDR control procedure.
8. Filter for $padj < 0.05$ and analyze effect sizes ($\log_2 \text{FC}$).

---

## 12. Limitations and Best Practices

- **Sample Size Constraints:** With fewer than 3 biological replicates per group, statistical power is severely restricted, and dispersion estimation becomes unreliable.
- **Low-Expression Noise:** Genes with minimal counts add multiple-testing burden without providing statistical power.
- **Confounding Variables:** Unaddressed batch effects can create false positives or mask true biological effects.
- **Single-Cell RNA-seq:** Standard DESeq2 bulk models are generally inappropriate for sparse single-cell matrices due to high zero-inflation. Pseudobulk aggregation is recommended before applying DESeq2 to scRNA-seq datasets.

---

## 13. References

1. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
2. Anders S, Huber W. (2010). Differential expression analysis for sequence count data. *Genome Biology*, 11, R106. [https://doi.org/10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106)
3. Robinson MD, McCarthy DJ, Smyth GK. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140. [https://doi.org/10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)
4. Benjamini Y, Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.
5. DESeq2 Bioconductor Documentation: [https://bioconductor.org/packages/DESeq2/](https://bioconductor.org/packages/DESeq2/)