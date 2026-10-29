# Gene Set Enrichment Analysis (GSEA)

A comprehensive guide to understanding pathway enrichment, the weighted running-sum statistic, leading-edge analysis, and implementation with GSEApy.

## Table of Contents

 1. [What Is This?](#1-what-is-this)

 2. [Why Is GSEA Useful?](#2-why-is-gsea-useful)

 3. [GSEA vs. Over-Representation Analysis (ORA)](#3-gsea-vs-over-representation-analysis-ora)

 4. [How GSEA Works](#4-how-gsea-works)

 5. [The Running-Sum Statistic in Detail](#5-the-running-sum-statistic-in-detail)

 6. [Worked Example](#6-worked-example)

 7. [Leading-Edge Analysis](#7-leading-edge-analysis)

 8. [Interpreting GSEA Results](#8-interpreting-gsea-results)

 9. [Practical Implementation with GSEApy](#9-practical-implementation-with-gseapy)

10. [Computational Complexity](#10-computational-complexity)

11. [Limitations and Best Practices](#11-limitations-and-best-practices)

12. [Bioinformatics Problem](#12-bioinformatics-problem)

13. [References](#13-references)

---

## 1. What Is This?

**Gene Set Enrichment Analysis (GSEA)** is a computational method used to determine whether groups of genes associated with a biological process, pathway, or molecular function show coordinated changes in an experiment.

After performing differential expression analysis, you may identify thousands of genes with different expression levels between treatment and control. Interpreting these genes individually can be difficult because biological processes typically involve many genes working together.

GSEA addresses this problem by analyzing **predefined gene sets** rather than evaluating each gene in isolation.

Examples of gene sets include:

* **Inflammatory response:** Genes involved in inflammation and immune signaling.
* **Oxidative phosphorylation:** Genes involved in mitochondrial ATP production.
* **Cell cycle:** Genes regulating cell division and proliferation.
* **Apoptosis:** Genes associated with programmed cell death.
* **DNA repair:** Genes involved in detecting and repairing DNA damage.

The central question is:

> **Are genes belonging to a particular biological process concentrated near the top or bottom of a ranked gene list more often than expected by chance?**

```
Ranked Gene List (Top = Up-regulated, Bottom = Down-regulated)
[ Gene 1 (+5.2) ] ───► Hit (+ P_hit)
[ Gene 2 (+4.6) ] ───► Miss (- P_miss)
[ Gene 3 (+3.9) ] ───► Hit (+ P_hit)  ───►  [ Calculate Running Sum ]
...                                                 │
[ Gene N (-4.4) ]                                   v
                                           Peak Deviation = Enrichment Score (ES)
```

For example, after a treatment, many inflammatory response genes might show modest increases in expression. Although individual genes may not pass a differential expression significance threshold, their coordinated behavior could indicate activation of the inflammatory response pathway.

GSEA can detect this type of coordinated signal without requiring every gene in the pathway to be individually significant.

---

## 2. Why Is GSEA Useful?

### 1. Biological Pathways Contain Multiple Genes
Biological processes rarely depend on a single gene. A pathway may involve dozens or hundreds of genes, each contributing a small part to the overall process. Examining the collective behavior of these genes can provide a more interpretable biological picture than examining a long list of individual differential expression results.

### 2. It Does Not Require a Significant-Gene Cutoff
A common alternative is over-representation analysis (ORA), which first selects significant differentially expressed genes and then tests whether particular pathways contain more selected genes than expected.

ORA depends on a threshold such as:

```text
adjusted_p_value < 0.05
```

Genes just above and below the threshold may have very similar evidence of differential expression, yet ORA treats them differently. GSEA instead uses a ranked list of genes, typically including all genes that were eligible for testing. It therefore preserves more information about the direction and strength of the expression changes.

### 3. It Can Detect Coordinated but Modest Changes
Suppose a pathway contains 100 genes, many of which show small positive changes in expression. Even if none passes the individual gene-level significance threshold, the pathway genes may be concentrated near the top of the ranked list. GSEA evaluates whether this concentration is stronger than expected under its enrichment null model.

### 4. It Provides Biological Interpretation
GSEA results can help connect gene-level statistical results to higher-level biological processes, pathways, and molecular functions.

> **Note:** Enrichment does not by itself prove that a pathway is activated, inhibited, or causally responsible for a phenotype. Interpretation depends on the ranking metric, direction of enrichment, gene-set definition, and experimental context.

---

## 3. GSEA vs. Over-Representation Analysis (ORA)

Both methods identify biological categories associated with experimental results, but they use different inputs and statistical approaches.

| Feature | GSEA | Over-Representation Analysis (ORA) |
|---|---|---|
| **Input** | Ranked gene list (all genes) | Selected gene list (e.g., DE cutoff) |
| **Requires Significant Genes First?** | No | Usually yes |
| **Uses Ranking and Direction?** | Yes, in standard ranked GSEA | Usually not in conventional ORA |
| **Can Detect Modest Coordinated Changes?** | Yes | May miss them if genes fail selection threshold |
| **Main Statistical Question** | Are gene-set members clustered at extremes? | Does selected list contain more members than expected? |
| **Common Statistical Framework** | Running-sum enrichment score & null distribution | Hypergeometric test or Fisher's exact test |
| **Background Set** | All genes eligible for ranking | All genes evaluated in the experiment |

Neither approach is universally superior. GSEA is particularly useful when the full ranked gene list is available and coordinated changes are of interest. ORA can be effective when a well-defined set of significant genes is the intended input.

---

## 4. How GSEA Works

A standard ranked GSEA analysis consists of several main steps.

### Step 1: Rank All Eligible Genes
First, perform a gene-level analysis comparing the experimental conditions using DESeq2, edgeR, or another suitable method. Each gene receives a ranking statistic representing its association with the comparison.

Common choices include:
* Signed Wald statistic from DESeq2.
* Signed test statistic from edgeR.
* Signal-to-noise ratio in appropriate expression-analysis settings.

For a treatment-versus-control comparison, a signed statistic places genes associated with increased treatment expression near the top and genes with decreased expression near the bottom:

| Rank | Gene | Ranking Statistic | Interpretation |
|---:|---|---:|---|
| 1 | Gene A | +5.2 | Strong positive treatment association |
| 2 | Gene B | +4.6 | Positive treatment association |
| 3 | Gene C | +3.9 | Positive treatment association |
| 4 | Gene D | +1.4 | Modest positive association |
| 5 | Gene E | +0.2 | Little evidence of a directional effect |
| ... | ... | ... | ... |
| $G-1$ | Gene Y | -4.4 | Strong negative treatment association |

> **Important:** Ranking by adjusted p-value alone is usually undesirable because it loses direction and creates ties. Use a signed test statistic instead. Include all tested genes, not just those passing a significance cutoff.

### Step 2: Define the Gene Sets
Gene sets are collections of genes associated with biological concepts, sourced from databases such as:
* **MSigDB:** Curated collections, including Hallmark, C2 (curated pathways), C5 (Gene Ontology), etc.
* **Gene Ontology (GO):** Biological processes, molecular functions, and cellular components.
* **KEGG / Reactome:** Curated biological pathway maps.

### Step 3: Calculate the Running Enrichment Score
For each gene set, GSEA walks down the ranked list and calculates a running sum:
* When encountering a gene in the set, the running score **increases** (weighted by its ranking statistic).
* When encountering a gene outside the set, the running score **decreases**.

### Step 4: Identify the Enrichment Score (ES)
The **Enrichment Score (ES)** is the running-sum deviation with the greatest absolute magnitude from zero:
* **Positive ES:** Indicates enrichment toward the top of the list (up-regulation).
* **Negative ES:** Indicates enrichment toward the bottom of the list (down-regulation).

### Step 5: Assess Statistical Significance
GSEA compares the observed ES against a null distribution generated via:
* **Phenotype permutation:** Shuffles sample labels and recomputes statistics (preferred when sample size permits, $N \ge 7$ per group).
* **Gene-set permutation:** Randomizes gene membership across the ranked list (used in preranked GSEA).

### Step 6: Normalize Enrichment Scores (NES)
To account for differences in gene-set sizes and null distribution variations, GSEA computes a **Normalized Enrichment Score (NES)**:

$$
\text{NES} = \frac{\text{Observed ES}}{\text{Mean(ES under Null for same direction)}}
$$

### Step 7: Correct for Multiple Testing
GSEA adjusts $p$-values across tested gene sets using the False Discovery Rate (FDR) $q$-value procedure.

---

## 5. The Running-Sum Statistic in Detail

Let:
* $L = (g_1, g_2, \ldots, g_G)$ be the ranked list of $G$ genes.
* $S$ be a gene set containing $N_H$ genes present in $L$.
* $r_i$ be the ranking statistic for gene $g_i$.
* $p$ be the weighting exponent (commonly $p = 1$).

### Hit Contribution
For a gene in the set ($g_i \in S$):

$$
P_{\text{hit}}(i) = \frac{|r_i|^p}{\sum_{g_j \in S} |r_j|^p}
$$

### Miss Penalty
For a gene outside the set ($g_i \notin S$):

$$
P_{\text{miss}}(i) = \frac{1}{G - N_H}
$$

### Running Sum
The running sum at position $i$ is:

$$
RS(i) = \sum_{\substack{j \le i \\ g_j \in S}} P_{\text{hit}}(j) - \sum_{\substack{j \le i \\ g_j \notin S}} P_{\text{miss}}(j)
$$

The **Enrichment Score (ES)** is the maximum absolute value of $RS(i)$ across all positions $i \in \{1, \ldots, G\}$.

---

## 6. Worked Example

Suppose a ranked list contains 10 genes, and the pathway of interest contains genes **A, C, D, and H**.

| Rank | Gene | In Pathway? | Illustrative Running Score |
|---:|---|---|---:|
| 1 | A | Yes | +0.40 |
| 2 | B | No | +0.23 |
| 3 | C | Yes | +0.53 |
| 4 | D | Yes | **+0.83** (Peak = ES) |
| 5 | E | No | +0.67 |
| 6 | F | No | +0.50 |
| 7 | G | No | +0.33 |
| 8 | H | Yes | +0.63 |
| 9 | I | No | +0.47 |
| 10 | J | No | 0.00 |

### Observations from Example:
1. The score increases at pathway hits (A, C, D, H) and decreases at misses.
2. The maximum peak occurs at **Rank 4 (Gene D)** with an $\text{ES} = +0.83$.
3. Because pathway genes cluster near the top, the pathway is positively enriched.

---

## 7. Leading-Edge Analysis

The **leading-edge subset** consists of the genes in the gene set encountered **before** the peak running sum score is reached.

```
Running Score
    ^        (Peak ES = +0.83)
0.83|            / \
    |   Hit     /   \
0.53|   /  \___/     \
0.40|  /              \_______  0.00
0.00+------------------------------> Gene Rank
    |  A  B  C  D  E  F  G  H  I  J
       [Leading-Edge]
```

In the worked example above, genes **A, C, and D** constitute the leading-edge subset.

### Applications of Leading-Edge Genes:
* Identify core driving genes responsible for pathway enrichment.
* Extract overlapping genes across multiple enriched pathways.
* Filter targets for experimental validation and functional follow-up.

---

## 8. Interpreting GSEA Results

Consider the following illustrative output from a treatment-vs-control comparison:

| Pathway | ES | NES | Nominal $p$-value | FDR $q$-value | Interpretation |
|---|---:|---:|---:|---:|---|
| **Inflammatory Response** | +0.72 | +2.10 | <0.001 | 0.010 | Strong up-regulation in treatment |
| **Oxidative Phosphorylation**| -0.65 | -1.80 | 0.002 | 0.030 | Down-regulated in treatment |
| **Cell Cycle** | +0.41 | +1.20 | 0.150 | 0.320 | Weak evidence; not significant |
| **DNA Repair** | -0.48 | -1.40 | 0.040 | 0.180 | Borderline significance |

### Key Parameters:
* **Positive NES:** Indicates general up-regulation of the pathway in the treatment group.
* **Negative NES:** Indicates general down-regulation in the treatment group.
* **FDR $q$-value:** Adjusted significance accounting for all tested gene sets. (Standard exploratory threshold is $q < 0.25$, though $q < 0.05$ is preferred for strict reporting).

---

## 9. Practical Implementation with GSEApy

### Installation

```bash
pip install gseapy
# Or via conda
conda install -c bioconda gseapy
```

### Python Script: Preranked GSEA

```python
import pandas as pd
import gseapy as gp

# 1. Load differential expression results (e.g., from DESeq2)
df = pd.read_csv("deseq2_results.csv")

# 2. Extract gene symbol and signed Wald statistic ('stat')
ranked_df = df[["gene_symbol", "stat"]].dropna().drop_duplicates(subset="gene_symbol")

# 3. Sort genes in descending order by statistic
ranked_df = ranked_df.sort_values(by="stat", ascending=False)

# 4. Run Preranked GSEA using MSigDB Hallmark collection
pre_res = gp.prerank(
    rnk=ranked_df,
    gene_sets="MSigDB_Hallmark_2020",
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42,
    outdir="gsea_output"
)

# 5. Access results dataframe
results = pre_res.res2d
print(results[["Term", "ES", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]].head())
```

---

## 10. Computational Complexity

Let:
* $G$: Number of genes in the ranked list.
* $P$: Number of gene sets tested.
* $R$: Number of permutations.

| Operation | Approximate Time Complexity | Description |
|---|---|---|
| **Gene List Sorting** | $O(G \log G)$ | Sort all genes by ranking statistic |
| **Single ES Calculation** | $O(G)$ | Calculate running sum for one gene set |
| **All Gene Sets (1 Pass)**| $O(P \cdot G)$ | Single-pass evaluation over $P$ sets |
| **Permutation Testing** | $O(P \cdot G \cdot R)$ | Compute ES across $R$ null iterations |

For a standard human dataset ($G \approx 20,000$, $P \approx 500$, $R = 1,000$), optimized libraries (such as GSEApy's Rust backend or `fastGSEA`) run in a few seconds to minutes.

---

## 11. Limitations and Best Practices

1. **Identifier Matching:** Ensure gene symbols or Ensembl IDs match the gene-set definition file (`.gmt`) exactly.
2. **Ranking Metric Selection:** Prefer signed statistics (e.g., Wald stat) over $p$-values or raw fold changes.
3. **Gene Set Overlap:** Redundant pathways often share leading-edge genes; group related pathways before drawing conclusions.
4. **Gene Set Size Bounds:** Filter out tiny sets ($<15$ genes) to avoid high variance, and huge sets ($>500$ genes) to maintain biological specificity.
5. **Correlation Assumptions:** Gene-set permutation in preranked GSEA assumes independent genes, which can inflate false positives if co-expressed genes cluster together. Use phenotype permutation when sample size allows ($N \ge 7$).

---

## 12. Bioinformatics Problem

**Scenario:** You have performed RNA-seq on 6 samples (3 Control vs. 3 Treated). Individual gene expression analysis yields 1,500 DE genes, but single-gene analysis does not clearly reveal the biological mechanism.

### Workflow Solution:
1. Obtain signed Wald statistics for all ~20,000 genes from DESeq2.
2. Construct a `.rnk` file sorted from highest positive to lowest negative statistic.
3. Execute GSEA against curated pathway databases (e.g., Reactome or KEGG).
4. Identify top enriched pathways with $\text{FDR } q < 0.05$.
5. Perform leading-edge analysis on top pathways to identify core driving genes for qPCR or functional knockdown.

---

## 13. References

1. Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545–15550. <https://doi.org/10.1073/pnas.0506580102>

2. Mootha VK, Lindgren CM, Eriksson KF, et al. (2003). PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nature Genetics*, 34, 267–273. <https://doi.org/10.1038/ng1180>

3. Fang Z, Liu X, Peltz G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*, 39(1), btac757. <https://doi.org/10.1093/bioinformatics/btac757>

4. Khatri P, Sirota M, Butte AJ. (2012). Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology*, 8(2), e1002375. <https://doi.org/10.1371/journal.pcbi.1002375>