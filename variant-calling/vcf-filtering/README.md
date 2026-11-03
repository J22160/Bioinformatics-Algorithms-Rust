# VCF Quality Filtering Guide

## What Is VCF Quality Filtering?

Variant calling identifies candidate genetic variants from sequencing data, but not every detected variant is reliable. Technical issues such as sequencing errors, poor read mapping, low coverage, strand bias, PCR artifacts, or regional paralogy can lead to false positives.

**VCF quality filtering** evaluates raw variant calls against confidence metrics and quality parameters. The primary goal is to minimize false-positive variant calls while retaining true biological variation.

Variant filtering is critical across diverse genomic domains:

* **Whole-Genome Sequencing (WGS)** & **Whole-Exome Sequencing (WES)**
* **Rare Variant Association Studies** & **Genome-Wide Association Studies (GWAS)**
* **Clinical & Diagnostic Genomics**
* **Population Genetics**

Filtering decisions heavily dictate the accuracy of downstream analyses. Retaining false-positive calls can skew allele-frequency distributions, inflate rare variant counts, and severely confound downstream statistical tests.

---

### Core Filtering Strategies

```
┌────────────────────────────────────────────────────────────────────────┐
│                        VCF Quality Filtering                           │
└────────────────────────────────────────────────────────────────────────┘
                                   │
         ┌─────────────────────────┴─────────────────────────┐
         ▼                                                   ▼
┌───────────────────────────────┐           ┌───────────────────────────────┐
│        Hard Filtering         │           │             VQSR              │
├───────────────────────────────┤           ├───────────────────────────────┤
│ • Predefined threshold cutoffs│           │ • Gaussian Mixture Model      │
│ • Rules evaluated independently│          │ • Multidimensional features   │
│ • Ideal for small datasets    │           │ • Learns true/false patterns  │
└───────────────────────────────┘           └───────────────────────────────┘
```

1. **Hard Filtering:** Applies predefined fixed thresholds to individual variant annotations.
2. **Variant Quality Score Recalibration (VQSR):** Employs a machine-learning model (Gaussian Mixture Model) to evaluate variants multidimensionally against high-confidence truth sets.

These approaches are not mutually exclusive; workflows select specific methodologies based on sample size, reference availability, and target applications.

---

## 1. Understanding VCF Quality Annotations

A standard Variant Call Format (VCF) file organizes data into key structural fields: `QUAL`, `FILTER`, `INFO`, and `FORMAT`.

### 1.1 `QUAL`: Variant Quality Score

The `QUAL` field represents the Phred-scaled probability that a variant exists at a specific genomic site, calculated under the caller's statistical model:

$$
\text{QUAL} = -10 \log_{10} \left( P(\text{variant is false}) \right)
$$

| `QUAL` Score | Error Probability ($P$) | Accuracy Confidence |
| :---: | :---: | :---: |
| **10** | $10\%$ | $90.0\%$ |
| **20** | $1\%$ | $99.0\%$ |
| **30** | $0.1\%$ | $99.9\%$ |
| **40** | $0.01\%$ | $99.99\%$ |
| **50** | $0.001\%$ | $99.999\%$ |
| **60** | $0.0001\%$ | $99.9999\%$ |

> **Note:** `QUAL` measures site-level quality rather than sample-level Genotype Quality (`GQ`). It should be interpreted according to the specific caller's calibration model.

---

### 1.2 `DP`: Read Depth

`DP` measures total sequencing coverage (number of reads) supporting observations at a position.

* **Low Depth ($DP \le 3$):** High sampling error and limited statistical confidence.
* **Moderate/High Depth ($DP \ge 20$):** Strong evidence for base calling.
* **Excessive Depth:** Potential indicator of repetitive regions, copy-number expansions, or paralogous alignment artifacts.

> **Important:** `INFO/DP` (total site-level depth across samples) must not be confused with `FORMAT/DP` (sample-specific read depth).

---

### 1.3 `QD`: Quality by Depth

`QD` normalizes the variant quality score relative to read coverage to prevent raw depth from inflating quality estimations:

$$
\text{QD} \approx \frac{\text{QUAL}}{\text{DP}_{\text{unfiltered alternate reads}}}
$$

Normalizing quality by depth resolves artifacts where low-confidence variants appear artificially reliable simply due to extremely high depth.

---

### 1.4 Strand Bias Metrics (`FS` & `SOR`)

Strand bias occurs when an alternate allele is disproportionately observed on either the forward or reverse sequencing strand, typically signaling a technical artifact.

```
Reference Reads:  Forward [══════════]  Reverse [══════════]  (Balanced)
Alternate Reads:  Forward [══════════]  Reverse [          ]  (Strand Biased Artifact)
```

#### `FS` (Fisher Strand Bias)
Calculates strand bias using Fisher's Exact Test. Reported as a Phred-scaled $p$-value:
* **Lower `FS` values:** Preferable (minimal bias).
* **High `FS` values ($>60$):** Strong evidence of technical strand bias.

#### `SOR` (Strand Odds Ratio)
Uses an odds-ratio formulation to evaluate strand bias. `SOR` avoids over-penalizing variants located near read ends or in high-coverage regions.
* **Lower `SOR` values ($<3.0$):** Preferable.
* **Elevated `SOR` values ($>3.0$ or $>4.0$):** Indicates potential strand-bias artifacts.

---

### 1.5 Mapping Quality Metrics (`MQ`, `MQRankSum`, `ReadPosRankSum`)

#### `MQ` (Root Mean Square Mapping Quality)
Reflects the average mapping quality of reads covering a site:

$$
\text{MQ} = -10 \log_{10} \left( P(\text{read is mismapped}) \right)
$$

* **Target Value:** $\text{MQ} \approx 60$ (Illumina standard). Low values ($\text{MQ} < 40$) signify mapping ambiguity (e.g., repetitive elements).

#### `MQRankSum` (Mapping Quality Rank Sum Test)
Uses a Wilcoxon rank-sum test to compare the mapping qualities of reads supporting the reference allele versus those supporting the alternate allele.
* **Near $0$:** Reference and alternate reads map with equal confidence.
* **Negative values ($<-12.5$):** Alternate reads have significantly lower mapping quality than reference reads.

#### `ReadPosRankSum` (Read Position Rank Sum Test)
Compares position distributions of reference vs. alternate alleles within sequencing reads.
* **Negative values ($<-8.0$):** Alternate alleles accumulate near read ends, where sequencing error rates naturally peak.

---

### Summary of Primary VCF Quality Annotations

| Metric | Level | Metric Description | Preferred Direction / Range |
| :--- | :--- | :--- | :--- |
| **`QUAL`** | Site | Phred-scaled variant confidence | Higher is better ($\ge 30$) |
| **`DP`** | Site / Sample | Read depth | Moderate depth (platform dependent) |
| **`QD`** | Site | Quality normalized by depth | Higher is better ($\ge 2.0$) |
| **`FS`** | Site | Fisher strand-bias statistic | Lower is better ($< 60.0$) |
| **`SOR`** | Site | Strand odds ratio | Lower is better ($< 3.0$) |
| **`MQ`** | Site | RMS mapping quality | Higher is better ($\ge 40.0$) |
| **`MQRankSum`** | Site | Ref vs Alt mapping quality rank sum | Close to zero ($> -12.5$) |
| **`ReadPosRankSum`** | Site | Ref vs Alt read position rank sum | Close to zero ($> -8.0$) |
| **`GQ`** | Sample | Genotype quality confidence | Higher is better ($\ge 20$) |
| **`AD`** | Sample | Allelic depth per allele | High alt depth relative to total depth |
| **`VAF`** | Sample | Variant Allele Fraction ($\frac{\text{Alt}}{\text{Total}}$) | Heterozygous $\approx 0.5$, Homozygous $\approx 1.0$ |

---

## 2. Hard Filtering

Hard filtering applies static, predefined conditional logic directly to individual VCF metrics. Sites failing any condition are tagged in the `FILTER` column.

### Standard Germline Threshold Example (SNPs)

```text
QUAL > 30.0
QD > 2.0
FS < 60.0
MQ > 40.0
MQRankSum > -12.5
ReadPosRankSum > -8.0
```

### Filtering Decision Matrix

| Variant | `QUAL` | `QD` | `FS` | `MQ` | Filter Status | Assigned Reason |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| `chr1:1000 A>G` | 65 | 8.0 | 2.0 | 60 | **PASS** | Passed all filters |
| `chr1:2000 C>T` | 18 | 1.2 | 5.0 | 55 | **FAIL** | Low `QUAL` and low `QD` |
| `chr1:3000 G>A` | 80 | 6.0 | 85.0 | 58 | **FAIL** | High `FS` (Strand Bias) |
| `chr1:4000 T>C` | 70 | 5.0 | 3.0 | 25 | **FAIL** | Low `MQ` (Mapping Quality) |

---

### Executing Hard Filtering with GATK & BCFtools

#### Step 1: Flag Variants with GATK `VariantFiltration`

```bash
gatk VariantFiltration \
    -R reference.fasta \
    -V input.vcf.gz \
    --filter-name "LowQD" --filter-expression "QD < 2.0" \
    --filter-name "HighFS" --filter-expression "FS > 60.0" \
    --filter-name "LowMQ" --filter-expression "MQ < 40.0" \
    --filter-name "LowMQRankSum" --filter-expression "MQRankSum < -12.5" \
    --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
    -O flagged.vcf.gz
```

#### Step 2: Extract Retained (`PASS`) Variants with BCFtools

```bash
# Retain only variants marked explicitly as PASS
bcftools view \
    -f PASS \
    flagged.vcf.gz \
    -Oz -o pass_only.vcf.gz

bcftools index -t pass_only.vcf.gz
```

> **Note:** A null filter entry (`.`) indicates an un-evaluated site, which is distinct from an explicit `PASS` call.

---

## 3. Variant Quality Score Recalibration (VQSR)

Hard filtering evaluates annotations independently using rigid boundaries. **VQSR** models annotation distributions holistically using a supervised machine-learning model (**Gaussian Mixture Model**).

```
                      ┌───────────────────────────────────────┐
                      │    Input Callset Annotations          │
                      │ (QD, FS, SOR, MQ, MQRankSum, etc.)    │
                      └───────────────────────────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │ High-Confidence Calibration Sets      │
                      │ (HapMap, Omni, 1000G, dbSNP)          │
                      └───────────────────────────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │   Gaussian Mixture Model (GMM)        │
                      │ Fits multidimensional distributions   │
                      └───────────────────────────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │ Calculate Log Odds Ratio (VQSLOD)     │
                      └───────────────────────────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │ Slice Variants into Truth Sensitivity │
                      │ Tranches (e.g., 99.0%, 99.5%)         │
                      └───────────────────────────────────────┘
```

### Gaussian Mixture Model (GMM) Formulation

VQSR models the multivariate annotation space as a combination of $K$ Gaussian distributions:

$$
P(\mathbf{x}) = \sum_{k=1}^{K} \pi_k \cdot \mathcal{N}\left(\mathbf{x} \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k\right)
$$

Where:
* $\mathbf{x}$ = Vector of combined variant annotations ($\text{QD}, \text{FS}, \text{SOR}, \text{MQ}, \text{MQRankSum}, \text{ReadPosRankSum}$).
* $K$ = Number of Gaussian components.
* $\pi_k$ = Mixture weight assigned to component $k$, where $\sum \pi_k = 1$.
* $\boldsymbol{\mu}_k$ = Mean vector of component $k$.
* $\boldsymbol{\Sigma}_k$ = Covariance matrix modeling inter-annotation dependencies.
* $\mathcal{N}$ = Multivariate normal probability density function.

The model computes a log-odds score (**VQSLOD**) evaluating whether a variant's annotation pattern aligns with true variants versus background noise.

---

### Executing VQSR with GATK (SNP Mode Example)

#### Step 1: Model Training (`VariantRecalibrator`)

```bash
gatk VariantRecalibrator \
    -R reference.fasta \
    -V input.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -O snps.recal \
    --tranches-file snps.tranches
```

#### Step 2: Apply Recalibration Cutoffs (`ApplyVQSR`)

```bash
gatk ApplyVQSR \
    -R reference.fasta \
    -V input.vcf.gz \
    --recal-file snps.recal \
    --tranches-file snps.tranches \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP \
    -O snps.recalibrated.vcf.gz
```

> **Important:** High-confidence training sets must match the genome build (e.g., GRCh38/hg19) and be indexed correctly. Indels must be recalibrated separately (`-mode INDEL`) using indel-specific truth resources.

---

## 4. Hard Filtering vs. VQSR Comparison

| Feature | Hard Filtering | VQSR |
| :--- | :--- | :--- |
| **Underlying Logic** | Static metric cutoffs | Statistical GMM clustering |
| **Annotation Interaction** | Independent criteria evaluation | Evaluates joint metric distributions |
| **Truth Resource Requirement** | None | Requires high-confidence truth sets |
| **Dataset Size Suitability** | Excellent for small/single samples | Requires large cohorts ($\ge 30$ WGS or $\ge 100$ WES samples) |
| **Computational Cost** | Minimal | High |
| **Interpretability** | Transparent and deterministic | Probabilistic score slices (VQSLOD) |
| **Non-Model Organisms** | Primary recommended strategy | Rarely applicable due to missing resources |

---

## 5. Multiple-Testing Correction: Benjamini-Hochberg (FDR)

While VCF quality filtering eliminates **technical false positives**, association studies require controlling **statistical false positives** across thousands or millions of concurrent hypothesis tests.

### The Multiple-Testing Problem

Performing $m = 10,000$ independent statistical tests at significance level $\alpha = 0.05$ under a true null hypothesis yields an expected number of false discoveries:

$$
E(\text{False Positives}) = m \cdot \alpha = 10,000 \times 0.05 = 500 \text{ false discoveries}
$$

---

### False Discovery Rate (FDR) Formulation

FDR represents the expected proportion of false positive discoveries ($V$) among all rejected null hypotheses ($R$):

$$
\text{FDR} = E \left[ \frac{V}{\max(R, 1)} \right]
$$

### The Benjamini-Hochberg (BH) Procedure

1. Sort raw $p$-values in ascending order: $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$.
2. Assign rank $i \in \{1, 2, \dots, m\}$ to each $p$-value.
3. Compute BH significance thresholds:

$$
p_{(i)} \le \frac{i}{m} \cdot q
$$

4. Identify the maximum rank $k$ satisfying the inequality.
5. Declare hypotheses with $p_{(1)}, \dots, p_{(k)}$ statistically significant at target FDR $q$.

---

### BH Worked Example ($m = 5$, Target $q = 0.05$)

| Rank ($i$) | Gene | Raw $p$-value ($p_{(i)}$) | BH Threshold ($\frac{i}{m} \cdot q$) | Condition ($p_{(i)} \le \text{Threshold}$) | Adjusted $p$-value ($q$-value) |
| :---: | :--- | :---: | :---: | :---: | :---: |
| **1** | Gene A | $0.001$ | $\frac{1}{5} \times 0.05 = 0.010$ | **Pass** | $0.005$ |
| **2** | Gene B | $0.008$ | $\frac{2}{5} \times 0.05 = 0.020$ | **Pass** | $0.020$ |
| **3** | Gene C | $0.012$ | $\frac{3}{5} \times 0.05 = 0.030$ | **Pass** | $0.020$ |
| **4** | Gene D | $0.030$ | $\frac{4}{5} \times 0.05 = 0.040$ | **Pass** | $0.0375$ |
| **5** | Gene E | $0.200$ | $\frac{5}{5} \times 0.05 = 0.050$ | Fail | $0.200$ |

*Result:* Tests 1 through 4 remain significant at $\text{FDR} \le 0.05$.

---

### Executing BH Correction in R

```r
# Input raw p-values from rare-variant burden testing or GWAS
p_values <- c(0.001, 0.008, 0.012, 0.030, 0.200)

# Compute Benjamini-Hochberg FDR adjusted p-values
fdr_adjusted <- p.adjust(p_values, method = "BH")

# Format output data frame
results <- data.frame(
  Raw_PValue = p_values,
  FDR_Adjusted = fdr_adjusted,
  Is_Significant = fdr_adjusted < 0.05
)

print(results)
```

---

### Technical QC vs. Statistical Significance

| Feature | Technical Quality Filtering (VCF QC) | Statistical Correction (BH FDR) |
| :--- | :--- | :--- |
| **Evaluates** | Was the base accurately sequenced and called? | Is the variant associated with a disease/trait beyond chance? |
| **Input** | Reads, alignments, Phred quality scores | Statistical test $p$-values across genes/variants |
| **Primary Risk** | Sequencing errors, misaligned reads | False positive discoveries from multiple hypothesis testing |
| **Methods** | `QUAL`, `QD`, `FS`, `MQ`, `VQSR` | Benjamini-Hochberg, Bonferroni |

---

## 6. End-to-End Quality Filtering Workflow

```
                   [ Raw FASTQ Files ]
                            │
                            ▼
               [ Read Alignment (BWA-MEM) ]
                            │
                            ▼
              [ Variant Calling (HaplotypeCaller) ]
                            │
                            ▼
                   [ Raw VCF File ]
                            │
                            ▼
        ┌──────────────────────────────────────┐
        │       Variant-Level Filtering        │
        ├──────────────────────────────────────┤
        │ Select Strategy based on Callset Size│
        │ • Hard Filtering (Small/Custom Set)  │
        │ • VQSR (Large WGS/WES Datasets)      │
        └──────────────────────────────────────┘
                            │
                            ▼
                [ PASS-Only Variant Subset ]
                            │
                            ▼
        ┌──────────────────────────────────────┐
        │      Genotype-Level Filtering        │
        ├──────────────────────────────────────┤
        │ Mask individual sample genotypes:    │
        │ • GQ < 20 (Genotype Quality)         │
        │ • DP < 10 (Sample Read Depth)        │
        │ • VAF deviation on Heterozygotes     │
        └──────────────────────────────────────┘
                            │
                            ▼
        ┌──────────────────────────────────────┐
        │        Sample-Level Filtering        │
        ├──────────────────────────────────────┤
        │ Filter out low-quality samples:      │
        │ • Excessive missingness rates        │
        │ • Abnormal Ti/Tv ratio outliers      │
        └──────────────────────────────────────┘
                            │
                            ▼
            [ Downstream Association Study ]
                            │
                            ▼
           [ Multiple-Testing Correction (BH) ]
                            │
                            ▼
             [ Validated Biological Hits ]
```

---

## 7. Additional Quality-Control Considerations

### 7.1 Variant Allele Fraction (VAF)
For diploid samples, heterozygous calls should exhibit $VAF \approx 0.5$:

$$
VAF = \frac{AD_{\text{ALT}}}{AD_{\text{REF}} + AD_{\text{ALT}}}
$$

Significant deviations ($VAF < 0.2$ or $VAF > 0.8$) in diploid germline calls strongly signal somatic mosaicism, allele dropout, or mapping bias.

---

### 7.2 Variant Class Differentiation
SNPs and indels exhibit distinct error profiles. Indels suffer disproportionately from local alignment shifts and homopolymer slippage. Filtering thresholds must be applied independently:

```bash
# Extract SNPs for separate filtering
gatk SelectVariants -V raw.vcf.gz -select-type SNP -O snps.vcf.gz

# Extract Indels for separate filtering
gatk SelectVariants -V raw.vcf.gz -select-type INDEL -O indels.vcf.gz
```

---

### 7.3 Missing Annotations
Annotations like `MQRankSum` or `ReadPosRankSum` require reads supporting both reference and alternate alleles. Sites lacking alternate read observations will report missing annotation values (`.`). **Missing values should not be automatically treated as failed filters.**

---

## 8. Command-Line Tool Cheatsheet

### 8.1 BCFtools Filtering Examples

```bash
# Filter variants by QUAL and total depth
bcftools filter -i 'QUAL > 30 && INFO/DP > 10' input.vcf.gz -Oz -o filtered.vcf.gz

# Filter individual sample genotypes lacking DP >= 10 or GQ >= 20
bcftools filter -S . -e 'FMT/DP < 10 || FMT/GQ < 20' input.vcf.gz -Oz -o genotype_masked.vcf.gz
```

### 8.2 VCFtools Filtering Examples

```bash
# Hard filter by site quality, depth, and genotype quality
vcftools \
    --gzvcf input.vcf.gz \
    --minQ 30 \
    --minDP 10 \
    --minGQ 20 \
    --recode --recode-INFO-all \
    --out vcftools_filtered
```

---

## 9. Best Practices Checklist

1. **Header Verification:** Always verify field definitions in the VCF header (`##INFO` and `##FORMAT`).
2. **Annotation Inspection:** Plot metric distributions (`QD`, `FS`, `MQ`) using R or Python prior to setting thresholds.
3. **Preserve Filtering Traceability:** Retain failing variants tagged with explicit reason flags rather than immediately deleting records.
4. **Separate QC Tiers:** Perform site-level filtering prior to sample genotype masking.
5. **Report Complete Parameters:** Document caller versions, reference genome build, training resource versions, and exact command flags for reproducible scientific reporting.