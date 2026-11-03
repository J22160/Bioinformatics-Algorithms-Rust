# Base Quality Score Recalibration (BQSR)

## What Is Base Quality Score Recalibration?

Sequencing technologies assign a quality score to every base in a sequencing read. This score estimates the probability that the observed nucleotide is incorrect.

However, the quality scores reported by sequencing instruments are not always perfectly calibrated. Systematic technical factors can cause the actual sequencing error rate to differ from the error rate predicted by the reported quality scores.

**Base Quality Score Recalibration (BQSR)** is a statistical procedure that adjusts reported base quality scores to better reflect their empirically observed error rates.

BQSR is commonly used in germline variant-calling workflows, particularly those based on the Genome Analysis Toolkit (GATK).

The main objectives are to:

- Identify systematic biases in reported base quality scores.
- Estimate empirical error rates using sequencing data.
- Account for technical factors associated with sequencing errors.
- Produce recalibrated quality scores that better represent base-level confidence.
- Improve the accuracy of downstream variant calling.

BQSR is generally performed after read alignment and before variant calling.

---

## 1. Understanding Base Quality Scores

A base quality score, or Phred quality score, represents the estimated probability that a particular base call is incorrect.

It is defined as:

$$Q = -10 \log_{10}(P_{\text{error}})$$

Where:

- $Q$ = Phred quality score.
- $P_{\text{error}}$ = probability that the base is incorrect.

Rearranging:

$$P_{\text{error}} = 10^{-Q / 10}$$

### Examples

| Phred score | Probability of error | Approximate accuracy |
| :--- | ---: | ---: |
| Q10 | 0.1 | 90% |
| Q20 | 0.01 | 99% |
| Q30 | 0.001 | 99.9% |
| Q40 | 0.0001 | 99.99% |
| Q50 | 0.00001 | 99.999% |

For example, a base with Q30 is expected to have an error probability of approximately 0.1%.

This does not guarantee that the base is correct. It means that, among bases assigned a comparable quality score under a well-calibrated model, approximately 0.1% are expected to be incorrect.

### Why do quality scores need recalibration?

Sequencing instruments estimate base quality scores using their internal error models.

These models may not fully account for all technical sources of error.

For example, base-calling errors can be influenced by:

- The position of a base within a sequencing read.
- The surrounding nucleotide sequence.
- Sequencing cycle effects.
- Read group or sequencing batch.
- Instrument-specific characteristics.
- Library preparation and sequencing chemistry.

Consequently, two bases assigned the same reported quality score may have different actual error probabilities.

Consider the following hypothetical example:

| Base group | Reported quality | Observed error rate |
| :--- | ---: | ---: |
| Group A | Q30 | 0.1% |
| Group B | Q30 | 0.5% |
| Group C | Q30 | 0.05% |

All three groups have the same reported quality score, but their empirical error rates differ.

BQSR attempts to identify these systematic differences and adjust the scores accordingly.

---

## 2. How BQSR Works

BQSR uses observed sequencing mismatches to estimate systematic deviations between reported and empirical base quality scores.

The general workflow consists of four major steps:

1. Collect covariates.
2. Identify likely sequencing errors.
3. Build recalibration tables.
4. Apply recalibrated quality scores.

### Step 1: Collect covariates

A covariate is a variable that may help explain systematic variation in sequencing error rates.

GATK BQSR primarily models the following covariates:

#### A. Reported quality score

The original Phred quality score assigned by the sequencing instrument.

For example:

```text
Q20, Q25, Q30, Q35, Q40
```

The model evaluates whether bases assigned the same reported quality score have similar empirical error rates.

#### B. Read group

A read group identifies the sequencing context associated with a read.

It can contain information such as:

- Sample identifier.
- Sequencing library.
- Sequencing platform.
- Platform unit or flowcell/lane identifier.

Read group information is important because different sequencing runs, libraries, or batches may have different error profiles.

For example, two sequencing runs may assign Q30 to bases but have different observed error rates.

BQSR can account for these differences through read-group-specific recalibration.

#### C. Read position

The position of a base within a sequencing read can influence its error probability.

For example, sequencing errors may increase toward the end of a read because of declining signal quality or other cycle-dependent effects.

A simplified example:

| Read position | Reported quality | Observed error rate |
| :--- | ---: | ---: |
| 1–20 | Q30 | 0.05% |
| 21–70 | Q30 | 0.10% |
| 71–100 | Q30 | 0.40% |

The original quality scores do not capture this positional difference, but recalibration can account for it.

#### D. Sequence context

The probability of a sequencing error can depend on the surrounding nucleotide sequence.

For example, errors may occur more frequently in certain repetitive or homopolymer contexts.

GATK uses a context-dependent covariate that considers the neighboring bases around the base being evaluated.

The default context covariate in GATK BQSR is based on the preceding nucleotide context, with a default context size of 2 bases.

This allows the model to identify sequence-dependent error patterns.

> **Note:** GATK's implementation uses specific covariate definitions and context encodings. The description above is conceptual rather than a complete representation of its internal encoding.

---

### Step 2: Identify likely sequencing errors

BQSR needs to estimate which observed bases are incorrect.

However, without a known correct sequence, it is impossible to directly determine whether every mismatch is a sequencing error or a genuine biological variant.

GATK addresses this problem using a known-sites VCF containing previously identified polymorphic positions.

These known variant positions are excluded from the mismatch-counting process.


**Why?**

Consider a genomic position where the reference genome contains A, but a sample genuinely carries G.

An aligned read containing G may appear to be a mismatch against the reference genome.

However, this mismatch is not necessarily a sequencing error. It may represent a genuine genetic variant.

If these mismatches were treated as sequencing errors, BQSR could incorrectly reduce the quality scores of genuine alternate alleles.

Therefore, known polymorphic positions are masked during recalibration.

#### How are errors counted?

At positions not masked by known variants, GATK compares the observed read bases with the reference genome.

A simplified example:

| Reference base | Observed base | Interpretation |
| :--- | :--- | :--- |
| A | A | Match |
| C | C | Match |
| G | T | Mismatch |
| T | T | Match |
| A | G | Mismatch |

The two mismatches contribute to the empirical error estimate.

The model also considers the reported quality scores and covariates associated with those bases.

> **Important:** A mismatch is only an estimate of a sequencing error. Some mismatches may still arise from unreported genuine variants, alignment artifacts, or differences between the sample and reference genome.

This is why a suitable known-sites resource and a well-matched reference genome are important.

---

### Step 3: Build recalibration tables

After collecting covariates and mismatch observations, BQSR estimates empirical error rates.

Conceptually, it compares:

- The number of observed sequencing errors.
- The number of bases examined.
- The quality scores originally reported by the sequencing instrument.

A simplified empirical error rate is:

$$P_{\text{empirical}} = \frac{\text{Number of observed errors}}{\text{Total number of observed bases}}$$

The corresponding empirical Phred score is:

$$Q_{\text{empirical}} = -10 \log_{10}(P_{\text{empirical}})$$

For example, suppose 100,000 bases assigned Q30 contain 300 observed errors.

Then:

$$P_{\text{empirical}} = \frac{300}{100000} = 0.003$$

The empirical quality score is:

$$Q_{\text{empirical}} = -10 \log_{10}(0.003) \approx 25.2$$

Although the instrument reported Q30, the observed error rate corresponds to approximately Q25.2.

This suggests that the original quality scores are overestimating base-call accuracy for this group.

In practice, GATK uses a more sophisticated statistical model than this simple calculation. It combines information across covariates and applies smoothing and statistical adjustments to avoid unreliable estimates from small groups.

#### Recalibration tables

GATK generates several recalibration tables.

Conceptually, these capture:

1. **Read-group table:** Systematic differences associated with read groups.
2. **Quality-score table:** Differences associated with reported quality scores within read groups.
3. **Covariate tables:** Additional adjustments associated with read position and sequence context.

The recalibration model uses these effects to estimate an adjusted quality score for each base.

---

### Step 4: Apply recalibrated quality scores

Once the recalibration model has been constructed, it can be applied to the original aligned reads.

The original base quality scores are replaced with recalibrated scores.

For example:

| Base | Original quality | Recalibrated quality |
| :--- | ---: | ---: |
| A | Q30 | Q27 |
| C | Q30 | Q29 |
| G | Q30 | Q24 |
| T | Q35 | Q33 |

These values are illustrative.

Some scores may decrease, while others may increase, depending on the empirical error patterns estimated by the model.

> **Important:** BQSR changes the quality score associated with a base. It does not change the nucleotide itself, remove the read, or directly correct the sequencing error.

---

## 3. Statistical Model Behind BQSR

BQSR estimates how the observed error rate differs from the error rate predicted by the original quality scores.

The original Phred score provides an expected error probability:

$$P_{\text{reported}} = 10^{-Q_{\text{reported}} / 10}$$

The observed mismatches provide an empirical estimate of the actual error probability.

A simplified conceptual recalibration model can be expressed as:

$$Q_{\text{recalibrated}} = -10 \log_{10}(P_{\text{recalibrated}})$$

Where the recalibrated error probability is estimated using the original quality score and relevant covariates.

GATK's model uses a hierarchical empirical-Bayes approach to estimate the effects of different covariates while accounting for limited observations.

### Example: Recalibrating a quality score

Suppose bases assigned Q30 have an empirical error rate of 0.003 in a particular read group and sequence context.

The empirical quality is:

$$Q_{\text{empirical}} = -10 \log_{10}(0.003) \approx 25.2$$

The reported score is Q30, but the empirical estimate is approximately Q25.

The recalibration model may therefore reduce the quality score for bases with similar characteristics.

However, it does not simply replace every Q30 score with Q25. It considers the combined covariate effects and the amount of supporting evidence.

---

## 4. GATK BQSR Workflow

In GATK, BQSR is commonly performed using two tools:

- `BaseRecalibrator`
- `ApplyBQSR`

The general workflow is:

```text
Raw sequencing reads
        |
        v
Read quality control
        |
        v
Read alignment
        |
        v
Aligned BAM/CRAM file
        |
        v
BaseRecalibrator
        |
        v
Recalibration table
        |
        v
ApplyBQSR
        |
        v
Recalibrated BAM/CRAM
        |
        v
Variant calling
        |
        v
VCF generation
```

BQSR is performed on aligned reads because the model needs to compare read bases with the reference genome and account for genomic positions.

### Step 1: Generate the recalibration table

GATK uses `BaseRecalibrator` to examine the aligned reads and estimate recalibration adjustments.

Example:

```bash
gatk BaseRecalibrator \
    -R reference.fasta \
    -I input.bam \
    --known-sites known_variants.vcf.gz \
    -O recalibration.table
```

Parameters:

- `-R`: Reference genome.
- `-I`: Input aligned BAM or CRAM file.
- `--known-sites`: Known variant sites excluded from mismatch counting.
- `-O`: Output recalibration table.

The known-sites file should be compatible with the reference genome assembly used for alignment.

For example, a GRCh38-aligned dataset should use an appropriately matched GRCh38 known-sites resource.

Depending on the dataset, multiple known-sites VCFs can be supplied.

### Step 2: Apply recalibration

The recalibration table is passed to `ApplyBQSR`.

```bash
gatk ApplyBQSR \
    -R reference.fasta \
    -I input.bam \
    --bqsr-recal-file recalibration.table \
    -O recalibrated.bam
```

This generates a BAM file containing recalibrated base quality scores.

The nucleotide sequences and alignments are not intentionally changed by this operation.

### Step 3: Index the recalibrated BAM

```bash
samtools index recalibrated.bam
```

The recalibrated BAM can now be used in downstream variant-calling workflows.

### Step 4: Optional recalibration diagnostics

GATK provides `AnalyzeCovariates` to visualize and compare reported versus empirical quality scores before and after recalibration.

```bash
gatk AnalyzeCovariates \
    -before recalibration.table \
    -after recalibration.table \
    -plots recalibration_plots.pdf
```

In a typical workflow, the `-before` and `-after` inputs should be separate recalibration tables generated before and after applying recalibration, respectively, using the same read dataset and compatible settings.

These plots help assess whether the recalibration model identified systematic biases and how the estimated quality scores changed.

---

## 5. Understanding BQSR Recalibration Plots

Recalibration plots help evaluate systematic differences between reported and empirical base quality scores.

Common plots include:

### A. Quality score plot

Compares reported quality scores with empirical quality estimates.

If a reported Q30 group has an empirical quality near Q25, the original scores may be overconfident.

After recalibration, the adjusted scores should better reflect the empirical error estimates.

### B. Read position plot

Shows how empirical quality varies with the position of a base within a read.

A systematic decrease in quality toward the end of reads may indicate cycle- or position-dependent errors.

### C. Sequence context plot

Shows how error rates vary across different nucleotide contexts.

Some contexts may exhibit higher mismatch rates than others.

### D. Read group plot

Compares recalibration behavior across sequencing read groups.

This can help identify differences between sequencing runs, libraries, or batches.

---

## 6. BQSR vs. VCF Quality Filtering

BQSR and variant quality filtering are both quality-control procedures used in variant-calling workflows, but they operate at different levels.

| Feature | BQSR | VCF Quality Filtering |
| :--- | :--- | :--- |
| **Main purpose** | Calibrate base quality scores | Identify unreliable variant calls |
| **Input** | Aligned sequencing reads (BAM/CRAM) | VCF |
| **Level** | Individual bases within reads | Variant sites or genotypes |
| **Main information** | Base qualities, read groups, sequence context, read position | QUAL, DP, QD, FS, MQ, and other annotations |
| **Output** | Recalibrated BAM/CRAM | Filtered or annotated VCF |
| **Common GATK tools** | BaseRecalibrator, ApplyBQSR | VariantFiltration, VariantRecalibrator, ApplyVQSR |
| **Timing** | Before variant calling | After variant calling |

### Why are both important?

Consider a sequencing read that contains an incorrect nucleotide but has been assigned an overly optimistic Q35 score.

Without recalibration, the variant caller may treat this base as stronger evidence than it deserves.

BQSR attempts to correct the quality score to better reflect the empirical error rate.

The variant caller can then use the recalibrated base qualities when calculating genotype and variant likelihoods.

After variant calling, VCF quality filtering evaluates the resulting calls using site-level and genotype-level evidence.

Therefore, BQSR can influence downstream variant calling, while VCF filtering evaluates the calls produced by that process.

---

## 7. When Should BQSR Be Used?

BQSR is commonly used in short-read germline variant-calling workflows, particularly GATK Best Practices pipelines.

It is useful when:

- The sequencing data contain systematic base-calling errors.
- Suitable known variant resources are available.
- The sequencing technology and data characteristics are compatible with the recalibration model.
- The downstream variant-calling workflow expects or benefits from recalibrated base qualities.

However, BQSR is not universally required for every sequencing dataset.

Considerations include:

### Known-sites resources

BQSR relies on known polymorphic sites to avoid treating genuine variants as sequencing errors.

If the known-sites resource is incomplete or poorly matched to the sample population, some genuine variants may remain unmasked and contribute to the error estimates.

This can be especially relevant for underrepresented populations or non-human organisms.

### Small datasets

When relatively few bases are available, some combinations of covariates may have limited observations.

GATK uses statistical adjustments to reduce the impact of sparse observations, but recalibration may still be less informative.

### Non-model organisms

Suitable known-sites resources may not exist for some organisms.

In such cases, the standard BQSR workflow may not be appropriate without a carefully validated alternative strategy.

### Sequencing technology

BQSR is most commonly associated with short-read sequencing workflows.

Its assumptions and usefulness should be evaluated for different sequencing technologies, including long-read platforms, rather than assuming that the same workflow applies universally.

### Already calibrated data

Some sequencing workflows may already produce quality scores with different calibration characteristics.

Avoid assuming that recalibration will always improve data quality. Its impact should be evaluated using the relevant sequencing and variant-calling workflow.

---

## 8. Limitations of BQSR

Although BQSR improves base quality calibration, it does not solve every source of sequencing error.

Important limitations include:

1. **Dependence on known-sites resources:** Incomplete known variant databases can result in genuine polymorphisms being counted as mismatches.
2. **Alignment-dependent error estimates:** Incorrect read alignments can produce apparent mismatches that are not sequencing errors.
3. **Limited observations:** Rare covariate combinations may not contain enough bases to estimate error rates reliably.
4. **Does not correct the nucleotide sequence:** BQSR modifies quality scores, not the observed base calls.
5. **Does not remove problematic reads:** Low-quality or misaligned reads require separate quality-control procedures.
6. **Does not eliminate biological variation:** Genuine variants should not be interpreted as sequencing errors merely because they differ from the reference genome.
7. **Technology dependence:** The usefulness of BQSR depends on the sequencing platform, base-calling process, and assumptions of the recalibration model.

---

## 9. Best Practices

1. Align reads to the correct reference genome before running BQSR.
2. Use known-sites resources compatible with the reference genome assembly.
3. Include appropriate known variant databases when available.
4. Preserve read group information in BAM/CRAM files.
5. Inspect recalibration plots to evaluate systematic quality-score biases.
6. Avoid assuming that every mismatch represents a sequencing error.
7. Keep the original BAM/CRAM and recalibration table for reproducibility.
8. Use consistent reference genome versions throughout alignment, recalibration, and variant calling.
9. Evaluate whether BQSR is appropriate for the sequencing technology and dataset.
10. Document GATK versions, input resources, and command-line parameters.

---

## References and Further Reading

1. **GATK Best Practices:** Data pre-processing for variant discovery.  
   [https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)

2. **GATK BaseRecalibrator documentation:**  
   [https://gatk.broadinstitute.org/hc/en-us/articles/360037426611](https://gatk.broadinstitute.org/hc/en-us/articles/360037426611)

3. **GATK ApplyBQSR documentation:**  
   [https://gatk.broadinstitute.org/hc/en-us/articles/360037430871](https://gatk.broadinstitute.org/hc/en-us/articles/360037430871)

4. **GATK AnalyzeCovariates documentation:**  
   [https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)

5. **DePristo, M. A. et al. (2011).** A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics*, 43, 491–498.

6. **McKenna, A. et al. (2010).** The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297–1303.

7. **Li, H. (2011).** A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*, 27(21), 2987–2993.