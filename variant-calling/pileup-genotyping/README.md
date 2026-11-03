# Pileup-Based Genotyping

## 1. What Is Pileup-Based Genotyping?

Pileup-based genotyping is a method for identifying genetic variants by examining the sequencing reads aligned to each position of a reference genome.

During whole-genome sequencing (WGS) or whole-exome sequencing (WES), millions of short DNA fragments are sequenced. These reads are aligned to a reference genome, allowing us to determine which genomic positions they cover and which nucleotide they contain at each position.

At a particular genomic position, multiple reads may overlap. If we stack the observed nucleotides from these reads vertically, we obtain a **pileup**.

For example, consider chromosome 7, position 100,000:

```text
Reference: A

Read 1:    A
Read 2:    A
Read 3:    G
Read 4:    A
Read 5:    G
Read 6:    A
Read 7:    G
Read 8:    A
```

The pileup contains:
- **A**: 5 reads
- **G**: 3 reads
- **Total depth**: 8 reads

The presence of G bases suggests a possible $A > G$ variant. However, observing an alternative nucleotide does not automatically mean that a true genetic variant exists.

Some alternative bases may result from:
- Sequencing errors.
- Poor-quality base calls.
- Alignment errors.
- PCR or library preparation artifacts.
- Contamination or sample mixture.

The central question is:

> **Does the observed pileup provide sufficient evidence that the individual carries a genuine A/G variant, or can the alternative G observations be explained by sequencing errors?**

Pileup-based genotyping addresses this question by calculating the probability of the observed sequencing data under each possible genotype.

---

## 2. Understanding the Pileup Representation

A pileup summarizes the sequencing evidence at a genomic position.

For each read covering a position, the pileup provides information such as:

| Information | Description |
|---|---|
| **Chromosome** | Chromosome containing the position |
| **Position** | Genomic coordinate |
| **Reference base** | Nucleotide in the reference genome |
| **Observed base** | Nucleotide observed in an aligned read |
| **Base quality** | Confidence in the observed nucleotide |
| **Mapping quality** | Confidence that the read is aligned to the correct genomic location |
| **Strand** | Whether the read originates from the forward or reverse strand |

Two important quality measures are used during variant calling:

### 2.1 Base Quality

Base quality measures the confidence of the sequencer in an individual nucleotide call. It is represented using a Phred quality score:

$$Q = -10 \log_{10}(e)$$

where $e$ is the probability that the base was called incorrectly.

Rearranging:

$$e = 10^{-Q/10}$$

**Phred Score Examples:**

| Phred score | Probability of error | Approximate accuracy |
|:---|:---|:---|
| **Q10** | 0.1 | 90% |
| **Q20** | 0.01 | 99% |
| **Q30** | 0.001 | 99.9% |
| **Q40** | 0.0001 | 99.99% |

A G observed with Q40 provides stronger evidence than a G observed with Q10.

### 2.2 Mapping Quality

Mapping quality measures confidence in the alignment of a sequencing read to the reference genome. A read may contain high-quality nucleotide calls but still be aligned to the wrong genomic location. For example, reads originating from repetitive genomic regions can sometimes map ambiguously.

Therefore, reliable genotyping should consider both:
- **Base quality:** Is the observed nucleotide trustworthy?
- **Mapping quality:** Is the read likely aligned to the correct position?

Variant callers may also account for strand bias, read position, local alignment characteristics, and other sources of systematic error.

---

## 3. How Pileup-Based Genotyping Works

Pileup-based genotyping evaluates the sequencing evidence at each genomic position and determines which genotype best explains the observed bases.

For a diploid organism, each individual normally carries two alleles at an autosomal genomic position.

At a biallelic A/G position, the possible genotypes are:
- **AA**: Homozygous reference
- **AG**: Heterozygous alternate
- **GG**: Homozygous alternate

For a general diploid position containing four possible nucleotides (A, C, G, T), there are ten unordered genotypes:

$$\frac{4(4+1)}{2} = 10$$

These ten genotypes are:

```text
AA, AC, AG, AT,
CC, CG, CT,
GG, GT,
TT
```

The genotyping process can be divided into four main stages:

### Step 1: Build the pileup
Collect all sequencing reads that overlap a genomic position. For every read, record:
- Observed nucleotide.
- Base quality score.
- Mapping quality.
- Read orientation and other relevant alignment information.

**Example:**

```text
Chromosome: 7
Position:   100000
Reference:  A

Observed bases:
A A A A A A A A G G G G G G G G G G G G G G G G G G G

Base qualities:
35 37 32 40 36 38 30 35 34 36 38 35
37 32 40 36 35 39 37 34 38 36 40 35 37
```

The pileup contains:
- **A**: 15 reads
- **G**: 12 reads
- **Total depth**: 27 reads

The observed alternate allele fraction is:

$$AF = \frac{12}{27} = 0.444$$

An allele fraction close to 0.5 is consistent with a heterozygous variant in a diploid sample. However, allele balance alone is not sufficient to establish a genotype. Sequencing quality, mapping quality, biological context, and other evidence must also be considered.

### Step 2: Calculate genotype likelihoods
A genotype likelihood measures how probable the observed sequencing data would be if a particular genotype were true.

For example, consider the three possible genotypes at an A/G position ($AA$, $AG$, $GG$). The caller calculates:

$$P(D \mid AA), \quad P(D \mid AG), \quad P(D \mid GG)$$

where $D$ represents all observed bases and their associated quality information.

A genotype that explains the observed pileup well will have a higher likelihood:
- Under **AA**, observing G bases requires sequencing errors.
- Under **AG**, both A and G bases are expected.
- Under **GG**, observing A bases requires sequencing errors.

If the observed A and G bases have high quality and are distributed approximately equally, $AG$ will have the highest likelihood.

### Step 3: Apply Bayesian inference
Genotype likelihoods describe how well each genotype explains the data. Bayesian inference combines these likelihoods with prior probabilities.

Using Bayes' theorem:

$$P(G \mid D) = \frac{P(D \mid G) \, P(G)}{P(D)}$$

where:
- $P(G \mid D)$: Posterior probability of genotype $G$.
- $P(D \mid G)$: Genotype likelihood.
- $P(G)$: Prior probability of genotype $G$.
- $P(D)$: Probability of observing the data across all possible genotypes.

For comparing genotypes:

$$P(G \mid D) \propto P(D \mid G) \, P(G)$$

The prior can incorporate information about expected genotype frequencies, population allele frequencies, or other assumptions.

### Step 4: Select the genotype and calculate confidence
The genotype with the highest posterior probability or genotype likelihood, depending on the caller, is selected as the best genotype.

```text
Position: 7:100000
Reference: A

Genotype: AG
Genotype likelihoods:
AA: 0.00001
AG: 0.99997
GG: 0.00002
```

*(Note: These values are illustrative and do not reflect specific software outputs.)*

The caller can also report confidence measures such as:
- **GQ (Genotype Quality):** Confidence in the selected genotype.
- **PL (Phred-scaled genotype likelihoods):** Relative likelihoods of possible genotypes.
- **GL (Genotype Likelihoods):** Log10-scaled genotype likelihoods.

Genotype quality is commonly defined using the posterior probability of the selected genotype ($G_{\mathrm{best}}$):

$$GQ = -10 \log_{10}\left(1 - P(G_{\mathrm{best}} \mid D)\right)$$

---

## 4. Mathematical Foundation

The central mathematical component of pileup-based genotyping is the genotype likelihood.

### 4.1 Converting Base Quality into Error Probability

For a sequencing read with base quality $Q_i$:

$$e_i = 10^{-Q_i/10}$$

For example, if a nucleotide has Q30:

$$e_i = 10^{-30/10} = 0.001$$

This indicates that the probability of an incorrect base call is approximately 0.1%.

### 4.2 Probability of Observing a Base Given a True Allele

Assume a simplified sequencing error model in which:
- The correct nucleotide is observed with probability $1 - e_i$.
- An incorrect nucleotide is observed with probability $e_i / 3$.
- Errors are equally distributed among the other three nucleotides.

$$P(b_i \mid a) = \begin{cases} 1 - e_i & \text{if } b_i = a \\ e_i / 3 & \text{if } b_i \neq a \end{cases}$$

where:
- $b_i$ is the nucleotide observed in read $i$.
- $a$ is the true allele.
- $e_i$ is the base-calling error probability.

### 4.3 Probability of an Observation Given a Diploid Genotype

For a diploid genotype containing alleles $a_1$ and $a_2$, assuming reads originate equally from both chromosomes:

$$P(b_i \mid G) = \frac{1}{2} P(b_i \mid a_1) + \frac{1}{2} P(b_i \mid a_2)$$

- **For a heterozygous AG genotype:**
  $$P(b_i \mid AG) = \frac{1}{2} P(b_i \mid A) + \frac{1}{2} P(b_i \mid G)$$

- **For a homozygous AA genotype:**
  $$P(b_i \mid AA) = P(b_i \mid A)$$

### 4.4 Calculating the Genotype Likelihood

Assuming observations are conditionally independent given the genotype:

$$L(G) = P(D \mid G) = \prod_{i=1}^{n} P(b_i \mid G)$$

To prevent numerical underflow when multiplying many small probabilities, variant callers work with log-likelihoods:

$$\log L(G) = \sum_{i=1}^{n} \log P(b_i \mid G)$$

### 4.5 Full Substitution Error Matrix

A more sophisticated error model uses a substitution probability matrix $M_{b,a}$:

$$M_{b,a} = P(\text{observed base } b \mid \text{true base } a)$$

| True Base | Observed A | Observed C | Observed G | Observed T |
|:---|:---|:---|:---|:---|
| **A** | 0.997 | 0.001 | 0.0015 | 0.0005 |
| **C** | 0.001 | 0.997 | 0.001 | 0.001 |
| **G** | 0.001 | 0.0015 | 0.997 | 0.0005 |
| **T** | 0.0005 | 0.001 | 0.001 | 0.9975 |

---

## 5. Worked Example: Is A/G a True Heterozygous Variant?

Consider WGS data at a suspected locus (chr7:100,000, Ref: A):

```text
Reference: A

A A A A A A A A A A A A A A A
G G G G G G G G G G G G

Total depth: 27
```

| Base | Read Count |
|:---|:---|
| **A** | 15 |
| **G** | 12 |
| **Total** | **27** |

The alternate allele fraction is:

$$AF = \frac{12}{27} = 0.444$$

### Evaluating Hypotheses:

1. **Hypothesis 1 (AA):** Every G observation must be explained by sequencing error. Observing 12 G bases out of 27 reads is extremely improbable under high read qualities.
2. **Hypothesis 2 (AG):** Both A and G are expected observations ($AF_A \approx 0.5, AF_G \approx 0.5$). An observed $AF$ of 0.444 strongly aligns with this hypothesis.
3. **Hypothesis 3 (GG):** The 15 observed A bases must be explained by sequencing error, which is also extremely unlikely.

**Conclusion:** The evidence strongly supports $AG$ as the most plausible genotype.

---

## 6. Genotype Likelihoods vs. Allele Frequency

- **Allele Frequency (AF):** Describes the fraction of observed reads supporting an allele:
  $$AF = \frac{\text{Alternate allele reads}}{\text{Total allele-supporting reads}}$$
- **Genotype Likelihood (GL):** Describes how probable the observed read data are under a candidate genotype, accounting for quality scores and confidence.

Two samples may both have $AF = 0.5$, but if Sample 1 features high-quality reads and Sample 2 features low-quality, ambiguously mapped reads, Sample 1 provides much stronger statistical support for a true variant.

---

## 7. Pileup-Based Variant Calling with BCFtools

**BCFtools** is a popular software suite for pileup-based variant calling:

```bash
# Step 1: Generate pileup likelihoods
bcftools mpileup \
    -f reference.fa \
    -Ou \
    sample.bam \
    -o sample.bcf

# Step 2: Call variants
bcftools call \
    -mv \
    -Ob \
    -o sample_variants.bcf
```

**Piped Single-Line Workflow:**

```bash
bcftools mpileup -f reference.fa -Ou sample.bam | \
bcftools call -mv -Ob -o sample_variants.bcf
```

| Option | Description |
|:---|:---|
| `-f` | Reference genome FASTA file |
| `-Ou` | Output uncompressed BCF |
| `-m` | Use multiallelic calling model |
| `-v` | Output variant sites only |
| `-Ob` | Output compressed BCF |

---

## 8. Pileup-Based Genotyping vs. Haplotype-Based Variant Calling

| Feature | Pileup-Based Approach | Haplotype-Based Approach |
|:---|:---|:---|
| **Primary Unit** | Individual genomic position | Local genomic region |
| **Evidence** | Bases and quality scores at a position | Read support for candidate haplotypes |
| **Nearby Variants** | Evaluated independently | Evaluated jointly |
| **Complex Variants** | Challenging to detect accurately | Better suited for indels and complex variants |
| **Computational Cost**| Lower | Higher |
| **Examples** | `BCFtools mpileup/call` | `GATK HaplotypeCaller` |

---

## 9. The Four-Base Mixture Model

A four-base mixture model evaluates the simultaneous contribution of A, C, G, and T at a given genomic coordinate:

```text
Reference: A

Observed pileup:
A = 40
C = 2
G = 35
T = 1

Total = 78
```

This model determines whether the observed distribution represents a single base with sequencing noise, a diploid heterozygous mixture, or a complex multiallelic/somatic mixture.

---

## 10. Computational Complexity

Let:
- $P$ = Number of genomic positions evaluated.
- $R_p$ = Number of reads covering position $p$.
- $G$ = Number of possible genotypes ($G = 10$ for diploid organisms).

### Time Complexity
- Per position: $O(R_p \times G)$
- Across all positions: $O\left(G \sum_{p=1}^{P} R_p\right) \approx O(P \times R \times G)$

Because $G = 10$ is fixed, computational time scales linearly with total read depth across all positions.

### Space Complexity
Streaming pileup algorithms require temporary per-position memory proportional to the depth: $O(R)$.

---

## 11. Applications in Bioinformatics

1. **Whole-Genome Sequencing (WGS):** Detection of genome-wide SNVs and small indels.
2. **Whole-Exome Sequencing (WES):** Identification of coding variants associated with Mendelian conditions.
3. **Rare Disease Diagnostics:** Classifying clinical variants (homozygous/heterozygous/artifact).
4. **Population Genetics:** Utilizing genotype likelihoods directly for low-coverage sequencing data to avoid bias from hard genotype calls.

---

## 12. Limitations and Important Considerations

- **Low Sequencing Depth:** High sampling variance can lead to missed heterozygous calls.
- **Mapping Errors:** Misaligned reads in repetitive regions create false-positive variants.
- **Strand Bias:** Alternate bases appearing on only one strand indicate technological artifacts.
- **Allelic Imbalance:** Reference bias, copy-number variations, or tumor purity can shift $AF$ away from expected ratios ($0.5$ or $1.0$).
- **Non-Diploid Systems:** Mitochondria, sex chromosomes, and polyploid organisms require specialized non-diploid models.

---

## 13. Connection to Core Computational Concepts

| Concept | Connection to Genotyping |
|:---|:---|
| **Statistical Inference** | Identifying the most probable underlying genotype from noisy read data |
| **Bayesian Statistics** | Incorporating prior allele frequencies with likelihood calculations |
| **Maximum Likelihood Estimation** | Identifying genotypes that maximize data observational probabilities |
| **Log-Space Computation** | Preventing floating-point underflow during likelihood multiplication |

---

## 14. Key Takeaways

- A **pileup** stacks aligned sequencing reads at a single genomic position.
- **Base Quality** and **Mapping Quality** convert read metadata into error probabilities.
- **Genotype Likelihoods** quantify how well different candidate genotypes explain the observed pileup.
- **BCFtools** separates likelihood generation (`mpileup`) from genotype calling (`call`).
- **Haplotype-based approaches** offer superior accuracy in complex or indel-heavy regions, whereas **pileup-based methods** offer rapid, linear-time evaluation.

---

## References and Further Reading

1. **BCFtools Documentation** – Pileup generation, genotype likelihoods, and variant calling.  
   [https://samtools.github.io/bcftools/bcftools.html](https://samtools.github.io/bcftools/bcftools.html)
2. **Li H. (2011).** A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*.  
   [https://doi.org/10.1093/bioinformatics/btr509](https://doi.org/10.1093/bioinformatics/btr509)
3. **Li H. (2011).** Improving SNP discovery by base alignment quality. *Bioinformatics*.  
   [https://doi.org/10.1093/bioinformatics/btr076](https://doi.org/10.1093/bioinformatics/btr076)
4. **McKenna et al. (2010).** The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*.  
   [https://doi.org/10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110)
5. **DePristo et al. (2011).** A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics*.  
   [https://doi.org/10.1038/ng.806](https://doi.org/10.1038/ng.806)
6. **Broad Institute GATK Documentation** – Genotype likelihoods and variant calling concepts.  
   [https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)