# Bayesian Genotyping

## 1. What Is Bayesian Genotyping?

Bayesian genotyping is a statistical approach for determining an individual's genotype by combining sequencing evidence with prior biological or population-level knowledge.

When sequencing reads are aligned to a reference genome, multiple reads may cover the same genomic position. These reads provide evidence about which nucleotide alleles are present in the individual's DNA.

However, sequencing is not perfect. Some observed bases may be incorrect because of sequencing errors, mapping errors, or other technical artifacts.

A basic genotype caller asks:

> **Which genotype most likely explains the observed sequencing reads?**

Bayesian genotyping asks a more complete question:

> **Which genotype is most probable given both the observed sequencing reads and our prior knowledge about genetic variation?**

This distinction is important when sequencing coverage is low or the observed evidence is ambiguous.

For example, suppose the reference genome contains A at a particular position, but some sequencing reads contain G.

Possible diploid genotypes include:

* **AA**: Homozygous reference.
* **AG**: Heterozygous.
* **GG**: Homozygous alternate.

If the reads contain both A and G, the heterozygous genotype may explain the observations well. However, if only a few low-quality reads contain G, the evidence may be insufficient to distinguish a true variant from sequencing noise.

Bayesian inference helps evaluate this uncertainty by combining three core elements:

1. **Genotype likelihood:** How well does the genotype explain the observed reads?
2. **Genotype prior:** How probable was the genotype before examining these reads?
3. **Posterior probability:** How probable is the genotype after considering both sources of information?

---

## 2. Bayes' Theorem in Genotype Calling

Bayesian inference provides a mathematical framework for updating prior beliefs after observing new evidence.

For genotype calling:

$$
P(G \mid D) = \frac{P(D \mid G) \, P(G)}{P(D)}
$$

where:

* $G$: A candidate genotype.
* $D$: Observed sequencing data.
* $P(G \mid D)$: Posterior probability of the genotype.
* $P(D \mid G)$: Likelihood of observing the sequencing data under that genotype.
* $P(G)$: Prior probability of the genotype.
* $P(D)$: Overall probability of observing the data.

The denominator normalizes the probabilities so that the posterior probabilities across all candidate genotypes sum to one.

For genotype comparison, the denominator is the same for every candidate genotype. Therefore:

$$
P(G \mid D) \propto P(D \mid G) \, P(G)
$$

The genotype with the highest posterior probability is the most probable genotype under the specified likelihood and prior models.

### 2.1 Genotype Likelihood

The likelihood measures how well a candidate genotype explains the observed sequencing reads.

For a diploid individual, consider three possible genotypes at an A/G position:

$$
G \in \{AA, AG, GG\}
$$

The caller evaluates:

$$
P(D \mid AA), \quad P(D \mid AG), \quad P(D \mid GG)
$$

These probabilities are calculated using the observed bases and their quality scores:

* If the true genotype is **AA**, G observations must be explained by sequencing errors.
* If the true genotype is **AG**, both A and G observations are expected.
* If the true genotype is **GG**, A observations must be explained by sequencing errors.

A genotype likelihood model therefore uses the sequencing evidence to determine which genotype best explains the data.

### 2.2 Genotype Prior

The prior represents the probability of a genotype before considering the sequencing reads. It can incorporate information from:

* Population allele frequencies.
* Expected genotype frequencies.
* Hardy-Weinberg equilibrium assumptions.
* Family or pedigree information.
* Known biological constraints.
* Ploidy and genomic context.

A common simplification is to use population allele frequencies to derive genotype priors.

For a biallelic position with reference allele A and alternate allele G, let $q = P(G)$ and $1-q = P(A)$. Under Hardy-Weinberg equilibrium, the genotype priors are:

$$
P(AA) = (1 - q)^2
$$

$$
P(AG) = 2q(1 - q)
$$

$$
P(GG) = q^2
$$

**Example Calculation:**

Suppose the alternate allele frequency is $q = 0.001$:

* $P(AA) = 0.999^2 = 0.998001$
* $P(AG) = 2(0.999)(0.001) = 0.001998$
* $P(GG) = 0.001^2 = 0.000001$

| Genotype | Prior Probability |
| :--- | ---: |
| **AA** | 0.998001 |
| **AG** | 0.001998 |
| **GG** | 0.000001 |

> **Important:** A heterozygosity rate of approximately 0.001 should not be treated as a universal prior for every genomic position or every genotype. Heterozygosity varies across genomic regions, populations, and loci. The genotype prior should be derived from an appropriate allele-frequency model or another explicitly specified prior.

### 2.3 Posterior Probability

The posterior combines the likelihood and prior.

Consider these illustrative values:

| Genotype | Likelihood | Prior |
| :--- | ---: | ---: |
| **AA** | 0.001 | 0.998001 |
| **AG** | 0.8 | 0.001998 |
| **GG** | 0.0005 | 0.000001 |

The unnormalized posterior scores are:

$$
S(G) = P(D \mid G) \, P(G)
$$

* $S(AA) = 0.001 \times 0.998001 = 0.000998$
* $S(AG) = 0.8 \times 0.001998 = 0.0015984$
* $S(GG) = 0.0005 \times 0.000001 = 0.0000000005$

The normalized posterior is:

$$
P(G \mid D) = \frac{S(G)}{S(AA) + S(AG) + S(GG)}
$$

Even though AA has a much larger prior probability, the strong likelihood for AG makes AG the most probable genotype after observing the sequencing data.

> **Key takeaway:** A strong sequencing signal can overcome a low prior probability, but weak or ambiguous evidence may leave the prior with substantial influence.

---

## 3. Calculating Genotype Likelihoods from Sequencing Reads

Bayesian genotyping relies on a likelihood model that describes how sequencing reads are generated under a candidate genotype.

Suppose a genomic position has $n$ reads:

$$
D = (r_1, r_2, \ldots, r_n)
$$

Each read contains an observed nucleotide $r_i$, a base quality score $Q_i$, and potentially additional alignment information.

### 3.1 Convert Phred Quality into Error Probability

The probability of an incorrect base call is:

$$
e_i = 10^{-Q_i/10}
$$

| Base Quality | Error Probability | Accuracy |
| :--- | ---: | ---: |
| **Q20** | 0.01 | 99% |
| **Q30** | 0.001 | 99.9% |
| **Q40** | 0.0001 | 99.99% |

### 3.2 Probability of Observing a Nucleotide

Under a simplified sequencing error model:

$$
P(r_i \mid b) = \begin{cases} 1 - e_i & \text{if } r_i = b \\ e_i / 3 & \text{if } r_i \neq b \end{cases}
$$

where $r_i$ is the observed nucleotide, $b$ is the true nucleotide, and $e_i$ is the error probability.

### 3.3 Probability of a Read Given a Diploid Genotype

For a diploid genotype containing two alleles $a_1$ and $a_2$, assuming equal contributions from both chromosomes:

$$
P(r_i \mid G) = \frac{1}{2} P(r_i \mid a_1) + \frac{1}{2} P(r_i \mid a_2)
$$

* **For an A/G heterozygote:**

$$
P(r_i \mid AG) = \frac{1}{2} P(r_i \mid A) + \frac{1}{2} P(r_i \mid G)
$$

* **For an A/A homozygote:**

$$
P(r_i \mid AA) = P(r_i \mid A)
$$

### 3.4 Combine Evidence Across Reads

Assuming conditional independence of read observations given the genotype:

$$
P(D \mid G) = \prod_{i=1}^{n} P(r_i \mid G)
$$

In log space to prevent numerical underflow:

$$
\log P(D \mid G) = \sum_{i=1}^{n} \log P(r_i \mid G)
$$

---

## 4. The Four-Base Mixture Model

A useful extension of the basic genotype model is the four-base mixture model ($B = \{A, C, G, T\}$).

Define the underlying nucleotide mixture as $f = (f_A, f_C, f_G, f_T)$ where $f_A + f_C + f_G + f_T = 1$.

* **Homozygous A/A:** $f = (1, 0, 0, 0)$
* **Heterozygous A/G:** $f = (0.5, 0, 0.5, 0)$
* **Heterozygous C/T:** $f = (0, 0.5, 0, 0.5)$

### 4.1 Substitution Error Matrix

A substitution error matrix describes the probability of observing one nucleotide when another is truly present:

$$
M_{r,b} = P(\text{observed base } r \mid \text{true base } b)
$$

| True Base | Observed A | Observed C | Observed G | Observed T |
| :--- | ---: | ---: | ---: | ---: |
| **A** | 0.997 | 0.001 | 0.0015 | 0.0005 |
| **C** | 0.001 | 0.997 | 0.001 | 0.001 |
| **G** | 0.001 | 0.0015 | 0.997 | 0.0005 |
| **T** | 0.0005 | 0.001 | 0.001 | 0.9975 |

### 4.2 Likelihood of a Read Under a Mixture

For an observed read base $r_i$, the probability under mixture $f$ is:

$$
P(r_i \mid f) = \sum_{b \in \{A,C,G,T\}} f_b M_{r_i,b}
$$

For multiple reads in log space:

$$
\log P(D \mid f) = \sum_{i=1}^{n} \log \left( \sum_{b} f_b M_{r_i,b} \right)
$$

### 4.3 Connection to Diploid Genotypes

A diploid genotype corresponds to discrete fixed mixtures (e.g., $AA = (1,0,0,0)$ or $AG = (0.5,0,0.5,0)$). Unrestricted mixture models can estimate continuous allele proportions for somatic variants, tumor samples, or pooled sequencing.

---

## 5. Maximum Likelihood Estimation of Allele Frequency

### 5.1 Simplified Binary Error Model

Consider $n$ total reads, $k$ alternate reads, true alternate allele fraction $p$, and error probability $e$. Assuming a symmetric binary error process:

$$
P(\text{observed alternate}) = (1 - p)e + p(1 - e) = e + p(1 - 2e)
$$

Setting the expected fraction equal to the observed fraction $\hat{q} = k/n$:

$$
\hat{q} = e + p(1 - 2e)
$$

Solving for $p$:

$$
p_{\mathrm{MLE}} = \frac{k/n - e}{1 - 2e}
$$

### 5.2 Interpreting the Estimate

If $n = 100, k = 30, e = 0.01$:

$$
p_{\mathrm{MLE}} = \frac{0.30 - 0.01}{1 - 2(0.01)} = \frac{0.29}{0.98} \approx 0.296 \quad (29.6\%)
$$

If $n = 100, k = 0, e = 0.01$:

$$
p_{\mathrm{MLE}} = \frac{0 - 0.01}{0.98} \approx -0.0102
$$

Because a negative allele frequency is non-physical, the constrained estimate is applied:

$$
\hat{p} = \max\left(0, \, p_{\mathrm{MLE}}\right) = 0
$$

### 5.3 Extension to a Uniform 4-Nucleotide Model

Under a model where errors are distributed uniformly across the 3 non-true bases ($e/3$ each):

$$
P(\text{observed alternate}) = (1 - p)\frac{e}{3} + p(1 - e) = \frac{e}{3} + p\left(1 - \frac{4e}{3}\right)
$$

Solving for $p$:

$$
p_{\mathrm{MLE}} = \frac{\hat{q} - e/3}{1 - 4e/3}
$$

---

## 6. Worked Example: Distinguishing a True Variant from Sequencing Errors

Suppose WGS data from a patient shows the following pileup at a genomic position (Reference: A):

```text
A A A A A A A A A A A A A A A
G G G G G G G G G G G G

Total depth: 27 (A: 15, G: 12)
```

The observed alternate allele fraction is $AF = 12/27 = 0.444$.

1. **Calculate Likelihoods:** $P(D \mid AA), P(D \mid AG), P(D \mid GG)$ are calculated using individual read base qualities. High-quality reads will yield a strong likelihood for $AG$.
2. **Apply Priors:** Assuming $q = 0.001$, Hardy-Weinberg priors give $P(AA) = 0.998001$, $P(AG) = 0.001998$, $P(GG) = 0.000001$.
3. **Calculate Posterior:** $P(G \mid D) \propto P(D \mid G) \, P(G)$.
4. **Conclusion:** Strong read support for $AG$ overcomes its lower prior probability ($0.001998$), resulting in a high posterior probability for $AG$.

---

## 7. Bayesian Genotyping vs. Maximum Likelihood Genotyping

| Feature | Maximum Likelihood (ML) | Bayesian (MAP) |
| :--- | :--- | :--- |
| **Uses sequencing reads** | Yes | Yes |
| **Uses genotype likelihoods** | Yes | Yes |
| **Uses prior information** | No | Yes |
| **Main result** | Most likely genotype | Posterior probability distribution |
| **Incorporates population AF** | No | Yes |
| **Quantifies posterior uncertainty** | No | Yes |
| **Objective formulation** | $G_{\mathrm{ML}} = \mathrm{argmax}_G P(D \mid G)$ | $G_{\mathrm{MAP}} = \mathrm{argmax}_G P(D \mid G)P(G)$ |

When all genotypes have equal prior probabilities, ML and MAP selections are mathematically identical.

---

## 8. Practical Considerations in Real Variant Calling

* **Population AF Priors:** Useful for low-coverage data, but must not override clear, high-quality variant calls in rare-disease contexts.
* **Mapping Quality:** Incorrectly aligned reads can mimic true variants; quality scores must adjust for alignment confidence.
* **Strand Bias:** Variants observed on only one strand often indicate systematic technology artifacts.
* **Allele Imbalance:** Real heterozygous sites can deviate from 0.5 AF due to copy number alterations, reference bias, or sample contamination.
* **Ploidy and Context:** Specialized models are required for non-diploid regions (e.g., mitochondria, male X/Y chromosomes, tumor samples).

---

## 9. Applications in Bioinformatics

* **Rare Disease Analysis:** Disambiguating true heterozygous/homozygous variants from noise in clinical panels and exomes.
* **Low-Coverage Sequencing:** Preserving statistical uncertainty instead of assigning premature hard genotype calls.
* **Population Genetics:** Directly utilizing genotype likelihoods for estimating allele frequency spectra.
* **Somatic & Tumor Analysis:** Adapting continuous mixture models for subclonal mutations and tumor purity estimation.

---

## 10. Computational Complexity

Let $P$ = positions evaluated, $R_p$ = read depth at position $p$, $G$ = candidate genotypes ($G=10$ for diploid), $B = 4$ bases.

* **Genotype Likelihood Time Complexity:** $O(R_p \times G)$ per site. Across the genome: $O\left(G \sum_{p=1}^{P} R_p\right)$.
* **Four-Base Mixture Time Complexity:** $O(R_p \times B)$ per site.
* **Bayesian Posterior Calculation:** $O(G)$ per site.
* **Space Complexity:** Streaming implementations require $O(R_p + G)$ memory per position.

---

## 11. Key Takeaways

* **Bayesian genotyping** combines sequencing evidence with prior genetic information.
* **Genotype likelihoods** quantify data probability given a genotype, while **priors** incorporate biological expectations.
* **Posterior probabilities** represent the updated confidence in each genotype.
* Unconstrained allele frequency MLE estimates that yield negative values are clipped to zero.
* Genotype calling is a rigorous statistical inference problem, not a simple base-counting exercise.

---

## 12. References and Further Reading

1. **Li H. (2011).** A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics*.  
   <https://doi.org/10.1093/bioinformatics/btr509>

2. **Li H. (2011).** Improving SNP discovery by base alignment quality. *Bioinformatics*.  
   <https://doi.org/10.1093/bioinformatics/btr076>

3. **DePristo et al. (2011).** A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics*.  
   <https://doi.org/10.1038/ng.806>

4. **McKenna et al. (2010).** The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*.  
   <https://doi.org/10.1101/gr.107524.110>

5. **BCFtools Documentation.** Pileup generation, genotype likelihoods, and variant calling.  
   <https://samtools.github.io/bcftools/bcftools.html>

6. **GATK Documentation.** Genotype likelihoods, Bayesian genotyping concepts, and variant calling.  
   <https://gatk.broadinstitute.org/>