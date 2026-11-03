# Variant Calling Algorithms

Variant calling is the process of identifying positions where a
sequenced genome differs from a reference genome. Every time
you run GATK HaplotypeCaller, bcftools mpileup, or DeepVariant,
the algorithms in this section are at work.

## Why Variant Calling Matters

Identifying genetic variants (single nucleotide polymorphisms (SNPs),
insertions, deletions (indels), and structural variants) is the
foundation of:
- **Clinical genomics**: Diagnosing genetic diseases
- **Cancer genomics**: Finding somatic mutations driving tumors
- **Population genetics**: Understanding human diversity
- **Pharmacogenomics**: Tailoring drug treatments to genotypes

## The Variant Calling Pipeline

```
Aligned Reads (BAM)
|
Pileup Construction <-- pileup-genotyping/
|
Genotype Likelihoods <-- bayesian-genotyping/
| pair-hmm/
Genotype Calling
|
Quality Filtering <-- vcf-filtering/
|
Base Quality Recalibration <-- quality-recalibration/
|
Final VCF
```

## Algorithms

| Algorithm | Purpose | Used By |
|-----------|---------|---------|
| Pileup Genotyping | Column-wise base counting + likelihood | bcftools mpileup |
| Bayesian Genotyping | Posterior P(genotype / data) | GATK, bcftools |
| Pair HMM | P(read / haplotype) with gaps + errors | GATK HaplotypeCaller |
| VCF Filtering | Quality-based variant filtration | GATK VQSR, bcftools |
| Base Quality Recalibration | Correcting systematic quality errors | GATK BQSR |

## Key Insight: Variant Calling Is Bayesian Inference

At its core, variant calling asks: given observed reads at a genomic
position, what is the most likely genotype? This is a direct
application of Bayes' theorem:

    P(genotype | reads) ∝ P(reads | genotype) × P(genotype)

The likelihood P(reads | genotype) is computed from base qualities
and alignment scores. The prior P(genotype) comes from population
genetics (e.g., expected heterozygosity).