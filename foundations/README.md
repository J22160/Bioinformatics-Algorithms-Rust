# Foundations

This section covers the fundamental algorithm design paradigms that underpin
every bioinformatics tool you will ever use.

Think of these as the grammar of computational biology — you need to understand
these patterns before you can read or write any bioinformatics "sentence."

## What You Will Learn

| Section | Key Idea | Bioinformatics Connection |
|---------|----------|--------------------------|
| **Complexity Analysis** | How to measure algorithm speed | Choosing BLAST vs. Smith-Waterman |
| **Sorting** | Ordering data efficiently | Sorting BAM files by coordinate |
| **Dynamic Programming** | Building solutions from subproblems | Sequence alignment (Needleman-Wunsch) |
| **Greedy Basics** | Making locally optimal choices | Exon chaining in gene prediction |

## Why These Matter for Biology

Every time you run `samtools sort`, it uses merge sort.
Every time you run BLAST, it uses heuristic pattern matching built on these foundations.
Every time GATK calls a variant, dynamic programming scores read-haplotype likelihoods.

Understanding these basics will make you a better bioinformatician, even if you
never implement them from scratch.