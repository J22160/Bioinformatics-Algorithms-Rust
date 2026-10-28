# Combinatorial Pattern Matching

This section covers the data structures and algorithms that make modern
genomic read mapping possible. Every time you run BWA, Bowtie2, or BLAST,
these algorithms are working behind the scenes.

## Why Pattern Matching Matters

A human genome is ~3 billion bases. A typical sequencing run produces
~30 million reads of 150 bases each. Finding where each read belongs
in the genome is like finding a specific sentence in a library of 1,000 books and doing it 30 million times.

## Algorithms

| Algorithm | Purpose | Time | Used By |
|-----------|---------|------|---------|
| Suffix Array | Sort all suffixes for fast lookup | O(n) build, O(m log n) search | Foundation for BWT/FM |
| BWT | Reversible text transform enabling compression + search | O(n) | BWA, Bowtie2, samtools |
| FM-Index | Compressed full-text index built on BWT | O(m) search | BWA, Bowtie2 |
| BLAST Heuristic | Seed-and-extend for fast approximate matching | Heuristic | BLAST, DIAMOND |

## The Pipeline

```
Suffix Array → BWT → FM-Index → Read Mapping
↑
BLAST uses a different approach: Seed → Extend → Score
```

These algorithms form the backbone of tools used in WGS, WES,
RNA-seq, snRNA-seq, and virtually every NGS workflow.
