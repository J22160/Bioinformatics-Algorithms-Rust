# Bioinformatics Algorithms in Rust

A comprehensive collection of algorithms commonly used in bioinformatics workflows,
implemented in **Rust** for performance and **Python** for accessibility, with
detailed explanations aimed at biologists learning programming.

## Who This Is For

- Biologists learning computational thinking for the first time
- Bioinformaticians wanting high-performance algorithm implementations
- Students studying transcriptomics, variant calling, WGS, WES, snRNASeq, or methylation analysis

## Why Rust?

Rust provides memory safety without a garbage collector, zero-cost abstractions,
and performance comparable to C/C++. Tools like `minimap2`, `noodles`, and `Sassy`
demonstrate that Rust is increasingly adopted in bioinformatics for its speed and safety.


## Repository Structure

```text
.
├── foundations/               # Core algorithm paradigms: sorting, DP, greedy, complexity
├── sequence-alignment/        # Needleman–Wunsch, Smith–Waterman, affine gaps
├── pattern-matching/          # BWT, FM-index, suffix arrays, BLAST heuristics
├── graph-algorithms/           # Eulerian paths, de Bruijn graphs
├── assembly/                  # De Bruijn graph assembly, overlap-layout-consensus (OLC)
├── read-mapping/              # BWA-MEM concepts, minimizers
├── variant-calling/           # Bayesian genotyping, Pair HMM
├── transcriptomics/            # EM quantification, pseudoalignment, DESeq2 models
├── single-cell/                # PCA, UMAP, Louvain clustering
├── methylation/                # Bisulfite alignment, differentially methylated regions (DMRs)
├── motif-finding/              # Gibbs sampling, random projections
├── hidden-markov-models/       # Viterbi, Baum–Welch, profile HMMs
├── clustering-and-trees/       # UPGMA, neighbor joining, k-means
├── genome-rearrangements/      # Sorting by reversals
├── data-structures/            # Bloom filters, interval trees
├── statistics/                 # Fisher's exact test, FDR correction
└── preprocessing/              # Adapter trimming, k-mer counting
```

## How to Build and Run

```bash
# Build all crates
cargo build --release

# Run a specific algorithm (example)
cargo run --release -p edit-distance

# Run tests
cargo test --workspace
```