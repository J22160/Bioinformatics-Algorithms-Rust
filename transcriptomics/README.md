# Transcriptomics Algorithms

This section covers the algorithms that power the entire RNA-seq
analysis pipeline — from raw read quantification to identifying
differentially expressed genes between conditions.

## Why These Algorithms Matter

Every time you run Salmon, kallisto, or DESeq2, these algorithms
are at work. Understanding them helps you interpret your results
correctly, choose the right tool, and troubleshoot when things
go wrong.

## The RNA-seq Analysis Pipeline

```
Raw Reads
↓
Pseudoalignment / Mapping ← pseudoalignment/
↓
Transcript Quantification ← em-quantification/
↓
Normalization ← normalization/
↓
Differential Expression ← differential-expression/
↓
Pathway Analysis ← gene-set-enrichment/
```


## Algorithms in This Section

| Algorithm | Purpose | Used By |
|-----------|---------|---------|
| EM for Quantification | Resolve multi-mapped reads | Salmon, RSEM, kallisto |
| Pseudoalignment | Fast read-transcript matching | kallisto, Salmon |
| Negative Binomial GLM | Model count data | DESeq2, edgeR |
| Normalization (TMM, MoR) | Remove technical bias | DESeq2, edgeR |
| Gene Set Enrichment | Pathway-level analysis | GSEA, GSEApy |

## Biological Context

- **WGS/WES**: Variant calling pipelines use similar statistical
  models for read assignment
- **snRNASeq**: Cell Ranger quantifies with EM; Seurat/Scanpy
  use similar normalization and DE methods
