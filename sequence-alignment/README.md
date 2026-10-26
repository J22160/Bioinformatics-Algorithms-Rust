# Sequence Alignment

Sequence alignment is the most fundamental operation in bioinformatics.
Every time you run BLAST, map reads with BWA, or call variants with GATK,
alignment algorithms are working behind the scenes.

## What Is Sequence Alignment?

Imagine you have two sentences:
- "THE QUICK BROWN FOX"
- "THE QUACK BROWN BOX"

Alignment finds the best way to line them up to see what changed:

```
THE QU-ICK BROWN FOX
THE QUACK- BROWN BOX
```

In biology, instead of English words, we align DNA sequences (A, C, G, T)
or protein sequences (20 amino acids) to find mutations, conserved regions,
and evolutionary relationships.

## Algorithms in This Section

| Algorithm | Time | Space | Use Case |
|-----------|------|-------|----------|
| Needleman-Wunsch | O(nm) | O(nm) | Global alignment of two sequences |
| Smith-Waterman | O(nm) | O(nm) | Local alignment (find best matching region) |
| Affine Gap | O(nm) | O(nm) | Realistic gap modeling (open vs extend) |
| Hirschberg | O(nm) | O(min(n,m)) | Global alignment with linear memory |
| Four Russians | O(n²/log n) | O(n²/log²n) | Subquadratic LCS/alignment |
| Multiple Alignment | O(n^k) | O(n^k) | Aligning k sequences simultaneously |

## Why These Matter for Your Workflows

- **WGS/WES**: BWA-MEM uses seed-and-extend with Smith-Waterman
- **Variant Calling**: GATK HaplotypeCaller uses a Pair HMM based on these DP ideas
- **Transcriptomics**: STAR aligner uses similar scoring for splice-aware alignment
- **snRNASeq**: Cell Ranger aligns reads using these foundational algorithms