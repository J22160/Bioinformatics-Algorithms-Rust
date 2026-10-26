# Longest Common Subsequence (LCS)

## What Is This?

You have two DNA sequences from different species and want to know how
much they share in common — not as a perfect match, but as a set of
characters that appear in the same order in both sequences. Those shared
characters in order are a **common subsequence**.

The **Longest Common Subsequence** is the longest such set.

## How It Works

Unlike a **substring** (which must be contiguous), a **subsequence** can
skip characters as long as the order is preserved.

Example:
- Sequence A: `A C G T A C`
- Sequence B: `A G C A T C`
- LCS: `A G A C` (length 4)

We build an (n+1) × (m+1) grid. Each cell (i,j) stores the length
of the LCS between the first i characters of A and first j characters of B.

**Rule for each cell**:
- If A[i] == B[j]: cell = diagonal + 1
- Otherwise: cell = max(cell above, cell to left)

## Complexity

- **Time**: O(n × m)
- **Space**: O(n × m) for full matrix, O(min(n,m)) if only length needed

## Bioinformatics Problem

**Scenario**: You are comparing conserved gene order (synteny) between
two bacterial genomes. Genome A has genes in order `[dnaA, gyrB, rpoB, recA, groEL]`
and Genome B has `[dnaA, recA, gyrB, groEL, rpoB]`.
What is the longest sequence of genes that appears in the same order in both?

LCS: `[dnaA, gyrB, groEL]` or `[dnaA, rpoB, groEL]` (length 3)

This tells you that at least 3 genes preserved their relative order
during evolution.

## Connection to Alignment

The LCS problem is a special case of global sequence alignment
where matches score +1 and indels score 0. The edit distance
between two sequences relates to LCS by: d(v,w) = n + m − 2·LCS(v,w)
when only insertions and deletions are allowed.
