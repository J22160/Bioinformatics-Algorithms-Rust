# Quick Sort

## What Is This?

Imagine you are a lab manager sorting tubes of samples by concentration.
You pick one tube at random and say: "All tubes with lower concentration
go on my left, all tubes with higher concentration go on my right."
Now that one tube is in its correct final position. You then repeat
the same process for the left pile and the right pile independently.

That is **Quick Sort** — a divide-and-conquer algorithm that
picks a "pivot" element, partitions the remaining elements around it,
and recursively sorts each partition.

## How It Works (Step by Step)

Given: `[38, 12, 27, 8, 41, 19]`

1. Choose pivot = `19` (last element)
2. Partition: elements < 19 go left, elements > 19 go right
   - Left: `[12, 8]`, Pivot: `19`, Right: `[38, 27, 41]`
3. Recursively sort left `[12, 8]` → `[8, 12]`
4. Recursively sort right `[38, 27, 41]` → `[27, 38, 41]`
5. Combined: `[8, 12, 19, 27, 38, 41]` ✓

## Complexity

| Case | Time | When? |
|------|------|-------|
| Best | O(n log n) | Pivot splits array roughly in half |
| Average | O(n log n) | Random data |
| Worst | O(n²) | Already sorted + bad pivot choice |

- **Space**: O(log n) average (recursion stack), O(n) worst case

## Bioinformatics Problem

**Scenario**: You are processing a VCF file from a whole-genome sequencing
experiment. You have variant quality scores (QUAL column) for 50,000 variants
and need to sort them to apply a quality threshold for filtering.

Given QUAL scores: `[934, 52, 1847, 300, 78, 2100, 415, 61]`

Sort them using quick sort to produce: `[52, 61, 78, 300, 415, 934, 1847, 2100]`

Then apply a QUAL threshold of 100 to keep only high-confidence variants.

## Why Quick Sort Matters

Quick Sort is the basis of most standard library sorting functions.
C++ `std::sort` uses a hybrid of Quick Sort (Introsort). When you sort
a VCF file by position, or sort read alignments by mapping quality,
Quick Sort is typically the engine doing the work.

Its average O(n log n) performance and in-place operation (no extra
memory allocation) make it ideal for large genomic datasets.

## Randomized Quick Sort

Choosing the pivot randomly avoids the O(n²) worst case on sorted data.
As Jones and Pevzner explain, if you choose m uniformly at random,
there is at least a 50% chance it will be a "good" splitter that
divides the array into reasonably balanced halves.
