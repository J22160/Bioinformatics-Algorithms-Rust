# Complexity Analysis

## What Is This?

Imagine you are a lab technician and you have two machines that can sequence DNA.
Machine A takes 1 minute per sample. Machine B takes 1 minute for the first sample,
but doubles the time for each additional sample (1 min, 2 min, 4 min, 8 min...).

For 10 samples: Machine A = 10 minutes, Machine B = 1023 minutes (~17 hours).

**Complexity analysis** is how computer scientists describe this difference formally.
It answers: "As my data grows, how much slower does my program get?"

## The Key Notation: Big-O

| Notation | Name | Example | 10 items | 1,000 items | 1,000,000 items |
|----------|------|---------|----------|-------------|-----------------|
| O(1) | Constant | Array lookup | instant | instant | instant |
| O(log n) | Logarithmic | Binary search | 3 steps | 10 steps | 20 steps |
| O(n) | Linear | Read a FASTQ file | 10 ops | 1,000 ops | 1,000,000 ops |
| O(n log n) | Linearithmic | Merge sort | 33 ops | 10,000 ops | 20,000,000 ops |
| O(n²) | Quadratic | Naive alignment | 100 ops | 1,000,000 ops | 10¹² ops |
| O(2ⁿ) | Exponential | Brute force motif | 1,024 ops | impossible | impossible |

## Problem Statement

**Given**: A list of gene expression values from an RNA-seq experiment.
**Task**: Demonstrate that different approaches to finding the maximum value
have different time complexities, and measure them empirically.

Specifically:
1. Find the maximum value by scanning the list once → O(n)
2. Sort the list and take the last element → O(n log n)
3. Compare every pair to find the max → O(n²) (deliberately wasteful, for teaching)

## How to Think About It

- O(n) means: if you double the data, the time doubles. This is like reading
  every page of a book — twice as many pages means twice the reading time.
- O(n²) means: if you double the data, the time quadruples. This is like
  comparing every student in a class with every other student — doubling the
  class size means four times as many comparisons.

## Why This Matters in Bioinformatics

When you align 30 million reads against a 3 billion base genome, the choice between
O(n²) and O(n log n) is the difference between finishing in minutes vs. never finishing.

## Files

- `complexity_demo.rs` — Rust implementation with timing benchmarks
- `complexity_demo.py` — Python implementation for comparison