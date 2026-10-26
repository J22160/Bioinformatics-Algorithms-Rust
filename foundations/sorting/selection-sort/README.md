# Selection Sort

## What Is This?

Imagine you have a stack of index cards, each with a patient's blood pressure reading.
You want to arrange them from lowest to highest. 

**Your strategy:** scan through all cards,
find the smallest value, put it first. Then scan the remaining cards, find the next
smallest, put it second. Repeat until done.

That is **selection sort**. It is simple but slow.

## How It Works (Step by Step)

Given: `[64, 25, 12, 22, 11]`

1. Find the minimum in the whole list → `11`. Swap with position 0.
   Result: `[11, 25, 12, 22, 64]`
2. Find the minimum from position 1 onward → `12`. Swap with position 1.
   Result: `[11, 12, 25, 22, 64]`
3. Find the minimum from position 2 onward → `22`. Swap with position 2.
   Result: `[11, 12, 22, 25, 64]`
4. Find the minimum from position 3 onward → `25`. Already in place.
   Result: `[11, 12, 22, 25, 64]` ✓

## Complexity

- **Time**: O(n²) always — even if the list is already sorted
- **Space**: O(1) — sorts in place, no extra memory needed

## Bioinformatics Problem

**Scenario**: You have quality scores (Phred scores) for a set of sequencing reads.
You want to sort them to identify the lowest-quality reads for filtering.

Given Phred scores: `[38, 12, 27, 35, 8, 41, 19]`

Sort them using selection sort to produce: `[8, 12, 19, 27, 35, 38, 41]`

Then filter out reads below quality threshold 20: reads with scores `8, 12, 19`
would be discarded.

## Why Selection Sort Is Not Used in Practice

For 30 million reads, selection sort would need ~9 × 10¹⁴ comparisons.
Merge sort would need ~7.5 × 10⁸ about a million times fewer.

We include it here because it is the simplest sorting algorithm to understand,
and the concepts (scanning, swapping, in-place operations) appear everywhere.