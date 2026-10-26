# Merge Sort

## What Is This?

Imagine you are organizing patient records. Instead of scanning through all of them
repeatedly (like selection sort), you try a different strategy:

1. **Split** the pile in half
2. **Sort each half** (by splitting again, recursively, until you have piles of 1)
3. **Merge** the two sorted halves by comparing the top card of each pile

This is **merge sort** — a divide-and-conquer algorithm invented by John von Neumann in 1945.

## How It Works

Merge sort follows three steps: **Divide, Sort, and Merge.**

1. **Divide:** Split the array into two halves. Continue splitting each half recursively until every subarray contains just one element. A single-element array is already sorted.

2. **Sort and Merge:** Merge adjacent subarrays by comparing their first elements. Move the smaller element into a new sorted array, then advance to the next element in that subarray. Repeat until one subarray is exhausted, then append the remaining elements from the other.

3. **Repeat:** Continue merging the sorted subarrays into progressively larger sorted arrays until a single array contains all the elements in ascending order.

**Example:** Sorting genomic coordinates:

```text
                 [155799, 5529, 87160, 24601]
                            |
                    Split into halves
                      /           \
            [155799, 5529]     [87160, 24601]
                 /   \              /   \
             [155799] [5529]    [87160] [24601]
                 \   /              \   /
              [5529, 155799]    [24601, 87160]
                      \           /
                     Merge sorted halves
                            |
             [5529, 24601, 87160, 155799]
```

During merging, the algorithm compares the smallest remaining element in each half and selects the smaller one. This guarantees that each merged subarray remains sorted.


## Complexity

- **Time**: O(n log n) — always, regardless of input order
- **Space**: O(n) — needs temporary storage for merging

## Bioinformatics Problem

**Scenario**: You have chromosome coordinates from a BED file representing
transcription factor binding sites identified by ChIP-seq.
You need to sort them by genomic position so they can be efficiently
intersected with gene annotations.

Given binding sites (as start positions on chr7):
`[155_799, 5_529, 87_160, 24_601, 117_482, 48_930, 72_215]`

Sort them to produce:
`[5_529, 24_601, 48_930, 72_215, 87_160, 117_482, 155_799]`

## Why Merge Sort Matters

`samtools sort` uses a variant of merge sort. When you sort a BAM file
by coordinate (required before variant calling with GATK), merge sort is
what makes it possible to handle billions of reads efficiently.

Merge sort guarantees O(n log n) even in the worst case, unlike quicksort
which can degrade to O(n²).