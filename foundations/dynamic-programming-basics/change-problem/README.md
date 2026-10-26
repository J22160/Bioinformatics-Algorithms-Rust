# The Change Problem (Dynamic Programming)

## What Is This?

Imagine you work at a DNA synthesis company. You need to build an
oligonucleotide of exactly 77 bases, and you have pre-made fragments
of lengths 1, 3, and 7 bases. What is the minimum number of fragments
you need to concatenate to reach exactly 77 bases?

This is the **Change Problem** — given a target amount and a set of
denominations (fragment lengths), find the minimum number of "coins"
(fragments) needed to make the exact amount.

## How It Works

The key insight is: the optimal solution for amount M depends on
optimal solutions for smaller amounts. If you have denominations
[1, 3, 7], then:

```
bestCoins(77) = 1 + min(bestCoins(76), bestCoins(74), bestCoins(70))
```


Instead of computing this recursively (which wastes enormous time
recalculating the same values), we build a table from 0 up to M,
filling in each entry exactly once.

### Example: Amount = 11, Coins = [1, 3, 5]

| Amount | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|--------|---|---|---|---|---|---|---|---|---|---|----|----|
| Min coins | 0 | 1 | 2 | 1 | 2 | 1 | 2 | 3 | 2 | 3 | 2 | 3 |

Answer: 3 coins (5 + 5 + 1, or 5 + 3 + 3)

## Complexity

- **Time**: O(M × d) where M is the target amount and d is the number of denominations
- **Space**: O(M) for the DP table

## Bioinformatics Problem

**Scenario**: You are designing a tiling path across a BAC clone library.
You have BAC inserts of sizes 50 kb, 120 kb, and 200 kb. What is the
minimum number of BAC clones needed to cover exactly 1000 kb of a
genomic region?

Denominations: [50, 120, 200]
Target: 1000

## Why This Matters

The Change Problem is the simplest dynamic programming example. The same
table-filling strategy powers Needleman-Wunsch sequence alignment,
the Viterbi HMM decoder, and exon chaining in gene prediction.
Understanding this problem is the gateway to understanding all of them.