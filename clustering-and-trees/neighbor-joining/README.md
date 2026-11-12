# Neighbor Joining

## What Is This?

Neighbor Joining (NJ) builds an unrooted phylogenetic tree
without assuming a molecular clock. It was developed by Saitou
and Nei in 1987 and finds pairs of taxa that are simultaneously
close to each other AND far from all other taxa [1][5].

## How It Differs from UPGMA

UPGMA assumes equal evolutionary rates (molecular clock).
NJ does not, making it more appropriate for real biological
data where lineages often evolve at different rates [1].

## The Key Insight

Instead of simply merging the closest pair, NJ computes a
corrected distance that accounts for each taxon's average
distance to all others. For each cluster C, define [1]:

    u(C) = [1/(L-2)] × Σ D(C, C')

Then merge the pair minimizing: D(C1, C2) - u(C1) - u(C2)

This ensures we select neighbors that are both close to each
other AND separated from the rest of the tree.

## Complexity
- **Time**: O(n³)
- **Space**: O(n²)

## Bioinformatics Problem

**Scenario**: You have protein sequence distances between
6 mammalian species and want to reconstruct their phylogeny
without assuming a molecular clock.
