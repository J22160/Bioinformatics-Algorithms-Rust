# Small Parsimony — Fitch's & Sankoff's Algorithms

## What Is This?

Given a phylogenetic tree with character labels on the leaves,
the Small Parsimony problem asks: what labels should the
internal nodes have to minimize the total number of character
changes across all edges? [1]

## Fitch's Algorithm (Unweighted)

A bottom-up dynamic programming approach [1]:
1. For each leaf, set S(v) = {character at v}
2. For each internal node v with children u and w:
   - If S(u) ∩ S(w) is non-empty: S(v) = S(u) ∩ S(w)
   - Otherwise: S(v) = S(u) ∪ S(w) and add 1 to the score
3. The parsimony score = number of union operations

## Sankoff's Algorithm (Weighted)

A more general version that allows different costs for
different character substitutions, using a scoring matrix δ [1]:

    s_t(v) = min_i{s_i(u) + δ(i,t)} + min_j{s_j(w) + δ(j,t)}

## Bioinformatics Problem

**Scenario**: You have a known phylogeny for 8 primate species
and an aligned nucleotide position. What was the ancestral
nucleotide at each internal node?
