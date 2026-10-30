# UMI Deduplication

## What Is This?

Before sequencing, the tiny amount of mRNA in a single cell is
amplified by PCR. This creates many copies of each molecule.
Without correction, a gene with 10 original molecules might
appear as 1000 reads after PCR — grossly inflating its count.

**Unique Molecular Identifiers (UMIs)** solve this. Each original
mRNA molecule gets tagged with a random short barcode (typically
10-12 bases) *before* PCR amplification. All PCR copies of the
same molecule share the same UMI. By counting unique UMIs per
gene rather than reads, we recover the true molecule count.

## The Problem

UMIs themselves can contain sequencing errors. If molecule A has
UMI `ACGTACGTAC` and a PCR copy gets sequenced as `ACGTACGTAT`
(1 error at position 10), naive counting would count these as
two separate molecules instead of one.

## How It Works (UMI-tools Directional Method)

The algorithm groups reads by (cell barcode, gene, alignment position),
then within each group clusters UMIs [7]:

1. **Build a graph**: Each unique UMI is a node. Connect two nodes
   if they differ by exactly 1 base (Hamming distance 1).

2. **Directional clustering**: The UMI with the most reads is the
   "seed" for the first cluster. Connect it to all UMIs within
   edit distance 1 that have *fewer* reads. Then connect those
   to their neighbors with still fewer reads. Repeat until no
   more UMIs can be added.

3. **Rationale**: PCR errors introduced in later cycles produce
   fewer copies. So if UMI-A has 100 reads and UMI-B has 2 reads,
   and they differ by 1 base, UMI-B is likely a PCR error of UMI-A.

4. **Count**: Each cluster = 1 original molecule.

## Complexity
- **Per group**: O(U² × L) where U = unique UMIs, L = UMI length
- **Total**: O(G × U² × L) for G gene-position groups

## Connection to Textbook

UMI deduplication is fundamentally a **graph clustering problem**
(Chapter 10) [1] combined with **approximate pattern matching**
(Chapter 9) [1]. Building the UMI adjacency graph uses the same
Hamming distance concept from the Motif Finding problem in
Chapter 4 [1]. The directional clustering resembles the
hierarchical clustering methods from Chapter 10 [1].
