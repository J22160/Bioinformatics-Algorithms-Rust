# BLAST Heuristic: Seed-and-Extend Algorithm


## What Is BLAST?

**BLAST (Basic Local Alignment Search Tool)** is a sequence similarity search algorithm used to identify regions of similarity between a query sequence and sequences in a biological database.

It is widely used in bioinformatics to identify homologous genes, predict protein functions, find conserved sequence regions, and investigate the evolutionary relationships between organisms.

The central challenge is computational efficiency. Given a query sequence and a large database, comparing the query against every database sequence using an optimal local alignment algorithm such as Smith–Waterman can be expensive.

BLAST addresses this problem using a **heuristic seed-and-extend strategy**:

1. Identify short sequence matches called **seeds**.
2. Extend promising seeds into longer regions of similarity.
3. Perform more computationally expensive gapped alignment on promising candidates.
4. Evaluate the statistical significance of the resulting alignments.

Rather than exhaustively evaluating every possible local alignment, BLAST concentrates computational effort on regions that are likely to produce meaningful matches.

> **Key idea:** Find promising short matches first, then spend computational resources extending only the most promising candidates.

BLAST is a heuristic algorithm. It is designed to find biologically meaningful similarities efficiently, but it does not guarantee that every possible optimal local alignment will be discovered.

---

## The Biological Problem

Suppose you have sequenced a previously uncharacterized gene and want to determine whether similar sequences already exist in a database containing 100,000 known genes.

Your input might be:

```text
Query:
ATGCTAGCTACGATCGTACGATCGATGCTAGC

```

A database sequence might contain a similar region:

```text
Database:
GGGATGCTAGCTACGATCGTTCGATGCTAGCTTT

```

The sequences are not identical, but they share a long region of similarity.

An exhaustive local alignment algorithm could examine many possible alignments between the query and every database sequence. BLAST instead searches for short matching words that suggest where a promising alignment may occur.

Once a seed is found, the algorithm investigates its surrounding sequence to determine whether the local similarity is strong enough to report.

This strategy is especially effective when the database is large and only a small fraction of its sequences are expected to match the query.

---

## How the BLAST Seed-and-Extend Algorithm Works

The general workflow consists of five stages.

### Step 1: Query processing and word generation

BLAST divides the query sequence into overlapping words of a specified length, commonly called the **word size**.

For example, consider this DNA query:

```text
Query: ATGCTAGC
Word size: W = 4

```

The resulting overlapping words are:

| Position | Word |
| --- | --- |
| 0 | ATGC |
| 1 | TGCT |
| 2 | GCTA |
| 3 | CTAG |
| 4 | TAGC |

For a query of length $m$ and word size $W$, the number of overlapping words is:

$$m - W + 1$$

assuming $m \ge W$.

For protein searches, BLAST commonly uses a word size of 3 with standard settings. Traditional nucleotide BLAST commonly uses a word size of 11, although the appropriate value depends on the BLAST program and search settings. For example, megablast uses a larger default word size.

> **Why use words?** Comparing complete sequences against every database position is expensive. Short words provide a fast way to identify candidate regions that may contain longer alignments.

---

### Step 2: Seed finding

BLAST searches for database locations that match the query words.

For exact matching, a seed is a short substring shared by the query and a database sequence:

```text
Query:      ATGCTAGC
Database:   GGATGCTAGCTT
              ||||||||
              ATGCTAGC

```

The shared region indicates a possible alignment.

However, protein sequences frequently contain substitutions even when the underlying proteins are evolutionarily related. Requiring every seed to match exactly could miss meaningful similarities.

BLAST therefore uses **scoring-based word neighborhoods** for many search modes. For a protein query word, it generates or considers similar words that score above a specified threshold under a substitution matrix such as BLOSUM62.

For example, a query word `MKT` may be associated with other words that achieve a sufficiently high substitution score, even if they are not identical.

The precise seeding strategy varies among BLAST programs and versions. Some nucleotide search modes emphasize exact matches, while protein BLAST uses neighborhood words. Modern BLAST implementations also use specialized strategies to improve sensitivity and speed.

> **Implementation concept:** A common implementation stores searchable words in a lookup structure, such as a hash table, or uses another indexed lookup strategy. The database is scanned or indexed to identify candidate seed locations efficiently, resulting in a collection of candidate matches called **hits**.

---

### Step 3: Ungapped seed extension

Finding a short seed does not necessarily mean that the surrounding region is biologically meaningful.

BLAST extends promising seed hits in both directions, initially without introducing gaps. This produces an ungapped local alignment.

Consider:

```text
Query:      ATGCTAGCTACG
Database:   ATGCTAGATACG
            ||||||| ||||

```

The sequences share a strong matching region but contain a mismatch.

During extension, BLAST evaluates the alignment using a scoring system that rewards matches and penalizes mismatches.

A simplified scoring scheme might be:

* **Match:** +2
* **Mismatch:** -3

For protein searches, substitution matrices assign different scores to different amino-acid substitutions.

BLAST uses an **X-drop-style stopping criterion** in extension procedures: extension stops when the score deteriorates sufficiently relative to the best score observed during the extension.

The exact extension thresholds and implementation depend on the BLAST program, preventing the system from extending every seed across an entire database sequence.

---

### Step 4: Gapped alignment

Ungapped extension is useful for identifying promising regions, but biological sequences frequently contain insertions and deletions, which create gaps in an alignment:

```text
Query:      ATGCTAGCTACG
Database:   ATGCT--CTACG

```

A gapped alignment can represent this relationship more accurately than an alignment that permits only matches and mismatches.

BLAST performs more computationally expensive gapped extension on promising candidate regions using dynamic-programming methods related to Smith–Waterman local alignment.

A typical alignment scoring system can be expressed as:

$$S = \sum_i s(q_i,d_i) - \sum_j g_j$$

where:

* $q_i$ is a query residue or nucleotide.
* $d_i$ is the aligned database residue or nucleotide.
* $s(q_i,d_i)$ is the substitution score.
* $g_j$ is the gap penalty for gap $j$.

With affine gap penalties, a gap of length $L$ is commonly penalized using a function of the form:

$$g(L) = g_{\mathrm{open}} + (L-1)g_{\mathrm{extend}}$$

The gap-opening penalty discourages introducing new gaps, while the gap-extension penalty controls the cost of lengthening existing gaps.

The resulting local alignment is often reported as a **high-scoring segment pair (HSP)**—a locally aligned region whose score is sufficiently high to merit further consideration.

Importantly, BLAST does not simply run unrestricted Smith–Waterman alignment across every position in every database sequence; its filtering and extension stages restrict the regions subjected to more expensive calculations.

---

### Step 5: Statistical evaluation

A high alignment score does not automatically mean that two sequences are biologically related. Short sequences, low-complexity regions, and repetitive motifs can produce matches by chance.

BLAST therefore estimates the statistical significance of alignments using methods based on Karlin–Altschul statistics for applicable search settings.

Two important measures are:

* **Bit score:** Normalizes the raw alignment score so that scores can be compared across searches using corresponding statistical parameters. A higher bit score indicates a stronger alignment.
* **E-value (expect value):** Estimates the number of alignments with a score at least as high as the observed score that would be expected by chance in a database search of the specified size, under the statistical model.

A commonly used approximation is:

$$E \approx Kmn e^{-\lambda S}$$

where:

* $E$ is the expected number of chance matches.
* $m$ and $n$ represent effective query and database lengths in the simplified model.
* $S$ is the raw alignment score.
* $K$ and $\lambda$ are statistical parameters determined by the scoring system and sequence background.

| E-value | General interpretation |
| --- | --- |
| **10** | Matches of this strength may occur relatively often by chance. |
| **0.01** | Approximately 0.01 chance-model matches are expected at this score or better under the model. |
| **$1\text{e}-10$** | A highly significant match under the model. |
| **$1\text{e}-50$** | An exceptionally strong statistical match under the model. |

> **Note:** An E-value is not the probability that a match is false, nor does it directly measure the probability that two sequences are homologous. Biological interpretation also depends on sequence coverage, alignment length, identity, composition, and the research question.

---

## The Filtration Principle and the Pigeonhole Principle

The mathematical intuition behind seed-based filtering is that a sufficiently long, similar alignment is likely to contain at least one exact matching segment.

Suppose two aligned sequences each contain $n$ positions and differ at no more than $k$ positions, assuming the alignment contains only substitutions with no insertions or deletions.

Partition the sequence into $k+1$ nonempty blocks of approximately equal length. Because there are at most $k$ mismatches, at most $k$ blocks can contain a mismatch if the blocks are disjoint. Therefore, at least one block must contain no mismatches—meaning that block is an exact match.

The approximate length of such a guaranteed matching block is:

$$\ell \approx \left\lfloor \frac{n}{k+1} \right\rfloor$$

This is the **pigeonhole principle** applied to sequence matching.

### Worked example

Consider two aligned sequences of length 12:

```text
Sequence A: ATGCTAGCTACG
Sequence B: ATGTTAGCTTCG

```

The sequences differ at two positions. If we partition the alignment into three blocks of four positions:

```text
Sequence A: ATGC | TAGC | TACG
Sequence B: ATGT | TAGC | TTCG

```

The first and third blocks contain mismatches, but the middle block matches exactly (`TAGC`) and can serve as a seed.

> **Important limitation:** This argument applies to aligned strings with a bounded number of mismatches and no indels within the alignment. It does **not** mean that ordinary BLAST always finds an exact seed of this length. Insertions, deletions, sequence rearrangements, low-complexity regions, and scoring-based seed rules complicate the guarantee.

---

## Worked Example: From Query to Database Hit

Suppose a researcher wants to search a novel DNA fragment against a database of known genes.

* **Query:** `ATGCTAGCTACGATCG`
* **Database Sequence A:** `GGGATGCTAGCTACGATCGTTT`
* **Database Sequence B:** `CCCTTTAGGCAATGCC`

Assume, for illustration, that the seed-finding step identifies the exact word `ATGC` in sequence A.

### Step 1: Find the seed

```text
Query:      ATGCTAGCTACGATCG
            ||||
Seed:       ATGC

Database A: GGGATGCTAGCTACGATCGTTT
                 ATGC

```

The seed identifies a promising region in database sequence A.

### Step 2: Extend the match

The algorithm examines the neighboring residues or nucleotides and discovers a long matching region:

```text
Query:      ATGCTAGCTACGATCG
            ||||||||||||||||
Database A: ATGCTAGCTACGATCG

```

The surrounding database sequence contains the complete query fragment, so this is a strong candidate alignment.

### Step 3: Evaluate other candidates

Sequence B does not contain the same strong matching region and is eliminated during seed finding or extension.

### Step 4: Report significant matches

BLAST reports alignments that meet the search's reporting criteria, including the configured E-value threshold.

---

## Why BLAST Is Faster Than Exhaustive Alignment

An optimal local alignment algorithm such as Smith–Waterman uses dynamic programming to find the best local alignment between two sequences. For sequences of lengths $m$ and $n$, its conventional time complexity is:

$$O(mn)$$

If a query of length 500 is compared independently against 100,000 database sequences with a total database length $N$, a simplified exhaustive comparison has a cost proportional to $O(mN)$.

BLAST reduces the amount of expensive alignment work by filtering candidates using seeds and extending only promising hits.

### Complexity Overview

| Operation | General computational consideration |
| --- | --- |
| **Query word generation** | $O(m)$ for fixed word size |
| **Seed lookup** | Depends on the lookup structure and number of candidate hits |
| **Seed extension** | Depends on the number and extent of extensions |
| **Gapped alignment** | Depends on candidate-region sizes and the alignment implementation |
| **Statistical evaluation** | Depends on the number of candidate alignments and associated calculations |

For a large database, BLAST can be dramatically faster than exhaustive alignment because most database regions never undergo expensive gapped extension. This is an empirical performance advantage, not a guarantee for every dataset.

---

## Sensitivity Versus Speed

BLAST parameters influence the trade-off between search speed and the ability to detect weak or distant similarities:

* **Word size:** A larger word size generally reduces candidate seeds, increasing speed but potentially missing weak/divergent similarities. A smaller word size improves sensitivity but increases candidate hits and runtime.
* **Seed threshold:** For protein searches, the neighborhood-word threshold controls which related words qualify for seed lookup.
* **Extension threshold:** Determines how much computation is spent following a seed.
* **E-value threshold:** Controls which results are reported.

---

## Common BLAST Programs

| Program | Query | Database | Typical purpose |
| --- | --- | --- | --- |
| `blastn` | Nucleotide | Nucleotide | Find similar DNA or RNA sequences |
| `blastp` | Protein | Protein | Find similar proteins |
| `blastx` | Translated nucleotide | Protein | Search translated query reading frames against a protein database |
| `tblastn` | Protein | Translated nucleotide | Find nucleotide regions that may encode a similar protein |
| `tblastx` | Translated nucleotide | Translated nucleotide | Compare translated reading frames of nucleotide sequences |

---

## BLAST and Other Seed-and-Extend Tools

| Tool | Main use | General strategy |
| --- | --- | --- |
| **BLAST** | General nucleotide and protein similarity searches | Word-based seed finding followed by extension |
| **BWA-MEM** | Short-read alignment to a reference genome | FM-index-based seed discovery and alignment extension |
| **minimap2** | Long-read mapping and general sequence alignment | Minimizer-based seeding followed by chaining and alignment |
| **Bowtie 2** | Short-read alignment | FM-index-based search with alignment and mismatch/gap handling |

---

## Applications in Bioinformatics

* **Gene and protein identification:** Compare novel genes against known sequences.
* **Genome annotation:** Identify candidate genes and conserved regions.
* **Pathogen identification:** Characterize microbial sequences.
* **Evolutionary analysis:** Suggest evolutionary relationships and identify conserved domains.
* **Variant interpretation:** Investigate sequence variants against reference sequences.

---

## Limitations and Caveats

* **Heuristic search:** Some alignments may be missed if they fail to generate qualifying seeds.
* **Low-complexity regions:** Repetitive or compositionally biased sequences can produce misleading matches.
* **Database dependence:** Missing or poorly annotated reference sequences limit conclusions.
* **Statistical vs. biological significance:** A low E-value supports statistical significance but does not automatically establish biological function.
* **Parameter dependence:** Word sizes, matrices, gap penalties, and thresholds heavily influence output.

---

## Why BLAST Matters

BLAST made large-scale sequence similarity searching practical by combining fast seed detection, selective alignment extension, and statistical evaluation. Its core lesson remains fundamental in algorithm design:

> **Avoid solving an expensive problem everywhere when a fast filtering stage can identify the regions worth investigating.**