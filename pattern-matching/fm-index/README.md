# FM-Index (Full-Text Index in Minute Space)

## What Is an FM-Index?

The **FM-index** is a compressed full-text index that supports efficient pattern searching while using substantially less memory than a conventional index that explicitly stores every suffix-array entry.

It was introduced by Paolo Ferragina and Giovanni Manzini and is based on the **Burrows–Wheeler Transform (BWT)**.

The central idea is to combine the BWT with character-counting and occurrence-rank operations so that a search can narrow down the possible locations of a pattern without scanning the entire text.

For example, suppose you want to find a 20-base DNA sequence in a human reference genome containing billions of bases. A naive approach could compare the pattern against many possible starting positions. An FM-index instead processes the pattern character by character, maintaining an interval of suffixes that match the pattern processed so far.

This makes FM-indexes especially useful in bioinformatics, where large reference genomes must be searched repeatedly for short sequencing reads.

### Why Is It Called a Compressed Index?

A conventional suffix array stores the starting position of every suffix of a text in lexicographically sorted order. Although suffix arrays enable fast searches, storing every position can require substantial memory.

An FM-index uses the BWT and compact auxiliary structures to represent the searchable text more efficiently. Implementations may also sample suffix-array positions to reduce memory further.

The name *FM-index* refers to its authors, Ferragina and Manzini, and the index's compressed full-text-search capabilities.

### Important Distinction

The FM-index is not simply the BWT. It is a data structure built around the BWT, usually including:

- The BWT of the indexed text.
- A cumulative character-count structure, commonly called the `C` array.
- An occurrence-rank data structure, often implemented using checkpoints, bitvectors, or other compressed representations.
- Optional sampled suffix-array entries for locating matching positions.

Together, these structures support efficient pattern counting and, with suitable suffix-array information, locating matches in the original text.

---

## The Core Components

Consider the text:

```text
T = ACAACG$

```

The `$` symbol is a unique end marker that does not occur elsewhere in the text. It sorts before the ordinary DNA characters.

The BWT of this text is:

```text
L = GC$AAAC

```

We use $L$ to construct the FM-index.

### 1. The C Array

The `C` array records how many characters in the indexed text are lexicographically smaller than a given character.

To calculate it:

1. Count how many times each character appears in the text.
2. Sort the alphabet.
3. Compute the cumulative counts of all characters preceding each character.

For `ACAACG$`, the character counts are:

| Character | Count |
| --- | --- |
| `$` | 1 |
| `A` | 3 |
| `C` | 2 |
| `G` | 1 |

The cumulative array is:

| Character $c$ | $C(c)$ |
| --- | --- |
| `$` | 0 |
| `A` | 1 |
| `C` | 4 |
| `G` | 6 |

For example:

$$C(\text{C}) = 1 + 3 = 4$$

There is one `$` and three `A` characters, so four characters sort before `C`.

Similarly:

$$C(\text{G}) = 1 + 3 + 2 = 6$$

Six characters sort before `G`.

**Why is this useful?**

The `C` array tells us where the block of suffix-array rows beginning with each character starts. It allows the search algorithm to move between the sorted first-column and last-column representations without explicitly constructing the entire BWT rotation matrix.

### 2. The Occurrence-Rank Structure

The occurrence function, often written as $\operatorname{Occ}(c,i)$, counts how many times character $c$ occurs in the BWT before position $i$.

For this document, we use zero-based indexing and a half-open interval:

$$\operatorname{Occ}(c,i) = \#\{j : 0 \le j < i,\ L[j] = c\}$$

In other words, it counts occurrences in $L[0:i]$, excluding position $i$.

For:

```text
L = G C $ A A A C
    0 1 2 3 4 5 6

```

some example rank queries are:

| Query | Result | Explanation |
| --- | --- | --- |
| $\operatorname{Occ}(\text{C},0)$ | 0 | No characters precede position 0 |
| $\operatorname{Occ}(\text{C},4)$ | 1 | One `C` occurs in `L[0:4]` |
| $\operatorname{Occ}(\text{C},7)$ | 2 | Two `C` characters occur in the entire BWT |
| $\operatorname{Occ}(\text{A},4)$ | 1 | One `A` occurs in `L[0:4]` |
| $\operatorname{Occ}(\text{A},6)$ | 3 | Three `A` characters occur in `L[0:6]` |

An implementation does not necessarily store a complete table containing every possible rank query. That would consume $O(n\vert{}\Sigma\vert{})$ space for a text of length $n$ and alphabet size $\vert{}\Sigma\vert{}$.

Instead, practical implementations often store checkpoint counts, bitvectors, or other succinct structures to calculate rank queries efficiently while reducing memory consumption.

---

## How Backward Search Works

**Backward search** is the main FM-index operation for exact pattern matching.

It starts with the last character of a query and processes the remaining characters from right to left. At each step, it updates an interval of suffix-array rows that match the pattern processed so far.

Instead of checking every possible starting position in the text, the algorithm progressively narrows the interval containing all matches.

### The interval representation

Use a half-open suffix-array interval:

$$[l,r)$$

This interval contains the rows with indices:

$$l, l+1, \ldots, r-1$$

The number of matching rows is:

$$\text{count} = r - l$$

Initially, every suffix is a candidate, so for a text of length $n$:

$$[l,r) = [0,n)$$

To extend the current pattern to the left with character $c$, update the interval as follows:

$$l' = C(c) + \operatorname{Occ}(c,l)$$

$$r' = C(c) + \operatorname{Occ}(c,r)$$

The new interval $[l',r')$ contains suffixes matching the extended pattern.

If $l' = r'$, the interval is empty and the pattern does not occur in the indexed text.

---

## Worked Example: Search for `AC`

Suppose the indexed text is:

```text
Text: ACAACG$
BWT:  GC$AAAC

```

We want to find the pattern:

```text
Pattern: AC

```

Because backward search proceeds from right to left, we process `C` first and then `A`.

### Step 1: Initialize the interval

The text contains seven characters, including the end marker.

```text
Initial interval: [0, 7)

```

All suffixes are candidates.

### Step 2: Search for `C`

The last character of the pattern is `C`.

Using:

$$l' = C(\text{C}) + \operatorname{Occ}(\text{C},l)$$

$$r' = C(\text{C}) + \operatorname{Occ}(\text{C},r)$$

and the initial interval $[0,7)$:

$$l' = 4 + \operatorname{Occ}(\text{C},0) = 4 + 0 = 4$$

$$r' = 4 + \operatorname{Occ}(\text{C},7) = 4 + 2 = 6$$

Therefore:

```text
New interval: [4, 6)

```

The interval contains two suffixes beginning with `C`.

### Step 3: Extend the match backward with `A`

Now process the preceding character, `A`.

The current interval is $[4,6)$.

$$l' = C(\text{A}) + \operatorname{Occ}(\text{A},4)$$

$$r' = C(\text{A}) + \operatorname{Occ}(\text{A},6)$$

Substitute the values:

$$l' = 1 + 1 = 2$$

$$r' = 1 + 3 = 4$$

The resulting interval is:

```text
Final interval: [2, 4)

```

Its size is:

$$\text{count} = 4 - 2 = 2$$

Therefore, `AC` occurs **twice** in the indexed text.

The two matching starting positions in `ACAACG` are:

```text
Position 0: AC AACG
Position 3: ACA ACG

```

Positions are zero-based, so the matches begin at positions 0 and 3.

The important point is that backward search determines the number of occurrences by narrowing the suffix-array interval. It does not need to inspect every starting position independently.

---

## Counting Matches vs. Locating Matches

An FM-index supports two related but distinct operations.

### Count

Count determines how many times a pattern occurs by calculating the size of the final suffix-array interval:

$$\text{count} = r - l$$

If the interval is empty, the pattern does not occur. Otherwise, its size gives the exact number of matches.

With efficient constant-time rank queries, counting a pattern of length $m$ takes $O(m)$ time.

### Locate

Locate retrieves the original text coordinates of the matches represented by the interval.

To locate matches, the index needs access to suffix-array information. Storing the complete suffix array allows direct access to each matching suffix's starting position. A memory-saving index may instead store only sampled suffix-array entries.

For example, if suffix-array entries are sampled at regular text-position intervals of size $d$, LF-mapping can move from an unsampled row to a sampled row. In a standard sampling scheme, locating one occurrence can require up to $O(d)$ LF steps, plus the cost of rank operations.

| Operation | Main result | Typical cost |
| --- | --- | --- |
| **Count** | Number of occurrences | $O(m)$, with constant-time rank |
| **Locate with full suffix array** | Coordinates of matches | $O(1)$ per reported match to retrieve its suffix-array value |
| **Locate with sampled suffix array** | Coordinates of matches | Up to $O(d)$ LF steps per match under regular sampling, plus rank-query costs |

These are typical bounds under the stated assumptions, not universal guarantees for every FM-index implementation.

If a pattern occurs $z$ times, reporting all coordinates necessarily incurs additional work proportional to the number of results.

---

## Time and Space Complexity

Let:

* $n$ = length of the indexed text, including the end marker.
* $m$ = length of the query pattern.
* $|\Sigma|$ = alphabet size.
* $z$ = number of occurrences of the query.
* $d$ = suffix-array sampling interval, when regular sampling is used.

### Construction

FM-index construction generally involves:

1. Constructing a suffix array or another suitable ordering structure.
2. Generating the BWT.
3. Building the cumulative character counts.
4. Constructing the rank data structure.
5. Optionally sampling suffix-array entries.

The BWT can be generated from a suffix array in $O(n)$ time, but the total construction complexity depends on the suffix-array algorithm and the rank representation.

For example, practical construction methods may use $O(n\log n)$ time, while suitable suffix-array algorithms can achieve linear-time construction under their assumptions.

Thus, **$O(n)$ construction time and space is possible with suitable algorithms and representations, but it should not be assumed for every implementation**.

### Pattern counting

With constant-time rank operations:

$$\boxed{T_{\text{count}} = O(m)}$$

The algorithm processes each query character once, updating the suffix-array interval.

### Locating occurrences

With a full suffix array, retrieving each matching position from the final interval takes constant time per occurrence after the search interval has been found.

With sampled suffix-array values, locating positions requires LF-mapping steps. For regular sampling at interval $d$, a typical bound is:

$$T_{\text{locate}} = O(zd)$$

assuming constant-time rank operations and a standard sampling scheme.

The exact bound depends on the sampling strategy and index implementation.

### Memory

A naive representation of the BWT alone uses $O(n)$ symbols. The complete FM-index also needs rank support and potentially sampled suffix-array values.

Its memory usage depends on:

* Alphabet size and encoding.
* Rank/checkpoint representation.
* Suffix-array sampling density.
* Whether the index is run-length compressed.
* Whether the complete suffix array is stored.

An FM-index is designed to use less memory than a conventional full-text index while preserving efficient search, but the precise space bound and compression ratio depend on its representation and the text being indexed.

---

## Bioinformatics Application: Mapping Sequencing Reads

### Scenario

Suppose you have millions of short DNA reads generated by a sequencing experiment and need to determine where they originate in a reference genome.

A naive method could compare each read against every possible starting position in the genome. Repeating this for millions of reads would be computationally expensive.

An FM-index makes reference searching more efficient.

### Typical workflow

1. **Prepare the reference genome:** Obtain the reference FASTA sequence and handle contigs or chromosomes consistently.
2. **Build the reference index:** Construct a BWT-based index and its auxiliary structures. This is generally done once for a given reference and reused across many reads.
3. **Process each read or seed:** Search for exact matches using backward search.
4. **Locate candidate positions:** Use suffix-array information to translate matching intervals into reference coordinates.
5. **Handle sequencing errors and biological variation:** Real reads may contain mismatches, insertions, or deletions. A read mapper must therefore extend beyond exact matching and explore candidate alignments under appropriate scoring and search constraints.
6. **Report alignments:** The mapper evaluates candidate positions and reports alignments according to its algorithm and parameters.

### How BWA uses the FM-index

**BWA (Burrows–Wheeler Aligner)** uses BWT-based indexing of the reference genome as a core component of its read-mapping algorithms.

For exact matching, backward search can efficiently identify the suffix-array interval corresponding to a query. For inexact matching, BWA uses additional search and alignment procedures to accommodate differences between reads and the reference.

The FM-index itself does not automatically solve approximate alignment. It provides the compressed search structure that makes the candidate-search process efficient.

Memory usage depends on the BWA version, index format, reference, and auxiliary data structures; there is no single universal RAM requirement for every human-genome FM-index.

---

## FM-Index in Rust

The Rust ecosystem includes the `bio` crate, commonly referred to as **Rust-Bio**, which provides bioinformatics data structures and algorithms, including FM-index functionality.

Rust-Bio is useful when you want to:

* Learn how a BWT-based index is represented in code.
* Explore backward search and occurrence-rank operations.
* Build a prototype sequence-search tool.
* Integrate sequence indexing into a larger Rust bioinformatics pipeline.

For implementation details and the supported API, consult the official project documentation:

* **Rust-Bio documentation:** [https://docs.rs/bio/](https://docs.rs/bio/)
* **Rust-Bio source code:** [https://github.com/rust-bio/rust-bio](https://github.com/rust-bio/rust-bio)

When benchmarking an implementation, measure index construction time, memory consumption, query-counting time, and match-location time separately. These operations can have substantially different costs.

---

## Improvements and Extensions

### 1. Rank checkpoints

Instead of storing a complete occurrence count for every character at every position, store counts at selected checkpoints. To answer a rank query, retrieve the nearest checkpoint and inspect the intervening symbols. This reduces auxiliary storage, at the cost of scanning a small part of the BWT for some queries.

### 2. Suffix-array sampling

Store only a subset of suffix-array entries. When a query returns a row that is not sampled, repeatedly apply LF-mapping until reaching a sampled row and recover the original position. This trades memory for additional work during coordinate recovery.

### 3. Run-length-compressed FM-indexes

When the BWT contains long runs of repeated characters, run-length encoding can reduce storage significantly. Specialized compressed indexes, such as the r-index, are designed for repetitive texts and genome collections with extensive shared sequence.

### 4. Bidirectional FM-index

A bidirectional index supports extending a query to the left or right. This can help search algorithms explore candidate matches in both directions and is useful in seed extension and approximate matching.

### 5. Approximate matching

Exact backward search only identifies occurrences of the query as written. Approximate matching algorithms can explore alternative branches involving substitutions, insertions, or deletions. Such searches require additional logic and can become more expensive as the allowed number of differences increases.

---

## Key Takeaway

The FM-index combines the Burrows–Wheeler Transform with cumulative character counts and occurrence-rank queries to search a text without storing every suffix-array entry explicitly. Backward search narrows a suffix-array interval one character at a time, allowing exact pattern counts in $O(m)$ time under constant-time rank operations. With suffix-array information, the index can also locate matches, making it a foundational component of efficient genome-search and short-read-mapping systems.

