# Suffix Arrays

## What Is a Suffix Array?

A **suffix array** is a data structure that stores the starting positions of all suffixes of a text in lexicographically sorted order.

A suffix is the substring that begins at a particular position and continues to the end of the text.

For example, consider:

```text
Text: BANANA$

```

The `$` symbol is a unique end-of-text marker that sorts before ordinary characters.

The text has seven suffixes:

| Starting position | Suffix |
| --- | --- |
| 0 | `BANANA$` |
| 1 | `ANANA$` |
| 2 | `NANA$` |
| 3 | `ANA$` |
| 4 | `NA$` |
| 5 | `A$` |
| 6 | `$` |

A suffix array sorts these suffixes alphabetically but stores their starting positions rather than storing the suffix strings themselves.

The sorted suffixes are:

| Sorted rank | Starting position | Suffix |
| --- | --- | --- |
| 0 | 6 | `$` |
| 1 | 5 | `A$` |
| 2 | 3 | `ANA$` |
| 3 | 1 | `ANANA$` |
| 4 | 0 | `BANANA$` |
| 5 | 4 | `NA$` |
| 6 | 2 | `NANA$` |

Therefore, the suffix array is:

```text
SA = [6, 5, 3, 1, 0, 4, 2]

```

Each entry identifies the starting position of a suffix in sorted order.

**Key idea:** Instead of storing every suffix as a separate string, a suffix array stores integer positions that refer back to the original text.

---

## Why Are Suffix Arrays Useful?

Suffix arrays provide a compact foundation for searching and analyzing long strings.

### 1. Efficient pattern matching

A suffix array sorts suffixes according to their contents. This ordering allows binary search to narrow down the range of suffixes that begin with a query pattern.

For a pattern of length $m$ in a text of length $n$, a straightforward binary-search implementation takes:

$$O(m\log n)$$

time, assuming each comparison can examine up to $m$ characters.

If the pattern occurs multiple times, all matching suffixes occupy a contiguous interval of the suffix array. The index can identify this interval and report the corresponding text positions.

### 2. More space-efficient than storing all suffixes

Explicitly storing every suffix would duplicate much of the text. A suffix array instead stores one integer per suffix.

For a text of length $n$, the array has $n$ entries.

If each entry uses a 32-bit integer, the array itself requires approximately:

$$4n \text{ bytes}$$

For example, a text containing 5 million characters requires approximately 20 MB for the suffix-array entries alone.

Actual memory consumption also includes the original text, temporary construction structures, and any auxiliary search data.

A 32-bit representation is sufficient only when the positions fit within the supported integer range. Larger texts may require 64-bit entries.

### 3. Foundation for other string indexes

Suffix arrays are closely related to several important bioinformatics data structures:

* **Burrows–Wheeler Transform (BWT):** Can be constructed from the suffix array by extracting the character preceding each suffix in sorted order, with a special case for the suffix beginning at position zero.
* **FM-index:** Uses the BWT, character counts, and occurrence-rank structures to support compressed pattern searching.
* **Longest Common Prefix (LCP) array:** Stores the lengths of common prefixes between neighboring suffixes in suffix-array order.
* **Repeated-substring analysis:** Suffix arrays and LCP arrays can identify repeated sequences, maximal matches, and other sequence patterns.

---

## How Is a Suffix Array Constructed?

The goal is to sort all suffixes of the input text without necessarily creating a separate copy of every suffix.

There are several construction algorithms, each with different time, space, and implementation trade-offs.

### 1. Naive construction

The simplest approach is:

1. Generate a suffix for every starting position.
2. Store each suffix alongside its starting position.
3. Sort the suffixes lexicographically.
4. Extract the starting positions in sorted order.

For a text of length $n$, there are $n$ suffixes, and their total length can be $O(n^2)$.

A straightforward comparison sort performs $O(n\log n)$ comparisons, but each comparison can take up to $O(n)$ time.

Consequently, a simple upper bound for this implementation is:

$$T(n) = O(n^2\log n)$$

The memory requirement can also be $O(n^2)$ if every suffix is explicitly materialized.

This approach is useful for teaching and testing, but it is unsuitable for large genomes.

### 2. Prefix-doubling algorithms

Prefix-doubling methods construct a suffix array by repeatedly sorting suffixes according to increasingly long prefixes.

The basic idea is:

1. Assign each character an initial rank.
2. Sort suffixes using ranks describing prefixes of length 1.
3. Combine rank information to represent prefixes of length 2, then 4, then 8, and so on.
4. Continue until the prefix lengths are sufficient to distinguish all suffixes.
5. Return the suffix positions in final sorted order.

At each iteration, the algorithm uses previously computed ranks instead of comparing entire suffix strings.

A common prefix-doubling implementation with efficient sorting achieves:

$$T(n) = O(n\log n)$$

time and typically $O(n)$ auxiliary space, depending on the sorting method and representation.

This approach is often easier to implement than more advanced linear-time algorithms.

### 3. SA-IS (Suffix Array by Induced Sorting)

**SA-IS** constructs suffix arrays using induced sorting.

Rather than sorting every suffix independently, it classifies suffix positions into types and uses relationships between selected suffixes to induce the ordering of the remaining suffixes.

At a high level, SA-IS:

1. Classifies suffix positions according to their relationships with neighboring characters.
2. Identifies special positions, commonly called LMS positions.
3. Sorts these selected suffixes.
4. Uses the sorted order to induce the order of other suffixes.
5. Recursively solves a reduced problem when necessary.
6. Produces the final suffix array.

Under its standard assumptions and with a suitable integer alphabet representation, SA-IS achieves:

$$T(n) = O(n)$$

time and $O(n)$ space.

Linear-time construction is theoretically attractive, but implementation complexity, memory allocation, alphabet handling, and cache behavior can affect practical performance.

### 4. DC3 / Skew algorithm

The DC3 algorithm is another linear-time suffix-array construction algorithm.

It partitions suffix positions into groups based on their indices modulo 3, recursively sorts a subset of suffixes, and combines the resulting orderings to produce the complete suffix array.

Under its standard assumptions, DC3 achieves:

$$T(n) = O(n)$$

time and linear auxiliary space.

It is an important theoretical and practical alternative to other suffix-array construction methods, although the best choice depends on the application and implementation.

---

## Worked Example: Find a DNA Motif

Suppose you have assembled a draft bacterial genome and want to locate every occurrence of a short regulatory motif.

Use this illustrative genome:

```text
Genome: ATCGATCGATCG$

```

The query motif is:

```text
Pattern: ATCG

```

For simplicity, use zero-based coordinates and include `$` as the unique end marker.

### Step 1: List the suffixes

The text contains 13 characters, including the end marker.

```text
Position 0:   ATCGATCGATCG$
Position 1:   TCGATCGATCG$
Position 2:   CGATCGATCG$
Position 3:   GATCGATCG$
Position 4:   ATCGATCG$
Position 5:   TCGATCG$
Position 6:   CGATCG$
Position 7:   GATCG$
Position 8:   ATCG$
Position 9:   TCG$
Position 10:  CG$
Position 11:  G$
Position 12:  $

```

### Step 2: Sort the suffixes

After lexicographic sorting, the suffix array is:

```text
SA = [12, 8, 4, 0, 9, 5, 1, 10, 6, 2, 11, 7, 3]

```

The suffixes beginning with `A` appear together near the beginning of the sorted suffix list. This ordering is what makes binary-search pattern matching possible.

### Step 3: Search for `ATCG`

Binary search compares the query with suffixes at selected positions in the suffix array.

At each step:

* If the suffix is lexicographically smaller than the pattern, search the upper half.
* If the suffix is lexicographically larger than the pattern, search the lower half.
* If the suffix begins with the complete pattern, record a match and identify the contiguous matching interval.

For this genome, the suffixes beginning with `ATCG` start at positions:

```text
0
4
8

```

Therefore, the motif occurs three times:

```text
Position 0: ATCGATCGATCG$
Position 4:     ATCGATCG$
Position 8:         ATCG$

```

The matching positions are:

```text
Occurrences = [0, 4, 8]

```

The positions are in the original genome coordinate system, not in suffix-array rank order.

**Important:** A complete implementation should use lower-bound and upper-bound binary searches over suffixes to find the entire interval matching a pattern prefix. A single successful binary-search comparison does not necessarily identify all occurrences.

---

## Pattern Matching with a Suffix Array

Let:

* $n$ = length of the indexed text.
* $m$ = length of the query pattern.
* $z$ = number of occurrences of the query.

A straightforward suffix-array search performs binary search over $n$ suffixes. Each comparison may inspect up to $m$ characters.

The search time is therefore:

$$T_{\text{search}} = O(m\log n)$$

Once the matching interval is identified, reporting the $z$ matching positions takes at least $O(z)$ time.

Thus, for a basic implementation:

$$T_{\text{total}} = O(m\log n + z)$$

This bound assumes ordinary suffix comparisons without additional longest-common-prefix acceleration.

### Can searching be faster?

Yes. An LCP-aware search can reuse information from previous comparisons and reduce repeated character comparisons. Other indexes, such as the FM-index, use different structures to achieve efficient pattern counting.

A suffix array alone does not automatically provide every optimization used by modern genome-search tools.

---

## Time and Space Complexity

Let $n$ be the length of the text.

| Operation | Typical time complexity | Auxiliary space |
| --- | --- | --- |
| **Naive construction by sorting explicit suffixes** | $O(n^2\log n)$ upper bound | Up to $O(n^2)$ |
| **Prefix doubling** | $O(n\log n)$ for common implementations | Typically $O(n)$ |
| **SA-IS construction** | $O(n)$ under standard assumptions | $O(n)$ |
| **DC3 / Skew construction** | $O(n)$ under standard assumptions | $O(n)$ |
| **Basic pattern search** | $O(m\log n)$ | $O(1)$ auxiliary search space |
| **Report all matches** | Additional $O(z)$ | Depends on whether results are materialized |

The suffix array itself requires $O(n)$ space. The original text must also remain available for suffix comparisons, and construction algorithms may require temporary arrays.

For a 32-bit suffix array:

$$S_{\text{SA}} \approx 4n \text{ bytes}$$

This is only the storage for the integer array. The total memory footprint depends on the text encoding, construction method, and auxiliary indexes.

---

## The Relationship Between Suffix Arrays, BWT, and FM-Index

These three data structures are closely related, but they serve different roles.

| Structure | Main purpose | Core idea |
| --- | --- | --- |
| **Suffix array** | Search and order suffixes | Stores starting positions in sorted suffix order |
| **BWT** | Reversible transformation and compressed indexing | Rearranges characters according to suffix/rotation order |
| **FM-index** | Compressed pattern searching | Uses the BWT and rank information to narrow matching intervals |

### From suffix array to BWT

Given a text $T$ and its suffix array $SA$, the BWT can be calculated as:

$$L[i] = \begin{cases} T[SA[i]-1], & SA[i] > 0 \\ T[n-1], & SA[i] = 0 \end{cases}$$

This extracts the character preceding each suffix, with wraparound for the suffix beginning at position zero.

The unique end marker makes the relationship between suffixes, rotations, and the original text well-defined.

### From BWT to FM-index

The FM-index builds on the BWT by adding:

* Cumulative character counts.
* Occurrence-rank support.
* Optional suffix-array samples for locating matches.

It can then perform backward search without requiring the complete suffix array to be stored explicitly.

### Why does this matter?

Suffix arrays are useful as a direct search index and as a building block for more advanced data structures. The BWT and FM-index extend these ideas toward compressed indexing, which is especially valuable when working with large genomic references.

---

## The kISS Algorithm and k-Ordered Suffix Arrays

A standard suffix array orders suffixes using their complete contents. Some specialized applications do not require the full suffix ordering.

A **k-ordered suffix array** orders suffixes according to their first $k$ characters, rather than comparing them all the way to the end of the text.

This can be useful when a downstream algorithm only needs to distinguish or group short prefixes, such as seeds used in certain read-mapping strategies.

The kISS algorithm explores this type of ordering for read-mapping applications.

However, a k-ordered suffix array is not generally interchangeable with a conventional suffix array:

* Suffixes that share the same first $k$ characters may remain tied under the truncated ordering.
* Searches longer than $k$ characters may need additional verification or refinement.
* The benefit depends on the application and the assumptions made about query lengths and matching behavior.

For a production implementation, consult the original kISS publication and compare its guarantees with those of a full suffix array. The exact performance advantage depends on the algorithm, data, and downstream search strategy.

---

## Practical Considerations for Genomic Data

### 1. Reference representation

For a genome containing multiple chromosomes or contigs, decide how each sequence is represented and whether artificial separators are needed. A query should not accidentally match across a boundary created by concatenating separate contigs.

### 2. Coordinate conventions

Suffix arrays usually store zero-based starting positions, but many downstream bioinformatics formats use one-based coordinates or half-open intervals. Convert coordinates explicitly when exporting results.

### 3. Repetitive sequences

Repeated motifs can appear at many genomic locations. The suffix array can identify their matching interval, but reporting and interpreting all occurrences may require substantial additional work.

### 4. Exact matching versus approximate matching

A basic suffix-array search identifies exact pattern occurrences. Sequencing reads often contain errors or biological variants, so practical read mapping may require seed-and-extend strategies, mismatch-aware search, or dynamic-programming alignment.

### 5. Construction memory

A suffix array is compact once built, but some construction algorithms require temporary arrays or substantial working memory. Measure construction memory separately from the final index size.

---

## Key Takeaway

A suffix array stores the starting positions of all suffixes in lexicographically sorted order, enabling efficient pattern searching without explicitly storing every suffix. It is useful for locating DNA motifs, analyzing repeated sequences, and building advanced indexes. It also forms a foundation for the BWT and FM-index, which support compressed searching of large genomic references.
