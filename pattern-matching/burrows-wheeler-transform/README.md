# Burrows–Wheeler Transform (BWT)

---

## What Is the Burrows–Wheeler Transform?

The **Burrows–Wheeler Transform (BWT)** is a reversible transformation that rearranges a string into an order that often makes it easier to compress and search.

It was originally developed for data compression and is used in compression algorithms such as `bzip2`. In bioinformatics, the BWT is a fundamental component of compressed text indexes used for searching large reference genomes.

Consider the string:

```text
BANANA
```

Its characters are not arranged in alphabetical order. The BWT reorganizes the string using the ordering of all its cyclic rotations. The resulting sequence is:

```text
ANNB$AA
```

Here, `$` is a special end-of-text marker that occurs exactly once and sorts before ordinary characters.

The transformed string may look unrelated to the original, but it retains enough information to reconstruct the original text exactly.

### Why Is BWT Useful?

The BWT has two important properties:

1. **It can improve compressibility:** Characters that occur in similar contexts often become adjacent in the transformed string. This can create runs of repeated characters that are suitable for run-length encoding and other compression techniques.
2. **It supports efficient searching:** When combined with rank information and a suffix-array representation, the BWT enables an **FM-index**, which can search for exact patterns without storing every suffix of the original text explicitly.

> **Important distinction:** BWT is not compression by itself. It rearranges the data in a reversible way; a compression algorithm can then encode the transformed data more efficiently.

---

## How Does BWT Work?

The classic definition uses cyclic rotations of a string containing a unique end marker.

Let the input text be:

$$T = \text{BANANA}\$$

The length of the text is $n = 7$.

### Step 1: Generate all cyclic rotations

Generate one rotation beginning at each character position:

```text
BANANA$
ANANA$B
NANA$BA
ANA$BAN
NA$BANA
A$BANAN
$BANANA
```

Each row contains the same characters in a different cyclic order.

For a string of length $n$, this process generates $n$ rotations. Although this is easy to understand, explicitly generating all rotations is inefficient for large strings because the total amount of stored text can grow to $O(n^2)$.

### Step 2: Sort the rotations lexicographically

Sort all rotations in alphabetical order, treating `$` as a unique character that sorts before `A`, `B`, and the other ordinary characters.

The sorted rotations are:

```text
$BANANA
A$BANAN
ANA$BAN
ANANA$B
BANANA$
NA$BANA
NANA$BA
```

The sorted order is the essential step: it establishes the relationship between the first and last characters of the rotations.

### Step 3: Extract the last column

Take the final character of every sorted row:

```text
$BANANA  -> A
A$BANAN  -> N
ANA$BAN  -> N
ANANA$B  -> B
BANANA$  -> $
NA$BANA  -> A
NANA$BA  -> A
```

Concatenating these characters gives:

```text
BWT(BANANA$) = ANNB$AA
```

The result is the **Burrows–Wheeler Transform** of $\text{BANANA}\$$.

### Why can characters cluster together?

The sorted rotations place strings with similar prefixes next to one another. Because each rotation is a cyclic rearrangement, the characters immediately preceding those prefixes appear in the last column.

This often brings identical characters together, creating runs that can be compressed efficiently.

For example, the transformed string contains:

```text
ANNB$AA
```

The two `N` characters and the two trailing `A` characters form runs.

However, BWT does not guarantee that every character will be grouped into one contiguous block. How much clustering occurs depends on the input and the contexts in which characters appear.

---

## The Mathematical Definition

Let $T$ be a string of length $n$ containing a unique end marker `$`.

Construct the matrix of all $n$ cyclic rotations of $T$, and sort its rows lexicographically.

The BWT is the string formed by taking the last character of every row in sorted order:

$$\operatorname{BWT}(T)[i] = M[i][n-1]$$

where $M$ is the sorted rotation matrix and $i$ is a row index.

An equivalent and more efficient formulation uses a **suffix array**.

Let $SA[i]$ denote the starting position of the $i$-th suffix in lexicographically sorted order. Then:

$$L[i] = \begin{cases} T[SA[i]-1], & SA[i] > 0 \\ T[n-1], & SA[i] = 0 \end{cases}$$

Here, $L$ is the BWT, and positions are zero-indexed.

This formulation makes it possible to construct the BWT from the suffix array without explicitly materializing the entire rotation matrix.

---

## The Key Property: BWT Is Reversible

The BWT is reversible because the transformed string preserves the ordering relationships needed to reconstruct the original sequence.

A central concept behind reconstruction is the **LF-mapping**.

### First and last columns

Let:
- $L$ be the last column, which is the BWT.
- $F$ be the first column of the sorted rotation matrix.

The first column is simply the BWT characters sorted lexicographically.

For $\text{BANANA}\$$:

```text
F = $AAABNN
L = ANNB$AA
```

The rows of the sorted rotation matrix are not needed to calculate the LF-mapping if the necessary character counts and occurrence ranks are available.

### What is LF-mapping?

The LF-mapping tells us how to move between corresponding rows of the sorted rotation matrix.

> The $r$-th occurrence of a character $c$ in the last column corresponds to the $r$-th occurrence of $c$ in the first column.

Occurrences are counted from left to right.

Define:
- $C(c)$: the number of characters in the text lexicographically smaller than $c$.
- $\operatorname{Occ}(c,i)$: the number of occurrences of $c$ in $L[0 \ldots i]$, including position $i$.

Then the LF-mapping is:

$$LF(i) = C(L[i]) + \operatorname{Occ}(L[i],i)-1$$

This version uses zero-based row indices and an inclusive occurrence count.

The mapping allows us to trace backward through the sorted rotations and recover the original text.

---

## Worked Example: Transform and Reconstruct `ACAACG`

Consider the genome fragment:

```text
ACAACG
```

Append the unique end marker:

```text
ACAACG$
```

### Part 1: Compute the BWT

The seven cyclic rotations are:

```text
ACAACG$
CAACG$A
AACG$AC
ACG$ACA
CG$ACAA
G$ACAAC
$ACAACG
```

After sorting lexicographically:

```text
$ACAACG
AACG$AC
ACAACG$
ACG$ACA
CAACG$A
CG$ACAA
G$ACAAC
```

Extract the last character of each row:

```text
$ACAACG  -> G
AACG$AC  -> C
ACAACG$  -> $
ACG$ACA  -> A
CAACG$A  -> A
CG$ACAA  -> A
G$ACAAC  -> C
```

Therefore:

```text
Original text: ACAACG$
BWT:           GC$AAAC
```

The transformed string contains a run of three consecutive `A` characters. The original text did not contain this exact run, illustrating how the transformation can create new runs.

### Part 2: Reconstruct the original text

Start with:

```text
L = GC$AAAC
F = $AAACCG
```

For each position in $L$, calculate the LF-mapping using character counts and occurrence ranks:

| Row $i$ | $L[i]$ | Occurrence rank | $LF(i)$ |
|---:|:---:|---:|---:|
| 0 | G | 1 | 6 |
| 1 | C | 1 | 4 |
| 2 | $ | 1 | 0 |
| 3 | A | 1 | 1 |
| 4 | A | 2 | 2 |
| 5 | A | 3 | 3 |
| 6 | C | 2 | 5 |

Start at the row containing `$` in the last column, which is row 2.

Repeatedly prepend the current last-column character and move to the row specified by LF-mapping:

```text
Start:  $
Next:   G$
Next:   CG$
Next:   ACG$
Next:   AACG$
Next:   CAACG$
Next:   ACAACG$
```

The original text is recovered:

```text
ACAACG$
```

Remove the special end marker to obtain the original input:

```text
ACAACG
```

The reconstruction is exact because BWT preserves the information needed to recover the original text.

> **Implementation note:** An inverse-BWT implementation can perform this operation in $O(n)$ time when it has efficient access to character counts and occurrence ranks. A simple educational implementation can also reconstruct the text by repeatedly sorting rows, but that approach uses substantially more time and memory.

---

## How BWT Enables Efficient Genome Searching

A BWT alone is not a complete genome-search index. Its real power for sequence searching comes from combining it with an **FM-index**, which stores the BWT along with character-count information, rank data, and often a sampled suffix array.

The FM-index supports exact pattern searches by processing the query from right to left.

### Backward search

Suppose the reference is:

```text
Reference: GATTACA
```

and the query is:

```text
Query: TAC
```

Backward search processes the pattern in reverse character order:

```text
C -> A -> T
```

At each step, the algorithm updates an interval of suffix-array rows that could match the pattern.

The LF-related operation is often written as:

$$[l',r'] = \left[ C(c)+\operatorname{Occ}(c,l-1), C(c)+\operatorname{Occ}(c,r)-1 \right]$$

This formula assumes an inclusive interval $[l,r]$, zero-based indexing, and an occurrence function that counts characters in the half-open prefix $L[0 \ldots i)$. In practice, implementations differ in their indexing conventions.

At a high level:
1. Initialize the interval to include all suffixes.
2. Process the query's last character.
3. Use the character counts and occurrence ranks to restrict the interval to matching suffixes.
4. Repeat for each preceding character.
5. If the interval becomes empty, the pattern does not occur in the indexed text.
6. Otherwise, the remaining interval identifies suffixes beginning with the query pattern.

The number of characters in the query is $m$. With constant-time rank operations, backward search takes $O(m)$ time to determine the matching interval, apart from implementation-specific overhead.

Locating the actual reference coordinates generally requires additional suffix-array information. With sampled suffix arrays, locating a position may involve multiple LF steps.

### Why is this useful for bioinformatics?

A bacterial or human reference genome can contain millions or billions of bases. Searching every possible starting position independently for each sequencing read is expensive.

An FM-index allows a read mapper to narrow the possible matching positions using the query's characters and the index, rather than scanning the entire reference for every read.

This is particularly useful for short-read sequencing, where millions of reads may need to be compared against a large reference genome.

---

## Applications in Bioinformatics

### 1. Short-read mapping
Tools such as **BWA** use BWT-based indexing and FM-index techniques to search large reference genomes efficiently. The index supports exact-match searching, while practical read mapping also needs to handle mismatches, insertions, deletions, and other sources of sequencing variation.

### 2. Reference genome indexing
A reference genome can be indexed once and reused to search many sequencing reads. This avoids rebuilding a complete search structure for every query and is especially valuable for high-throughput sequencing workflows.

### 3. Data compression
The transformed text may contain repeated character runs that are suitable for run-length encoding or other compression techniques. The actual compression ratio depends on the input sequence and the compression method used after the transform.

### 4. Repetitive sequence analysis
The BWT and related compressed indexes can represent and search large repetitive texts efficiently. Specialized structures, including run-length-compressed BWT indexes and r-indexes, are useful when the transformed sequence contains relatively few runs.

---

## Time and Space Complexity

Let $n$ be the length of the input text, including the end marker.

### BWT Construction & Operations

| Method | Time Complexity | Space Considerations |
| :--- | :--- | :--- |
| **Explicit rotation matrix** | At least $O(n^2)$ to materialize all rotations, plus sorting | $O(n^2)$ for the matrix |
| **Suffix-array-based construction** | $O(n)$ to derive the BWT once the suffix array is available | $O(n)$ additional output space, plus suffix-array construction storage |
| **Practical suffix-array construction** | Depends on the algorithm; commonly $O(n\log n)$ or linear-time variants | Usually linear or near-linear storage |
| **Inverse BWT** (with counts/ranks) | $O(n)$ | $O(n)$ for transformed text and reconstruction output |
| **FM-index backward search** (query length $m$) | $O(m)$ assuming constant-time rank queries | Auxiliary rank data structures and sampled suffix arrays |

> **Important distinction:** The often-quoted $O(m)$ search complexity does not mean that every genome-mapping operation takes exactly $O(m)$ time. Practical read mapping includes locating matches and evaluating approximate alignments, which may require additional computation.

---

## Improvements and Extensions

### 1. Run-length encoding
Run-length encoding stores repeated characters as a character and a run length (e.g., `AAAAACCCGG` becomes `A5C3G2`). The BWT often creates runs that can make this type of compression effective.

### 2. Suffix-array sampling
Storing every suffix-array value can consume substantial memory. Sampling suffix-array entries reduces storage, while LF-mapping can be used to reach sampled entries at the trade-off of extra calculation steps.

### 3. Run-length-compressed BWT and r-index
A run-length-compressed BWT stores runs rather than every character independently. The **r-index** is a compressed full-text indexing approach designed to support searching and locating matches in repetitive texts.

### 4. Bidirectional FM-index
A bidirectional FM-index supports extension of a query in either direction, which is useful when exploring alternative matches or extending seeds during sequence alignment.

### 5. External-memory construction
For very large references or collections of genomes, external-memory algorithms use disk-based intermediate data to handle datasets that do not fit comfortably in RAM, at the cost of additional I/O overhead.

---

## Practical Considerations for Genomic Data

* **Use a unique end marker:** The sentinel must not appear as an ordinary character in the indexed text. For DNA, `$` is a convenient illustrative choice.
* **Be explicit about indexing:** LF-mapping formulas differ depending on whether indices are zero- or one-based and whether occurrence counts include the current position.
* **Distinguish exact search from read mapping:** An FM-index finds exact pattern matches efficiently, but real reads often contain errors or biological variation.
* **Account for reference representation:** Concatenating multiple chromosomes or contigs requires care so that searches do not report matches crossing artificial sequence boundaries.
* **Consider repetitive regions:** Highly repetitive genomes can produce many matching locations, and locating or resolving them can be more expensive than finding the search interval.
* **Benchmark memory and speed:** Index construction, index storage, query processing, and coordinate recovery have different performance costs.

---

## Key Takeaway

The Burrows–Wheeler Transform rearranges a string into a reversible representation that often improves compressibility and supports efficient indexing. In bioinformatics, its combination with suffix arrays, character-rank structures, and FM-index techniques makes it possible to search large reference genomes efficiently. Understanding the BWT, LF-mapping, and backward search provides a foundation for understanding how modern sequencing tools index and query genomic sequences.