# Shortest Superstring Problem

## What Is the Shortest Superstring Problem?

The **Shortest Superstring Problem (SSP)** is a classic string optimization problem in which we are given a collection of strings and want to construct the shortest possible string that contains every input string as a substring.

Formally, given a set of strings

$$S = \{s_1, s_2, \ldots, s_n\}$$

we want to find a string $T$ such that:

1. Every input string $s_i$ occurs as a contiguous substring of $T$.
2. The length of $T$, denoted by $|T|$, is minimized.

The objective is:

$$\boxed{T^* = \arg\min_T |T|}$$

subject to:

$$s_i \text{ is a substring of } T \quad \forall s_i \in S$$

The shortest superstring problem is closely related to genome assembly. Sequencing technologies produce DNA reads sampled from a genome, and genome assembly attempts to reconstruct the underlying sequence from these reads. Overlapping reads can provide evidence about how different portions of the genome fit together.

However, the shortest superstring formulation is a simplified model of genome assembly. Real sequencing data contain errors, repeated sequences, uneven coverage, and sometimes multiple possible reconstructions.

---

## Understanding String Overlaps

The central idea behind the shortest superstring problem is that two strings may share a suffix-prefix overlap.

An **overlap** occurs when a suffix of one string matches a prefix of another string.

For example:

```text
String A: ATGCCA
String B:     CAGT
```

The suffix `CA` of String A matches the prefix `CA` of String B.

Therefore, the two strings can be merged without repeating the shared characters:

```text
ATGCCA
    CAGT
--------
ATGCCAGT
```

The merged string is:

```text
ATGCCAGT
```

The overlap length is 2, so the merged length is:

$$|A \oplus B| = |A| + |B| - 2$$

where $A \oplus B$ denotes merging the strings using their overlap.

More generally, if the maximum suffix-prefix overlap between $A$ and $B$ has length $k$, then:

$$|A \oplus B| = |A| + |B| - k$$

**Important:** For the standard greedy shortest-superstring algorithm, we consider suffix-prefix overlaps between different strings. An arbitrary shared substring is not sufficient to merge two strings this way.

---

## The Greedy Strategy

Finding the globally shortest superstring is computationally difficult. A common heuristic is the **greedy shortest-superstring algorithm**.

The algorithm repeatedly chooses the pair of strings with the largest available suffix-prefix overlap, merges them, and continues until no strings remain to be merged.

### Algorithm

1. Start with the original collection of strings.
2. Calculate the suffix-prefix overlap for each eligible ordered pair of strings.
3. Identify the pair with the largest overlap.
4. Merge the pair, keeping the shared region only once.
5. Remove the two original strings and insert the merged string.
6. Recalculate overlaps involving the updated string.
7. Repeat until one string remains or no positive overlaps are available.

The greedy algorithm makes a locally optimal choice at each iteration. It does not guarantee that the final superstring is globally optimal.

### Worked Example

Consider the following DNA reads:

```text
R1: ATGCCA
R2: GCCATG
R3: CATGCC
R4: TGCCAT
```

**Step 1: Find the largest overlap.**

The suffix `GCCA` of R1 matches the prefix `GCCA` of R3 only if the strings are oriented accordingly; here, the actual suffix-prefix overlaps should be calculated directly.

For example:

```text
R1: ATGCCA
R2:     GCCATG
```

The suffix `GCC` of R1 matches the prefix `GCC` of R2, producing:

```text
ATGCCATG
```

The overlap length is 3.

**Step 2: Merge the selected reads.**

Replace R1 and R2 with the merged string:

```text
M1: ATGCCATG
R3: CATGCC
R4: TGCCAT
```

**Step 3: Recalculate overlaps.**

The merged string is treated as a new string. Its overlaps with the remaining reads may differ from the original pairwise overlaps, so they must be recalculated.

Continue selecting the largest available overlap, merging strings, and updating the collection.

**Step 4: Stop when no useful overlaps remain.**

If positive overlaps are exhausted before all strings have been merged, the remaining strings can be concatenated with no overlap, producing a valid superstring. This may not be the shortest possible result.

> The greedy method is a heuristic, not an exact solution. The order of merges can influence the final sequence, especially when several pairs have equal or similar overlap lengths.

---

## Python Implementation

The following implementation demonstrates the basic greedy strategy using exact suffix-prefix overlaps.

```python
def overlap_length(a: str, b: str) -> int:
    """Return the longest suffix of a matching a prefix of b."""
    max_len = min(len(a), len(b))

    for length in range(max_len, 0, -1):
        if a.endswith(b[:length]):
            return length

    return 0


def merge_strings(a: str, b: str, overlap: int) -> str:
    """Merge a and b, retaining their shared overlap once."""
    return a + b[overlap:]


def greedy_shortest_superstring(strings: list[str]) -> str:
    """Construct a superstring using greedy maximum-overlap merges."""
    strings = list(dict.fromkeys(strings))

    if not strings:
        return ""

    # Remove strings already contained in another input string.
    strings = [
        s for s in strings
        if not any(
            s != t and s in t
            for t in strings
        )
    ]

    while len(strings) > 1:
        best_i, best_j = -1, -1
        best_overlap = 0

        for i in range(len(strings)):
            for j in range(len(strings)):
                if i == j:
                    continue

                overlap = overlap_length(strings[i], strings[j])

                if overlap > best_overlap:
                    best_overlap = overlap
                    best_i, best_j = i, j

        if best_i == -1:
            # No positive overlaps remain.
            strings[0] += "".join(strings[1:])
            break

        merged = merge_strings(
            strings[best_i],
            strings[best_j],
            best_overlap
        )

        # Remove the merged strings, then add their replacement.
        strings = [
            s for k, s in enumerate(strings)
            if k not in (best_i, best_j)
        ]
        strings.append(merged)

    return strings[0]


reads = ["ATG", "TGC", "GCA", "CAT"]

result = greedy_shortest_superstring(reads)
print(result)
```

One possible result is:

```text
ATGCAT
```

The result contains every input read as a substring:

- `ATG` occurs at positions 1–3.
- `TGC` occurs at positions 2–4.
- `GCA` occurs at positions 3–5.
- `CAT` occurs at positions 4–6.

The implementation uses exact overlaps and does not attempt to correct sequencing errors or reverse-complement reads.

---

## Time Complexity

Let:

- $n$ be the number of input strings.
- $L$ be the maximum length of an input string.

### Pairwise overlap calculation

For a pair of strings, the simple implementation tests possible overlap lengths from largest to smallest. Each suffix-prefix comparison can take up to $O(L)$ time, and there are up to $O(L)$ candidate lengths.

Therefore, a straightforward overlap calculation can take:

$$O(L^2)$$

### Greedy merging

At each iteration, the algorithm examines $O(n^2)$ ordered pairs. There are $O(n)$ merge iterations, and the straightforward overlap computation can cost $O(L^2)$ per pair.

A conservative bound for this implementation is therefore:

$$\boxed{O(n^3L^2)}$$

This bound assumes string lengths remain $O(L)$. In practice, merged strings can grow to length $O(nL)$, which can increase the cost further if overlap calculations repeatedly scan entire merged strings.

The commonly quoted $O(n^2L)$ complexity is achievable only under different assumptions or with more efficient overlap computation and data structures; it should not be treated as a universal bound for a naive implementation.

### Space complexity

The algorithm stores the current strings and intermediate merged strings. If the final superstring has length $M$, the output itself requires $O(M)$ space, with additional temporary storage depending on the implementation.

---

## Is the Greedy Algorithm Optimal?

No. The greedy strategy can produce a valid superstring that is longer than the globally shortest one.

Why?

A large overlap selected early may prevent the algorithm from using a different combination of overlaps that would yield a shorter final sequence. Each merge changes the strings available for later merges.

### Approximation guarantee

The classical greedy shortest-superstring algorithm has a proven approximation guarantee: its output length is at most twice the optimum length, under the standard formulation.

If $G$ is the length of the greedy result and $OPT$ is the length of an optimal superstring, then:

$$OPT \leq G \leq 2OPT$$

Equivalently, the greedy algorithm achieves at least half the optimal solution's value when expressed as a maximization of overlap savings under the standard approximation analysis.

This distinction matters: the shortest-superstring problem minimizes the **length of the final string**, while the phrase "half the optimal compression" refers to a related overlap-saving objective. It is clearer to state the approximation guarantee directly in terms of superstring length.

The factor of 2 is a worst-case bound, not a prediction that the greedy result will typically be twice as long as the optimum.

---

## Why Is the Shortest Superstring Problem Difficult?

The shortest superstring problem is **NP-hard** in its general form. Its decision version—whether a superstring of length at most a specified bound exists—is NP-complete.

An exact algorithm must consider how strings can be ordered and overlapped. The number of possible orderings grows factorially with the number of strings.

For example, considering all possible orders of $n$ strings can involve:

$$n!$$

candidate orderings, before accounting for other details.

For small inputs, exact approaches such as dynamic programming over subsets can find an optimal solution. For large datasets, heuristics and specialized algorithms are generally more practical.

---

## Applications in Bioinformatics

### Genome assembly

In sequencing, reads are short DNA fragments sampled from a longer genome. If two reads overlap, their shared sequence can help determine their relative positions.

For example:

```text
Read 1: ATGCCA
Read 2:     CCATGC
Read 3:         TGCGAA
```

Overlap information can help combine compatible reads into longer sequences called **contigs**.

However, real genome assembly is more complicated than the basic shortest-superstring model because reads may contain sequencing errors, repeats, and ambiguous overlaps.

### Long-read assembly

Long-read technologies, such as PacBio HiFi and Oxford Nanopore sequencing, can generate reads that span larger genomic regions. Overlap-based methods can use these reads to construct assembly graphs and infer contigs.

Real assemblers generally use additional steps for error correction, repeat handling, graph simplification, and consensus generation.

### Sequence reconstruction

The shortest-superstring formulation is also useful for understanding how overlapping fragments can reconstruct a longer DNA, RNA, or text sequence.

---

## Limitations of the Greedy Approach

The basic greedy algorithm assumes that overlaps can be detected reliably and that strings can be merged without introducing ambiguity.

Its main limitations include:

- **Local decisions:** The largest immediate overlap may not lead to the globally shortest superstring.
- **Sequencing errors:** Even a single mismatch can disrupt an exact overlap.
- **Repeated sequences:** Different reads may overlap equally well in multiple locations.
- **Read orientation:** DNA reads may need to be considered in both forward and reverse-complement orientations.
- **Coverage gaps:** If no reads bridge a region, overlap information alone may not determine the missing sequence.
- **Biological ambiguity:** Multiple reconstructions may be consistent with the same reads.

For these reasons, the greedy shortest-superstring algorithm is primarily an educational model and a useful introduction to overlap-based sequence reconstruction, rather than a complete modern genome assembler.

---

## Shortest Superstring vs. de Bruijn Graph Assembly

| Feature | Shortest Superstring / Overlap-Based | de Bruijn Graph |
|---|---|---|
| Main representation | Reads or merged strings | k-mers and graph edges |
| Core operation | Find suffix-prefix overlaps | Connect overlapping k-mers |
| Common strength | Naturally uses long-read overlaps | Efficiently handles large collections of short reads |
| Major challenge | Many pairwise overlaps and repeat ambiguity | Graph complexity, errors, and repeat resolution |
| Common applications | Long-read assembly and overlap-layout-consensus methods | Short-read assembly and many genome assembly pipelines |

The shortest-superstring problem provides a useful mathematical foundation for understanding overlap-based assembly. However, modern assemblers do not simply apply the basic greedy algorithm to all reads. They use specialized data structures, graph algorithms, statistical models, and error-handling strategies to make assembly practical and reliable.

---

## Key Takeaways

- The shortest superstring problem asks for the shortest string containing every input string.
- A suffix-prefix overlap allows two strings to be merged without repeating their shared sequence.
- The greedy algorithm repeatedly merges the pair with the largest available overlap.
- Greedy merging is a heuristic: it does not guarantee the globally shortest result.
- The general shortest-superstring optimization problem is NP-hard.
- Runtime depends on both the number and length of the strings, as well as the overlap implementation.
- The problem is conceptually important in genome assembly, but real assembly pipelines must also handle sequencing errors, repeats, orientation, and coverage.