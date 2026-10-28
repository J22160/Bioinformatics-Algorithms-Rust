# Overlap Graphs for Genome Assembly

## What Is an Overlap Graph?

An **overlap graph** is a directed graph used to represent relationships between sequencing reads. Each vertex represents a read, and a directed edge connects two reads when the suffix of one read overlaps the prefix of another.

Overlap graphs are an important concept in genome assembly, particularly in the **overlap-layout-consensus (OLC)** approach, which is widely used in long-read sequencing workflows.

The central idea is straightforward: if two reads originate from overlapping regions of the same DNA molecule, their shared sequence can be used to determine their relative order.

For example, consider these reads:

```text
Read A: ATGCTAG
Read B:     CTAGGCA
```

The suffix of Read A, `CTAG`, matches the prefix of Read B.

```text
Read A: ATGCTAG
              ||||
Read B:     CTAGGCA
```

This creates a directed edge from Read A to Read B:

```text
Read A ──overlap: CTAG──> Read B
```

The reads can then be merged into a longer sequence:

```text
ATGCTAGGCA
```

The overlapping bases are included only once.

**Key idea:** An overlap graph converts read-to-read sequence similarities into a graph whose paths can help reconstruct the original DNA molecule.

---

## Why Are Overlap Graphs Useful?

DNA sequencing produces many shorter reads from different regions of a genome. The assembler must determine which reads overlap, how they should be ordered, and how their sequence information can be combined.

For example, imagine receiving these reads:

```text
Read 1: ATGCCA
Read 2: GCCATG
Read 3: CATGCC
```

The reads share sequence regions, but their correct order is not explicitly provided.

By identifying suffix–prefix overlaps, we can infer connections between the reads and construct a graph that represents possible assembly arrangements.

Overlap graphs are particularly useful when reads are long enough to provide substantial overlap information. Long reads can span repetitive regions and connect distant genomic regions that short reads may not resolve.

The overlap-layout-consensus approach separates the problem into three broad stages:

1. **Overlap:** Identify reads that share sufficiently strong sequence overlaps.
2. **Layout:** Determine the order and orientation of reads using the overlap evidence.
3. **Consensus:** Combine the aligned reads to estimate the underlying DNA sequence.

Real assemblers often iterate between these stages and use additional evidence to resolve ambiguity.

---

## How Does an Overlap Graph Work?

### Step 1: Identify suffix–prefix overlaps

Suppose we have two reads:

```text
Read A: ATGCCA
Read B:     CCATGC
```

The suffix of Read A is `CCA`, which matches the prefix of Read B.

```text
Read A: ATGCCA
              |||
Read B:     CCATGC
```

The overlap length is 3 bases.

We therefore create a directed edge:

```text
Read A ──overlap length: 3──> Read B
```

An overlap graph usually considers overlaps that meet a minimum length and quality or alignment-score threshold. The exact criteria depend on the assembler.

In real sequencing data, reads may contain errors, so overlaps do not always match perfectly. Overlap detection commonly uses sequence alignment or specialized algorithms to identify approximate matches.

### Step 2: Create graph vertices

Each read becomes a vertex.

For example:

```text
R1: ATGCCA
R2: GCCATG
R3: CATGCC
```

The graph initially contains three vertices:

```text
R1       R2       R3
```

### Step 3: Create directed edges

For every suitable pair of reads, determine whether the suffix of one read overlaps the prefix of the other.

If the overlap is valid, add a directed edge.

```text
R1 ──> R2
```

An edge may store information such as:

- Overlap length.
- Alignment score.
- Number of mismatches.
- Number of inserted or deleted bases.
- Estimated relative offset between reads.
- Read orientation.

These attributes help the assembler distinguish strong overlaps from unreliable ones.

### Step 4: Determine the read layout

The assembler uses the overlap graph to infer a consistent ordering of reads.

A simplified example is:

```text
R1 ──> R2 ──> R3 ──> R4
```

If each consecutive pair overlaps correctly, the reads can be merged into a longer sequence.

However, real graphs often contain branches caused by repetitive DNA, sequencing errors, heterozygosity, or alternative genomic paths. Consequently, the assembler must determine which connections are supported by the data.

### Step 5: Generate a consensus sequence

Once the read layout has been established, overlapping reads are aligned and combined.

In ideal error-free data, merging matching overlaps may be sufficient. In real data, multiple reads often disagree at individual positions because of sequencing errors or biological variation.

Consensus algorithms use the aligned evidence to estimate the most likely underlying sequence.

For example:

```text
Read A: ATGCTAG
Read B:     CTAGGCA
Read C:     CTAGGCA
```

The shared region supports the overlap, while the additional bases extend the reconstructed sequence.

A production assembler may use quality scores, read-specific error models, and multiple alignments to generate a more accurate consensus.

---

## Worked Example: Building an Overlap Graph

Consider these five reads from a hypothetical circular DNA molecule:

```text
R1: ATGCCA
R2: GCCATG
R3: CATGCC
R4: TGCCAT
R5: CCATGC
```

For this simplified example, assume that:

- The reads are error-free.
- Reads are already in the correct orientation.
- An overlap of at least 3 matching bases is accepted.
- We are looking for exact suffix–prefix overlaps.

These assumptions make the graph easier to understand than a real sequencing dataset.

### Step 1: Find the overlaps

Compare the suffix of each read with the prefix of the other reads.

| From read | To read | Suffix–prefix overlap | Length |
|---|---|---|---:|
| R1: ATGCCA | R5: CCATGC | CCA | 3 |
| R2: GCCATG | R1: ATGCCA | ATG | 3 |
| R3: CATGCC | R2: GCCATG | GCC | 3 |
| R4: TGCCAT | R3: CATGCC | CAT | 3 |
| R5: CCATGC | R4: TGCCAT | TGC | 3 |

For example, the suffix `GCC` of R3 matches the prefix `GCC` of R2. This creates an edge from R3 to R2.

### Step 2: Construct the graph

The corresponding directed edges are:

```text
R3 ──GCC──> R2
R2 ──ATG──> R1
R1 ──CCA──> R5
R5 ──TGC──> R4
R4 ──CAT──> R3
```

The graph forms a directed cycle:

```text
        ┌────────────────────────────┐
        │                            ▼
       R3 ──> R2 ──> R1 ──> R5 ──> R4
        ▲                            │
        └────────────────────────────┘
```

Each read has one incoming and one outgoing edge in this simplified graph.

### Step 3: Determine the assembly order

One valid cyclic order is:

```text
R3 → R2 → R1 → R5 → R4 → R3
```

Because the final read connects back to the first read, the graph supports a circular arrangement under the assumptions of this example.

### Step 4: Merge the reads

Begin with R3:

```text
CATGCC
```

R2 overlaps the last three bases, `GCC`, and contributes the remaining bases `ATG`:

```text
CATGCCATG
```

R1 overlaps the suffix `ATG` and contributes `CCA`:

```text
CATGCCATGCCA
```

R5 overlaps the suffix `CCA` and contributes `TGC`:

```text
CATGCCATGCCATGC
```

R4 overlaps the suffix `TGC` and contributes `CAT`, which crosses the proposed circular junction.

The cyclic sequence can therefore be represented as:

```text
CATGCCATGCCATGC
```

This 15-base circular sequence produces the five 6-base reads when sampled at offsets 0, 3, 6, 9, and 12, with reads allowed to wrap around the circular junction.

**Result:** The reads are consistent with the circular sequence `CATGCCATGCCATGC`.

This is a toy reconstruction, not proof that an unknown molecule is circular. In real assembly, circularity must be supported by reliable overlaps and evidence spanning the proposed junction. The same reads could otherwise be consistent with a longer linear sequence.

---

## What Is the Hamiltonian Path Connection?

A Hamiltonian path is a path that visits every vertex exactly once.

In a simplified overlap-graph formulation of genome assembly:

- Each vertex represents a read.
- Each edge represents a compatible overlap.
- A candidate assembly order corresponds to a path through the reads.

If a path visits every read exactly once, it is Hamiltonian.

This gives an intuitive connection between overlap graphs and genome assembly: the assembler is trying to arrange reads into a consistent sequence of overlaps.

However, **real OLC assemblers do not generally solve the unrestricted Hamiltonian path problem directly**. Sequencing reads may contain errors, repeated regions, redundant coverage, and alternative alignments. Practical assemblers use filtering, graph simplification, repeat handling, layout heuristics, and consensus methods.

The Hamiltonian path formulation is useful for understanding the theoretical problem, but it is not a complete description of modern assembly algorithms.

---

## Overlap Graphs Versus de Bruijn Graphs

Overlap graphs and de Bruijn graphs both represent sequence relationships, but they use different graph structures.

| Property | Overlap graph | de Bruijn graph |
|---|---|---|
| Vertices | Sequencing reads | Usually $(k-1)$-mers |
| Edges | Overlaps between reads | K-mers connecting prefix and suffix |
| Typical assembly formulation | Read layout using overlap evidence | Sequence reconstruction using k-mer graph traversal |
| Main computational challenge | Finding and resolving read overlaps | Building, simplifying, and resolving the k-mer graph |
| Repeats | Can cause ambiguous read connections | Can create shared paths and branching |
| Common association | OLC assembly | de Bruijn graph assembly |

### Why does the distinction matter?

In an overlap graph, the graph size is closely related to the number of reads and the candidate overlaps between them.

In a de Bruijn graph, the graph is built from k-mer observations. Identical k-mers can be combined into shared graph elements, which can be advantageous when processing large volumes of short-read data.

Neither representation eliminates repeat ambiguity. Both require careful handling of errors, coverage, and alternative paths.

---

## Computational Complexity

Let:

- $n$ be the number of reads.
- $L$ be a representative read length.
- $O$ be the number of candidate overlap edges.

### Overlap computation

A naive algorithm compares every read against every other read.

The number of ordered read pairs is:

$$n(n-1) = O(n^2)$$

If each comparison uses an alignment method whose cost is $O(L^2)$ for reads of length $L$, a naive all-pairs alignment approach can cost:

$$O(n^2L^2)$$

The exact cost depends on the alignment algorithm, read lengths, and whether the algorithm can exploit the suffix–prefix structure.

For a simplified fixed-length exact-overlap comparison that examines at most $O(L)$ positions per pair, the cost is:

$$O(n^2L)$$

Therefore, $O(n^2L)$ is a reasonable bound for a naive pairwise overlap procedure with linear-time comparisons, but it is **not a universal bound for all overlap detection algorithms**.

Modern assemblers use techniques such as minimizers, k-mer indexing, hashing, filtering, and specialized alignment algorithms to avoid performing expensive comparisons for every read pair.

### Graph construction

Once the overlaps have been identified, creating the corresponding vertices and edges requires approximately:

$$O(n+O)$$

time, assuming efficient data structures.

Here, $O$ denotes the number of retained overlap edges.

### Hamiltonian path search

Finding a Hamiltonian path in a general graph is NP-complete as a decision problem. An exhaustive search can require exponential time in the number of vertices.

However, practical genome assemblers do not generally enumerate all Hamiltonian paths. They use biological constraints and heuristics to reduce the search space.

### Consensus generation

Consensus computation depends on the number of reads in each layout, their lengths, the alignment representation, and the consensus algorithm. It cannot be summarized by a single complexity bound for all assemblers.

### Summary

| Stage | Simplified computational cost |
|---|---|
| Naive overlap detection | $O(n^2L)$ with linear-time comparison per pair |
| More expensive pairwise alignment | Potentially $O(n^2L^2)$ under conventional dynamic programming |
| Graph construction after overlaps | $O(n+O)$ |
| General Hamiltonian path search | NP-complete in general |
| Practical assembly | Depends on indexing, filtering, graph structure, heuristics, and consensus methods |

---

## Challenges in Real Overlap Graphs

Real sequencing data are much more complicated than the error-free reads in the worked example.

### Sequencing errors

Errors can prevent true overlaps from matching exactly or create false overlaps between unrelated reads.

Assemblers therefore use approximate alignment, quality information, and overlap scoring to distinguish likely true overlaps from noise.

### Repetitive sequences

If a genome contains repeated regions, a read may overlap several different reads equally well.

This creates branching paths in the graph and makes it difficult to determine which regions belong next to one another.

Long reads can help resolve these ambiguities when they span the repeated region and reach unique sequence on either side.

### Redundant coverage

Many reads may originate from the same genomic region. The assembler must avoid interpreting redundant reads as distinct consecutive regions of the genome.

Read-layout algorithms account for overlapping evidence and coverage rather than simply concatenating every read.

### Reverse-complement orientation

DNA is double-stranded, and sequencing reads may represent either orientation of a genomic region.

A real overlap detection pipeline must account for reverse-complement relationships when appropriate, rather than assuming every read is already in the correct orientation.

### Chimeric reads and contamination

A read that incorrectly joins unrelated genomic regions can create a misleading edge. Contamination can also produce overlap structures unrelated to the target genome.

Assemblers use filtering and validation to identify suspicious connections.

---

## Modern Overlap-Based Assemblers

Several widely used assemblers employ overlap-based or related long-read assembly strategies.

| Assembler | General use |
|---|---|
| Canu | Long-read assembly, with correction and overlap-based processing |
| Hifiasm | High-quality PacBio HiFi genome assembly, including approaches for resolving haplotypes |
| Flye | Long-read genome assembly using repeat graphs and specialized graph-processing strategies |
| Celera Assembler | Historically important whole-genome assembly system based on overlap-layout-consensus concepts |

These assemblers are not identical implementations of a simple overlap graph. They incorporate specialized overlap detection, error handling, graph construction, repeat resolution, and consensus algorithms.

The Celera project was one of the major efforts in producing an early draft of the human genome, alongside the public Human Genome Project. It is more accurate to describe it as one of the landmark human-genome assembly efforts than to attribute the first human genome assembly to Celera alone.

---

## When Should You Use an Overlap-Based Approach?

Overlap-based methods are especially relevant when:

- Reads are long enough to span repetitive regions.
- Read-to-read alignment provides useful information about genomic layout.
- The sequencing technology has an error profile that the assembler can model.
- The goal is to reconstruct long contiguous sequences or resolve genome structure.

For large short-read datasets, de Bruijn graph approaches have historically been popular because k-mer graphs can represent shared sequence information efficiently.

For long-read datasets, overlap-based methods and related repeat-graph approaches are often more suitable because individual reads can provide longer-range sequence context.

The best choice depends on read length, error rate, sequencing coverage, genome complexity, and the intended assembly output.

---

## Why Overlap Graphs Matter

Overlap graphs demonstrate how local sequence similarities can be converted into a global reconstruction problem.

Their central contributions are:

- Representing reads as vertices and overlaps as directed edges.
- Inferring the relative order and orientation of reads.
- Supporting layout construction and consensus generation.
- Providing a foundation for overlap-layout-consensus assembly methods.
- Illustrating the computational challenges caused by repeats and ambiguous connections.

The most important lesson is that **finding overlaps is only the first step**. A successful assembler must also decide which overlaps are reliable, resolve competing paths, account for sequencing errors, and construct a consensus sequence.

For bioinformatics developers, overlap graphs are a useful foundation for learning sequence alignment, graph representations, indexing, approximate string matching, and assembly algorithms. They also provide an instructive contrast with de Bruijn graphs and Eulerian traversal.