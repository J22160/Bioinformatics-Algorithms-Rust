# String Graphs for Genome Assembly

## What Is a String Graph?

A **string graph** is a graph-based representation of overlaps between sequencing reads, developed by Gene Myers in 2005 as a foundation for overlap-based genome assembly.

It is a refined representation of an overlap graph that removes redundant, transitive overlaps while preserving the essential relationships needed to reconstruct sequences.

To understand the motivation, imagine having thousands of overlapping newspaper clippings and wanting to reconstruct the original article. Many clippings may overlap directly, while other overlaps can already be inferred through intermediate clippings. A string graph removes these redundant connections, making the graph more compact and easier to analyze.

In genome assembly, the same principle applies to DNA reads:

- **Vertices** represent sequencing reads or their oriented representations.
- **Edges** represent valid suffix-prefix overlaps between reads.
- **Edge labels or overlap information** describe how the reads align and how their sequences can be combined.
- **Transitive reduction** removes overlaps that can be inferred through alternative paths.
- **Graph simplification** can subsequently collapse unambiguous paths into longer sequences.

The objective is not simply to create the smallest possible graph. It is to eliminate redundant overlap information without losing the sequence relationships needed for assembly.

### Why are string graphs useful?

Sequencing reads often overlap multiple other reads. If every detectable overlap becomes an edge, the graph can become dense, particularly in highly covered regions.

String graphs reduce this redundancy while retaining read-level information. They are especially useful for understanding overlap-based assembly with long reads, including PacBio and Oxford Nanopore sequencing.

However, string graphs do not automatically resolve every repeat or guarantee a unique genome reconstruction. Repeated sequences can still create branching paths and ambiguous connections.

---

## Overlap Graph vs. String Graph vs. de Bruijn Graph

These three graph representations describe sequence relationships differently.

| Property | Overlap Graph | String Graph | de Bruijn Graph |
|---|---|---|---|
| Basic representation | Sequencing reads | Sequencing reads with redundant overlaps removed | k-mers represented as edges or vertices, depending on convention |
| Edges | Detected read overlaps | Irreducible read overlaps | Connections between adjacent k-mers or their equivalent |
| Redundancy | Can be high | Reduced through transitive reduction | Depends on graph construction and simplification |
| Read length | Preserved | Preserved | Reads are decomposed into k-mers |
| Repeat handling | Can create many ambiguous overlaps | Can still contain repeat-induced branches | Can contain repeat-induced branches and collapsed repeated paths |
| Typical use | Conceptual overlap representation | Overlap-based assembly, particularly long reads | Widely used for short-read assembly and other sequencing applications |

**Key distinction:** A string graph is not simply an overlap graph with fewer nodes. Its main refinement is the removal of redundant overlap edges. Nodes still represent reads or oriented reads rather than individual k-mers.

de Bruijn graphs use a different representation. They decompose sequences into k-mers, which can reduce the need to compare every pair of complete reads.

Neither representation automatically solves repeat ambiguity. Their effectiveness depends on the read data, graph construction, and subsequent assembly algorithms.

---

## How String Graph Construction Works

A simplified string-graph construction workflow consists of four main stages.

### Step 1: Detect read overlaps

Given a set of sequencing reads, identify suffix-prefix overlaps between compatible read orientations.

Suppose we have two reads:

```text
Read X: ATGCTAGCTAGC
Read Y:       GCTAGCATCG
```

Their suffix-prefix overlap is:

```text
Read X: ATGCTAGCTAGC
Read Y:       GCTAGCATCG
                 ^^^^^^
```

The shared sequence is:

```text
GCTAGC
```

The overlap length is 6 bases.

The reads can therefore be merged into:

```text
ATGCTAGCTAGCATCG
```

The shared six bases are represented only once in the merged sequence.

For reads A and B, let $ov(A,B)$ denote the length of the selected suffix-prefix overlap. The length of the merged sequence is:

$$|A \oplus B| = |A| + |B| - ov(A,B)$$

In practice, overlap detection must account for read orientation, sequencing errors, overlap quality, and minimum overlap thresholds.

For DNA assembly, reads may need to be compared in both forward and reverse-complement orientations. Real overlap algorithms also commonly permit mismatches and insertions or deletions.

#### Efficient overlap detection

Naively comparing every pair of reads is expensive for large sequencing datasets. Indexing methods, such as suffix arrays and FM-indexes, can help identify candidate overlaps without performing exhaustive character-by-character comparisons for every possible pair.

The FM-index is based on the Burrows–Wheeler transform and supports efficient substring searches. However, the runtime of a complete overlap-detection pipeline depends on the index, search strategy, overlap criteria, read characteristics, and number of reported overlaps.

---

### Step 2: Remove contained reads

A read is **contained** if its entire sequence occurs within another read under the relevant orientation and alignment rules.

For example:

```text
Read A: ATGCCATGCCA
Read B:     GCCATGC
```

Read B is entirely contained within Read A.

If Read A already provides the sequence represented by Read B, retaining Read B as an independent sequence may be unnecessary for basic sequence reconstruction.

Removing contained reads can reduce the number of vertices and simplify the overlap structure.

However, containment handling in real assemblers requires care. Read quality, coverage, repeat structure, and the evidence supplied by individual reads can affect whether a read should be discarded or retained for downstream analysis.

---

### Step 3: Remove transitive edges

This is the defining operation of string-graph construction.

An overlap edge from X to Z is transitive when the overlap and sequence relationships implied by an alternative path through intermediate reads make the direct edge redundant.

Consider the following graph:

```text
Before transitive reduction:

X ───────→ Y ───────→ Z
 \___________________→ Z
        direct edge
```

Suppose the overlap between X and Z can be completely explained by the path from X through Y to Z. The direct X → Z edge does not contribute a new sequence relationship.

After removing that edge:

```text
After transitive reduction:

X ───────→ Y ───────→ Z
```

The path still connects X to Z, and the sequence implied by the path is consistent with the direct overlap.

#### Why is the edge redundant?

Assume the overlap between X and Y has length $o_{XY}$, and the overlap between Y and Z has length $o_{YZ}$.

The two-step path spells:

$$X \oplus Y \oplus Z$$

If the suffix-prefix relationship between X and Z is completely explained by this path, the direct X → Z edge adds no independent overlap information.

**Important:** Transitivity is not determined solely by the existence of any two-edge path. The overlap lengths, alignments, orientations, and sequence consistency must satisfy the string-graph definition. An arbitrary graph path is not enough to justify removing an edge.

Transitive reduction removes redundant edges; it does not remove the underlying reads simply because they are connected by a path.

#### How is transitive reduction performed?

A traditional approach is to construct an overlap graph and identify edges that can be explained through alternative paths. Efficient string-graph algorithms exploit overlap structure to avoid checking every possible path explicitly.

Some implementations identify irreducible overlaps directly during overlap search, reducing the need to construct all redundant edges first.

The precise algorithm and its complexity depend on the representation and implementation.

---

### Step 4: Simplify the reduced graph

After transitive reduction, some paths may consist of vertices with unambiguous connections.

For example:

```text
A → B → C → D
```

If the intermediate vertices have no competing connections and satisfy the assembler's criteria for simplification, the path can be represented more compactly as:

```text
A ─────────→ D
```

The simplified edge or path representation must preserve the sequence contributed by the intermediate reads.

This process is related to **unitig construction**: a unitig is a sequence represented by a non-branching path under a graph's particular definition and simplification rules.

A string graph and a unitig representation are not identical concepts. Transitive reduction produces a graph with redundant overlaps removed; unitig construction is a subsequent graph-simplification step that may collapse suitable chains.

Branching regions generally require additional reasoning because they can indicate repeats, alternative haplotypes, sequencing errors, or other ambiguities.

---

## Worked Example: Constructing a String Graph

Consider eight sequencing reads obtained from a hypothetical small plasmid.

```text
R0: ATGCCATGCC
R1: TGCCATGCCA
R2: GCCATGCCAG
R3: CCATGCCAGT
R4: CATGCCAGTC
R5: ATGCCAGTCG
R6: TGCCAGTCGA
R7: GCCAGTCGAT
```

These reads are deliberately constructed to illustrate transitive overlaps. They are a simplified teaching example, not a realistic noisy sequencing dataset.

### Step 1: Identify overlaps

Consider R0 and R1:

```text
R0: ATGCCATGCC
R1:  TGCCATGCCA
```

R0 has a suffix of length 9 matching the prefix of R1:

```text
TGCCATGCC
```

Therefore:

$$ov(R0,R1)=9$$

The same pattern continues for neighboring reads:

| Read pair | Overlap | Length |
|---|---|---:|
| R0 → R1 | TGCCATGCC | 9 |
| R1 → R2 | GCCATGCCA | 9 |
| R2 → R3 | CCATGCCAG | 9 |
| R3 → R4 | CATGCCAGT | 9 |
| R4 → R5 | ATGCCAGTC | 9 |
| R5 → R6 | TGCCAGTCG | 9 |
| R6 → R7 | GCCAGTCGA | 9 |

These overlaps form a chain:

```text
R0 → R1 → R2 → R3 → R4 → R5 → R6 → R7
```

### Step 2: Identify longer-range overlaps

The reads also have longer-range connections with smaller overlaps.

For example:

```text
R0: ATGCCATGCC
R2:   GCCATGCCAG
```

The suffix of R0 matches the prefix of R2:

```text
GCCATGCC
```

This is an overlap of length 8:

$$ov(R0,R2)=8$$

But R0 already overlaps R1 by 9 bases, and R1 overlaps R2 by 9 bases. The sequence implied by the path R0 → R1 → R2 is consistent with the direct R0 → R2 overlap.

Thus, under the exact-overlap assumptions of this example, R0 → R2 is transitive.

The same reasoning applies to longer-range connections such as R0 → R3, provided the selected overlap threshold and alignment rules admit those edges.

### Step 3: Build a simplified overlap graph

For clarity, consider only forward overlaps of at least 7 bases. This threshold excludes short accidental overlaps in this toy dataset.

The resulting overlap graph contains the adjacent edges as well as longer-range edges:

```text
R0 ──→ R1 ──→ R2 ──→ R3 ──→ R4 ──→ R5 ──→ R6 ──→ R7
│       │       │       │       │       │
├──────→R2      ├──────→R4      ├──────→R6
│               │               │
└──────────────→R3              └──────→R7
```

This is a schematic illustration rather than a complete rendering of every edge. In particular, every forward pair separated by one, two, or three positions meets the chosen overlap threshold.

### Step 4: Remove transitive edges

The direct R0 → R2 edge can be inferred from R0 → R1 → R2, so it is removed.

Likewise, the direct R0 → R3 edge can be inferred through R1 and R2. The same principle removes other longer-range edges whose overlaps are explained by intermediate paths.

After transitive reduction, the simplified graph is:

```text
R0 → R1 → R2 → R3 → R4 → R5 → R6 → R7
```

The remaining chain retains the intended read order while eliminating redundant long-range connections.

### Step 5: Reconstruct the sequence

Merge the reads along the chain, counting each 9-base overlap only once.

The resulting sequence is:

```text
ATGCCATGCCAGTCGAT
```

Its length is:

$$10 + 7(10-9)=17$$

The first read contributes 10 bases, and each of the remaining seven reads contributes one new base.

This reconstruction illustrates how a string graph can preserve the sequence relationships needed for assembly while removing redundant overlap edges.

**Limitations of this example:** The reads are intentionally designed with a clear ordering and exact overlaps. Real datasets may contain reverse-complement overlaps, sequencing errors, repeated regions, and competing paths. A transitive reduction alone does not guarantee that a real plasmid or genome can be reconstructed unambiguously.

---

## Complexity and Scalability

Let:

- $N$ be the total number of bases across all reads.
- $n$ be the number of reads.
- $E$ be the number of candidate overlap edges.
- $C$ be the number of reported overlaps or overlap candidates, depending on the algorithm's definition.

### Pairwise overlap detection

A naive method compares all pairs of reads. With $n$ reads of typical length $L$, there are $O(n^2)$ read pairs. The cost per pair depends on the overlap algorithm.

This approach becomes expensive for large long-read datasets.

### Indexed overlap detection

FM-index-based methods can accelerate overlap discovery by searching an index instead of comparing every read pair independently.

A frequently quoted bound for particular algorithms is $O(N+C)$, where the cost of reporting overlaps is included. Other methods can have different bounds depending on the indexing scheme, search strategy, and overlap model.

It is not correct to assume that every FM-index-based overlap detector, or every complete string-graph construction pipeline, runs in $O(N)$ time regardless of its output.

### Transitive reduction

The cost of transitive reduction depends on the graph structure and algorithm. Specialized string-graph algorithms can exploit the properties of read overlaps to perform this operation efficiently.

The general graph problem and a specialized string-graph implementation should not be assigned the same complexity without specifying the algorithm and its assumptions.

### Practical takeaway

Efficient genome assemblers use indexed overlap detection, careful filtering, and specialized graph algorithms to control computational costs. The number of reads, read lengths, sequencing error rates, coverage, and repeat content all influence runtime and memory use.

---

## String Graphs in Modern Genome Assembly

String graphs are particularly relevant to overlap-layout-consensus (OLC) assembly and long-read assembly.

The traditional OLC workflow consists of:

1. **Overlap:** Identify overlaps between reads.
2. **Layout:** Infer a consistent arrangement of reads, often using an overlap graph or string graph.
3. **Consensus:** Generate a more accurate sequence from the aligned reads.

Modern assemblers often combine these ideas with additional algorithms for error correction, repeat resolution, haplotype separation, and graph simplification.

### Examples of assemblers

| Assembler | Typical data | Relevant characteristics |
|---|---|---|
| SGA | Short reads | Uses FM-index-based methods and string-graph concepts |
| FALCON | PacBio CLR long reads | Uses overlap-based assembly and string-graph methods |
| hifiasm | PacBio HiFi reads | Uses string-graph-related representations and advanced graph processing for high-quality, often haplotype-resolved assembly |
| Flye | Long reads, including noisy reads | Uses a repeat graph rather than being simply classified as a conventional string-graph assembler |

These tools differ substantially in their algorithms. It is more accurate to say that they use related overlap or graph-based approaches than to treat them as interchangeable implementations of the same string-graph algorithm.

---

## Applications in Bioinformatics

### Whole-genome sequencing

Long reads can span genomic regions that are difficult to reconstruct from short reads alone. Overlap-based graphs help establish the order and orientation of reads and can support assembly across some repeated regions.

### Plasmid assembly

For bacterial sequencing, overlap graphs can help connect reads across plasmid regions and identify candidate plasmid contigs. Circularity must be supported by appropriate sequence and graph evidence; a graph path or cycle alone is not proof that a molecule is circular.

### Metagenomics

Microbial communities contain DNA from multiple organisms and often closely related strains. Assembly graphs can expose alternative paths caused by strain variation, repeats, and shared genomic regions.

However, string graphs do not inherently resolve strain mixtures. Additional coverage, sequence variation, and graph-analysis methods may be required.

### Structural variation and haplotype analysis

Graph structures can reveal alternative paths corresponding to insertions, deletions, inversions, or other genomic differences. Haplotype-aware assembly methods may preserve biologically meaningful alternatives rather than collapsing them into a single consensus.

The graph alone does not establish the biological interpretation of a branch; that requires evidence from the reads and downstream analysis.

---

## Potential Improvements to a String-Graph Pipeline

### 1. Accelerate overlap discovery

Use suitable indexing and filtering methods to reduce unnecessary pairwise comparisons. FM-indexes, suffix arrays, minimizers, and other sequence-indexing approaches can help identify candidate overlaps.

### 2. Detect irreducible overlaps efficiently

Where supported by the algorithm, identify overlaps that cannot be explained by intermediate reads during overlap discovery. This can reduce the cost of constructing a dense graph and then removing transitive edges.

### 3. Handle sequencing errors

Exact string matching is insufficient for most real sequencing datasets. Error-tolerant overlap detection uses alignment scores, identity thresholds, and minimum overlap lengths to distinguish plausible overlaps from weak matches.

Thresholds should reflect the sequencing technology and expected error profile.

### 4. Preserve haplotype information

Diploid genomes may contain alternative sequences inherited from different parental chromosomes. Haplotype-aware assembly methods aim to preserve these alternatives when the reads provide enough evidence, rather than incorrectly collapsing them into one sequence.

### 5. Simplify graphs carefully

Remove unsupported or redundant connections while preserving meaningful branches. Aggressive simplification can erase real variation or create incorrect assemblies, particularly in repetitive regions.

---

## Key Takeaways

- A string graph represents overlaps between sequencing reads while removing redundant, transitive overlap edges.
- Transitive reduction is the defining distinction between a conventional overlap graph and a string graph.
- Contained-read removal and unitig construction are related operations, but they are not the same as transitive reduction.
- Efficient overlap detection can use indexing methods, including FM-index-based approaches, but complexity claims depend on the implementation and output size.
- String graphs preserve read-level sequence information, which is useful for overlap-based long-read assembly.
- String graphs can still contain repeat-induced ambiguity and do not automatically resolve repeats or determine genome circularity.
- SGA and FALCON are examples of assemblers associated with string-graph approaches; hifiasm uses related graph methods for high-quality and haplotype-aware assembly.
- Understanding string graphs provides a foundation for learning how sequencing reads are organized into contigs and how modern genome assemblers represent uncertainty.