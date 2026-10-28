# Eulerian Paths and Cycles: Graph Algorithms for Bioinformatics

## What Are Eulerian Paths and Cycles?

Imagine a city containing islands connected by bridges. You want to take a walk that crosses every bridge exactly once.

Two questions arise:

- Can you cross every bridge exactly once?
- Can you do so and return to your starting point?

These questions lead to two fundamental graph concepts: **Eulerian paths** and **Eulerian cycles**.

An **Eulerian path** is a walk through a graph that traverses every edge exactly once.

An **Eulerian cycle** is an Eulerian path that returns to its starting vertex.

The important distinction is that Eulerian traversal focuses on edges, not vertices. A walk may visit the same vertex multiple times, provided that it uses every edge exactly once.

These concepts are named after the mathematician Leonhard Euler, whose work on the Seven Bridges of Königsberg problem helped establish graph theory.

### Eulerian Path Versus Eulerian Cycle

| Property | Eulerian Path | Eulerian Cycle |
|---|---|---|
| Uses every edge exactly once | Yes | Yes |
| Must return to starting vertex | No | Yes |
| Can revisit vertices | Yes | Yes |
| Can start and end at different vertices | Yes | No |

For example, consider this directed graph:

```text
A ──> B ──> C ──> D
```

The walk `A → B → C → D` is an Eulerian path because it traverses all three edges exactly once.

It is not an Eulerian cycle because it ends at `D`, not `A`.

Now consider:

```text
A ──> B
^     |
|     v
D <── C
```

The walk `A → B → C → D → A` is an Eulerian cycle because it traverses every edge exactly once and returns to its starting vertex.

---

## Understanding Graph Terminology

Before studying Eulerian traversal, it helps to understand a few basic graph concepts.

### Vertices and Edges

A graph consists of vertices connected by edges.

```text
A ──> B ──> C
```

Here:

- `A`, `B`, and `C` are vertices.
- `A → B` and `B → C` are directed edges.

### Directed and Undirected Graphs

In an undirected graph, an edge can be traversed in either direction:

```text
A ─── B
```

In a directed graph, an edge has a specified direction:

```text
A ──> B
```

The direction matters when determining whether an Eulerian path or cycle exists.

### In-degree and Out-degree

For a directed graph:

- **In-degree** is the number of edges entering a vertex.
- **Out-degree** is the number of edges leaving a vertex.

Consider:

```text
A ──> B ──> C
      ^
      |
      D
```

The vertex `B` has two incoming edges, from `A` and `D`, and one outgoing edge, to `C`.

Therefore:

$$\deg^-(B) = 2$$

$$\deg^+(B) = 1$$

where $\deg^-(v)$ denotes in-degree and $\deg^+(v)$ denotes out-degree.

These degree values determine whether a directed graph can have an Eulerian path or cycle.

---

## Conditions for an Eulerian Path or Cycle

Not every graph admits an Eulerian traversal. The graph must satisfy specific degree and connectivity conditions.

### Undirected Graphs

For a connected undirected graph, ignoring isolated vertices:

**Eulerian Cycle**

Every vertex must have an even degree.

$$\deg(v) \equiv 0 \pmod 2$$

**Eulerian Path**

The graph must have either:

- Zero vertices of odd degree, in which case an Eulerian cycle exists; or
- Exactly two vertices of odd degree, which can serve as the start and end of an open Eulerian path.

If more than two vertices have odd degree, no Eulerian path exists.

### Directed Graphs

For a directed graph, the conditions are based on in-degree and out-degree.

**Eulerian Cycle**

Every vertex must be balanced:

$$\deg^+(v) = \deg^-(v)$$

In addition, all vertices incident to edges must satisfy the appropriate connectivity condition. For a finite directed graph, balance at every vertex together with connectivity of the underlying undirected graph over the non-isolated vertices is sufficient for an Eulerian cycle.

**Eulerian Path**

For a directed Eulerian path that does not form a cycle:

- One vertex has out-degree exactly one greater than its in-degree.
- One vertex has in-degree exactly one greater than its out-degree.
- All other vertices have equal in-degree and out-degree.
- The graph must have the appropriate connectivity.

The start vertex satisfies:

$$\deg^+(s) = \deg^-(s) + 1$$

The end vertex satisfies:

$$\deg^-(t) = \deg^+(t) + 1$$

All other vertices satisfy:

$$\deg^+(v) = \deg^-(v)$$

These conditions are necessary and, with the appropriate connectivity, sufficient.

### Example

Consider:

```text
A ──> B ──> C ──> D
```

The degrees are:

| Vertex | In-degree | Out-degree | Difference |
|---|---:|---:|---|
| A | 0 | 1 | Out-degree exceeds in-degree by 1 |
| B | 1 | 1 | Balanced |
| C | 1 | 1 | Balanced |
| D | 1 | 0 | In-degree exceeds out-degree by 1 |

The graph has an Eulerian path from `A` to `D`, but it does not have an Eulerian cycle.

---

## Hierholzer's Algorithm

**Hierholzer's algorithm** finds an Eulerian cycle in a graph that satisfies the necessary conditions. A small modification also handles Eulerian paths.

The algorithm works by constructing a cycle, identifying vertices with unused edges, and joining additional cycles into the existing traversal.

### Algorithm Intuition

Suppose we have the following directed graph:

```text
A ──> B ──> C
^           |
|           v
└───────────D
```

A possible traversal is:

```text
A → B → C → D → A
```

This traversal uses every edge exactly once and returns to `A`.

For more complicated graphs, there may be several cycles that share vertices. Hierholzer's algorithm combines them into one Eulerian cycle.

### Step-by-Step Procedure

1. Start at a valid vertex with an unused outgoing edge.
2. Follow unused outgoing edges until the current walk cannot continue.
3. If the walk has returned to a valid completed cycle and unused edges remain, find a vertex on the current traversal with an unused outgoing edge.
4. Start another walk from that vertex, using only unused edges.
5. Insert the new walk into the existing traversal at that vertex.
6. Repeat until every edge has been used exactly once.

For a directed Eulerian cycle, degree balance ensures that a walk cannot become stranded at a different vertex before returning to its starting point, provided the graph satisfies the required conditions.

For an open Eulerian path, begin at the vertex whose out-degree exceeds its in-degree by one. The traversal then ends at the vertex whose in-degree exceeds its out-degree by one.

### Pseudocode

```text
function Hierholzer(graph, start):
    stack = [start]
    circuit = []

    while stack is not empty:
        v = top(stack)

        if v has an unused outgoing edge:
            edge = choose an unused outgoing edge from v
            mark edge as used
            push destination(edge) onto stack
        else:
            append pop(stack) to circuit

    reverse(circuit)
    return circuit
```

For an implementation that returns the actual edge sequence, it is useful to store edge identifiers and record which edges are consumed. This matters when multiple edges connect the same vertices.

### Why Does the Algorithm Work?

The algorithm exploits the degree structure of an Eulerian graph.

When every vertex is balanced, a walk that follows unused edges cannot terminate at a different vertex while leaving a valid closed traversal incomplete. If unused edges remain, connectivity ensures that another unused portion can be reached from a vertex already in the current circuit.

By joining these portions together, the algorithm eventually constructs a traversal containing every edge exactly once.

---

## Worked Example: Finding an Eulerian Path

Consider this directed graph:

```text
A ──> B ──> C ──> A
      |
      v
      D ──> B
```

Its edges are:

```text
A → B
B → C
C → A
B → D
D → B
```

### Step 1: Calculate the Degrees

| Vertex | In-degree | Out-degree |
|---|---:|---:|
| A | 1 | 1 |
| B | 2 | 2 |
| C | 1 | 1 |
| D | 1 | 1 |

Every vertex is balanced. The graph is connected in the required sense, so an Eulerian cycle exists.

### Step 2: Find a Cycle

Start at `A`:

```text
A → B → C → A
```

This uses three edges.

Two edges remain:

```text
B → D
D → B
```

The vertex `B` is part of the existing cycle and has unused outgoing edges.

### Step 3: Find Another Cycle

Starting from `B`, follow the unused edges:

```text
B → D → B
```

This is another cycle.

### Step 4: Combine the Cycles

Insert the second cycle into the first at vertex `B`:

```text
A → B → D → B → C → A
```

Check the edges:

```text
A → B
B → D
D → B
B → C
C → A
```

Every edge appears exactly once.

Therefore, the result is an Eulerian cycle.

---

## Why Eulerian Paths Matter in Genome Assembly

One of the most important applications of Eulerian traversal is **de Bruijn graph-based genome assembly**.

DNA sequencing does not necessarily produce one continuous genome-length sequence. Instead, it generates many shorter reads sampled from different regions of the genome.

An assembler must reconstruct a longer sequence from these fragments.

A de Bruijn graph represents k-mers as edges and their $(k-1)$-mer prefixes and suffixes as vertices. Under the classical formulation, reconstructing a sequence consistent with the k-mer spectrum can be expressed as finding an Eulerian traversal through the graph.

### Example

Suppose the observed DNA k-mers are:

```text
ATG
TGC
GCA
```

For $k=3$, each k-mer becomes an edge between its 2-base prefix and 2-base suffix.

| K-mer | Prefix | Suffix | Edge |
|---|---|---|---|
| ATG | AT | TG | AT → TG |
| TGC | TG | GC | TG → GC |
| GCA | GC | CA | GC → CA |

The graph is:

```text
AT ──ATG──> TG ──TGC──> GC ──GCA──> CA
```

An Eulerian path traverses all three edges:

```text
AT → TG → GC → CA
```

To reconstruct the sequence, begin with the first vertex `AT` and append the final base of each edge:

```text
AT + G + C + A = ATGCA
```

The reconstructed sequence is:

```text
ATGCA
```

Its 3-mers are:

```text
ATG, TGC, GCA
```

These match the supplied spectrum.

### Why Use Edges Instead of Vertices?

In the de Bruijn graph, each edge represents a k-mer observation. The assembly problem therefore requires using the k-mer edges in a consistent order.

Eulerian traversal is a natural fit because it uses every edge once in the classical k-mer-spectrum formulation.

In an overlap graph, by contrast, vertices commonly represent reads, and the task may be formulated as finding a path that visits each read vertex once. This is related to the Hamiltonian path problem.

---

## Worked Bioinformatics Problem: Reconstruct a Sequence from a 3-mer Spectrum

Suppose sequencing-by-hybridization produces the following 3-mer spectrum:

```text
ATG, TGC, GCA, CAG, AGT, GTG, TGT
```

The goal is to find a DNA sequence consistent with these k-mers.

For this example, assume the spectrum contains each k-mer occurrence exactly once. In real data, multiplicity and sequencing errors must be considered.

### Step 1: Build the De Bruijn Graph

Each 3-mer creates a directed edge between its 2-base prefix and 2-base suffix.

| 3-mer | Prefix Node | Suffix Node |
|---|---|---|
| ATG | AT | TG |
| TGC | TG | GC |
| GCA | GC | CA |
| CAG | CA | AG |
| AGT | AG | GT |
| GTG | GT | TG |
| TGT | TG | GT |

The graph is:

```text
AT ──ATG──> TG ──TGC──> GC ──GCA──> CA
             ^                       |
             |                       v
             GT <──AGT── AG <──CAG───┘
             |
             └──GTG──> TG

TG ──TGT──> GT
```

The diagram is schematic; the edge table above defines the exact graph.

### Step 2: Check the Degree Conditions

| Node | In-degree | Out-degree |
|---|---:|---:|
| AT | 0 | 1 |
| TG | 2 | 2 |
| GC | 1 | 1 |
| CA | 1 | 1 |
| AG | 1 | 1 |
| GT | 2 | 1 |

The vertex `AT` has out-degree one greater than in-degree, so it is the start.

The vertex `GT` has in-degree one greater than out-degree, so it is the end.

All other vertices are balanced. The graph has an Eulerian path.

### Step 3: Find an Eulerian Path

One valid edge order is:

```text
ATG → TGC → GCA → CAG → AGT → GTG → TGT
```

Check the overlaps:

```text
ATG
 TGC
  GCA
   CAG
    AGT
     GTG
      TGT
```

Each consecutive pair of k-mers shares the required 2-base overlap.

### Step 4: Reconstruct the Sequence

Start with the first k-mer:

```text
ATG
```

Append the final base from each subsequent k-mer:

```text
ATG + C + A + G + T + G + T
```

The reconstructed sequence is:

```text
ATGCAGTGT
```

### Step 5: Validate the Reconstruction

Extract the 3-mers from `ATGCAGTGT`:

```text
ATG
TGC
GCA
CAG
AGT
GTG
TGT
```

These match the input spectrum exactly.

**Result:** `ATGCAGTGT` is a valid reconstruction for the given 3-mer spectrum.

This example demonstrates the mathematical relationship between Eulerian paths and sequence reconstruction. Real genome assembly is more complicated because observed k-mers can be missing, repeated, or erroneous, and multiple valid reconstructions may exist.

---

## Eulerian Path Versus Hamiltonian Path

Eulerian and Hamiltonian paths solve different graph traversal problems.

An Eulerian path traverses every edge exactly once. A Hamiltonian path visits every vertex exactly once.

| Property | Eulerian Path | Hamiltonian Path |
|---|---|---|
| Main objective | Use every edge exactly once | Visit every vertex exactly once |
| May revisit vertices? | Yes | No |
| May leave edges unused? | No | Yes, edges are not the traversal target |
| General computational difficulty | Polynomial-time algorithm | NP-complete decision problem in general |
| Typical genome-assembly connection | de Bruijn graph formulation | Classical overlap-graph formulation |

### Why Is This Distinction Important?

Finding an Eulerian path can be solved in $O(V+E)$ time using standard graph representations and algorithms.

The general Hamiltonian path problem is NP-complete: no polynomial-time algorithm is known for all possible inputs, and finding a Hamiltonian path in a general graph is computationally difficult.

This difference helped motivate de Bruijn graph approaches to sequence assembly.

However, the historical development of genome assembly is more nuanced than a simple replacement of Hamiltonian paths with Eulerian paths. Practical assemblers must still resolve repeats, correct errors, manage coverage, and choose among ambiguous paths. Efficient Eulerian traversal solves the graph traversal subproblem, not every assembly challenge.

---

## Computational Complexity

Let:

- $V$ be the number of vertices.
- $E$ be the number of edges.

Hierholzer's algorithm runs in:

$$O(V+E)$$

time, assuming an adjacency-list representation and efficient tracking of unused edges.

For a graph containing many edges but relatively few vertices, the runtime is often described informally as linear in the number of edges:

$$O(E)$$

This simplification assumes the relevant vertices are already represented and that $V$ is bounded by the graph size.

### Space Complexity

The graph itself requires:

$$O(V+E)$$

space.

The traversal stack and output circuit can require up to $O(E)$ additional space.

Thus, the overall auxiliary storage including the graph and output is:

$$O(V+E)$$

### Important Distinction: Traversal Versus Assembly

For genome assembly, total runtime includes more than Eulerian traversal:

1. K-mer generation.
2. Graph construction.
3. Graph cleaning and simplification.
4. Repeat resolution.
5. Graph traversal and contig generation.
6. Validation and scaffolding.

Eulerian traversal is efficient, but the complete assembly pipeline can be substantially more expensive and may require heuristics or additional biological evidence.

---

## Common Implementation Considerations

When implementing Eulerian traversal in a bioinformatics tool, several details matter.

### Track Individual Edges

Two edges may connect the same pair of vertices. Do not assume that a pair of vertex labels uniquely identifies an edge.

Assign edge identifiers or maintain multiplicities so that each edge occurrence is consumed exactly once.

### Preserve K-mer Multiplicity

In sequencing data, repeated k-mers may occur multiple times. If the goal is to reconstruct a spectrum that contains repeated k-mer observations, a set of distinct k-mers is insufficient.

Represent multiplicity explicitly, such as with an edge count.

### Check Graph Connectivity

Degree conditions alone are not enough if the graph contains disconnected components with edges. Validate that the relevant graph structure is connected in the sense required for an Eulerian traversal.

### Handle Cycles and Paths Separately

For a cycle, any suitable starting vertex may be used. For an open directed Eulerian path, start at the vertex whose out-degree exceeds its in-degree by one.

### Validate the Output

After constructing a traversal, verify that:

- Every required edge occurrence was used exactly once.
- Consecutive edges connect correctly.
- The path begins and ends at valid vertices.
- The reconstructed sequence produces the expected k-mer spectrum.

These checks are especially important when the graph contains parallel edges, repeated k-mers, or noisy sequencing data.

---

## Applications Beyond Genome Assembly

Eulerian paths and cycles are useful wherever a problem requires using every connection exactly once.

### Genome Assembly

Reconstruct sequences from k-mer spectra using de Bruijn graphs.

### Sequencing-by-Hybridization

Reconstruct a sequence from observed short oligonucleotides, subject to ambiguity and measurement limitations.

### DNA Sequence Analysis

Represent overlapping k-mers as graph edges and search for a sequence consistent with the observed k-mer collection.

### Network Routing and Circuit Design

Find routes that traverse each connection exactly once, when the graph meets the required conditions.

### Algorithm Design

Study graph connectivity, degree balance, stack-based traversal, and efficient processing of large sparse graphs.

---

## Why Eulerian Paths Matter

Eulerian traversal is a fundamental example of how the right mathematical formulation can make a computational problem tractable.

In sequence assembly, a de Bruijn graph converts overlapping k-mers into directed edges. An Eulerian path or cycle then provides a way to reconstruct a sequence consistent with those edges.

The algorithm is efficient because it processes each edge a bounded number of times instead of exploring all possible edge orderings.

The central lessons are:

- **Eulerian paths use every edge exactly once.**
- **Eulerian cycles use every edge exactly once and return to the starting vertex.**
- **Degree balance and connectivity determine whether a traversal exists.**
- **Hierholzer's algorithm finds an Eulerian traversal in linear time in the graph size.**
- **de Bruijn graphs connect Eulerian traversal to DNA sequence reconstruction.**
- **Efficient traversal does not eliminate the biological ambiguities caused by repeats, errors, and incomplete sequencing.**

For bioinformatics developers, Eulerian algorithms are an excellent foundation for learning graph representations, adjacency lists, edge multiplicity, stack-based algorithms, and the design of efficient sequence-analysis tools in languages such as Rust.