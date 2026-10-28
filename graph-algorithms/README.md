# Graph Algorithms for Genome Assembly

Graph algorithms form the computational backbone of DNA sequence assembly.
Every time you assemble a genome from sequencing reads, whether using
SPAdes, Velvet, or the Verkko assembler, graph algorithms are at work.

## Why Graphs for Genome Assembly?

Imagine shredding a book into millions of tiny overlapping snippets.
Reassembling the book requires figuring out which snippets overlap
and in what order they should be joined. Graphs provide the
mathematical framework to represent these overlaps and find
the correct arrangement.

## Key Insight: Hamiltonian vs. Eulerian

Two fundamentally different graph formulations for genome assembly:

| Approach | Graph Type | Nodes | Edges | Path Type | Complexity |
|----------|-----------|-------|-------|-----------|------------|
| Overlap-Layout-Consensus | Overlap graph | Reads | Overlaps | Hamiltonian | NP-complete |
| de Bruijn graph | de Bruijn graph | (k-1)-mers | k-mers | Eulerian | Linear time |

Finding a Hamiltonian path (visiting every NODE once) is NP-complete.
Finding an Eulerian path (visiting every EDGE once) can be done in
linear time. This is why modern assemblers use de Bruijn graphs.

## Algorithms in This Section

| Algorithm | Purpose | Used By |
|-----------|---------|---------|
| Eulerian Path/Cycle | Traverse every edge exactly once | All de Bruijn assemblers |
| de Bruijn Graph | k-mer-based assembly graph | SPAdes, Velvet, SOAPdenovo |
| Overlap Graph | Read-overlap assembly graph | Celera Assembler, Canu |
| String Graph | Compressed overlap graph | SGA, Hifiasm |
| Shortest Superstring | Greedy assembly approximation | Conceptual foundation |

## Biological Context

- **WGS**: Short-read assemblers (SPAdes) build de Bruijn graphs from
  Illumina k-mers, then find Eulerian paths
- **WES**: Targeted assembly of exome regions uses the same graph framework
- **Transcriptomics**: Trinity builds de Bruijn graphs for de novo
  transcriptome assembly from RNA-seq reads
- **snRNASeq**: Cell Ranger uses graph-based assembly for barcode correction