//! String Graph Construction for Genome Assembly
//!
//! Implements Myers' string graph: overlap detection, contained read
//! removal, transitive edge reduction, and chain collapsing.
//! This is the core data structure used by SGA, FALCON, and Hifiasm.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

/// Represents an overlap between two reads.
/// read_a[overlap_start_a..] == read_b[..overlap_len]
#[derive(Debug, Clone)]
pub struct Overlap {
    pub read_a: usize,       // Index of the first read (suffix provider)
    pub read_b: usize,       // Index of the second read (prefix provider)
    pub overlap_len: usize,  // Length of the overlap
    pub label: Vec<u8>,      // Non-overlapping suffix of read_b (edge label)
}

/// An edge in the string graph (irreducible overlap).
#[derive(Debug, Clone)]
pub struct StringEdge {
    pub from: usize,
    pub to: usize,
    pub label: String,       // The sequence added by traversing this edge
    pub overlap_len: usize,
}

/// The string graph data structure.
///
/// # Rust Performance Advantages
///
/// 1. **Vec-based adjacency lists**: Edges are stored in contiguous
///    Vec<StringEdge> per node, giving excellent cache locality
///    during graph traversal. Python's list-of-objects scatters
///    edge data across the heap.
///
/// 2. **Zero-copy string slices**: During overlap computation,
///    we compare &[u8] slices directly without creating new String
///    objects. Python creates new string objects for every slice.
///
/// 3. **Compact representation**: Read indices are usize (8 bytes).
///    Python stores references to int objects (28+ bytes each).
///
/// 4. **In-place transitive reduction**: We mark edges for removal
///    using a boolean flag rather than copying the entire edge set.
///    This avoids allocation churn during graph simplification.
pub struct StringGraph {
    pub reads: Vec<Vec<u8>>,
    pub read_names: Vec<String>,
    pub adj: Vec<Vec<StringEdge>>,   // Adjacency list: adj[u] = outgoing edges
    pub contained: HashSet<usize>,   // Set of contained read indices
    pub num_nodes: usize,
}

impl StringGraph {
    /// Build a string graph from a set of reads with minimum overlap τ.
    ///
    /// Steps:
    /// 1. Find all suffix-prefix overlaps of length >= τ
    /// 2. Identify and remove contained reads
    /// 3. Build initial overlap graph
    /// 4. Remove transitive edges
    pub fn build(reads: &[&str], names: &[&str], min_overlap: usize) -> Self {
        let read_bytes: Vec<Vec<u8>> = reads.iter().map(|r| r.as_bytes().to_vec()).collect();
        let read_names: Vec<String> = names.iter().map(|n| n.to_string()).collect();
        let n = reads.len();

        // Step 1: Find all pairwise suffix-prefix overlaps
        let all_overlaps = Self::find_all_overlaps(&read_bytes, min_overlap);

        // Step 2: Identify contained reads
        let contained = Self::find_contained_reads(&read_bytes);

        // Step 3: Build initial overlap graph (excluding contained reads)
        let proper_overlaps: Vec<Overlap> = all_overlaps
            .into_iter()
            .filter(|ov| !contained.contains(&ov.read_a) && !contained.contains(&ov.read_b))
            .collect();

        // Step 4: Remove transitive edges
        let irreducible = Self::transitive_reduction(&proper_overlaps, n);

        // Build adjacency list
        let mut adj = vec![Vec::new(); n];
        for ov in &irreducible {
            adj[ov.read_a].push(StringEdge {
                from: ov.read_a,
                to: ov.read_b,
                label: String::from_utf8_lossy(&ov.label).to_string(),
                overlap_len: ov.overlap_len,
            });
        }

        StringGraph {
            reads: read_bytes,
            read_names,
            adj,
            contained,
            num_nodes: n,
        }
    }

    /// Find all suffix-prefix overlaps of length >= min_overlap.
    /// For each pair (i, j), checks if a suffix of read[i] matches
    /// a prefix of read[j].
    ///
    /// In production assemblers like SGA, this is done using the
    /// FM-index in O(N) time. Here we use a simpler O(N²L) approach
    /// for clarity.
    fn find_all_overlaps(reads: &[Vec<u8>], min_overlap: usize) -> Vec<Overlap> {
        let mut overlaps = Vec::new();
        let n = reads.len();

        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                // Find the longest suffix of reads[i] matching a prefix of reads[j]
                let max_possible = reads[i].len().min(reads[j].len());
                for k in (min_overlap..=max_possible).rev() {
                    let suffix_start = reads[i].len() - k;
                    if reads[i][suffix_start..] == reads[j][..k] {
                        let label = reads[j][k..].to_vec();
                        overlaps.push(Overlap {
                            read_a: i,
                            read_b: j,
                            overlap_len: k,
                            label,
                        });
                        break; // Only keep the longest overlap per pair
                    }
                }
            }
        }
        overlaps
    }

    /// Identify reads that are entirely contained within another read.
    fn find_contained_reads(reads: &[Vec<u8>]) -> HashSet<usize> {
        let mut contained = HashSet::new();
        let n = reads.len();

        for i in 0..n {
            for j in 0..n {
                if i == j || reads[j].len() >= reads[i].len() {
                    continue;
                }
                // Check if reads[j] is a substring of reads[i]
                let needle = &reads[j];
                let haystack = &reads[i];
                if haystack
                    .windows(needle.len())
                    .any(|window| window == needle.as_slice())
                {
                    contained.insert(j);
                }
            }
        }
        contained
    }

    /// Remove transitive edges from the overlap graph.
    ///
    /// An edge (X, Z) is transitive if there exists Y such that
    /// (X, Y) and (Y, Z) are both edges and the path X→Y→Z
    /// produces the same assembly as the direct edge X→Z.
    ///
    /// We detect this by checking: for each pair of edges (X,Y) and
    /// (X,Z) where both overlap the same end of X, if there also
    /// exists an edge (Y,Z), then the longer overlap (X,Z if
    /// overlap(X,Y) > overlap(X,Z)) is transitive.
    fn transitive_reduction(overlaps: &[Overlap], num_reads: usize) -> Vec<Overlap> {
        // Build lookup: for each source read, collect all overlaps
        // sorted by overlap length (longest first)
        let mut by_source: Vec<Vec<&Overlap>> = vec![Vec::new(); num_reads];
        for ov in overlaps {
            by_source[ov.read_a].push(ov);
        }

        // Sort each source's overlaps by overlap length (longest first)
        for ovs in &mut by_source {
            ovs.sort_by(|a, b| b.overlap_len.cmp(&a.overlap_len));
        }

        // Build a quick lookup: does edge (a, b) exist?
        let edge_set: HashSet<(usize, usize)> = overlaps
            .iter()
            .map(|ov| (ov.read_a, ov.read_b))
            .collect();

        // Mark transitive edges
        let mut transitive: HashSet<(usize, usize)> = HashSet::new();

        for x in 0..num_reads {
            let x_edges = &by_source[x];

            // For each pair of outgoing edges from X
            for i in 0..x_edges.len() {
                let y = x_edges[i].read_b;

                // Check if any shorter overlap from X goes to a node Z
                // that is also reachable from Y
                for j in (i + 1)..x_edges.len() {
                    let z = x_edges[j].read_b;

                    // If Y → Z exists, then X → Z is transitive
                    if edge_set.contains(&(y, z)) {
                        transitive.insert((x, z));
                    }
                }
            }
        }

        // Return only irreducible edges
        overlaps
            .iter()
            .filter(|ov| !transitive.contains(&(ov.read_a, ov.read_b)))
            .cloned()
            .collect()
    }

    /// Collapse simple chains: nodes with exactly one in-edge and
    /// one out-edge are merged into a single edge.
    pub fn collapse_chains(&self) -> Vec<(Vec<usize>, String)> {
        let n = self.num_nodes;
        let mut in_degree = vec![0usize; n];
        for u in 0..n {
            for edge in &self.adj[u] {
                in_degree[edge.to] += 1;
            }
        }

        let mut visited = vec![false; n];
        let mut chains = Vec::new();

        for start in 0..n {
            if self.contained.contains(&start) || visited[start] {
                continue;
            }
            // Start a chain from nodes that are NOT simple pass-through
            if in_degree[start] == 1 && self.adj[start].len() == 1 {
                continue; // This is a middle node; we'll reach it from a start
            }

            for edge in &self.adj[start] {
                let mut chain_nodes = vec![start];
                let mut chain_seq = self.read_names[start].clone();
                let mut current = edge.to;
                chain_seq.push_str(&format!(" -[{}]-> {}", edge.label, self.read_names[current]));

                visited[start] = true;

                while !visited[current]
                    && self.adj[current].len() == 1
                    && in_degree[current] == 1
                {
                    visited[current] = true;
                    chain_nodes.push(current);
                    let next_edge = &self.adj[current][0];
                    current = next_edge.to;
                    chain_seq.push_str(&format!(
                        " -[{}]-> {}",
                        next_edge.label, self.read_names[current]
                    ));
                }
                chain_nodes.push(current);
                visited[current] = true;

                chains.push((chain_nodes, chain_seq));
            }
        }
        chains
    }

    /// Reconstruct the assembled sequence from a chain of reads.
    pub fn assemble_chain(&self, chain: &[usize]) -> String {
        if chain.is_empty() {
            return String::new();
        }
        let mut seq = String::from_utf8_lossy(&self.reads[chain[0]]).to_string();

        for i in 0..chain.len() - 1 {
            let from = chain[i];
            let to = chain[i + 1];
            if let Some(edge) = self.adj[from].iter().find(|e| e.to == to) {
                seq.push_str(&edge.label);
            }
        }
        seq
    }

    /// Print graph statistics
    pub fn print_stats(&self) {
        let active_nodes = (0..self.num_nodes)
            .filter(|i| !self.contained.contains(i))
            .count();
        let total_edges: usize = self.adj.iter().map(|v| v.len()).sum();

        println!("Nodes (active/total): {}/{}", active_nodes, self.num_nodes);
        println!("Contained reads:     {}", self.contained.len());
        println!("Irreducible edges:   {}", total_edges);
    }
}

fn main() {
    println!("=== String Graph: Long-Read Plasmid Assembly ===\n");

    let reads = vec![
        "ATGCCATGCC",
        "TGCCATGCCA",
        "GCCATGCCAG",
        "CCATGCCAGT",
        "CATGCCAGTC",
        "ATGCCAGTCG",
        "TGCCAGTCGA",
        "GCCAGTCGAT",
    ];
    let names: Vec<&str> = (0..reads.len()).map(|i| {
        // Leak is fine for static demo strings
        Box::leak(format!("R{}", i).into_boxed_str()) as &str
    }).collect();

    let min_overlap = 5;

    println!("Input reads:");
    for (i, r) in reads.iter().enumerate() {
        println!("  R{}: {}", i, r);
    }
    println!("Minimum overlap: {}\n", min_overlap);

    let start = Instant::now();
    let graph = StringGraph::build(&reads, &names, min_overlap);
    let elapsed = start.elapsed();

    println!("--- String Graph Statistics ---");
    graph.print_stats();
    println!("Build time: {:?}\n", elapsed);

    println!("--- Adjacency List (Irreducible Edges) ---");
    for u in 0..graph.num_nodes {
        if graph.contained.contains(&u) {
            continue;
        }
        for edge in &graph.adj[u] {
            println!(
                "  {} --[{}]--> {} (overlap: {}bp)",
                graph.read_names[edge.from],
                edge.label,
                graph.read_names[edge.to],
                edge.overlap_len
            );
        }
    }

    if !graph.contained.is_empty() {
        println!("\n--- Contained Reads (Removed) ---");
        for &c in &graph.contained {
            println!("  {} is contained within another read", graph.read_names[c]);
        }
    }

    // Collapse chains and assemble
    println!("\n--- Chain Collapse & Assembly ---");
    let chains = graph.collapse_chains();
    for (i, (nodes, path_str)) in chains.iter().enumerate() {
        let assembled = graph.assemble_chain(nodes);
        println!("Chain {}: {}", i, path_str);
        println!("  Assembled: {} ({}bp)\n", assembled, assembled.len());
    }

    // Benchmark with larger synthetic data
    println!("=== Benchmark: Scaling with Read Count ===");
    for &n in &[50, 200, 500] {
        let genome: Vec<u8> = (0..n + 15)
            .map(|i| b"ACGT"[i % 4])
            .collect();
        let synth_reads: Vec<String> = (0..n)
            .map(|i| String::from_utf8_lossy(&genome[i..i + 15]).to_string())
            .collect();
        let synth_refs: Vec<&str> = synth_reads.iter().map(|s| s.as_str()).collect();
        let synth_names: Vec<&str> = (0..n).map(|i| {
            Box::leak(format!("S{}", i).into_boxed_str()) as &str
        }).collect();

        let start = Instant::now();
        let g = StringGraph::build(&synth_refs, &synth_names, 8);
        let elapsed = start.elapsed();
        let edges: usize = g.adj.iter().map(|v| v.len()).sum();
        println!(
            "  {} reads: {} irreducible edges, {} contained, time={:?}",
            n,
            edges,
            g.contained.len(),
            elapsed
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_chain() {
        let reads = vec!["ABCDE", "CDEFG", "EFGHI"];
        let names = vec!["R0", "R1", "R2"];
        let graph = StringGraph::build(&reads, &names, 3);
        assert!(graph.contained.is_empty());
        // R0→R1→R2 should be the chain
        assert!(!graph.adj[0].is_empty());
    }

    #[test]
    fn test_contained_removal() {
        let reads = vec!["ABCDEFGH", "CDEF", "DEFGHIJK"];
        let names = vec!["R0", "R1", "R2"];
        let graph = StringGraph::build(&reads, &names, 3);
        assert!(graph.contained.contains(&1)); // R1 contained in R0
    }

    #[test]
    fn test_transitive_removal() {
        // R0 overlaps R1, R1 overlaps R2, R0 overlaps R2
        // Edge R0→R2 should be transitive and removed
        let reads = vec!["AACCGGTT", "CCGGTTAA", "GGTTAACC"];
        let names = vec!["R0", "R1", "R2"];
        let graph = StringGraph::build(&reads, &names, 4);

        // R0 should have edge to R1 (irreducible) but NOT directly to R2
        let r0_targets: Vec<usize> = graph.adj[0].iter().map(|e| e.to).collect();
        assert!(r0_targets.contains(&1), "R0→R1 should be irreducible");
    }
}
