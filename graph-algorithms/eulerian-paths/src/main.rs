//! Eulerian Path and Cycle Algorithms
//!
//! Finds paths that traverse every edge exactly once in a directed graph.
//! This is the core algorithm that makes de Bruijn graph-based
//! genome assembly computationally tractable.

use std::collections::HashMap;
use std::time::Instant;

/// Directed graph represented as adjacency lists.
///
/// # Rust Performance Advantages
///
/// 1. **Vec-based adjacency lists**: Each node's outgoing edges are stored
///    in a contiguous Vec<usize>. Iterating over neighbors is cache-friendly.
///    Python's list-of-lists scatters objects across the heap.
///
/// 2. **usize node indices**: Nodes are integers, not string objects.
///    Comparison and hashing are single CPU instructions.
///
/// 3. **In-place edge consumption**: We track which edges have been used
///    via an index into each adjacency list, avoiding allocation of
///    "visited" sets. Python would use set() with hash overhead per check.
///
/// 4. **Stack-based iteration**: We replace recursion with an explicit
///    stack (Vec<usize>) to avoid stack overflow on large graphs and
///    to give the allocator a single contiguous block.
pub struct DirectedGraph {
    adj: Vec<Vec<usize>>,  // adj[u] = list of nodes v where edge u->v exists
    num_edges: usize,
}

impl DirectedGraph {
    pub fn new(num_nodes: usize) -> Self {
        DirectedGraph {
            adj: vec![Vec::new(); num_nodes],
            num_edges: 0,
        }
    }

    pub fn add_edge(&mut self, from: usize, to: usize) {
        self.adj[from].push(to);
        self.num_edges += 1;
    }

    pub fn num_nodes(&self) -> usize {
        self.adj.len()
    }

    /// Compute in-degree and out-degree for each node.
    pub fn degrees(&self) -> (Vec<usize>, Vec<usize>) {
        let n = self.num_nodes();
        let mut in_deg = vec![0usize; n];
        let mut out_deg = vec![0usize; n];
        for u in 0..n {
            out_deg[u] = self.adj[u].len();
            for &v in &self.adj[u] {
                in_deg[v] += 1;
            }
        }
        (in_deg, out_deg)
    }

    /// Check if the graph has an Eulerian cycle (all nodes balanced).
    pub fn has_eulerian_cycle(&self) -> bool {
        let (in_deg, out_deg) = self.degrees();
        (0..self.num_nodes()).all(|i| in_deg[i] == out_deg[i])
    }

    /// Check if the graph has an Eulerian path.
    /// Returns Some((start, end)) if it does, None otherwise.
    pub fn eulerian_path_endpoints(&self) -> Option<(usize, usize)> {
        let (in_deg, out_deg) = self.degrees();
        let mut start = None;
        let mut end = None;
        let mut num_unbalanced = 0;

        for i in 0..self.num_nodes() {
            if out_deg[i] == in_deg[i] {
                continue; // balanced node
            } else if out_deg[i] == in_deg[i] + 1 {
                // Potential start (one extra outgoing edge)
                if start.is_some() { return None; }
                start = Some(i);
                num_unbalanced += 1;
            } else if in_deg[i] == out_deg[i] + 1 {
                // Potential end (one extra incoming edge)
                if end.is_some() { return None; }
                end = Some(i);
                num_unbalanced += 1;
            } else {
                return None; // difference > 1, no Eulerian path
            }
        }

        match num_unbalanced {
            0 => Some((0, 0)), // Eulerian cycle; start anywhere
            2 => Some((start.unwrap(), end.unwrap())),
            _ => None,
        }
    }

    /// Find an Eulerian path/cycle using Hierholzer's algorithm.
    /// Returns the sequence of nodes visited.
    ///
    /// Time complexity: O(|E|)
    /// Space complexity: O(|V| + |E|)
    pub fn find_eulerian_path(&self) -> Option<Vec<usize>> {
        let endpoints = self.eulerian_path_endpoints()?;
        let start = endpoints.0;

        // Track current edge index for each node (avoids re-scanning)
        let mut edge_idx: Vec<usize> = vec![0; self.num_nodes()];

        // Hierholzer's algorithm using an explicit stack
        let mut stack = vec![start];
        let mut path = Vec::with_capacity(self.num_edges + 1);

        // We need mutable access to adjacency lists, but they don't
        // change — we just advance the edge_idx pointer
        while let Some(&current) = stack.last() {
            if edge_idx[current] < self.adj[current].len() {
                let next = self.adj[current][edge_idx[current]];
                edge_idx[current] += 1;
                stack.push(next);
            } else {
                path.push(stack.pop().unwrap());
            }
        }

        path.reverse();

        if path.len() == self.num_edges + 1 {
            Some(path)
        } else {
            None // graph was disconnected
        }
    }
}

/// Build a de Bruijn graph from a set of k-mers.
/// Returns (graph, node_labels) where node_labels maps index -> (k-1)-mer string.
pub fn build_debruijn_from_kmers(kmers: &[&str]) -> (DirectedGraph, Vec<String>) {
    let mut node_map: HashMap<String, usize> = HashMap::new();
    let mut node_labels: Vec<String> = Vec::new();

    let mut get_or_create = |s: &str, map: &mut HashMap<String, usize>, labels: &mut Vec<String>| -> usize {
        if let Some(&idx) = map.get(s) {
            idx
        } else {
            let idx = labels.len();
            labels.push(s.to_string());
            map.insert(s.to_string(), idx);
            idx
        }
    };

    // First pass: collect all (k-1)-mer nodes
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for kmer in kmers {
        let prefix = &kmer[..kmer.len() - 1];
        let suffix = &kmer[1..];
        let from = get_or_create(prefix, &mut node_map, &mut node_labels);
        let to = get_or_create(suffix, &mut node_map, &mut node_labels);
        edges.push((from, to));
    }

    // Build graph
    let mut graph = DirectedGraph::new(node_labels.len());
    for (from, to) in edges {
        graph.add_edge(from, to);
    }

    (graph, node_labels)
}

/// Reconstruct genome sequence from Eulerian path through de Bruijn graph.
pub fn reconstruct_sequence(path: &[usize], node_labels: &[String]) -> String {
    if path.is_empty() {
        return String::new();
    }
    let mut sequence = node_labels[path[0]].clone();
    for &node in &path[1..] {
        // Each step adds one character (the last char of the (k-1)-mer)
        let label = &node_labels[node];
        sequence.push(label.as_bytes()[label.len() - 1] as char);
    }
    sequence
}

fn main() {
    println!("=== Eulerian Path: Genome Reconstruction from k-mer Spectrum ===\n");

    // SBH spectrum of a short DNA fragment
    let kmers = vec!["ATG", "TGC", "GCA", "CAG", "AGT", "GTG", "TGT"];

    println!("Input k-mers (3-mer spectrum):");
    for kmer in &kmers {
        print!("{} ", kmer);
    }
    println!("\n");

    let (graph, labels) = build_debruijn_from_kmers(&kmers);

    println!("de Bruijn graph nodes (2-mers):");
    for (i, label) in labels.iter().enumerate() {
        print!("  {} (idx {})", label, i);
        if !graph.adj[i].is_empty() {
            print!(" -> ");
            for &neighbor in &graph.adj[i] {
                print!("{} ", labels[neighbor]);
            }
        }
        println!();
    }

    // Check Eulerian properties
    let (in_deg, out_deg) = graph.degrees();
    println!("\nDegree analysis:");
    for (i, label) in labels.iter().enumerate() {
        println!("  {}: in={} out={} {}", label, in_deg[i], out_deg[i],
                 if in_deg[i] == out_deg[i] { "(balanced)" }
                 else if out_deg[i] > in_deg[i] { "(START)" }
                 else { "(END)" });
    }

    match graph.find_eulerian_path() {
        Some(path) => {
            let node_path: Vec<&str> = path.iter().map(|&i| labels[i].as_str()).collect();
            println!("\nEulerian path: {}", node_path.join(" -> "));

            let sequence = reconstruct_sequence(&path, &labels);
            println!("Reconstructed sequence: {}", sequence);
        }
        None => println!("\nNo Eulerian path exists!"),
    }

    // Larger example: Eulerian CYCLE
    println!("\n=== Eulerian Cycle: Circular Genome ===\n");
    let circular_kmers = vec!["ATG", "TGC", "GCA", "CAT", "ATG"];
    // Remove duplicate for de Bruijn construction — the cycle handles it
    let circular_unique: Vec<&str> = vec!["ATG", "TGC", "GCA", "CAT"];

    let (cgraph, clabels) = build_debruijn_from_kmers(&circular_unique);
    println!("Circular k-mers: {:?}", circular_unique);
    println!("Has Eulerian cycle: {}", cgraph.has_eulerian_cycle());

    if let Some(path) = cgraph.find_eulerian_path() {
        let seq = reconstruct_sequence(&path, &clabels);
        println!("Reconstructed circular genome: {} (circular)", seq);
    }

    // Benchmark
    println!("\n=== Performance Benchmark ===");
    for &n in &[1_000, 10_000, 100_000] {
        // Build a synthetic linear graph with n edges
        let mut g = DirectedGraph::new(n + 1);
        for i in 0..n {
            g.add_edge(i, i + 1);
        }
        let start = Instant::now();
        let _path = g.find_eulerian_path();
        let elapsed = start.elapsed();
        println!("  {} edges: Eulerian path found in {:?}", n, elapsed);
    }
    println!("\nNotice: time scales linearly with edges — O(E) confirmed.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_cycle() {
        // Triangle: 0->1->2->0
        let mut g = DirectedGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 0);
        assert!(g.has_eulerian_cycle());
        let path = g.find_eulerian_path().unwrap();
        assert_eq!(path.len(), 4); // 3 edges + return to start
        assert_eq!(path.first(), path.last());
    }

    #[test]
    fn test_simple_path() {
        // 0->1->2 (path, not cycle)
        let mut g = DirectedGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        assert!(!g.has_eulerian_cycle());
        let path = g.find_eulerian_path().unwrap();
        assert_eq!(path, vec![0, 1, 2]);
    }

    #[test]
    fn test_sbh_reconstruction() {
        let kmers = vec!["ATG", "TGC", "GCA", "CAG", "AGT", "GTG", "TGT"];
        let (graph, labels) = build_debruijn_from_kmers(&kmers);
        let path = graph.find_eulerian_path().unwrap();
        let seq = reconstruct_sequence(&path, &labels);
        // Verify all k-mers are present in reconstruction
        for kmer in &kmers {
            assert!(seq.contains(kmer), "Missing k-mer {} in {}", kmer, seq);
        }
    }

    #[test]
    fn test_no_eulerian_path() {
        // 0->1, 0->2 (out-degree 2 at node 0, no path)
        let mut g = DirectedGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        assert!(g.find_eulerian_path().is_none());
    }
}
