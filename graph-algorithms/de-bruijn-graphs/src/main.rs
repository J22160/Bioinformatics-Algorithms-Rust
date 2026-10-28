//! de Bruijn Graph Construction and Genome Assembly
//!
//! Builds a de Bruijn graph from sequencing reads, then finds
//! an Eulerian path to reconstruct the genome.

use std::collections::HashMap;
use std::time::Instant;

/// A de Bruijn graph with edge coverage tracking.
pub struct DeBruijnGraph {
    /// Node label -> node index
    node_map: HashMap<String, usize>,
    /// Node index -> label
    node_labels: Vec<String>,
    /// Adjacency list: adj[u] contains (v, kmer_label, coverage)
    adj: Vec<Vec<(usize, String, u32)>>,
    /// Total number of edges
    num_edges: usize,
}

impl DeBruijnGraph {
    pub fn new() -> Self {
        DeBruijnGraph {
            node_map: HashMap::new(),
            node_labels: Vec::new(),
            adj: Vec::new(),
            num_edges: 0,
        }
    }

    fn get_or_create_node(&mut self, label: &str) -> usize {
        if let Some(&idx) = self.node_map.get(label) {
            idx
        } else {
            let idx = self.node_labels.len();
            self.node_labels.push(label.to_string());
            self.node_map.insert(label.to_string(), idx);
            self.adj.push(Vec::new());
            idx
        }
    }

    /// Add a k-mer as an edge from its prefix to its suffix.
    pub fn add_kmer(&mut self, kmer: &str) {
        let prefix = &kmer[..kmer.len() - 1];
        let suffix = &kmer[1..];
        let from = self.get_or_create_node(prefix);
        let to = self.get_or_create_node(suffix);

        // Check if edge already exists (increment coverage)
        let existing = self.adj[from].iter_mut()
            .find(|(dest, label, _)| *dest == to && label == kmer);

        if let Some((_, _, cov)) = existing {
            *cov += 1;
        } else {
            self.adj[from].push((to, kmer.to_string(), 1));
            self.num_edges += 1;
        }
    }

    /// Build de Bruijn graph from reads with given k.
    pub fn from_reads(reads: &[&str], k: usize) -> Self {
        let mut graph = DeBruijnGraph::new();
        for read in reads {
            if read.len() < k {
                continue;
            }
            for i in 0..=(read.len() - k) {
                let kmer = &read[i..i + k];
                graph.add_kmer(kmer);
            }
        }
        graph
    }

    /// Get degree information for each node.
    pub fn degrees(&self) -> Vec<(usize, usize)> {
        let n = self.node_labels.len();
        let mut in_deg = vec![0usize; n];
        let mut out_deg = vec![0usize; n];

        for u in 0..n {
            for (v, _, cov) in &self.adj[u] {
                out_deg[u] += *cov as usize;
                in_deg[*v] += *cov as usize;
            }
        }

        in_deg.into_iter().zip(out_deg.into_iter()).collect()
    }

    /// Find Eulerian path using Hierholzer's algorithm.
    /// Handles multigraph edges by respecting coverage counts.
    pub fn find_eulerian_path(&self) -> Option<Vec<usize>> {
        let degs = self.degrees();
        let n = self.node_labels.len();

        // Find start node
        let mut start = 0;
        for i in 0..n {
            let (in_d, out_d) = degs[i];
            if out_d > in_d {
                start = i;
                break;
            }
        }

        // Expand multigraph: create flat adjacency lists
        let mut flat_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut total_edges = 0;
        for u in 0..n {
            for (v, _, cov) in &self.adj[u] {
                for _ in 0..*cov {
                    flat_adj[u].push(*v);
                    total_edges += 1;
                }
            }
        }

        // Hierholzer's
        let mut edge_idx = vec![0usize; n];
        let mut stack = vec![start];
        let mut path = Vec::with_capacity(total_edges + 1);

        while let Some(&current) = stack.last() {
            if edge_idx[current] < flat_adj[current].len() {
                let next = flat_adj[current][edge_idx[current]];
                edge_idx[current] += 1;
                stack.push(next);
            } else {
                path.push(stack.pop().unwrap());
            }
        }

        path.reverse();

        if path.len() == total_edges + 1 {
            Some(path)
        } else {
            None
        }
    }

    /// Reconstruct sequence from Eulerian path.
    pub fn path_to_sequence(&self, path: &[usize]) -> String {
        if path.is_empty() {
            return String::new();
        }
        let mut seq = self.node_labels[path[0]].clone();
        for &node in &path[1..] {
            let label = &self.node_labels[node];
            seq.push(label.as_bytes()[label.len() - 1] as char);
        }
        seq
    }

    /// Print graph statistics.
    pub fn print_stats(&self) {
        println!("Nodes: {}", self.node_labels.len());
        println!("Distinct edges: {}", self.num_edges);
        let total: u32 = self.adj.iter()
            .flat_map(|v| v.iter().map(|(_, _, c)| c))
            .sum();
        println!("Total edges (with multiplicity): {}", total);
    }
}

fn main() {
    println!("=== de Bruijn Graph: Genome Assembly from Reads ===\n");

    let reads = vec!["AAGCTGCA", "GCTGCAAG", "TGCAAGCT"];
    let k = 4;

    println!("Input reads:");
    for r in &reads {
        println!("  {}", r);
    }
    println!("k = {}\n", k);

    let graph = DeBruijnGraph::from_reads(&reads, k);

    println!("de Bruijn graph:");
    graph.print_stats();

    println!("\nEdges (k-mers):");
    for u in 0..graph.node_labels.len() {
        for (v, label, cov) in &graph.adj[u] {
            println!("  {} --[{}]--> {} (coverage: {})",
                     graph.node_labels[u], label, graph.node_labels[*v], cov);
        }
    }

    // Degree analysis
    println!("\nNode degrees:");
    let degs = graph.degrees();
    for (i, (in_d, out_d)) in degs.iter().enumerate() {
        let status = if in_d == out_d { "balanced" }
        else if out_d > in_d { "START" }
        else { "END" };
        println!("  {}: in={} out={} ({})",
                 graph.node_labels[i], in_d, out_d, status);
    }

    match graph.find_eulerian_path() {
        Some(path) => {
            let node_path: Vec<&str> = path.iter()
                .map(|&i| graph.node_labels[i].as_str())
                .collect();
            println!("\nEulerian path:");
            println!("  {}", node_path.join(" -> "));
            let sequence = graph.path_to_sequence(&path);
            println!("\nAssembled sequence: {}", sequence);
            println!("Length: {} bp", sequence.len());

            // Verify reads are substrings
            println!("\nVerification:");
            for read in &reads {
                let found = sequence.contains(read)
                    || format!("{}{}", &sequence, &sequence[..sequence.len()-1]).contains(read);
                println!("  {} in assembly: {}", read, if found { "✓" } else { "✗" });
            }
        }
        None => println!("\nNo Eulerian path found — graph may be disconnected."),
    }

    // Benchmark: k-mer extraction speed
    println!("\n=== Benchmark: k-mer Extraction + Graph Build ===");
    for &n in &[10_000, 100_000, 1_000_000] {
        let read: String = (0..150).map(|i| b"ACGT"[i % 4] as char).collect();
        let reads_vec: Vec<&str> = std::iter::repeat(read.as_str()).take(n / 150).collect();
        let start = Instant::now();
        let g = DeBruijnGraph::from_reads(&reads_vec, 31);
        let elapsed = start.elapsed();
        println!("  {} bases: {} nodes, {} edges, built in {:?}",
                 n, g.node_labels.len(), g.num_edges, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_assembly() {
        let reads = vec!["ATGCGT", "GCGTAC"];
        let graph = DeBruijnGraph::from_reads(&reads, 4);
        let path = graph.find_eulerian_path().unwrap();
        let seq = graph.path_to_sequence(&path);
        assert!(seq.contains("ATGCGTAC"));
    }

    #[test]
    fn test_repeated_kmer() {
        // Genome: ATGATG (circular with repeat ATG)
        let kmers = vec!["ATG", "TGA", "GAT", "ATG"];
        let mut graph = DeBruijnGraph::new();
        for kmer in &kmers {
            graph.add_kmer(kmer);
        }
        let path = graph.find_eulerian_path();
        assert!(path.is_some());
    }
}
