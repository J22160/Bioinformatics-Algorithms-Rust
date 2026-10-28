//! Overlap Graph Construction for Genome Assembly

use std::time::Instant;

/// Compute the length of the longest suffix of `a` that matches
/// a prefix of `b`.
pub fn overlap_length(a: &[u8], b: &[u8], min_overlap: usize) -> usize {
    let max_possible = a.len().min(b.len());
    for k in (min_overlap..=max_possible).rev() {
        if a[a.len() - k..] == b[..k] {
            return k;
        }
    }
    0
}

/// Build overlap graph from reads.
/// Returns adjacency matrix of overlap lengths.
pub fn build_overlap_graph(reads: &[&[u8]], min_overlap: usize) -> Vec<Vec<usize>> {
    let n = reads.len();
    let mut graph = vec![vec![0usize; n]; n];

    for i in 0..n {
        for j in 0..n {
            if i != j {
                graph[i][j] = overlap_length(reads[i], reads[j], min_overlap);
            }
        }
    }
    graph
}

/// Greedy Hamiltonian path: repeatedly pick the best unused overlap.
pub fn greedy_hamiltonian(graph: &[Vec<usize>], reads: &[&[u8]]) -> Vec<usize> {
    let n = reads.len();
    let mut used = vec![false; n];
    let mut path = Vec::with_capacity(n);

    // Start with node 0 (arbitrary)
    let mut current = 0;
    used[current] = true;
    path.push(current);

    for _ in 1..n {
        let mut best_next = None;
        let mut best_overlap = 0;
        for j in 0..n {
            if !used[j] && graph[current][j] > best_overlap {
                best_overlap = graph[current][j];
                best_next = Some(j);
            }
        }
        if let Some(next) = best_next {
            used[next] = true;
            path.push(next);
            current = next;
        }
    }
    path
}

/// Reconstruct sequence from ordered reads using overlaps.
pub fn reconstruct_from_path(
    path: &[usize],
    graph: &[Vec<usize>],
    reads: &[&[u8]],
) -> String {
    let mut seq = Vec::from(reads[path[0]]);
    for i in 1..path.len() {
        let prev = path[i - 1];
        let curr = path[i];
        let ol = graph[prev][curr];
        seq.extend_from_slice(&reads[curr][ol..]);
    }
    String::from_utf8(seq).unwrap()
}

fn main() {
    println!("=== Overlap Graph: Genome Assembly ===\n");

    let reads: Vec<&[u8]> = vec![
        b"ATGCCA", b"GCCATG", b"CATGCC", b"TGCCAT", b"CCATGC",
    ];
    let read_strs: Vec<&str> = vec!["ATGCCA", "GCCATG", "CATGCC", "TGCCAT", "CCATGC"];
    let min_overlap = 3;

    println!("Reads:");
    for (i, r) in read_strs.iter().enumerate() {
        println!("  R{}: {}", i, r);
    }

    let graph = build_overlap_graph(&reads, min_overlap);

    println!("\nOverlap matrix (overlap lengths):");
    print!("     ");
    for i in 0..reads.len() { print!("R{}  ", i); }
    println!();
    for i in 0..reads.len() {
        print!("R{}:  ", i);
        for j in 0..reads.len() {
            print!("{:<4}", graph[i][j]);
        }
        println!();
    }

    let path = greedy_hamiltonian(&graph, &reads);
    let path_str: Vec<String> = path.iter().map(|&i| format!("R{}", i)).collect();
    println!("\nGreedy Hamiltonian path: {}", path_str.join(" -> "));

    let assembled = reconstruct_from_path(&path, &graph, &reads);
    println!("Assembled sequence: {}", assembled);

    // Benchmark
    println!("\n=== Benchmark ===");
    for &n in &[100, 500, 1000] {
        let synth_reads: Vec<Vec<u8>> = (0..n)
            .map(|i| {
                (0..50).map(|j| b"ACGT"[((i * 7 + j * 3) % 4)]).collect()
            })
            .collect();
        let refs: Vec<&[u8]> = synth_reads.iter().map(|v| v.as_slice()).collect();
        let start = Instant::now();
        let _g = build_overlap_graph(&refs, 10);
        let elapsed = start.elapsed();
        println!("  {} reads: overlap graph built in {:?}", n, elapsed);
    }
    println!("\nNote: O(n²) scaling — this is why de Bruijn graphs are preferred.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlap_length() {
        assert_eq!(overlap_length(b"ABCDEF", b"DEFGHI", 3), 3);
        assert_eq!(overlap_length(b"ABCDEF", b"CDEFGH", 3), 4);
        assert_eq!(overlap_length(b"ABCDEF", b"GHIJKL", 3), 0);
    }
}
