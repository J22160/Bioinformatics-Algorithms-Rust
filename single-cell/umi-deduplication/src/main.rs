//! UMI Deduplication — Removing PCR Duplicates in scRNA-seq
//!
//! Implements the directional adjacency method from UMI-tools
//! (Smith et al. 2017) to collapse PCR-duplicated UMIs.

use std::collections::{HashMap, HashSet, VecDeque};
use std::time::Instant;

/// Compute Hamming distance between two byte slices.
fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(&x, &y)| x != y).count()
}

/// A UMI with its read count.
#[derive(Debug, Clone)]
struct UmiEntry {
    sequence: Vec<u8>,
    read_count: u32,
    cluster_id: Option<usize>,
}

/// Deduplicate UMIs using the directional adjacency method.
///
/// Returns the number of unique molecules (clusters).
pub fn directional_dedup(umis: &[(Vec<u8>, u32)]) -> Vec<Vec<usize>> {
    let n = umis.len();
    if n == 0 { return vec![]; }

    // Build adjacency list: connect UMIs within Hamming distance 1
    // where the higher-count UMI points to the lower-count UMI
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if hamming_distance(&umis[i].0, &umis[j].0) == 1 {
                // Directional: higher count -> lower count
                if umis[i].1 >= umis[j].1 {
                    adj[i].push(j);
                }
                if umis[j].1 >= umis[i].1 {
                    adj[j].push(i);
                }
            }
        }
    }

    // Directional BFS clustering
    let mut visited = vec![false; n];
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    // Process UMIs in order of decreasing read count
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| umis[b].1.cmp(&umis[a].1));

    for &seed in &order {
        if visited[seed] { continue; }

        let mut cluster = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(seed);
        visited[seed] = true;

        while let Some(current) = queue.pop_front() {
            cluster.push(current);
            for &neighbor in &adj[current] {
                if !visited[neighbor] && umis[neighbor].1 <= umis[current].1 {
                    visited[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        clusters.push(cluster);
    }
    clusters
}

fn main() {
    println!("=== UMI Deduplication: PCR Duplicate Removal ===\n");

    // Simulated UMIs for one gene in one cell
    // Format: (UMI sequence, read count)
    let umis: Vec<(Vec<u8>, u32)> = vec![
        (b"ACGTACGTAC".to_vec(), 150), // True molecule 1
        (b"ACGTACGTAT".to_vec(), 3),   // PCR error of molecule 1
        (b"ACGTACGTAG".to_vec(), 2),   // PCR error of molecule 1
        (b"TGCATGCATG".to_vec(), 200), // True molecule 2
        (b"TGCATGCTTG".to_vec(), 5),   // PCR error of molecule 2
        (b"GGCCAATTGG".to_vec(), 80),  // True molecule 3
        (b"AACCGGTTAA".to_vec(), 120), // True molecule 4 (distant from others)
    ];

    println!("Input UMIs:");
    for (seq, count) in &umis {
        println!("  {} : {} reads", std::str::from_utf8(seq).unwrap(), count);
    }

    let start = Instant::now();
    let clusters = directional_dedup(&umis);
    let elapsed = start.elapsed();

    println!("\nDeduplication result: {} unique molecules\n", clusters.len());
    for (i, cluster) in clusters.iter().enumerate() {
        let representative = &umis[cluster[0]];
        let total_reads: u32 = cluster.iter().map(|&idx| umis[idx].1).sum();
        println!("  Molecule {}: {} (representative)", i + 1,
                 std::str::from_utf8(&representative.0).unwrap());
        println!("    {} UMIs collapsed, {} total reads",
                 cluster.len(), total_reads);
    }
    println!("\nTime: {:?}", elapsed);
    println!("\nExpected: ~4 true molecules (7 UMIs collapsed to 4)");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_umis() {
        let umis = vec![
            (b"ACGT".to_vec(), 100),
        ];
        assert_eq!(directional_dedup(&umis).len(), 1);
    }

    #[test]
    fn test_one_error_collapses() {
        let umis = vec![
            (b"ACGT".to_vec(), 100),
            (b"ACTT".to_vec(), 2), // 1 error, fewer reads
        ];
        assert_eq!(directional_dedup(&umis).len(), 1);
    }

    #[test]
    fn test_distant_umis_separate() {
        let umis = vec![
            (b"ACGT".to_vec(), 100),
            (b"TGCA".to_vec(), 100), // Hamming distance 4
        ];
        assert_eq!(directional_dedup(&umis).len(), 2);
    }
}
