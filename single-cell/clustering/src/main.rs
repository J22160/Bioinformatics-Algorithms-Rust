//! Graph-Based Clustering — Simplified Louvain Community Detection
//!
//! Implements a simplified version of the Louvain algorithm for
//! community detection on a k-nearest-neighbor graph.

use std::collections::HashMap;
use std::time::Instant;

/// Compute Euclidean distance between two points.
fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Build a KNN graph from data points.
/// Returns adjacency list with (neighbor_index, weight) pairs.
pub fn build_knn_graph(data: &[Vec<f64>], k: usize) -> Vec<Vec<(usize, f64)>> {
    let n = data.len();
    let mut graph = vec![Vec::new(); n];

    for i in 0..n {
        let mut dists: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, euclidean_dist(&data[i], &data[j])))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        // Convert distance to similarity (weight)
        for &(j, d) in dists.iter().take(k) {
            let w = 1.0 / (1.0 + d);
            graph[i].push((j, w));
        }
    }
    graph
}

/// Compute modularity of a partition.
pub fn modularity(graph: &[Vec<(usize, f64)>], communities: &[usize]) -> f64 {
    let n = graph.len();
    let mut total_weight = 0.0;
    for edges in graph {
        for &(_, w) in edges {
            total_weight += w;
        }
    }
    if total_weight == 0.0 { return 0.0; }

    let mut q = 0.0;
    for i in 0..n {
        let ki: f64 = graph[i].iter().map(|&(_, w)| w).sum();
        for &(j, w) in &graph[i] {
            if communities[i] == communities[j] {
                let kj: f64 = graph[j].iter().map(|&(_, w)| w).sum();
                q += w - (ki * kj) / total_weight;
            }
        }
    }
    q / total_weight
}

/// Simplified Louvain community detection.
pub fn louvain_clustering(graph: &[Vec<(usize, f64)>], max_iter: usize) -> Vec<usize> {
    let n = graph.len();
    let mut communities: Vec<usize> = (0..n).collect();
    let mut improved = true;
    let mut iteration = 0;

    while improved && iteration < max_iter {
        improved = false;
        iteration += 1;

        for i in 0..n {
            let current_comm = communities[i];
            let mut best_comm = current_comm;
            let mut best_gain = 0.0f64;

            // Try moving node i to each neighbor's community
            let neighbor_comms: Vec<usize> = graph[i]
                .iter()
                .map(|&(j, _)| communities[j])
                .collect();

            for &candidate_comm in &neighbor_comms {
                if candidate_comm == current_comm { continue; }

                // Compute modularity gain of moving i to candidate_comm
                communities[i] = candidate_comm;
                let new_q = modularity(graph, &communities);
                communities[i] = current_comm;
                let old_q = modularity(graph, &communities);
                let gain = new_q - old_q;

                if gain > best_gain {
                    best_gain = gain;
                    best_comm = candidate_comm;
                }
            }

            if best_comm != current_comm {
                communities[i] = best_comm;
                improved = true;
            }
        }
    }

    // Relabel communities to be contiguous 0, 1, 2, ...
    let mut label_map: HashMap<usize, usize> = HashMap::new();
    let mut next_label = 0;
    for c in &mut communities {
        let new_label = *label_map.entry(*c).or_insert_with(|| {
            let l = next_label;
            next_label += 1;
            l
        });
        *c = new_label;
    }
    communities
}

fn main() {
    println!("=== Louvain Clustering: Cell Type Discovery ===\n");

    // Simulated 2D cell embeddings (e.g., from UMAP)
    // Three clusters: T cells, B cells, Monocytes
    let data: Vec<Vec<f64>> = vec![
        // T cells (cluster around (2, 2))
        vec![1.8, 2.1], vec![2.2, 1.9], vec![2.0, 2.3], vec![1.9, 1.8], vec![2.1, 2.2],
        // B cells (cluster around (8, 8))
        vec![7.9, 8.1], vec![8.2, 7.8], vec![8.0, 8.3], vec![7.8, 8.0], vec![8.1, 7.9],
        // Monocytes (cluster around (2, 8))
        vec![1.9, 7.8], vec![2.1, 8.2], vec![2.0, 8.0], vec![1.8, 8.1], vec![2.2, 7.9],
    ];
    let labels = ["T","T","T","T","T","B","B","B","B","B","M","M","M","M","M"];

    let k = 4;
    println!("Building KNN graph (k={})...", k);
    let graph = build_knn_graph(&data, k);

    let start = Instant::now();
    let clusters = louvain_clustering(&graph, 50);
    let elapsed = start.elapsed();

    let q = modularity(&graph, &clusters);
    println!("Clustering complete: {} clusters, modularity={:.4}, time={:?}\n",
             clusters.iter().collect::<std::collections::HashSet<_>>().len(), q, elapsed);

    for (i, (&cluster, &label)) in clusters.iter().zip(labels.iter()).enumerate() {
        println!("  Cell {:2}: cluster {} (true: {})", i, cluster, label);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_clusters() {
        let data = vec![
            vec![0.0, 0.0], vec![0.1, 0.1], vec![0.2, 0.0],
            vec![10.0, 10.0], vec![10.1, 10.1], vec![9.9, 10.0],
        ];
        let graph = build_knn_graph(&data, 2);
        let clusters = louvain_clustering(&graph, 20);
        assert_eq!(clusters[0], clusters[1]);
        assert_eq!(clusters[3], clusters[4]);
        assert_ne!(clusters[0], clusters[3]);
    }
}
