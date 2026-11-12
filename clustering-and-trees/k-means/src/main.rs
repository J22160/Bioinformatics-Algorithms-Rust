//! K-Means Clustering — Gene Expression Grouping

use std::time::Instant;

/// Compute Euclidean distance between two points.
fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// K-Means clustering using Lloyd's algorithm.
///
/// Returns (assignments, centers, iterations, wcss).
pub fn kmeans(
    data: &[Vec<f64>],
    k: usize,
    max_iter: usize,
    seed: u64,
) -> (Vec<usize>, Vec<Vec<f64>>, usize, f64) {
    let n = data.len();
    let d = data[0].len();

    // Initialize centers using a simple deterministic spread
    let mut centers: Vec<Vec<f64>> = Vec::with_capacity(k);
    let mut state = seed;
    for _ in 0..k {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let idx = (state >> 33) as usize % n;
        centers.push(data[idx].clone());
    }

    let mut assignments = vec![0usize; n];
    let mut iterations = 0;

    for iter in 0..max_iter {
        iterations = iter + 1;
        let mut changed = false;

        // Assignment step
        for i in 0..n {
            let mut best_cluster = 0;
            let mut best_dist = f64::INFINITY;
            for c in 0..k {
                let dist = euclidean_dist(&data[i], &centers[c]);
                if dist < best_dist {
                    best_dist = dist;
                    best_cluster = c;
                }
            }
            if assignments[i] != best_cluster {
                assignments[i] = best_cluster;
                changed = true;
            }
        }

        if !changed { break; }

        // Update step: recompute centers
        let mut new_centers = vec![vec![0.0; d]; k];
        let mut counts = vec![0usize; k];

        for i in 0..n {
            let c = assignments[i];
            counts[c] += 1;
            for j in 0..d {
                new_centers[c][j] += data[i][j];
            }
        }

        for c in 0..k {
            if counts[c] > 0 {
                for j in 0..d {
                    new_centers[c][j] /= counts[c] as f64;
                }
            }
        }
        centers = new_centers;
    }

    // Compute WCSS
    let wcss: f64 = (0..n)
        .map(|i| euclidean_dist(&data[i], &centers[assignments[i]]).powi(2))
        .sum();

    (assignments, centers, iterations, wcss)
}

fn main() {
    println!("=== K-Means: Gene Expression Clustering ===\n");

    // Simulated gene expression: 3 groups of genes
    let data: Vec<Vec<f64>> = vec![
        // Group 1: high in condition 1
        vec![10.0, 1.0, 1.0], vec![9.5, 1.5, 0.5], vec![10.5, 0.5, 1.5],
        vec![9.0, 2.0, 1.0], vec![10.2, 1.2, 0.8],
        // Group 2: high in condition 2
        vec![1.0, 10.0, 1.0], vec![0.5, 9.5, 1.5], vec![1.5, 10.5, 0.5],
        vec![2.0, 9.0, 1.0], vec![1.2, 10.2, 0.8],
        // Group 3: high in condition 3
        vec![1.0, 1.0, 10.0], vec![1.5, 0.5, 9.5], vec![0.5, 1.5, 10.5],
        vec![1.0, 2.0, 9.0], vec![0.8, 1.2, 10.2],
    ];

    let gene_names: Vec<&str> = (0..15).map(|i| {
        if i < 5 { "GroupA" } else if i < 10 { "GroupB" } else { "GroupC" }
    }).collect();

    let start = Instant::now();
    let (assignments, centers, iters, wcss) = kmeans(&data, 3, 100, 42);
    let elapsed = start.elapsed();

    println!("Converged in {} iterations ({:?})", iters, elapsed);
    println!("WCSS: {:.2}\n", wcss);

    for c in 0..3 {
        let members: Vec<usize> = (0..data.len())
            .filter(|&i| assignments[i] == c)
            .collect();
        println!("Cluster {}: {} genes", c, members.len());
        for &m in &members {
            println!("  Gene {} ({}): {:?}", m, gene_names[m],
                     data[m].iter().map(|v| format!("{:.1}", v)).collect::<Vec<_>>());
        }
    }

    // Elbow method
    println!("\n=== Elbow Method ===");
    for k in 1..=6 {
        let (_, _, _, wcss) = kmeans(&data, k, 100, 42);
        println!("  k={}: WCSS={:.2}", k, wcss);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_clusters() {
        let data = vec![
            vec![0.0, 0.0], vec![0.1, 0.1],
            vec![10.0, 10.0], vec![10.1, 10.1],
        ];
        let (asgn, _, _, _) = kmeans(&data, 2, 50, 42);
        assert_eq!(asgn[0], asgn[1]);
        assert_eq!(asgn[2], asgn[3]);
        assert_ne!(asgn[0], asgn[2]);
    }
}
