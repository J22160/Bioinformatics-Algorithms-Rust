//! Neighbor Joining — Unrooted Phylogenetic Tree Construction

use std::time::Instant;

/// Edge in the NJ tree
#[derive(Debug)]
pub struct NJEdge {
    pub from: usize,
    pub to: usize,
    pub length: f64,
}

/// Build a neighbor-joining tree. Returns a list of edges.
pub fn neighbor_joining(
    dist_matrix: &[Vec<f64>],
    labels: &[&str],
) -> (Vec<NJEdge>, Vec<String>) {
    let n = labels.len();
    let mut d: Vec<Vec<f64>> = dist_matrix.to_vec();
    let mut active: Vec<usize> = (0..n).collect();
    let mut all_labels: Vec<String> = labels.iter().map(|s| s.to_string()).collect();
    let mut edges: Vec<NJEdge> = Vec::new();
    let mut next_id = n;

    while active.len() > 2 {
        let m = active.len();

        // Compute u(i) for each active node
        let mut u: Vec<f64> = vec![0.0; d.len()];
        for &i in &active {
            let sum: f64 = active.iter()
                .filter(|&&j| j != i)
                .map(|&j| d[i][j])
                .sum();
            u[i] = sum / (m as f64 - 2.0).max(1.0);
        }

        // Find pair minimizing D(i,j) - u(i) - u(j)
        let mut best_val = f64::INFINITY;
        let mut best_i = 0;
        let mut best_j = 0;
        for idx_a in 0..m {
            for idx_b in (idx_a + 1)..m {
                let i = active[idx_a];
                let j = active[idx_b];
                let val = d[i][j] - u[i] - u[j];
                if val < best_val {
                    best_val = val;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        // Compute branch lengths
        let dist_ij = d[best_i][best_j];
        let len_i = 0.5 * dist_ij + 0.5 * (u[best_i] - u[best_j]);
        let len_j = 0.5 * dist_ij + 0.5 * (u[best_j] - u[best_i]);

        // Create new node
        let new_id = next_id;
        next_id += 1;
        all_labels.push(format!("Int_{}", new_id));

        // Expand distance matrix
        let new_size = d.len().max(new_id + 1);
        while d.len() <= new_id {
            d.push(vec![0.0; new_size]);
        }
        for row in &mut d {
            while row.len() <= new_id {
                row.push(0.0);
            }
        }

        // Compute distances from new node to all other active nodes
        for &k in &active {
            if k == best_i || k == best_j { continue; }
            let new_dist = (d[best_i][k] + d[best_j][k] - dist_ij) / 2.0;
            d[new_id][k] = new_dist;
            d[k][new_id] = new_dist;
        }

        edges.push(NJEdge { from: new_id, to: best_i, length: len_i });
        edges.push(NJEdge { from: new_id, to: best_j, length: len_j });

        // Update active set
        active.retain(|&x| x != best_i && x != best_j);
        active.push(new_id);
    }

    // Connect the final two nodes
    if active.len() == 2 {
        let i = active[0];
        let j = active[1];
        edges.push(NJEdge { from: i, to: j, length: d[i][j] });
    }

    (edges, all_labels)
}

fn main() {
    println!("=== Neighbor Joining: Mammalian Phylogeny ===\n");

    let labels = vec!["Human", "Chimp", "Gorilla", "Mouse", "Rat", "Dog"];
    let dist = vec![
        vec![0.0, 1.2, 1.8, 6.0, 6.2, 5.0],
        vec![1.2, 0.0, 1.6, 5.8, 6.0, 4.8],
        vec![1.8, 1.6, 0.0, 6.2, 6.4, 5.2],
        vec![6.0, 5.8, 6.2, 0.0, 1.0, 5.5],
        vec![6.2, 6.0, 6.4, 1.0, 0.0, 5.7],
        vec![5.0, 4.8, 5.2, 5.5, 5.7, 0.0],
    ];

    let start = Instant::now();
    let (edges, all_labels) = neighbor_joining(&dist, &labels.iter().map(|s| *s).collect::<Vec<_>>());
    let elapsed = start.elapsed();

    println!("Edges:");
    for e in &edges {
        println!("  {} -- {} (length {:.4})",
                 all_labels[e.from], all_labels[e.to], e.length);
    }
    println!("\nTime: {:?}", elapsed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_three_taxa() {
        let dist = vec![
            vec![0.0, 2.0, 4.0],
            vec![2.0, 0.0, 4.0],
            vec![4.0, 4.0, 0.0],
        ];
        let labels = vec!["A", "B", "C"];
        let (edges, _) = neighbor_joining(&dist, &labels);
        assert!(!edges.is_empty());
    }
}
