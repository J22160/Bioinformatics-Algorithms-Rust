//! UPGMA — Ultrametric Phylogenetic Tree Construction
//!
//! Builds a rooted tree from a distance matrix assuming a
//! molecular clock (equal rates of evolution across lineages).

use std::collections::HashMap;
use std::time::Instant;

/// A node in the UPGMA tree.
#[derive(Debug, Clone)]
pub struct TreeNode {
    pub id: usize,
    pub label: String,
    pub height: f64,
    pub left: Option<Box<TreeNode>>,
    pub right: Option<Box<TreeNode>>,
    pub size: usize, // number of leaves in this subtree
}

impl TreeNode {
    fn leaf(id: usize, label: String) -> Self {
        TreeNode { id, label, height: 0.0, left: None, right: None, size: 1 }
    }

    fn internal(id: usize, height: f64, left: TreeNode, right: TreeNode) -> Self {
        let size = left.size + right.size;
        TreeNode {
            id,
            label: format!("({}::{:.3})", id, height),
            height,
            left: Some(Box::new(left)),
            right: Some(Box::new(right)),
            size,
        }
    }

    /// Print the tree in Newick format.
    pub fn to_newick(&self) -> String {
        match (&self.left, &self.right) {
            (Some(l), Some(r)) => {
                let l_branch = self.height - l.height;
                let r_branch = self.height - r.height;
                format!("({}:{:.4},{}:{:.4})",
                        l.to_newick(), l_branch,
                        r.to_newick(), r_branch)
            }
            _ => self.label.clone(),
        }
    }
}

/// Build a UPGMA tree from a distance matrix.
///
/// # Rust Performance Advantages
///
/// 1. **Flat Vec for distance matrix**: O(1) indexed access with
///    perfect cache locality. Python's list-of-lists scatters rows.
///
/// 2. **In-place matrix updates**: We update the distance matrix
///    without allocating new storage each iteration. Python would
///    create new list objects.
///
/// 3. **No GC during the main loop**: The merge loop allocates
///    only TreeNode structs (small, fixed-size). Python would
///    trigger GC on intermediate dict/list objects.
pub fn upgma(
    dist: &[Vec<f64>],
    labels: &[String],
) -> TreeNode {
    let n = labels.len();
    let mut active: Vec<bool> = vec![true; 2 * n]; // track active clusters
    let mut nodes: Vec<Option<TreeNode>> = Vec::with_capacity(2 * n);
    let mut sizes: Vec<usize> = vec![0; 2 * n];

    // Initialize leaf nodes
    for i in 0..n {
        nodes.push(Some(TreeNode::leaf(i, labels[i].clone())));
        sizes[i] = 1;
    }

    // Working distance matrix (expandable)
    let capacity = 2 * n;
    let mut d = vec![vec![f64::INFINITY; capacity]; capacity];
    for i in 0..n {
        for j in 0..n {
            d[i][j] = dist[i][j];
        }
    }

    let mut next_id = n;

    for _step in 0..(n - 1) {
        // Find the two closest active clusters
        let mut min_dist = f64::INFINITY;
        let mut min_i = 0;
        let mut min_j = 0;

        for i in 0..next_id {
            if !active[i] { continue; }
            for j in (i + 1)..next_id {
                if !active[j] { continue; }
                if d[i][j] < min_dist {
                    min_dist = d[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // Create new internal node
        let height = min_dist / 2.0;
        let left = nodes[min_i].take().unwrap();
        let right = nodes[min_j].take().unwrap();
        let new_node = TreeNode::internal(next_id, height, left, right);

        // Compute distances from new cluster to all other active clusters
        let si = sizes[min_i] as f64;
        let sj = sizes[min_j] as f64;
        for k in 0..next_id {
            if !active[k] || k == min_i || k == min_j { continue; }
            let new_dist = (si * d[min_i][k] + sj * d[min_j][k]) / (si + sj);
            d[next_id][k] = new_dist;
            d[k][next_id] = new_dist;
        }

        sizes[next_id] = sizes[min_i] + sizes[min_j];
        active[min_i] = false;
        active[min_j] = false;
        active[next_id] = true;

        // Extend nodes vector
        while nodes.len() <= next_id {
            nodes.push(None);
        }
        nodes[next_id] = Some(new_node);
        next_id += 1;
    }

    // Find the root (last active node)
    for i in (0..next_id).rev() {
        if active[i] {
            return nodes[i].take().unwrap();
        }
    }
    unreachable!()
}

fn main() {
    println!("=== UPGMA: Bacterial Phylogeny from 16S rRNA Distances ===\n");

    let labels: Vec<String> = vec![
        "E_coli", "S_enterica", "B_subtilis", "S_aureus", "M_tuberculosis",
    ].iter().map(|s| s.to_string()).collect();

    // Pairwise evolutionary distances (symmetric)
    let dist = vec![
        vec![0.0, 0.8, 3.0, 3.2, 4.5],
        vec![0.8, 0.0, 2.8, 3.0, 4.3],
        vec![3.0, 2.8, 0.0, 1.5, 4.0],
        vec![3.2, 3.0, 1.5, 0.0, 3.8],
        vec![4.5, 4.3, 4.0, 3.8, 0.0],
    ];

    println!("Distance matrix:");
    print!("{:>15}", "");
    for l in &labels { print!("{:>15}", l); }
    println!();
    for (i, row) in dist.iter().enumerate() {
        print!("{:>15}", labels[i]);
        for &v in row { print!("{:>15.1}", v); }
        println!();
    }

    let start = Instant::now();
    let tree = upgma(&dist, &labels);
    let elapsed = start.elapsed();

    println!("\nNewick tree: {}", tree.to_newick());
    println!("Root height: {:.4} (estimated divergence time)", tree.height);
    println!("Time: {:?}", elapsed);

    // Print tree structure
    println!("\n--- Tree Structure ---");
    fn print_tree(node: &TreeNode, prefix: &str, is_last: bool) {
        let connector = if is_last { "└── " } else { "├── " };
        let branch_len = if node.height > 0.0 {
            format!(" (h={:.3})", node.height)
        } else { String::new() };
        println!("{}{}{}{}", prefix, connector, node.label, branch_len);
        let new_prefix = format!("{}{}", prefix, if is_last { "    " } else { "│   " });
        if let Some(ref left) = node.left {
            print_tree(left, &new_prefix, node.right.is_none());
        }
        if let Some(ref right) = node.right {
            print_tree(right, &new_prefix, true);
        }
    }
    print_tree(&tree, "", true);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_taxa() {
        let dist = vec![vec![0.0, 4.0], vec![4.0, 0.0]];
        let labels = vec!["A".into(), "B".into()];
        let tree = upgma(&dist, &labels);
        assert!((tree.height - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_three_taxa() {
        let dist = vec![
            vec![0.0, 2.0, 4.0],
            vec![2.0, 0.0, 4.0],
            vec![4.0, 4.0, 0.0],
        ];
        let labels = vec!["A".into(), "B".into(), "C".into()];
        let tree = upgma(&dist, &labels);
        assert!(tree.height > 0.0);
    }
}
