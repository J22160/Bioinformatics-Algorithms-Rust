//! Dimensionality Reduction: PCA for scRNA-seq
//!
//! Implements PCA via power iteration for the top principal components.
//! This is the first step in most scRNA-seq analysis pipelines.

use std::time::Instant;

/// Center the data (subtract column means).
fn center_data(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = data.len();
    let d = data[0].len();
    let mut means = vec![0.0; d];
    for row in data {
        for (j, &val) in row.iter().enumerate() {
            means[j] += val;
        }
    }
    for m in &mut means { *m /= n as f64; }
    data.iter()
        .map(|row| row.iter().zip(means.iter()).map(|(&v, &m)| v - m).collect())
        .collect()
}

/// Compute dot product of two vectors.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Normalize a vector to unit length.
fn normalize(v: &mut Vec<f64>) {
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for x in v.iter_mut() { *x /= norm; }
    }
}

/// Compute top-k principal components via power iteration.
///
/// Returns (components, explained_variance).
pub fn pca(data: &[Vec<f64>], k: usize, max_iter: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let centered = center_data(data);
    let n = centered.len();
    let d = centered[0].len();

    let mut components = Vec::new();
    let mut variances = Vec::new();
    let mut residual = centered.clone();

    for _comp in 0..k {
        // Power iteration to find dominant eigenvector
        let mut pc = vec![1.0 / (d as f64).sqrt(); d];
        normalize(&mut pc);

        for _ in 0..max_iter {
            // Compute X^T * X * pc
            let mut projections: Vec<f64> = residual.iter()
                .map(|row| dot(row, &pc))
                .collect();
            let mut new_pc = vec![0.0; d];
            for (i, row) in residual.iter().enumerate() {
                for (j, &val) in row.iter().enumerate() {
                    new_pc[j] += val * projections[i];
                }
            }
            normalize(&mut new_pc);
            pc = new_pc;
        }

        // Compute variance explained by this component
        let projections: Vec<f64> = residual.iter().map(|row| dot(row, &pc)).collect();
        let variance: f64 = projections.iter().map(|p| p * p).sum::<f64>() / n as f64;
        variances.push(variance);

        // Deflate: remove this component from the data
        for (i, row) in residual.iter_mut().enumerate() {
            for j in 0..d {
                row[j] -= projections[i] * pc[j];
            }
        }

        components.push(pc);
    }

    (components, variances)
}

/// Project data onto principal components.
pub fn project(data: &[Vec<f64>], components: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let centered = center_data(data);
    centered.iter()
        .map(|row| components.iter().map(|pc| dot(row, pc)).collect())
        .collect()
}

fn main() {
    println!("=== PCA: Dimensionality Reduction for scRNA-seq ===\n");

    // Simulated gene expression: 15 cells, 8 genes
    // 3 cell types with distinct expression patterns
    let data: Vec<Vec<f64>> = vec![
        // T cells: high CD3D, CD3E; low CD14, CD19
        vec![10.0, 9.5, 0.5, 0.3, 8.0, 7.5, 1.0, 0.8],
        vec![9.8, 10.2, 0.3, 0.4, 7.8, 8.2, 0.9, 0.7],
        vec![10.5, 9.0, 0.6, 0.2, 8.5, 7.0, 1.1, 0.6],
        vec![9.2, 9.8, 0.4, 0.5, 7.2, 7.8, 0.8, 0.9],
        vec![10.1, 10.0, 0.5, 0.3, 8.1, 8.0, 1.0, 0.7],
        // B cells: high CD19, MS4A1; low CD3D, CD14
        vec![0.4, 0.3, 9.0, 10.0, 0.5, 0.4, 8.5, 9.5],
        vec![0.5, 0.2, 9.5, 9.8, 0.4, 0.3, 9.0, 10.0],
        vec![0.3, 0.4, 8.8, 10.2, 0.6, 0.5, 8.2, 9.8],
        vec![0.6, 0.3, 9.2, 9.5, 0.5, 0.4, 8.8, 9.2],
        vec![0.4, 0.5, 9.0, 10.1, 0.4, 0.3, 8.5, 9.7],
        // Monocytes: high CD14, LYZ; low CD3D, CD19
        vec![0.3, 0.5, 0.4, 0.3, 10.0, 9.5, 0.5, 0.4],
        vec![0.4, 0.3, 0.5, 0.4, 9.8, 10.2, 0.6, 0.3],
        vec![0.2, 0.4, 0.3, 0.5, 10.5, 9.0, 0.4, 0.5],
        vec![0.5, 0.6, 0.4, 0.3, 9.2, 9.8, 0.5, 0.4],
        vec![0.3, 0.4, 0.5, 0.4, 10.1, 10.0, 0.6, 0.3],
    ];

    let gene_names = ["CD3D", "CD3E", "CD19", "MS4A1", "CD14", "LYZ", "IGHM", "CD79A"];
    let cell_types = ["T","T","T","T","T","B","B","B","B","B","M","M","M","M","M"];

    let start = Instant::now();
    let (components, variances) = pca(&data, 3, 100);
    let projected = project(&data, &components);
    let elapsed = start.elapsed();

    let total_var: f64 = variances.iter().sum();
    println!("Top 3 PCs (time: {:?}):", elapsed);
    for (i, var) in variances.iter().enumerate() {
        println!("  PC{}: {:.1}% variance explained", i + 1, var / total_var * 100.0);
        let mut loadings: Vec<(usize, f64)> = components[i].iter()
            .enumerate().map(|(j, &v)| (j, v.abs())).collect();
        loadings.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        print!("    Top genes: ");
        for &(j, _) in loadings.iter().take(3) {
            print!("{} ", gene_names[j]);
        }
        println!();
    }

    println!("\nProjected cells (PC1, PC2):");
    for (i, proj) in projected.iter().enumerate() {
        println!("  Cell {:2} ({}): ({:+.2}, {:+.2})",
                 i, cell_types[i], proj[0], proj[1]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pca_reduces_dimensions() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0], vec![7.0, 8.0, 9.0]];
        let (components, _) = pca(&data, 2, 50);
        assert_eq!(components.len(), 2);
        assert_eq!(components[0].len(), 3);
    }
}
