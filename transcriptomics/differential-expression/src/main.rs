//! Negative Binomial GLM for Differential Expression
//!
//! Simplified implementation of the core DESeq2 statistical model.

use std::time::Instant;

/// Compute size factors using the median-of-ratios method (DESeq2).
///
/// For each gene, compute a "pseudo-reference" as the geometric mean
/// across samples. Then for each sample, the size factor is the
/// median of (count / pseudo-reference) across all genes.
pub fn estimate_size_factors(counts: &[Vec<u32>]) -> Vec<f64> {
    let n_genes = counts.len();
    let n_samples = counts[0].len();

    // Compute geometric mean per gene (pseudo-reference)
    let mut geo_means = vec![0.0f64; n_genes];
    for i in 0..n_genes {
        let log_sum: f64 = counts[i].iter()
            .map(|&c| if c > 0 { (c as f64).ln() } else { 0.0 })
            .sum();
        let n_nonzero = counts[i].iter().filter(|&&c| c > 0).count();
        if n_nonzero == n_samples {
            geo_means[i] = (log_sum / n_samples as f64).exp();
        }
    }

    // For each sample, compute median of (count / geo_mean)
    let mut size_factors = vec![1.0; n_samples];
    for j in 0..n_samples {
        let mut ratios: Vec<f64> = Vec::new();
        for i in 0..n_genes {
            if geo_means[i] > 0.0 && counts[i][j] > 0 {
                ratios.push(counts[i][j] as f64 / geo_means[i]);
            }
        }
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !ratios.is_empty() {
            size_factors[j] = ratios[ratios.len() / 2];
        }
    }
    size_factors
}

/// Compute normalized counts.
pub fn normalize_counts(counts: &[Vec<u32>], size_factors: &[f64]) -> Vec<Vec<f64>> {
    counts.iter().map(|gene_counts| {
        gene_counts.iter().zip(size_factors.iter())
            .map(|(&c, &sf)| c as f64 / sf)
            .collect()
    }).collect()
}

/// Compute log2 fold change between two groups using normalized counts.
pub fn compute_log2fc(
    norm_counts: &[Vec<f64>],
    group_a: &[usize],  // sample indices for condition A
    group_b: &[usize],  // sample indices for condition B
) -> Vec<f64> {
    norm_counts.iter().map(|gene| {
        let mean_a: f64 = group_a.iter().map(|&j| gene[j]).sum::<f64>()
            / group_a.len() as f64;
        let mean_b: f64 = group_b.iter().map(|&j| gene[j]).sum::<f64>()
            / group_b.len() as f64;
        if mean_a > 0.0 && mean_b > 0.0 {
            (mean_b / mean_a).log2()
        } else {
            0.0
        }
    }).collect()
}

/// Estimate dispersion for each gene using method of moments.
/// This is a simplified version; DESeq2 uses Cox-Reid adjusted
/// profile likelihood with shrinkage toward a fitted trend.
pub fn estimate_dispersions(
    counts: &[Vec<u32>],
    size_factors: &[f64],
    group_a: &[usize],
    group_b: &[usize],
) -> Vec<f64> {
    let n_genes = counts.len();
    let mut dispersions = vec![0.1; n_genes];

    for i in 0..n_genes {
        let mut all_mu = Vec::new();
        let mut all_var = Vec::new();

        for group in [group_a, group_b] {
            let norm: Vec<f64> = group.iter()
                .map(|&j| counts[i][j] as f64 / size_factors[j])
                .collect();
            let mu = norm.iter().sum::<f64>() / norm.len() as f64;
            if norm.len() > 1 && mu > 0.0 {
                let var = norm.iter()
                    .map(|&x| (x - mu).powi(2))
                    .sum::<f64>() / (norm.len() - 1) as f64;
                all_mu.push(mu);
                all_var.push(var);
            }
        }

        if !all_mu.is_empty() {
            let avg_mu: f64 = all_mu.iter().sum::<f64>() / all_mu.len() as f64;
            let avg_var: f64 = all_var.iter().sum::<f64>() / all_var.len() as f64;
            // NB variance = mu + alpha * mu^2
            // So alpha = (var - mu) / mu^2
            if avg_mu > 0.0 {
                let alpha = (avg_var - avg_mu) / (avg_mu * avg_mu);
                dispersions[i] = alpha.max(0.01); // floor at 0.01
            }
        }
    }
    dispersions
}

/// Simple Wald-like test for differential expression.
/// Returns (log2FC, p-value) for each gene.
pub fn wald_test(
    norm_counts: &[Vec<f64>],
    group_a: &[usize],
    group_b: &[usize],
) -> Vec<(f64, f64)> {
    let n_genes = norm_counts.len();
    let na = group_a.len() as f64;
    let nb = group_b.len() as f64;

    (0..n_genes).map(|i| {
        let vals_a: Vec<f64> = group_a.iter().map(|&j| norm_counts[i][j]).collect();
        let vals_b: Vec<f64> = group_b.iter().map(|&j| norm_counts[i][j]).collect();

        let mean_a: f64 = vals_a.iter().sum::<f64>() / na;
        let mean_b: f64 = vals_b.iter().sum::<f64>() / nb;

        let log2fc = if mean_a > 0.0 && mean_b > 0.0 {
            (mean_b / mean_a).log2()
        } else { 0.0 };

        let var_a = vals_a.iter().map(|&x| (x - mean_a).powi(2)).sum::<f64>()
            / (na - 1.0).max(1.0);
        let var_b = vals_b.iter().map(|&x| (x - mean_b).powi(2)).sum::<f64>()
            / (nb - 1.0).max(1.0);

        let se = ((var_a / na) + (var_b / nb)).sqrt();
        let stat = if se > 0.0 { (mean_b - mean_a) / se } else { 0.0 };

        // Approximate p-value using normal distribution
        let p_value = 2.0 * normal_cdf(-stat.abs());
        (log2fc, p_value)
    }).collect()
}

/// Standard normal CDF approximation (Abramowitz & Stegun).
fn normal_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let d = 0.3989422804014327; // 1/sqrt(2*pi)
    let p = d * (-x * x / 2.0).exp()
        * (t * (0.319381530
            + t * (-0.356563782
                + t * (1.781477937
                    + t * (-1.821255978
                        + t * 1.330274429)))));
    if x >= 0.0 { 1.0 - p } else { p }
}

/// Benjamini-Hochberg FDR correction.
pub fn bh_adjust(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    let mut indexed: Vec<(usize, f64)> = p_values.iter()
        .enumerate()
        .map(|(i, &p)| (i, p))
        .collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut adjusted = vec![1.0; n];
    let mut min_so_far = 1.0f64;

    for rank in (0..n).rev() {
        let (orig_idx, p) = indexed[rank];
        let adj = (p * n as f64 / (rank + 1) as f64).min(1.0);
        min_so_far = min_so_far.min(adj);
        adjusted[orig_idx] = min_so_far;
    }
    adjusted
}

fn main() {
    println!("=== Differential Expression: Negative Binomial Model ===\n");

    // Simulated count matrix: 10 genes × 6 samples (3 ctrl, 3 treated)
    let gene_names = vec![
        "BRCA1", "TP53", "MYC", "GAPDH", "ACTB",
        "IL6", "TNF", "EGFR", "VEGFA", "KRAS",
    ];

    // Some genes DE (BRCA1, MYC, IL6, TNF), others not
    let counts: Vec<Vec<u32>> = vec![
        vec![150, 130, 160, 450, 420, 480],   // BRCA1 - UP
        vec![500, 520, 490, 510, 530, 505],   // TP53 - no change
        vec![300, 280, 310, 80,  95,  70],     // MYC - DOWN
        vec![5000, 4800, 5200, 4900, 5100, 5000], // GAPDH - housekeeping
        vec![3000, 3100, 2900, 3050, 2950, 3100], // ACTB - housekeeping
        vec![20,  25,  18,  200, 180, 220],   // IL6 - UP (strong)
        vec![50,  45,  55,  250, 230, 270],   // TNF - UP (strong)
        vec![800, 780, 820, 790, 810, 800],   // EGFR - no change
        vec![400, 380, 420, 410, 390, 430],   // VEGFA - no change
        vec![600, 620, 580, 610, 590, 605],   // KRAS - no change
    ];

    let group_ctrl = vec![0usize, 1, 2];
    let group_treat = vec![3usize, 4, 5];

    // Step 1: Size factors
    let start = Instant::now();
    let sf = estimate_size_factors(&counts);
    println!("Size factors: {:?}", sf.iter().map(|s| format!("{:.3}", s)).collect::<Vec<_>>());

    // Step 2: Normalize
    let norm = normalize_counts(&counts, &sf);

    // Step 3: Test
    let results = wald_test(&norm, &group_ctrl, &group_treat);
    let p_values: Vec<f64> = results.iter().map(|r| r.1).collect();
    let padj = bh_adjust(&p_values);
    let elapsed = start.elapsed();

    println!("\nDifferential Expression Results ({:?}):\n", elapsed);
    println!("{:<8} {:>10} {:>10} {:>10} {:>6}",
             "Gene", "log2FC", "pvalue", "padj", "Sig?");
    println!("{}", "-".repeat(54));

    for (i, name) in gene_names.iter().enumerate() {
        let (lfc, pval) = results[i];
        let sig = if padj[i] < 0.05 { "***" }
        else if padj[i] < 0.1 { "*" }
        else { "" };
        println!("{:<8} {:>10.3} {:>10.2e} {:>10.2e} {:>6}",
                 name, lfc, pval, padj[i], sig);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_factors() {
        let counts = vec![
            vec![100, 200],
            vec![50, 100],
        ];
        let sf = estimate_size_factors(&counts);
        assert!(sf[1] / sf[0] > 1.5);
    }

    #[test]
    fn test_bh_correction() {
        let p = vec![0.01, 0.04, 0.03, 0.20];
        let adj = bh_adjust(&p);
        assert!(adj[0] <= adj[2]);
        for &a in &adj { assert!(a <= 1.0); }
    }
}