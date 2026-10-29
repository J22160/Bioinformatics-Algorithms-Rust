//! Gene Set Enrichment Analysis (GSEA)
//!
//! Tests whether predefined sets of genes show statistically
//! significant enrichment at the top or bottom of a ranked list.

/// Compute the enrichment score for a gene set against a ranked list.
///
/// gene_set_mask: boolean mask where true = gene is in the set
/// ranked_scores: the ranking metric (e.g., log2FC × -log10(p))
/// Returns: (enrichment_score, running_scores for plotting)
pub fn enrichment_score(
    gene_set_mask: &[bool],
    ranked_scores: &[f64],
) -> (f64, Vec<f64>) {
    let n = ranked_scores.len();
    let n_hit: usize = gene_set_mask.iter().filter(|&&b| b).count();
    let n_miss = n - n_hit;

    if n_hit == 0 || n_miss == 0 {
        return (0.0, vec![0.0; n]);
    }

    // Weighted increment for hits
    let nr: f64 = gene_set_mask.iter()
        .zip(ranked_scores.iter())
        .filter(|(&in_set, _)| in_set)
        .map(|(_, &score)| score.abs())
        .sum();

    let miss_penalty = 1.0 / n_miss as f64;
    let mut running = Vec::with_capacity(n);
    let mut current = 0.0f64;
    let mut max_dev = 0.0f64;

    for i in 0..n {
        if gene_set_mask[i] {
            current += ranked_scores[i].abs() / nr;
        } else {
            current -= miss_penalty;
        }
        running.push(current);
        if current.abs() > max_dev.abs() {
            max_dev = current;
        }
    }

    (max_dev, running)
}

/// Simple permutation test for enrichment significance.
pub fn permutation_test(
    gene_set_mask: &[bool],
    ranked_scores: &[f64],
    n_perms: usize,
    seed: u64,
) -> (f64, f64) {
    let (observed_es, _) = enrichment_score(gene_set_mask, ranked_scores);

    let n = gene_set_mask.len();
    let mut state = seed;
    let mut count_extreme = 0usize;

    for _ in 0..n_perms {
        // Fisher-Yates shuffle of mask
        let mut perm_mask = gene_set_mask.to_vec();
        for i in (1..n).rev() {
            state = state.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let j = (state >> 33) as usize % (i + 1);
            perm_mask.swap(i, j);
        }
        let (perm_es, _) = enrichment_score(&perm_mask, ranked_scores);
        if perm_es.abs() >= observed_es.abs() {
            count_extreme += 1;
        }
    }

    let p_value = (count_extreme + 1) as f64 / (n_perms + 1) as f64;
    (observed_es, p_value)
}

fn main() {
    println!("=== GSEA: Gene Set Enrichment Analysis ===\n");

    let gene_names = vec![
        "IL6", "TNF", "IFNG", "CXCL8", "IL1B",     // inflammatory
        "GAPDH", "ACTB", "TUBB", "HSP90", "RPS6",   // housekeeping
        "MYC", "KRAS", "TP53", "BRCA1", "EGFR",      // cancer
    ];

    // Ranked scores (e.g., log2FC × sign)
    let ranked_scores: Vec<f64> = vec![
        3.2, 2.8, 2.5, 2.1, 1.9,    // inflammatory genes ranked high
        0.1, -0.05, 0.02, -0.1, 0.0, // housekeeping near zero
        -0.5, -0.3, 0.4, 0.6, -0.2,  // mixed cancer genes
    ];

    // Gene set: "Inflammatory Response"
    let inflammatory_set: Vec<bool> = vec![
        true, true, true, true, true,
        false, false, false, false, false,
        false, false, false, false, false,
    ];

    let (es, running) = enrichment_score(&inflammatory_set, &ranked_scores);
    println!("Inflammatory Response pathway:");
    println!("  Enrichment Score: {:.4}", es);

    let (es_perm, p_val) = permutation_test(&inflammatory_set, &ranked_scores, 10000, 42);
    println!("  Permutation p-value: {:.4} (n_perm=10000)", p_val);
    println!("  Significant: {}", if p_val < 0.05 { "YES" } else { "NO" });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_enrichment() {
        let mask = vec![true, true, true, false, false, false];
        let scores = vec![3.0, 2.0, 1.0, -1.0, -2.0, -3.0];
        let (es, _) = enrichment_score(&mask, &scores);
        assert!(es > 0.0);
    }

    #[test]
    fn test_no_enrichment() {
        let mask = vec![true, false, true, false, true, false];
        let scores = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let (es, _) = enrichment_score(&mask, &scores);
        assert!(es.abs() < 0.5);
    }
}
