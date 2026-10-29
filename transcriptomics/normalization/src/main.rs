//! RNA-seq Normalization Methods

/// Median of Ratios normalization (DESeq2 method).
pub fn median_of_ratios(counts: &[Vec<u32>]) -> Vec<f64> {
    let n_genes = counts.len();
    let n_samples = counts[0].len();

    let mut geo_means = vec![0.0f64; n_genes];
    for i in 0..n_genes {
        if counts[i].iter().all(|&c| c > 0) {
            let log_sum: f64 = counts[i].iter()
                .map(|&c| (c as f64).ln()).sum();
            geo_means[i] = (log_sum / n_samples as f64).exp();
        }
    }

    let mut size_factors = vec![1.0; n_samples];
    for j in 0..n_samples {
        let mut ratios: Vec<f64> = (0..n_genes)
            .filter(|&i| geo_means[i] > 0.0 && counts[i][j] > 0)
            .map(|i| counts[i][j] as f64 / geo_means[i])
            .collect();
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !ratios.is_empty() {
            size_factors[j] = ratios[ratios.len() / 2];
        }
    }
    size_factors
}

/// TPM normalization.
pub fn compute_tpm(counts: &[Vec<u32>], lengths: &[f64]) -> Vec<Vec<f64>> {
    let n_samples = counts[0].len();
    let mut tpm = vec![vec![0.0; n_samples]; counts.len()];

    for j in 0..n_samples {
        let rpk: Vec<f64> = counts.iter().zip(lengths.iter())
            .map(|(c, &l)| c[j] as f64 / (l / 1000.0))
            .collect();
        let scale = rpk.iter().sum::<f64>() / 1e6;
        for (i, &r) in rpk.iter().enumerate() {
            tpm[i][j] = if scale > 0.0 { r / scale } else { 0.0 };
        }
    }
    tpm
}

fn main() {
    let counts = vec![
        vec![100u32, 200, 150],
        vec![500, 1000, 750],
        vec![50, 100, 75],
    ];
    let lengths = vec![1000.0, 2000.0, 500.0];

    let sf = median_of_ratios(&counts);
    println!("Size factors: {:?}", sf);

    let tpm = compute_tpm(&counts, &lengths);
    for (i, row) in tpm.iter().enumerate() {
        println!("Gene {}: TPM = {:?}", i, row.iter().map(|v| format!("{:.1}", v)).collect::<Vec<_>>());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tpm_sums_to_million() {
        let counts = vec![vec![100u32, 200], vec![300, 600]];
        let lengths = vec![1000.0, 2000.0];
        let tpm = compute_tpm(&counts, &lengths);
        let sum: f64 = tpm.iter().map(|r| r[0]).sum();
        assert!((sum - 1e6).abs() < 1.0);
    }
}
