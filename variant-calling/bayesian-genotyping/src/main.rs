//! Bayesian Genotyping — Full Posterior Genotype Inference
//!
//! Computes posterior genotype probabilities using Bayes' theorem
//! with population-informed priors and a substitution error matrix.

use std::time::Instant;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Substitution error matrix: P(observed | truth)
/// Rows = truth base (A,C,G,T), Cols = observed base (A,C,G,T)
/// Diagonal = correct base probability, off-diagonal = error
fn default_substitution_matrix() -> [[f64; 4]; 4] {
    [
        [0.995, 0.002, 0.001, 0.002], // truth=A
        [0.001, 0.995, 0.002, 0.002], // truth=C
        [0.002, 0.002, 0.995, 0.001], // truth=G
        [0.002, 0.001, 0.002, 0.995], // truth=T
    ]
}

/// Diploid genotype as allele frequency mixture
fn genotype_to_mixture(a1: usize, a2: usize) -> [f64; 4] {
    let mut f = [0.0; 4];
    if a1 == a2 {
        f[a1] = 1.0;
    } else {
        f[a1] = 0.5;
        f[a2] = 0.5;
    }
    f
}

/// Compute log-likelihood of observed base counts given a mixture
/// using the substitution error matrix approach.
pub fn log_likelihood_mixture(
    counts: &[u32; 4],
    mixture: &[f64; 4],
    sub_matrix: &[[f64; 4]; 4],
) -> f64 {
    let mut ll = 0.0f64;
    for obs_base in 0..4 {
        if counts[obs_base] == 0 { continue; }
        // P(obs | mixture) = Σ_b f_b × P(obs | b)
        let p_obs: f64 = (0..4)
            .map(|truth| mixture[truth] * sub_matrix[truth][obs_base])
            .sum();
        ll += counts[obs_base] as f64 * p_obs.max(1e-300).ln();
    }
    ll
}

/// MLE of allele frequency given binary model
pub fn allele_frequency_mle(k: u32, n: u32, e: f64) -> f64 {
    let ratio = k as f64 / n as f64;
    (ratio - e) / (1.0 - 2.0 * e)
}

/// Full Bayesian genotyping using substitution matrix
pub fn bayesian_genotype(
    counts: &[u32; 4],
    ref_base: usize,
    prior_het: f64,
    sub_matrix: &[[f64; 4]; 4],
) -> (usize, usize, f64, Vec<f64>) {
    let mut genotypes = Vec::new();
    let mut log_posteriors = Vec::new();

    for a1 in 0..4usize {
        for a2 in a1..4 {
            genotypes.push((a1, a2));
            let mixture = genotype_to_mixture(a1, a2);
            let ll = log_likelihood_mixture(counts, &mixture, sub_matrix);

            // Prior
            let prior = if a1 == a2 {
                if a1 == ref_base { 1.0 - 3.0 * prior_het }
                else { prior_het * 0.05 }
            } else {
                prior_het
            };
            log_posteriors.push(ll + prior.max(1e-15).ln());
        }
    }

    // Normalize
    let max_lp = log_posteriors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_sum: f64 = log_posteriors.iter()
        .map(|&lp| (lp - max_lp).exp()).sum::<f64>().ln() + max_lp;
    let posteriors: Vec<f64> = log_posteriors.iter()
        .map(|&lp| (lp - log_sum).exp()).collect();

    let best = posteriors.iter().enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap().0;
    let gq = -10.0 * (1.0 - posteriors[best]).max(1e-300).log10();

    (genotypes[best].0, genotypes[best].1, gq, posteriors)
}

fn main() {
    println!("=== Bayesian Genotyping with Substitution Matrix ===\n");

    let sub_matrix = default_substitution_matrix();
    let counts: [u32; 4] = [15, 0, 12, 1]; // A=15, C=0, G=12, T=1

    println!("Observed counts: A={}, C={}, G={}, T={}", counts[0], counts[1], counts[2], counts[3]);

    let (a1, a2, gq, posts) = bayesian_genotype(&counts, 0, 0.001, &sub_matrix);
    println!("Called genotype: {}/{} (GQ={:.1})\n",
             BASES[a1] as char, BASES[a2] as char, gq);

    // MLE of allele frequency
    let total = counts.iter().sum::<u32>();
    let p_mle = allele_frequency_mle(counts[2], total, 0.005);
    println!("MLE allele frequency (G): {:.3}", p_mle);

    println!("\nPosterior probabilities:");
    let mut idx = 0;
    for a1 in 0..4 {
        for a2 in a1..4 {
            if posts[idx] > 0.001 {
                println!("  {}/{}: {:.4}", BASES[a1] as char, BASES[a2] as char, posts[idx]);
            }
            idx += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hom_ref_call() {
        let sub = default_substitution_matrix();
        let counts = [30, 0, 1, 0];
        let (a1, a2, _, _) = bayesian_genotype(&counts, 0, 0.001, &sub);
        assert_eq!((a1, a2), (0, 0)); // A/A
    }

    #[test]
    fn test_mle_no_error() {
        let p = allele_frequency_mle(15, 30, 0.0);
        assert!((p - 0.5).abs() < 0.01);
    }
}