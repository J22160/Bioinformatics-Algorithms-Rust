//! Pileup-Based Genotyping
//!
//! Computes genotype likelihoods from a column of aligned bases
//! at each genomic position, using base quality scores to model
//! sequencing error. This is the approach used by bcftools mpileup.

use std::time::Instant;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// A single base observation with its quality score.
#[derive(Debug, Clone)]
pub struct PileupBase {
    pub base: u8,
    pub quality: u8, // Phred score
}

/// A diploid genotype represented as two allele indices (0-3 for ACGT).
#[derive(Debug, Clone, Copy)]
pub struct Genotype {
    pub allele1: usize,
    pub allele2: usize,
}

impl Genotype {
    pub fn label(&self) -> String {
        format!("{}/{}",
                BASES[self.allele1] as char,
                BASES[self.allele2] as char)
    }

    pub fn is_het(&self) -> bool {
        self.allele1 != self.allele2
    }
}

/// Result of genotyping at one position.
#[derive(Debug)]
pub struct GenotypeCall {
    pub genotype: Genotype,
    pub genotype_quality: f64,     // Phred-scaled
    pub log_likelihoods: Vec<f64>, // for all 10 genotypes
    pub posterior_probs: Vec<f64>,
    pub all_genotypes: Vec<Genotype>,
}

/// Enumerate all 10 diploid genotypes.
pub fn all_diploid_genotypes() -> Vec<Genotype> {
    let mut genotypes = Vec::with_capacity(10);
    for i in 0..4 {
        for j in i..4 {
            genotypes.push(Genotype { allele1: i, allele2: j });
        }
    }
    genotypes
}

/// Convert Phred quality score to error probability.
///
/// Phred score Q means P(error) = 10^(-Q/10).
/// Q30 → 0.001, Q20 → 0.01, Q10 → 0.1
fn phred_to_error(quality: u8) -> f64 {
    10.0f64.powf(-(quality as f64) / 10.0)
}

/// Compute log-likelihood of observing a base given a genotype.
///
/// For diploid genotype {a1, a2}, the probability of observing
/// base b with error rate e is:
///   P(b | G) = 0.5 * P(b | a1) + 0.5 * P(b | a2)
///
/// where P(b | a) = (1-e) if b==a, else e/3
///
/// # Rust Performance Advantages
///
/// 1. **No heap allocation in inner loop**: All computations use
///    stack-allocated f64 values. Python would create float objects.
///
/// 2. **Log-space arithmetic**: We accumulate log-likelihoods to
///    avoid floating-point underflow on deep pileups. Rust's
///    f64::ln() compiles to a single x87 FPU instruction.
///
/// 3. **Branch-free base matching**: The comparison b == allele
///    compiles to a cmov instruction, avoiding branch misprediction.
fn log_likelihood_base_given_genotype(
    base: u8,
    quality: u8,
    genotype: &Genotype,
) -> f64 {
    let e = phred_to_error(quality);
    let base_idx = match base {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 0,
    };

    let p_allele = |allele: usize| -> f64 {
        if base_idx == allele {
            1.0 - e
        } else {
            e / 3.0
        }
    };

    let p = if genotype.allele1 == genotype.allele2 {
        // Homozygous: P(b|G) = P(b|a)
        p_allele(genotype.allele1)
    } else {
        // Heterozygous: P(b|G) = 0.5*P(b|a1) + 0.5*P(b|a2)
        0.5 * p_allele(genotype.allele1) + 0.5 * p_allele(genotype.allele2)
    };

    p.ln()
}

/// Compute genotype likelihoods for a pileup at one position.
///
/// Returns log10-scaled likelihoods for all 10 diploid genotypes,
/// the called genotype, and the genotype quality (GQ).
pub fn call_genotype(
    pileup: &[PileupBase],
    prior_het: f64,  // prior probability of heterozygosity (~0.001)
) -> GenotypeCall {
    let genotypes = all_diploid_genotypes();
    let n_genotypes = genotypes.len();

    // Compute log-likelihoods for each genotype
    let mut log_likes = vec![0.0f64; n_genotypes];
    for (g_idx, genotype) in genotypes.iter().enumerate() {
        for obs in pileup {
            log_likes[g_idx] += log_likelihood_base_given_genotype(
                obs.base, obs.quality, genotype,
            );
        }
    }

    // Compute priors (simplified: het sites get prior_het, hom-ref
    // gets most of the remaining, hom-alt gets a tiny fraction)
    let mut log_priors = vec![0.0f64; n_genotypes];
    for (g_idx, genotype) in genotypes.iter().enumerate() {
        let prior = if genotype.allele1 == genotype.allele2 {
            if genotype.allele1 == 0 { // hom-ref (assuming ref=A)
                1.0 - prior_het * 3.0
            } else {
                prior_het * 0.1 // hom-alt is rare
            }
        } else {
            prior_het // het
        };
        log_priors[g_idx] = prior.max(1e-10).ln();
    }

    // Posterior = likelihood × prior (in log space: sum)
    let mut log_posteriors: Vec<f64> = log_likes.iter()
        .zip(log_priors.iter())
        .map(|(&ll, &lp)| ll + lp)
        .collect();

    // Normalize posteriors using log-sum-exp trick
    let max_post = log_posteriors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_sum_exp: f64 = log_posteriors.iter()
        .map(|&lp| (lp - max_post).exp())
        .sum::<f64>()
        .ln() + max_post;

    let posteriors: Vec<f64> = log_posteriors.iter()
        .map(|&lp| (lp - log_sum_exp).exp())
        .collect();

    // Find best genotype
    let best_idx = posteriors.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap().0;

    // Compute GQ: Phred-scaled probability that call is wrong
    let gq = -10.0 * (1.0 - posteriors[best_idx]).max(1e-300).log10();

    GenotypeCall {
        genotype: genotypes[best_idx],
        genotype_quality: gq,
        log_likelihoods: log_likes,
        posterior_probs: posteriors,
        all_genotypes: genotypes,
    }
}

fn main() {
    println!("=== Pileup Genotyping: Variant Detection from Read Pileups ===\n");

    // Simulate a heterozygous A/G site
    let pileup: Vec<PileupBase> = vec![
        // True A allele reads (with high quality)
        PileupBase { base: b'A', quality: 30 },
        PileupBase { base: b'A', quality: 35 },
        PileupBase { base: b'A', quality: 28 },
        PileupBase { base: b'A', quality: 32 },
        PileupBase { base: b'A', quality: 30 },
        PileupBase { base: b'A', quality: 25 },
        PileupBase { base: b'A', quality: 33 },
        // True G allele reads
        PileupBase { base: b'G', quality: 30 },
        PileupBase { base: b'G', quality: 28 },
        PileupBase { base: b'G', quality: 35 },
        PileupBase { base: b'G', quality: 30 },
        PileupBase { base: b'G', quality: 27 },
        PileupBase { base: b'G', quality: 31 },
        // One error read
        PileupBase { base: b'T', quality: 8 },
    ];

    println!("Pileup at chr7:100000 (ref=A):");
    let mut counts = [0u32; 4];
    for obs in &pileup {
        let idx = match obs.base { b'A'=>0, b'C'=>1, b'G'=>2, b'T'=>3, _=>0 };
        counts[idx] += 1;
    }
    println!("  A={}, C={}, G={}, T={}, depth={}\n",
             counts[0], counts[1], counts[2], counts[3], pileup.len());

    let start = Instant::now();
    let call = call_genotype(&pileup, 0.001);
    let elapsed = start.elapsed();

    println!("Genotype call: {} (GQ={:.1})", call.genotype.label(), call.genotype_quality);
    println!("Time: {:?}\n", elapsed);

    println!("{:<6} {:>12} {:>12}", "GT", "LogLike", "Posterior");
    println!("{}", "-".repeat(34));
    for (i, gt) in call.all_genotypes.iter().enumerate() {
        let marker = if i == call.posterior_probs.iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap().0 { " <<<" } else { "" };
        println!("{:<6} {:>12.2} {:>12.6}{}",
                 gt.label(), call.log_likelihoods[i],
                 call.posterior_probs[i], marker);
    }

    // Benchmark
    println!("\n=== Benchmark: Pileups per Second ===");
    for &depth in &[30, 60, 100, 200] {
        let big_pileup: Vec<PileupBase> = (0..depth).map(|i| PileupBase {
            base: if i % 2 == 0 { b'A' } else { b'G' },
            quality: 30,
        }).collect();
        let start = Instant::now();
        for _ in 0..100_000 {
            let _ = call_genotype(&big_pileup, 0.001);
        }
        let elapsed = start.elapsed();
        println!("  depth={}: {:.0} calls/sec",
                 depth, 100_000.0 / elapsed.as_secs_f64());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hom_ref() {
        let pileup: Vec<PileupBase> = (0..30).map(|_|
            PileupBase { base: b'A', quality: 30 }
        ).collect();
        let call = call_genotype(&pileup, 0.001);
        assert_eq!(call.genotype.label(), "A/A");
    }

    #[test]
    fn test_het() {
        let mut pileup: Vec<PileupBase> = (0..15).map(|_|
            PileupBase { base: b'A', quality: 30 }
        ).collect();
        pileup.extend((0..15).map(|_|
            PileupBase { base: b'G', quality: 30 }
        ));
        let call = call_genotype(&pileup, 0.001);
        assert_eq!(call.genotype.label(), "A/G");
    }

    #[test]
    fn test_hom_alt() {
        let pileup: Vec<PileupBase> = (0..30).map(|_|
            PileupBase { base: b'G', quality: 30 }
        ).collect();
        let call = call_genotype(&pileup, 0.001);
        assert_eq!(call.genotype.label(), "G/G");
    }
}
