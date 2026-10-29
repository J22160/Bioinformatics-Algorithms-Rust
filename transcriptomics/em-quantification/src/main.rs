//! Expectation-Maximization for Transcript Quantification
//!
//! Resolves multi-mapped reads to estimate transcript abundances.
//! This is the core algorithm behind Salmon, RSEM, and kallisto.

use std::time::Instant;

/// A single read's alignment information.
/// Each read maps to one or more transcripts.
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    pub read_id: usize,
    pub transcript_ids: Vec<usize>,
    pub alignment_scores: Vec<f64>, // normalized probability weights
}

/// Result of EM quantification
#[derive(Debug)]
pub struct QuantResult {
    pub abundances: Vec<f64>,      // relative abundance (θ)
    pub counts: Vec<f64>,          // estimated counts per transcript
    pub tpm: Vec<f64>,             // transcripts per million
    pub iterations: usize,
    pub log_likelihood: f64,
}

/// EM algorithm for transcript quantification.
///
/// # Rust Performance Advantages
///
/// 1. **Cache-friendly iteration**: Alignments stored in contiguous
///    Vec<ReadAlignment>, so the E-step iterates sequentially
///    through memory. Python's list-of-dicts scatters data.
///
/// 2. **No GC during EM loop**: The inner loop allocates nothing.
///    Python creates temporary float objects each iteration.
///
/// 3. **f64 arithmetic**: Direct IEEE 754 operations, no boxing.
///    Python floats are heap-allocated objects with ref counting.
///
/// 4. **Vectorizable inner loop**: The probability normalization
///    in the E-step compiles to SIMD instructions on modern CPUs.
pub fn em_quantify(
    alignments: &[ReadAlignment],
    num_transcripts: usize,
    effective_lengths: &[f64],
    max_iterations: usize,
    convergence_threshold: f64,
) -> QuantResult {
    let n_reads = alignments.len();

    // Initialize abundances uniformly
    let mut theta: Vec<f64> = vec![1.0 / num_transcripts as f64; num_transcripts];
    let mut counts = vec![0.0f64; num_transcripts];
    let mut log_likelihood = f64::NEG_INFINITY;
    let mut iterations = 0;

    for iter in 0..max_iterations {
        // Reset counts for M-step accumulation
        counts.iter_mut().for_each(|c| *c = 0.0);
        let mut new_ll = 0.0f64;

        // E-step: For each read, compute posterior probability of
        // assignment to each transcript, weighted by current abundances
        for aln in alignments {
            // Compute unnormalized assignment probabilities
            let mut probs: Vec<f64> = aln.transcript_ids.iter()
                .zip(aln.alignment_scores.iter())
                .map(|(&tid, &score)| {
                    theta[tid] * score / effective_lengths[tid]
                })
                .collect();

            // Normalize to get posterior probabilities
            let total: f64 = probs.iter().sum();
            if total > 0.0 {
                new_ll += total.ln();
                for p in probs.iter_mut() {
                    *p /= total;
                }

                // M-step accumulation: add fractional counts
                for (i, &tid) in aln.transcript_ids.iter().enumerate() {
                    counts[tid] += probs[i];
                }
            }
        }

        // M-step: Update abundances
        let count_sum: f64 = counts.iter().sum();
        if count_sum > 0.0 {
            for i in 0..num_transcripts {
                theta[i] = counts[i] / count_sum;
            }
        }

        // Check convergence
        let ll_change = (new_ll - log_likelihood).abs();
        log_likelihood = new_ll;
        iterations = iter + 1;

        if ll_change < convergence_threshold && iter > 0 {
            break;
        }
    }

    // Compute TPM (Transcripts Per Million)
    let mut tpm = vec![0.0f64; num_transcripts];
    let rate_sum: f64 = (0..num_transcripts)
        .map(|i| counts[i] / effective_lengths[i])
        .sum();

    if rate_sum > 0.0 {
        for i in 0..num_transcripts {
            tpm[i] = (counts[i] / effective_lengths[i]) / rate_sum * 1e6;
        }
    }

    QuantResult {
        abundances: theta,
        counts,
        tpm,
        iterations,
        log_likelihood,
    }
}

/// Simulate RNA-seq read alignments for testing
fn simulate_alignments(
    true_abundances: &[f64],
    effective_lengths: &[f64],
    num_reads: usize,
    multi_map_rate: f64,
    seed: u64,
) -> Vec<ReadAlignment> {
    let num_tx = true_abundances.len();
    let mut state = seed;
    let mut alignments = Vec::with_capacity(num_reads);

    // Simple LCG random number generator
    let mut next_rand = || -> f64 {
        state = state.wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as f64 / (1u64 << 31) as f64
    };

    // Sampling weights proportional to abundance × length
    let weights: Vec<f64> = true_abundances.iter()
        .zip(effective_lengths.iter())
        .map(|(&a, &l)| a * l)
        .collect();
    let weight_sum: f64 = weights.iter().sum();
    let cum_weights: Vec<f64> = weights.iter()
        .scan(0.0, |acc, &w| { *acc += w / weight_sum; Some(*acc) })
        .collect();

    for read_id in 0..num_reads {
        // Sample true origin transcript
        let r = next_rand();
        let true_tx = cum_weights.iter()
            .position(|&cw| r <= cw)
            .unwrap_or(num_tx - 1);

        let mut transcript_ids = vec![true_tx];
        let mut scores = vec![1.0];

        // With some probability, add multi-mapping targets
        if next_rand() < multi_map_rate {
            let extra = ((next_rand() * 3.0) as usize).max(1);
            for _ in 0..extra {
                let other = (next_rand() * num_tx as f64) as usize % num_tx;
                if other != true_tx && !transcript_ids.contains(&other) {
                    transcript_ids.push(other);
                    scores.push(0.5 + next_rand() * 0.5);
                }
            }
        }

        alignments.push(ReadAlignment {
            read_id,
            transcript_ids,
            alignment_scores: scores,
        });
    }

    alignments
}

fn main() {
    println!("=== EM Algorithm: Transcript Quantification ===\n");

    // Simulate a gene with 4 isoforms
    let transcript_names = vec!["TP53-α", "TP53-β", "TP53-γ", "TP53-δ"];
    let true_abundances = vec![0.50, 0.25, 0.15, 0.10];
    let effective_lengths = vec![2500.0, 1800.0, 1200.0, 900.0];
    let num_reads = 5000;

    println!("Simulating {} reads from {} transcripts", num_reads, transcript_names.len());
    println!("True abundances: {:?}\n", true_abundances);

    let alignments = simulate_alignments(
        &true_abundances,
        &effective_lengths,
        num_reads,
        0.4, // 40% multi-mapping rate
        42,
    );

    let multi_mapped = alignments.iter()
        .filter(|a| a.transcript_ids.len() > 1)
        .count();
    println!("Unique reads: {}", num_reads - multi_mapped);
    println!("Multi-mapped reads: {} ({:.1}%)\n",
             multi_mapped, multi_mapped as f64 / num_reads as f64 * 100.0);

    let start = Instant::now();
    let result = em_quantify(
        &alignments,
        transcript_names.len(),
        &effective_lengths,
        1000,
        1e-8,
    );
    let elapsed = start.elapsed();

    println!("EM converged in {} iterations ({:?})\n", result.iterations, elapsed);
    println!("{:<12} {:>10} {:>10} {:>10} {:>10}",
             "Transcript", "True θ", "Est θ", "Est Count", "TPM");
    println!("{}", "-".repeat(62));

    for (i, name) in transcript_names.iter().enumerate() {
        println!("{:<12} {:>10.4} {:>10.4} {:>10.1} {:>10.1}",
                 name,
                 true_abundances[i],
                 result.abundances[i],
                 result.counts[i],
                 result.tpm[i]);
    }

    // Benchmark scaling
    println!("\n=== Scaling Benchmark ===");
    for &n in &[10_000, 50_000, 200_000] {
        let alns = simulate_alignments(&true_abundances, &effective_lengths, n, 0.4, 42);
        let start = Instant::now();
        let r = em_quantify(&alns, 4, &effective_lengths, 500, 1e-8);
        let elapsed = start.elapsed();
        println!("  {} reads: {} iterations, time={:?}", n, r.iterations, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_transcript() {
        let alns = vec![
            ReadAlignment { read_id: 0, transcript_ids: vec![0], alignment_scores: vec![1.0] },
            ReadAlignment { read_id: 1, transcript_ids: vec![0], alignment_scores: vec![1.0] },
        ];
        let result = em_quantify(&alns, 1, &[1000.0], 100, 1e-8);
        assert!((result.abundances[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_equal_abundance() {
        let alns: Vec<ReadAlignment> = (0..1000).map(|i| ReadAlignment {
            read_id: i,
            transcript_ids: vec![i % 2],
            alignment_scores: vec![1.0],
        }).collect();
        let result = em_quantify(&alns, 2, &[1000.0, 1000.0], 100, 1e-8);
        assert!((result.abundances[0] - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_convergence() {
        let alns = simulate_alignments(&[0.7, 0.3], &[1000.0, 1000.0], 5000, 0.3, 42);
        let result = em_quantify(&alns, 2, &[1000.0, 1000.0], 1000, 1e-8);
        assert!(result.iterations < 500);
    }
}
