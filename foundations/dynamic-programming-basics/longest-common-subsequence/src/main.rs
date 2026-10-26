//! Longest Common Subsequence — Conserved Gene Order Analysis

use std::time::Instant;

/// Compute the LCS length and the actual subsequence via traceback.
pub fn lcs<T: PartialEq + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let n = a.len();
    let m = b.len();
    let idx = |i: usize, j: usize| i * (m + 1) + j;

    // Build DP table
    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    for i in 1..=n {
        for j in 1..=m {
            if a[i - 1] == b[j - 1] {
                dp[idx(i, j)] = dp[idx(i - 1, j - 1)] + 1;
            } else {
                dp[idx(i, j)] = dp[idx(i - 1, j)].max(dp[idx(i, j - 1)]);
            }
        }
    }

    // Traceback
    let mut result = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        if a[i - 1] == b[j - 1] {
            result.push(a[i - 1].clone());
            i -= 1;
            j -= 1;
        } else if dp[idx(i - 1, j)] >= dp[idx(i, j - 1)] {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    result.reverse();
    result
}

/// Space-efficient LCS length only (two rows).
pub fn lcs_length(a: &[u8], b: &[u8]) -> u32 {
    let m = b.len();
    let mut prev = vec![0u32; m + 1];
    let mut curr = vec![0u32; m + 1];

    for i in 1..=a.len() {
        for j in 1..=m {
            if a[i - 1] == b[j - 1] {
                curr[j] = prev[j - 1] + 1;
            } else {
                curr[j] = prev[j].max(curr[j - 1]);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        curr.iter_mut().for_each(|x| *x = 0);
    }
    prev[m]
}

fn main() {
    println!("=== LCS: Conserved Gene Order Analysis ===\n");

    let genome_a = vec!["dnaA", "gyrB", "rpoB", "recA", "groEL"];
    let genome_b = vec!["dnaA", "recA", "gyrB", "groEL", "rpoB"];

    println!("Genome A gene order: {:?}", genome_a);
    println!("Genome B gene order: {:?}", genome_b);

    let conserved = lcs(&genome_a, &genome_b);
    println!("LCS (conserved order): {:?} (length {})\n", conserved, conserved.len());

    // DNA sequence example
    let seq_a = b"ACGTACGTACGT";
    let seq_b = b"AGCATCGATCG";

    let result = lcs(seq_a.as_ref(), seq_b.as_ref());
    let result_str: String = result.iter().map(|&b| b as char).collect();
    println!("Sequence A: {}", std::str::from_utf8(seq_a).unwrap());
    println!("Sequence B: {}", std::str::from_utf8(seq_b).unwrap());
    println!("LCS: {} (length {})\n", result_str, result.len());

    // Benchmark
    println!("=== Benchmark ===");
    for &n in &[100, 500, 1000, 2000] {
        let a: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..n).map(|i| b"ACGT"[(i + 1) % 4]).collect();
        let start = Instant::now();
        let len = lcs_length(&a, &b);
        let elapsed = start.elapsed();
        println!("  {}x{}: LCS length={}, time={:?}", n, n, len, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        assert_eq!(lcs(b"ACGT".as_ref(), b"ACGT".as_ref()), b"ACGT".to_vec());
    }

    #[test]
    fn test_no_common() {
        assert_eq!(lcs(b"AAA".as_ref(), b"TTT".as_ref()), Vec::<u8>::new());
    }

    #[test]
    fn test_empty() {
        assert_eq!(lcs(b"".as_ref(), b"ACGT".as_ref()), Vec::<u8>::new());
    }

    #[test]
    fn test_gene_order() {
        let a = vec!["dnaA", "gyrB", "rpoB", "recA", "groEL"];
        let b = vec!["dnaA", "recA", "gyrB", "groEL", "rpoB"];
        let result = lcs(&a, &b);
        assert_eq!(result.len(), 3);
    }
}
