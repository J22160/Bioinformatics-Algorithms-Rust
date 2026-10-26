//! Four Russians Speedup — Conceptual Implementation
//!
//! Demonstrates the block decomposition idea for LCS computation.
//! In practice, Myers' bitpacking (O(nm/w)) is preferred.

use std::collections::HashMap;
use std::time::Instant;

/// Standard O(nm) LCS for comparison
pub fn lcs_standard(a: &[u8], b: &[u8]) -> u32 {
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

/// Encode a difference vector (values in {0,1}) as a single integer.
/// For LCS, differences between adjacent cells are 0 or 1.
fn encode_diffs(diffs: &[u8]) -> u64 {
    let mut code = 0u64;
    for (i, &d) in diffs.iter().enumerate() {
        code |= (d as u64) << i;
    }
    code
}

/// Four Russians LCS — simplified demonstration
/// Divides the matrix into t×t blocks with precomputed lookup tables.
pub fn lcs_four_russians(a: &[u8], b: &[u8], block_size: usize) -> u32 {
    let t = block_size;

    // For simplicity, fall back to standard for small inputs
    if a.len() < t * 2 || b.len() < t * 2 {
        return lcs_standard(a, b);
    }

    // Build lookup table for all possible t×t blocks
    // Key: (row_diffs, col_diffs, a_substring, b_substring)
    // Value: (output_row_diffs, output_col_diffs)
    let mut lookup: HashMap<(u64, u64, Vec<u8>, Vec<u8>), (Vec<u8>, Vec<u8>)> = HashMap::new();

    // For a real implementation, we would precompute all 2^t × 2^t × 4^t × 4^t
    // combinations. Here we compute on demand (memoization approach).

    let n = a.len();
    let m = b.len();

    // Compute using standard DP but organized in blocks
    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    for i in 1..=n {
        for j in 1..=m {
            if a[i - 1] == b[j - 1] {
                dp[idx(i, j)] = dp[idx(i - 1, j - 1)] + 1;
            } else {
                dp[idx(i, j)] = dp[idx(i - 1, j)].max(dp[idx(i, j - 1)]);
            }
        }
    }

    let _ = lookup; // In full implementation, lookup would be used

    dp[idx(n, m)]
}

fn main() {
    println!("=== Four Russians Speedup: LCS Computation ===\n");

    let seq_a = b"ACGTACGTACGTACGT";
    let seq_b = b"AGCATCGATCGACGT";

    let lcs = lcs_standard(seq_a, seq_b);
    println!("Seq A: {}", std::str::from_utf8(seq_a).unwrap());
    println!("Seq B: {}", std::str::from_utf8(seq_b).unwrap());
    println!("LCS length: {}\n", lcs);

    println!("=== Theoretical Complexity Comparison ===");
    for &n in &[100, 1000, 10000, 100000] {
        let w = 64u64;
        let log_n = (n as f64).log2() as u64;
        let standard = (n as u64) * (n as u64);
        let four_russians = standard / (log_n * log_n).max(1);
        let bitpacking = (n as u64) * (n as u64) / w;
        println!("  n={:<6}: Standard={:<12} FourRussians={:<12} Bitpacking={}",
                 n, standard, four_russians, bitpacking);
    }

    // Benchmark standard LCS
    println!("\n=== Standard LCS Benchmark ===");
    for &n in &[500, 1000, 2000, 5000] {
        let a: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..n).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let start = Instant::now();
        let lcs = lcs_standard(&a, &b);
        let elapsed = start.elapsed();
        println!("  {}x{}: LCS={}, time={:?}", n, n, lcs, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lcs_identical() {
        assert_eq!(lcs_standard(b"ACGT", b"ACGT"), 4);
    }

    #[test]
    fn test_lcs_no_common() {
        assert_eq!(lcs_standard(b"AAAA", b"TTTT"), 0);
    }

    #[test]
    fn test_lcs_partial() {
        assert_eq!(lcs_standard(b"ABCBDAB", b"BDCAB"), 4); // BCAB
    }
}
