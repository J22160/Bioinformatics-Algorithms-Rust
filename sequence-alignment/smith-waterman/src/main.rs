//! Smith-Waterman Local Alignment
//!
//! Finds the highest-scoring local alignment between two sequences.
//! Used to find conserved domains, motifs, or genes within longer sequences.

use std::time::Instant;

#[derive(Debug)]
pub struct LocalAlignmentResult {
    pub score: i32,
    pub aligned_a: String,
    pub aligned_b: String,
    pub start_a: usize,
    pub start_b: usize,
    pub end_a: usize,
    pub end_b: usize,
}

/// Smith-Waterman local alignment.
///
/// Key difference from Needleman-Wunsch: cells never go below 0,
/// and the best alignment can start/end anywhere in the matrix.
///
/// # Rust Performance Advantages
///
/// 1. **Tracking max during fill**: We track the maximum score and its
///    position while filling the matrix. No second pass needed.
///    Rust compiles this to branch-free cmov instructions.
///
/// 2. **Flat array layout**: Same cache benefits as NW implementation.
///    The sequential access pattern during fill is perfectly predictable
///    by the CPU hardware prefetcher.
///
/// 3. **Zero-copy traceback**: We read directly from the DP matrix
///    without copying it. Rust's borrow checker guarantees this is safe.
pub fn smith_waterman(
    a: &[u8],
    b: &[u8],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> LocalAlignmentResult {
    let n = a.len();
    let m = b.len();
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    let mut dp = vec![0i32; (n + 1) * (m + 1)];
    // No initialization needed — cells default to 0

    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    for i in 1..=n {
        for j in 1..=m {
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };

            let val = 0i32
                .max(dp[idx(i - 1, j - 1)] + score)
                .max(dp[idx(i - 1, j)] + gap_penalty)
                .max(dp[idx(i, j - 1)] + gap_penalty);

            dp[idx(i, j)] = val;

            if val > max_score {
                max_score = val;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback from max_score position until we hit 0
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let (end_a, end_b) = (max_i, max_j);
    let mut i = max_i;
    let mut j = max_j;

    while i > 0 && j > 0 && dp[idx(i, j)] > 0 {
        let score = if a[i - 1] == b[j - 1] {
            match_score
        } else {
            mismatch_penalty
        };

        if dp[idx(i, j)] == dp[idx(i - 1, j - 1)] + score {
            aligned_a.push(a[i - 1]);
            aligned_b.push(b[j - 1]);
            i -= 1;
            j -= 1;
        } else if dp[idx(i, j)] == dp[idx(i - 1, j)] + gap_penalty {
            aligned_a.push(a[i - 1]);
            aligned_b.push(b'-');
            i -= 1;
        } else {
            aligned_a.push(b'-');
            aligned_b.push(b[j - 1]);
            j -= 1;
        }
    }

    aligned_a.reverse();
    aligned_b.reverse();

    LocalAlignmentResult {
        score: max_score,
        aligned_a: String::from_utf8(aligned_a).unwrap(),
        aligned_b: String::from_utf8(aligned_b).unwrap(),
        start_a: i,
        start_b: j,
        end_a,
        end_b,
    }
}

fn main() {
    println!("=== Smith-Waterman: Finding a Resistance Gene in a Contig ===\n");

    // Simulated contig containing an embedded resistance gene
    let contig   = b"TTGCAATGCAATCGATCGAATGGCATGCTAGCTAGCATCGATATCGCAATGCAATGCAAGC";
    let res_gene = b"GCATGCTAGCTAGCATCGA";

    let result = smith_waterman(contig, res_gene, 2, -1, -2);

    println!("Contig ({}bp): {}...{}",
             contig.len(),
             std::str::from_utf8(&contig[..20]).unwrap(),
             std::str::from_utf8(&contig[contig.len()-10..]).unwrap());
    println!("Gene   ({}bp): {}\n", res_gene.len(), std::str::from_utf8(res_gene).unwrap());
    println!("Local alignment (score = {}):", result.score);
    println!("  Position in contig: {}-{}", result.start_a, result.end_a);
    println!("  Position in gene:   {}-{}", result.start_b, result.end_b);

    let match_line: String = result.aligned_a.chars()
        .zip(result.aligned_b.chars())
        .map(|(a, b)| if a == b { '|' } else if a == '-' || b == '-' { ' ' } else { '.' })
        .collect();
    println!("  {}", result.aligned_a);
    println!("  {}", match_line);
    println!("  {}", result.aligned_b);

    // Benchmark
    println!("\n=== Benchmark ===");
    for &n in &[500, 1000, 2000] {
        let a: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..n).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let start = Instant::now();
        let r = smith_waterman(&a, &b, 2, -1, -2);
        let elapsed = start.elapsed();
        println!("  {}x{}: score={}, time={:?}", n, n, r.score, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_substring() {
        let result = smith_waterman(b"XXXACGTXXX", b"ACGT", 1, -1, -2);
        assert_eq!(result.score, 4);
    }

    #[test]
    fn test_no_similarity() {
        let result = smith_waterman(b"AAAA", b"TTTT", 1, -1, -2);
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_partial_match() {
        let result = smith_waterman(b"AACGT", b"CGTAA", 2, -1, -2);
        assert!(result.score > 0);
    }
}
