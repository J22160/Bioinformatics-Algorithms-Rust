//! Needleman-Wunsch Global Alignment
//!
//! Aligns two sequences end-to-end, finding the minimum-cost edit path.
//! This is the foundational algorithm behind all sequence alignment in bioinformatics.

use std::time::Instant;

/// Scoring parameters for alignment
#[derive(Debug, Clone, Copy)]
pub struct ScoringParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_penalty: i32,
}

impl Default for ScoringParams {
    fn default() -> Self {
        ScoringParams {
            match_score: 1,
            mismatch_penalty: -1,
            gap_penalty: -2,
        }
    }
}

/// Result of a global alignment
#[derive(Debug)]
pub struct AlignmentResult {
    pub score: i32,
    pub aligned_a: String,
    pub aligned_b: String,
    pub match_line: String,
    pub identity: f64,
    pub gaps: usize,
}

/// Needleman-Wunsch global alignment with full traceback.
///
/// # Rust Performance Advantages
///
/// 1. **Flat 1D Vec for the matrix**: We use `vec![0i32; (n+1)*(m+1)]` with
///    `idx(i,j) = i*(m+1)+j`. This gives contiguous memory layout.
///    Python's list-of-lists scatters rows across the heap.
///
/// 2. **i32 cells (4 bytes)**: Python int objects are 28+ bytes each.
///    For a 1000×1000 matrix that's 4 MB vs 28+ MB.
///
/// 3. **No bounds checking in release**: The compiler proves safety
///    at compile time, then emits unchecked indexing.
///
/// 4. **Cache-friendly traversal**: Row-major filling matches memory layout,
///    so the CPU prefetcher efficiently loads next cache lines.
pub fn needleman_wunsch(a: &[u8], b: &[u8], params: &ScoringParams) -> AlignmentResult {
    let n = a.len();
    let m = b.len();
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    // Allocate DP matrix as flat 1D array
    let mut dp = vec![0i32; (n + 1) * (m + 1)];

    // Initialize first column (gaps in b)
    for i in 1..=n {
        dp[idx(i, 0)] = i as i32 * params.gap_penalty;
    }
    // Initialize first row (gaps in a)
    for j in 1..=m {
        dp[idx(0, j)] = j as i32 * params.gap_penalty;
    }

    // Fill the matrix
    for i in 1..=n {
        for j in 1..=m {
            let match_score = if a[i - 1] == b[j - 1] {
                params.match_score
            } else {
                params.mismatch_penalty
            };

            let diag = dp[idx(i - 1, j - 1)] + match_score;
            let up = dp[idx(i - 1, j)] + params.gap_penalty;
            let left = dp[idx(i, j - 1)] + params.gap_penalty;

            dp[idx(i, j)] = diag.max(up).max(left);
        }
    }

    // Traceback
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let match_score = if a[i - 1] == b[j - 1] {
                params.match_score
            } else {
                params.mismatch_penalty
            };
            if dp[idx(i, j)] == dp[idx(i - 1, j - 1)] + match_score {
                aligned_a.push(a[i - 1]);
                aligned_b.push(b[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[idx(i, j)] == dp[idx(i - 1, j)] + params.gap_penalty {
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

    // Build match line and compute statistics
    let mut matches = 0usize;
    let mut gaps = 0usize;
    let match_line: String = aligned_a
        .iter()
        .zip(aligned_b.iter())
        .map(|(&ca, &cb)| {
            if ca == b'-' || cb == b'-' {
                gaps += 1;
                ' '
            } else if ca == cb {
                matches += 1;
                '|'
            } else {
                '.'
            }
        })
        .collect();

    let aln_len = aligned_a.len();
    let identity = if aln_len > 0 {
        matches as f64 / aln_len as f64 * 100.0
    } else {
        0.0
    };

    AlignmentResult {
        score: dp[idx(n, m)],
        aligned_a: String::from_utf8(aligned_a).unwrap(),
        aligned_b: String::from_utf8(aligned_b).unwrap(),
        match_line,
        identity,
        gaps,
    }
}

/// Score-only version using O(m) space (two rows).
/// Useful when you only need the alignment score, not the traceback.
pub fn nw_score_only(a: &[u8], b: &[u8], params: &ScoringParams) -> i32 {
    let m = b.len();
    let mut prev = vec![0i32; m + 1];
    let mut curr = vec![0i32; m + 1];

    for j in 0..=m {
        prev[j] = j as i32 * params.gap_penalty;
    }

    for i in 1..=a.len() {
        curr[0] = i as i32 * params.gap_penalty;
        for j in 1..=m {
            let match_score = if a[i - 1] == b[j - 1] {
                params.match_score
            } else {
                params.mismatch_penalty
            };
            curr[j] = (prev[j - 1] + match_score)
                .max(prev[j] + params.gap_penalty)
                .max(curr[j - 1] + params.gap_penalty);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[m]
}

fn main() {
    println!("=== Needleman-Wunsch: Hemoglobin Beta Chain Comparison ===\n");

    // Human vs Gorilla HBB gene fragment (protein sequence)
    let human  = b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDL";
    let gorilla = b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFDSFGDL";

    let params = ScoringParams {
        match_score: 2,
        mismatch_penalty: -1,
        gap_penalty: -3,
    };

    let start = Instant::now();
    let result = needleman_wunsch(human, gorilla, &params);
    let elapsed = start.elapsed();

    println!("Human HBB:   {}", std::str::from_utf8(human).unwrap());
    println!("Gorilla HBB: {}\n", std::str::from_utf8(gorilla).unwrap());
    println!("Alignment:");
    println!("  {}", result.aligned_a);
    println!("  {}", result.match_line);
    println!("  {}", result.aligned_b);
    println!("\nScore: {}", result.score);
    println!("Identity: {:.1}%", result.identity);
    println!("Gaps: {}", result.gaps);
    println!("Time: {:?}", elapsed);

    // DNA alignment example
    println!("\n=== DNA Alignment: Promoter Region Comparison ===\n");
    let seq_a = b"ATGCTAGCTAGCATCGATCGATCG";
    let seq_b = b"ATGCAAGCTAGATCGATCGATCG";

    let dna_params = ScoringParams::default();
    let dna_result = needleman_wunsch(seq_a, seq_b, &dna_params);

    println!("Seq A: {}", std::str::from_utf8(seq_a).unwrap());
    println!("Seq B: {}\n", std::str::from_utf8(seq_b).unwrap());
    println!("  {}", dna_result.aligned_a);
    println!("  {}", dna_result.match_line);
    println!("  {}", dna_result.aligned_b);
    println!("\nScore: {} | Identity: {:.1}%", dna_result.score, dna_result.identity);

    // Benchmark
    println!("\n=== Performance Benchmark ===");
    for &n in &[100, 500, 1000, 2000] {
        let a: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..n).map(|i| b"ACGT"[(i + 1) % 4]).collect();
        let start = Instant::now();
        let score = nw_score_only(&a, &b, &dna_params);
        let elapsed = start.elapsed();
        println!("  {}x{}: score={}, time={:?}", n, n, score, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_sequences() {
        let params = ScoringParams::default();
        let result = needleman_wunsch(b"ACGT", b"ACGT", &params);
        assert_eq!(result.score, 4); // 4 matches × +1
        assert_eq!(result.gaps, 0);
    }

    #[test]
    fn test_one_substitution() {
        let params = ScoringParams::default();
        let result = needleman_wunsch(b"ACGT", b"ACTT", &params);
        assert_eq!(result.score, 2); // 3 matches - 1 mismatch
    }

    #[test]
    fn test_gap_insertion() {
        let params = ScoringParams::default();
        let result = needleman_wunsch(b"AGT", b"ACGT", &params);
        assert!(result.gaps >= 1);
    }

    #[test]
    fn test_empty_sequences() {
        let params = ScoringParams::default();
        let result = needleman_wunsch(b"", b"ACGT", &params);
        assert_eq!(result.score, -8); // 4 gaps × -2
    }

    #[test]
    fn test_score_only_matches_full() {
        let params = ScoringParams::default();
        let a = b"AGTACGCA";
        let b_seq = b"TATGC";
        let full = needleman_wunsch(a, b_seq, &params);
        let score = nw_score_only(a, b_seq, &params);
        assert_eq!(full.score, score);
    }
}
