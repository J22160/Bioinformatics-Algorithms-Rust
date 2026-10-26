//! Edit Distance (Levenshtein Distance) — CRISPR Off-Target Analysis
//!
//! Computes the minimum number of insertions, deletions, and substitutions
//! needed to transform one DNA sequence into another.
//!
//! This is the foundational algorithm behind sequence alignment, which is
//! described in Chapter 6 of Jones & Pevzner's bioinformatics textbook.

use std::time::Instant;

/// Compute edit distance using full O(n*m) dynamic programming matrix.
///
/// # Rust Performance Advantages
///
/// 1. **Flat 1D vector for the matrix**: We use a single `Vec<u32>` indexed
///    by `i * (m+1) + j` instead of `Vec<Vec<u32>>`. This gives perfect
///    memory locality — the CPU prefetcher can predict access patterns.
///    Python's list-of-lists scatters data across the heap.
///
/// 2. **u32 instead of arbitrary-precision int**: Each cell is exactly 4 bytes.
///    Python integers are 28+ byte objects on the heap.
///
/// 3. **No bounds checking in release mode**: The compiler proves array
///    accesses are safe, then removes the checks.
///
/// 4. **Branch prediction friendly**: The min-of-three pattern compiles to
///    conditional moves (cmov) which avoid branch mispredictions.
pub fn edit_distance(a: &[u8], b: &[u8]) -> u32 {
    let n = a.len();
    let m = b.len();

    // Flat 2D matrix stored as 1D vector for cache efficiency
    let mut dp = vec![0u32; (n + 1) * (m + 1)];

    // Helper to index into flat array
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    // Initialize first row and column
    for i in 0..=n {
        dp[idx(i, 0)] = i as u32;
    }
    for j in 0..=m {
        dp[idx(0, j)] = j as u32;
    }

    // Fill the matrix
    for i in 1..=n {
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0u32 } else { 1u32 };

            let delete = dp[idx(i - 1, j)] + 1;
            let insert = dp[idx(i, j - 1)] + 1;
            let substitute = dp[idx(i - 1, j - 1)] + cost;

            dp[idx(i, j)] = delete.min(insert).min(substitute);
        }
    }

    dp[idx(n, m)]
}

/// Space-optimized version using only two rows: O(min(n,m)) space.
/// This is critical when comparing long genomic sequences where
/// storing the full matrix would require gigabytes of RAM.
pub fn edit_distance_linear_space(a: &[u8], b: &[u8]) -> u32 {
    let n = a.len();
    let m = b.len();

    // Ensure we iterate over the shorter dimension for space efficiency
    let (a, b, n, m) = if n < m { (b, a, m, n) } else { (a, b, n, m) };

    let mut prev = (0..=m as u32).collect::<Vec<_>>();
    let mut curr = vec![0u32; m + 1];

    for i in 1..=n {
        curr[0] = i as u32;
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j] + 1)
                .min(curr[j - 1] + 1)
                .min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[m]
}

/// Reconstruct the alignment (traceback) from the DP matrix.
/// Returns a pair of aligned sequences with '-' for gaps.
pub fn align(a: &[u8], b: &[u8]) -> (String, String) {
    let n = a.len();
    let m = b.len();
    let mut dp = vec![0u32; (n + 1) * (m + 1)];
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    for i in 0..=n { dp[idx(i, 0)] = i as u32; }
    for j in 0..=m { dp[idx(0, j)] = j as u32; }

    for i in 1..=n {
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            dp[idx(i, j)] = (dp[idx(i - 1, j)] + 1)
                .min(dp[idx(i, j - 1)] + 1)
                .min(dp[idx(i - 1, j - 1)] + cost);
        }
    }

    // Traceback
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            if dp[idx(i, j)] == dp[idx(i - 1, j - 1)] + cost {
                aligned_a.push(a[i - 1]);
                aligned_b.push(b[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[idx(i, j)] == dp[idx(i - 1, j)] + 1 {
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

    (
        String::from_utf8(aligned_a).unwrap(),
        String::from_utf8(aligned_b).unwrap(),
    )
}

fn main() {
    println!("=== Edit Distance: CRISPR Off-Target Analysis ===\n");

    // CRISPR guide RNA vs candidate off-target site
    let guide    = b"GAGTCCGAGCAGAAGAAGAA";
    let offtarget = b"GAGTCCAAGCAGAAGATGAA";

    let dist = edit_distance(guide, offtarget);
    let dist_linear = edit_distance_linear_space(guide, offtarget);
    assert_eq!(dist, dist_linear);

    println!("Guide RNA:        {}", std::str::from_utf8(guide).unwrap());
    println!("Off-target site:  {}", std::str::from_utf8(offtarget).unwrap());
    println!("Edit distance:    {}", dist);

    let (aligned_a, aligned_b) = align(guide, offtarget);
    let match_line: String = aligned_a
        .chars()
        .zip(aligned_b.chars())
        .map(|(a, b)| if a == b { '|' } else { ' ' })
        .collect();

    println!("\nAlignment:");
    println!("  {}", aligned_a);
    println!("  {}", match_line);
    println!("  {}", aligned_b);

    let threshold = 5;
    if dist <= threshold {
        println!("\n⚠  WARNING: This site is within {} edits — potential off-target!", threshold);
    } else {
        println!("\n✓  This site exceeds the {} edit threshold — likely safe.", threshold);
    }

    // Performance benchmark
    println!("\n=== Performance Benchmark ===");
    let seq_a: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
    let seq_b: Vec<u8> = (0..1000).map(|i| b"ACGT"[(i + 1) % 4]).collect();

    let start = Instant::now();
    let d = edit_distance(&seq_a, &seq_b);
    let elapsed = start.elapsed();
    println!("  1000×1000 matrix: distance={}, time={:?}", d, elapsed);

    let start = Instant::now();
    let d2 = edit_distance_linear_space(&seq_a, &seq_b);
    let elapsed2 = start.elapsed();
    println!("  1000×1000 linear: distance={}, time={:?}", d2, elapsed2);
    assert_eq!(d, d2);

    println!("\n=== Memory Usage ===");
    println!("  Full matrix:    {} bytes ({} KB)", (1001 * 1001 * 4), (1001 * 1001 * 4) / 1024);
    println!("  Linear space:   {} bytes",  1001 * 4 * 2);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        assert_eq!(edit_distance(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn test_empty() {
        assert_eq!(edit_distance(b"", b"ACGT"), 4);
        assert_eq!(edit_distance(b"ACGT", b""), 4);
    }

    #[test]
    fn test_single_substitution() {
        assert_eq!(edit_distance(b"ACGT", b"ACTT"), 1);
    }

    #[test]
    fn test_single_insertion() {
        assert_eq!(edit_distance(b"ACT", b"ACGT"), 1);
    }

    #[test]
    fn test_single_deletion() {
        assert_eq!(edit_distance(b"ACGT", b"ACT"), 1);
    }

    #[test]
    fn test_completely_different() {
        assert_eq!(edit_distance(b"AAAA", b"TTTT"), 4);
    }

    #[test]
    fn test_crispr_example() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let offtarget = b"GAGTCCAAGCAGAAGATGAA";
        let d = edit_distance(guide, offtarget);
        assert!(d <= 5, "Expected ≤5 edits, got {}", d);
    }

    #[test]
    fn test_linear_space_agrees() {
        let pairs = vec![
            (b"ACGT".as_ref(), b"ACGT".as_ref()),
            (b"KITTEN", b"SITTING"),
            (b"GATTACA", b"GCATGCU"),
            (b"", b"ABC"),
        ];
        for (a, b) in pairs {
            assert_eq!(
                edit_distance(a, b),
                edit_distance_linear_space(a, b),
                "Mismatch for {:?} vs {:?}",
                std::str::from_utf8(a),
                std::str::from_utf8(b)
            );
        }
    }
}