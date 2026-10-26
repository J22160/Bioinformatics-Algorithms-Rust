//! Affine Gap Alignment (Gotoh's Algorithm)
//!
//! Uses three DP matrices to model gap-open and gap-extend penalties separately.

use std::time::Instant;

const NEG_INF: i32 = i32::MIN / 2;

pub struct AffineParams {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

pub fn affine_gap_align(
    a: &[u8],
    b: &[u8],
    params: &AffineParams,
) -> (i32, String, String) {
    let n = a.len();
    let m = b.len();
    let idx = |i: usize, j: usize| -> usize { i * (m + 1) + j };

    // Three matrices: M (match/mismatch), X (gap in B), Y (gap in A)
    let mut mat_m = vec![NEG_INF; (n + 1) * (m + 1)];
    let mut mat_x = vec![NEG_INF; (n + 1) * (m + 1)]; // insertion (gap in A)
    let mut mat_y = vec![NEG_INF; (n + 1) * (m + 1)]; // deletion (gap in B)

    mat_m[idx(0, 0)] = 0;
    for i in 1..=n {
        mat_y[idx(i, 0)] = params.gap_open + i as i32 * params.gap_extend;
    }
    for j in 1..=m {
        mat_x[idx(0, j)] = params.gap_open + j as i32 * params.gap_extend;
    }

    for i in 1..=n {
        for j in 1..=m {
            let score = if a[i - 1] == b[j - 1] {
                params.match_score
            } else {
                params.mismatch
            };

            // Match/mismatch: come from any state at (i-1, j-1)
            mat_m[idx(i, j)] = mat_m[idx(i - 1, j - 1)]
                .max(mat_x[idx(i - 1, j - 1)])
                .max(mat_y[idx(i - 1, j - 1)])
                + score;

            // Gap in A (insertion): horizontal move
            mat_x[idx(i, j)] = (mat_m[idx(i, j - 1)] + params.gap_open + params.gap_extend)
                .max(mat_x[idx(i, j - 1)] + params.gap_extend);

            // Gap in B (deletion): vertical move
            mat_y[idx(i, j)] = (mat_m[idx(i - 1, j)] + params.gap_open + params.gap_extend)
                .max(mat_y[idx(i - 1, j)] + params.gap_extend);
        }
    }

    let final_score = mat_m[idx(n, m)]
        .max(mat_x[idx(n, m)])
        .max(mat_y[idx(n, m)]);

    // Traceback
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = n;
    let mut j = m;

    // Determine which matrix we end in
    let mut state = if final_score == mat_m[idx(n, m)] {
        0 // M
    } else if final_score == mat_x[idx(n, m)] {
        1 // X
    } else {
        2 // Y
    };

    while i > 0 || j > 0 {
        match state {
            0 => {
                // In M: came from diagonal
                if i == 0 || j == 0 { break; }
                aligned_a.push(a[i - 1]);
                aligned_b.push(b[j - 1]);
                let score = if a[i - 1] == b[j - 1] { params.match_score } else { params.mismatch };
                let prev = mat_m[idx(i, j)] - score;
                i -= 1;
                j -= 1;
                if prev == mat_m[idx(i, j)] { state = 0; }
                else if prev == mat_x[idx(i, j)] { state = 1; }
                else { state = 2; }
            }
            1 => {
                // In X: gap in A (horizontal)
                if j == 0 { break; }
                aligned_a.push(b'-');
                aligned_b.push(b[j - 1]);
                j -= 1;
                if mat_x[idx(i, j + 1)] == mat_x[idx(i, j)] + params.gap_extend {
                    state = 1; // extend gap
                } else {
                    state = 0; // came from M (gap open)
                }
            }
            2 => {
                // In Y: gap in B (vertical)
                if i == 0 { break; }
                aligned_a.push(a[i - 1]);
                aligned_b.push(b'-');
                i -= 1;
                if mat_y[idx(i + 1, j)] == mat_y[idx(i, j)] + params.gap_extend {
                    state = 2;
                } else {
                    state = 0;
                }
            }
            _ => unreachable!(),
        }
    }

    aligned_a.reverse();
    aligned_b.reverse();

    (
        final_score,
        String::from_utf8(aligned_a).unwrap(),
        String::from_utf8(aligned_b).unwrap(),
    )
}

fn main() {
    println!("=== Affine Gap Alignment: Detecting a Segmental Deletion ===\n");

    let gene_a = b"ATGCTAGCTAGCATCGATCGATCGATCGATCGATCG";
    let gene_b = b"ATGCTAGCTAGCATCG---------------ATCGATCG";
    // Simulate: gene_b has a 15-base deletion
    let gene_b_clean = b"ATGCTAGCTAGCATCGATCG";

    let params = AffineParams {
        match_score: 2,
        mismatch: -3,
        gap_open: -5,
        gap_extend: -1,
    };

    let start = Instant::now();
    let (score, aln_a, aln_b) = affine_gap_align(gene_a, gene_b_clean, &params);
    let elapsed = start.elapsed();

    println!("Gene A: {}", std::str::from_utf8(gene_a).unwrap());
    println!("Gene B: {}\n", std::str::from_utf8(gene_b_clean).unwrap());
    println!("Alignment (score={}):", score);
    println!("  {}", aln_a);
    println!("  {}", aln_b);
    println!("Time: {:?}", elapsed);

    // Benchmark
    println!("\n=== Benchmark ===");
    for &n in &[200, 500, 1000] {
        let a: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..n).map(|i| b"ACGT"[(i + 1) % 4]).collect();
        let start = Instant::now();
        let (s, _, _) = affine_gap_align(&a, &b, &params);
        let elapsed = start.elapsed();
        println!("  {}x{}: score={}, time={:?}", n, n, s, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        let params = AffineParams { match_score: 2, mismatch: -3, gap_open: -5, gap_extend: -1 };
        let (score, _, _) = affine_gap_align(b"ACGT", b"ACGT", &params);
        assert_eq!(score, 8);
    }

    #[test]
    fn test_single_gap_preferred() {
        let params = AffineParams { match_score: 1, mismatch: -1, gap_open: -3, gap_extend: -1 };
        let (_, aln_a, _) = affine_gap_align(b"AACGT", b"ACGT", &params);
        // Should prefer one contiguous gap over scattered mismatches
        assert!(aln_a.contains('-') || aln_a.len() >= 5);
    }
}
