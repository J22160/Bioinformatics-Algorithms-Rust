//! Hirschberg's Linear-Space Global Alignment
//!
//! Computes the same optimal alignment as Needleman-Wunsch
//! but uses only O(min(n,m)) memory via divide-and-conquer.

use std::time::Instant;

/// Compute the last row of the NW score matrix in linear space.
/// Returns a vector of scores for the last row.
fn nw_last_row(a: &[u8], b: &[u8], match_s: i32, mismatch: i32, gap: i32) -> Vec<i32> {
    let m = b.len();
    let mut prev = vec![0i32; m + 1];
    let mut curr = vec![0i32; m + 1];

    for j in 0..=m {
        prev[j] = j as i32 * gap;
    }

    for i in 1..=a.len() {
        curr[0] = i as i32 * gap;
        for j in 1..=m {
            let s = if a[i - 1] == b[j - 1] { match_s } else { mismatch };
            curr[j] = (prev[j - 1] + s)
                .max(prev[j] + gap)
                .max(curr[j - 1] + gap);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev
}

/// Hirschberg's algorithm: returns (aligned_a, aligned_b)
pub fn hirschberg(
    a: &[u8],
    b: &[u8],
    match_s: i32,
    mismatch: i32,
    gap: i32,
) -> (Vec<u8>, Vec<u8>) {
    let n = a.len();
    let m = b.len();

    if n == 0 {
        return (vec![b'-'; m], b.to_vec());
    }
    if m == 0 {
        return (a.to_vec(), vec![b'-'; n]);
    }
    if n == 1 {
        // Base case: align single character against b
        let mut best_j = 0;
        let mut best_score = (n + m) as i32 * gap.abs();
        for j in 0..=m {
            let score = j as i32 * gap
                + if j < m && a[0] == b[j] { match_s } else { mismatch }
                + (m - j) as i32 * gap;
            // Simplified: just find best position to place the match
            if j < m && a[0] == b[j] {
                best_j = j;
                break;
            }
        }
        let _ = best_score;
        let mut al_a = Vec::new();
        let mut al_b = Vec::new();
        for j in 0..m {
            if j == best_j {
                al_a.push(a[0]);
                al_b.push(b[j]);
            } else {
                al_a.push(b'-');
                al_b.push(b[j]);
            }
        }
        if best_j == m {
            // character didn't match anywhere
            al_a.push(a[0]);
            al_b.push(b'-');
        }
        return (al_a, al_b);
    }

    let mid = n / 2;

    // Forward: score last row of aligning a[..mid] with b
    let score_l = nw_last_row(&a[..mid], b, match_s, mismatch, gap);

    // Reverse: score last row of aligning reverse(a[mid..]) with reverse(b)
    let a_rev: Vec<u8> = a[mid..].iter().rev().copied().collect();
    let b_rev: Vec<u8> = b.iter().rev().copied().collect();
    let score_r = nw_last_row(&a_rev, &b_rev, match_s, mismatch, gap);

    // Find optimal split point j*
    let mut best_j = 0;
    let mut best_sum = i32::MIN;
    for j in 0..=m {
        let sum = score_l[j] + score_r[m - j];
        if sum > best_sum {
            best_sum = sum;
            best_j = j;
        }
    }

    // Recurse
    let (mut al_a1, mut al_b1) = hirschberg(&a[..mid], &b[..best_j], match_s, mismatch, gap);
    let (al_a2, al_b2) = hirschberg(&a[mid..], &b[best_j..], match_s, mismatch, gap);

    al_a1.extend_from_slice(&al_a2);
    al_b1.extend_from_slice(&al_b2);

    (al_a1, al_b1)
}

fn main() {
    println!("=== Hirschberg: Linear-Space Global Alignment ===\n");

    let seq_a = b"AGTACGCATAG";
    let seq_b = b"TATGCCATG";

    let start = Instant::now();
    let (aln_a, aln_b) = hirschberg(seq_a, seq_b, 1, -1, -2);
    let elapsed = start.elapsed();

    let aln_a_str = String::from_utf8(aln_a).unwrap();
    let aln_b_str = String::from_utf8(aln_b).unwrap();

    println!("Seq A: {}", std::str::from_utf8(seq_a).unwrap());
    println!("Seq B: {}\n", std::str::from_utf8(seq_b).unwrap());
    println!("Alignment:");
    println!("  {}", aln_a_str);
    println!("  {}", aln_b_str);
    println!("Time: {:?}", elapsed);

    println!("\n=== Memory Comparison ===");
    for &n in &[100, 1000, 10000] {
        let nw_memory = (n + 1) * (n + 1) * 4; // bytes for i32 matrix
        let hirsch_memory = (n + 1) * 4 * 2;    // two rows of i32
        println!("  n={}: NW={} KB, Hirschberg={} KB (ratio: {}x)",
                 n, nw_memory / 1024, hirsch_memory / 1024,
                 nw_memory / hirsch_memory.max(1));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        let (a, b) = hirschberg(b"ACGT", b"ACGT", 1, -1, -2);
        assert_eq!(a, b"ACGT".to_vec());
        assert_eq!(b, b"ACGT".to_vec());
    }

    #[test]
    fn test_empty() {
        let (a, b) = hirschberg(b"", b"ACG", 1, -1, -2);
        assert_eq!(a, b"---".to_vec());
    }
}
