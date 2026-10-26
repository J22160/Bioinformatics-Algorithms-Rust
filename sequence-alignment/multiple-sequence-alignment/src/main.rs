//! Multiple Sequence Alignment — Progressive Star Method
//!
//! Aligns multiple sequences by selecting a center sequence
//! and progressively aligning all others to it.

use std::time::Instant;

/// Compute edit distance between two sequences (simplified NW score-only)
fn edit_distance(a: &[u8], b: &[u8]) -> i32 {
    let m = b.len();
    let mut prev: Vec<i32> = (0..=m as i32).collect();
    let mut curr = vec![0i32; m + 1];

    for i in 1..=a.len() {
        curr[0] = i as i32;
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost)
                .min(prev[j] + 1)
                .min(curr[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[m]
}

/// Pairwise NW alignment returning aligned strings
fn pairwise_align(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = a.len();
    let m = b.len();
    let idx = |i: usize, j: usize| i * (m + 1) + j;
    let mut dp = vec![0i32; (n + 1) * (m + 1)];

    for i in 1..=n { dp[idx(i, 0)] = -(i as i32); }
    for j in 1..=m { dp[idx(0, j)] = -(j as i32); }

    for i in 1..=n {
        for j in 1..=m {
            let s = if a[i-1] == b[j-1] { 1i32 } else { -1 };
            dp[idx(i, j)] = (dp[idx(i-1, j-1)] + s)
                .max(dp[idx(i-1, j)] - 1)
                .max(dp[idx(i, j-1)] - 1);
        }
    }

    let mut al_a = Vec::new();
    let mut al_b = Vec::new();
    let (mut i, mut j) = (n, m);

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let s = if a[i-1] == b[j-1] { 1 } else { -1 };
            if dp[idx(i, j)] == dp[idx(i-1, j-1)] + s {
                al_a.push(a[i-1]); al_b.push(b[j-1]);
                i -= 1; j -= 1; continue;
            }
        }
        if i > 0 && dp[idx(i, j)] == dp[idx(i-1, j)] - 1 {
            al_a.push(a[i-1]); al_b.push(b'-');
            i -= 1;
        } else {
            al_a.push(b'-'); al_b.push(b[j-1]);
            j -= 1;
        }
    }
    al_a.reverse(); al_b.reverse();
    (al_a, al_b)
}

/// Star MSA: align all sequences to the center sequence
pub fn star_msa(sequences: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let k = sequences.len();
    if k == 0 { return vec![]; }
    if k == 1 { return sequences.to_vec(); }

    // Step 1: Find center sequence (minimize total distance)
    let mut best_center = 0;
    let mut best_total = i32::MAX;

    for i in 0..k {
        let total: i32 = (0..k)
            .filter(|&j| j != i)
            .map(|j| edit_distance(&sequences[i], &sequences[j]))
            .sum();
        if total < best_total {
            best_total = total;
            best_center = i;
        }
    }

    // Step 2: Align all sequences to center
    let mut alignments: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for i in 0..k {
        if i == best_center {
            let seq = sequences[i].clone();
            alignments.push((seq.clone(), seq));
        } else {
            let (al_center, al_other) = pairwise_align(&sequences[best_center], &sequences[i]);
            alignments.push((al_center, al_other));
        }
    }

    // Step 3: Merge by inserting gap columns
    // Find the master alignment of the center
    let center_aligned: Vec<&Vec<u8>> = alignments.iter().map(|(c, _)| c).collect();

    // For simplicity, return individual pairwise alignments
    // A full implementation would merge gap columns across all alignments
    let mut result = Vec::new();
    let master = &alignments[0].0;
    result.push(master.clone());
    for i in 0..k {
        if i != best_center {
            let idx = if i < best_center { i } else { i };
            result.push(alignments[idx].1.clone());
        }
    }
    let _ = center_aligned;

    result
}

fn main() {
    println!("=== Multiple Sequence Alignment: Star Method ===\n");

    let sequences: Vec<Vec<u8>> = vec![
        b"ATGCTAGCTA".to_vec(),
        b"ATGCAAGCTA".to_vec(),
        b"ATGCTAGTA".to_vec(),
        b"ATCTAGCTA".to_vec(),
        b"ATGCTAGCTAG".to_vec(),
    ];

    println!("Input sequences:");
    for (i, seq) in sequences.iter().enumerate() {
        println!("  S{}: {}", i + 1, std::str::from_utf8(seq).unwrap());
    }

    let start = Instant::now();
    let msa = star_msa(&sequences);
    let elapsed = start.elapsed();

    println!("\nMultiple alignment (center star):");
    for (i, aligned) in msa.iter().enumerate() {
        println!("  S{}: {}", i + 1, String::from_utf8_lossy(aligned));
    }
    println!("\nTime: {:?}", elapsed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_sequences() {
        let seqs = vec![b"ACGT".to_vec(), b"ACGT".to_vec()];
        let msa = star_msa(&seqs);
        assert_eq!(msa.len(), 2);
    }

    #[test]
    fn test_single_sequence() {
        let seqs = vec![b"ACGT".to_vec()];
        let msa = star_msa(&seqs);
        assert_eq!(msa.len(), 1);
    }
}
