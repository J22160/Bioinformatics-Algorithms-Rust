//! Burrows-Wheeler Transform (BWT)
//!
//! Reversible text transform that groups similar characters together,
//! enabling both compression and fast pattern searching.
//! This is the foundation of BWA and Bowtie2 read aligners.

use std::time::Instant;

/// Compute the BWT of a text using suffix array construction.
///
/// # Rust Performance Advantages
///
/// 1. **In-place suffix sorting**: We sort integer indices rather than
///    copying strings. Rust's sort_unstable avoids allocation.
///
/// 2. **Byte-level operations**: Text is &[u8], giving direct memory
///    access without Unicode overhead. Python str requires UTF-8 decoding.
///
/// 3. **Zero-cost abstractions**: The iterator chain in `bwt_from_sa`
///    compiles to a single loop with no intermediate allocations.
///
/// 4. **Cache-friendly access**: Suffix array is a contiguous Vec<usize>,
///    so sequential access during BWT extraction is prefetcher-friendly.
pub fn suffix_array_naive(text: &[u8]) -> Vec<usize> {
    let n = text.len();
    let mut sa: Vec<usize> = (0..n).collect();
    sa.sort_unstable_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

/// Build BWT from a precomputed suffix array.
pub fn bwt_from_sa(text: &[u8], sa: &[usize]) -> Vec<u8> {
    let n = text.len();
    sa.iter()
        .map(|&i| {
            if i == 0 {
                text[n - 1]
            } else {
                text[i - 1]
            }
        })
        .collect()
}

/// Inverse BWT: reconstruct original text from BWT string.
///
/// Uses the LF-mapping property: the i-th occurrence of character c
/// in the last column (BWT) corresponds to the i-th occurrence
/// of c in the first column (sorted BWT).
pub fn inverse_bwt(bwt: &[u8]) -> Vec<u8> {
    let n = bwt.len();
    if n == 0 {
        return vec![];
    }

    // Build the first column (sorted BWT)
    let mut first_col = bwt.to_vec();
    first_col.sort_unstable();

    // Build LF-mapping using counting sort logic
    // C[c] = number of characters in BWT lexicographically smaller than c
    let mut counts = [0usize; 256];
    for &b in bwt {
        counts[b as usize] += 1;
    }
    let mut c_array = [0usize; 256];
    let mut total = 0;
    for i in 0..256 {
        c_array[i] = total;
        total += counts[i];
    }

    // Occ[i] = how many times bwt[i] appears in bwt[0..i]
    let mut occ = vec![0usize; n];
    let mut running = [0usize; 256];
    for i in 0..n {
        let c = bwt[i] as usize;
        occ[i] = running[c];
        running[c] += 1;
    }

    // LF-mapping: LF[i] = C[bwt[i]] + Occ[i]
    let lf: Vec<usize> = (0..n)
        .map(|i| c_array[bwt[i] as usize] + occ[i])
        .collect();

    // Find the position of '$' (sentinel) in BWT
    let start = bwt.iter().position(|&b| b == b'$').unwrap_or(0);

    // Reconstruct by following LF-mapping from '$'
    let mut result = vec![0u8; n];
    let mut idx = start;
    for i in (0..n).rev() {
        result[i] = bwt[idx];
        idx = lf[idx];
    }

    result
}

/// Count character occurrences — useful for demonstrating BWT clustering
fn count_runs(data: &[u8]) -> usize {
    if data.is_empty() {
        return 0;
    }
    let mut runs = 1;
    for i in 1..data.len() {
        if data[i] != data[i - 1] {
            runs += 1;
        }
    }
    runs
}

fn main() {
    println!("=== Burrows-Wheeler Transform: Genome Compression ===\n");

    // Example: small genome fragment
    let text = b"ACAACG$";
    println!("Original text: {}", std::str::from_utf8(text).unwrap());

    let sa = suffix_array_naive(text);
    let bwt = bwt_from_sa(text, &sa);
    let bwt_str = String::from_utf8_lossy(&bwt);

    println!("Suffix Array:  {:?}", sa);
    println!("BWT:           {}", bwt_str);

    // Show clustering effect
    let original_runs = count_runs(text);
    let bwt_runs = count_runs(&bwt);
    println!("\nCharacter runs in original: {}", original_runs);
    println!("Character runs in BWT:      {}", bwt_runs);
    println!("(Fewer runs = more compressible)\n");

    // Verify inverse
    let reconstructed = inverse_bwt(&bwt);
    println!("Reconstructed: {}", String::from_utf8_lossy(&reconstructed));
    assert_eq!(text.as_ref(), reconstructed.as_slice());
    println!("✓ Inverse BWT perfectly reconstructs original\n");

    // Larger example with repetitive DNA
    println!("=== Compression Demo on Repetitive Sequence ===\n");
    let repetitive = b"ATGATGATGATGATGATGATGATGATGATG$";
    let sa2 = suffix_array_naive(repetitive);
    let bwt2 = bwt_from_sa(repetitive, &sa2);
    println!(
        "Original ({} bp): {} runs",
        repetitive.len(),
        count_runs(repetitive)
    );
    println!(
        "BWT      ({} bp): {} runs",
        bwt2.len(),
        count_runs(&bwt2)
    );
    println!("BWT: {}", String::from_utf8_lossy(&bwt2));

    // Benchmark
    println!("\n=== Performance Benchmark ===");
    for &n in &[1_000, 10_000, 100_000] {
        let mut seq: Vec<u8> = (0..n)
            .map(|i| b"ACGT"[i % 4])
            .collect();
        seq.push(b'$');

        let start = Instant::now();
        let sa = suffix_array_naive(&seq);
        let _bwt = bwt_from_sa(&seq, &sa);
        let elapsed = start.elapsed();
        println!("  n={}: BWT construction time={:?}", n, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_banana() {
        let text = b"BANANA$";
        let sa = suffix_array_naive(text);
        let bwt = bwt_from_sa(text, &sa);
        assert_eq!(bwt, b"ANNB$AA");
    }

    #[test]
    fn test_inverse_bwt() {
        let text = b"ACAACG$";
        let sa = suffix_array_naive(text);
        let bwt = bwt_from_sa(text, &sa);
        let reconstructed = inverse_bwt(&bwt);
        assert_eq!(text.as_ref(), reconstructed.as_slice());
    }

    #[test]
    fn test_dna_sequence() {
        let text = b"GATTACA$";
        let sa = suffix_array_naive(text);
        let bwt = bwt_from_sa(text, &sa);
        let reconstructed = inverse_bwt(&bwt);
        assert_eq!(text.as_ref(), reconstructed.as_slice());
    }
}
