//! FM-Index — Compressed Full-Text Index
//!
//! Enables pattern matching in O(m) time independent of text length.
//! Core data structure of BWA and Bowtie2 read aligners.

use std::time::Instant;

const ALPHABET: [u8; 5] = [b'$', b'A', b'C', b'G', b'T'];

/// Simple FM-Index with full occurrence table
pub struct FMIndex {
    bwt: Vec<u8>,
    sa: Vec<usize>,        // suffix array for locating
    c_array: [usize; 256], // C[c] = count of chars < c
    occ: Vec<[usize; 256]>, // Occ[i][c] = count of c in bwt[0..=i]
    n: usize,
}

impl FMIndex {
    /// Build FM-Index from text.
    ///
    /// # Rust Performance Advantages
    ///
    /// 1. **Fixed-size arrays for C and Occ**: The C array and per-position
    ///    occurrence counts use [usize; 256], which is stack-allocated
    ///    and provides O(1) lookup by character value.
    ///
    /// 2. **Contiguous Occ table**: Vec<[usize; 256]> stores all occurrence
    ///    data in a single allocation, giving excellent cache locality
    ///    during backward search.
    ///
    /// 3. **No hash map overhead**: Python would use dict for character
    ///    counting, adding hash computation + collision handling per lookup.
    ///    Rust uses direct array indexing.
    pub fn new(text: &[u8]) -> Self {
        let n = text.len();

        // Build suffix array
        let mut sa: Vec<usize> = (0..n).collect();
        sa.sort_unstable_by(|&a, &b| text[a..].cmp(&text[b..]));

        // Build BWT from suffix array
        let bwt: Vec<u8> = sa.iter()
            .map(|&i| if i == 0 { text[n - 1] } else { text[i - 1] })
            .collect();

        // Build C array
        let mut counts = [0usize; 256];
        for &b in &bwt {
            counts[b as usize] += 1;
        }
        let mut c_array = [0usize; 256];
        let mut total = 0;
        for i in 0..256 {
            c_array[i] = total;
            total += counts[i];
        }

        // Build full occurrence table
        let mut occ = vec![[0usize; 256]; n + 1];
        for i in 0..n {
            occ[i + 1] = occ[i];
            occ[i + 1][bwt[i] as usize] += 1;
        }

        FMIndex { bwt, sa, c_array, occ, n }
    }

    /// Count occurrences of pattern using backward search.
    /// Returns the SA range [lo, hi) where all matching suffixes reside.
    pub fn backward_search(&self, pattern: &[u8]) -> (usize, usize) {
        let mut lo = 0;
        let mut hi = self.n;

        for &c in pattern.iter().rev() {
            let cu = c as usize;
            lo = self.c_array[cu] + self.occ[lo][cu];
            hi = self.c_array[cu] + self.occ[hi][cu];
            if lo >= hi {
                return (0, 0); // pattern not found
            }
        }
        (lo, hi)
    }

    /// Count how many times pattern occurs in text.
    pub fn count(&self, pattern: &[u8]) -> usize {
        let (lo, hi) = self.backward_search(pattern);
        hi - lo
    }

    /// Locate all positions where pattern occurs in text.
    pub fn locate(&self, pattern: &[u8]) -> Vec<usize> {
        let (lo, hi) = self.backward_search(pattern);
        (lo..hi).map(|i| self.sa[i]).collect()
    }
}

fn main() {
    println!("=== FM-Index: Genomic Pattern Search ===\n");

    let text = b"ATGCTAGCTAGCATCGATCGATCG$";
    let fm = FMIndex::new(text);

    println!("Text: {}", std::str::from_utf8(text).unwrap());
    println!("Text length: {}\n", text.len());

    let patterns = [b"ATCG" as &[u8], b"TAG", b"CGA", b"ZZZZZ"];
    for pat in &patterns {
        let pat_str = std::str::from_utf8(pat).unwrap();
        let count = fm.count(pat);
        let positions = fm.locate(pat);
        println!(
            "Pattern '{}': {} occurrences at positions {:?}",
            pat_str, count, positions
        );
    }

    // Demonstrate O(m) search independence from text length
    println!("\n=== Search Time vs Text Length ===");
    let pattern = b"ACGT";
    for &n in &[10_000, 50_000, 100_000] {
        let mut seq: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        seq.push(b'$');

        let build_start = Instant::now();
        let fm = FMIndex::new(&seq);
        let build_time = build_start.elapsed();

        let search_start = Instant::now();
        for _ in 0..10_000 {
            let _ = fm.count(pattern);
        }
        let search_time = search_start.elapsed();

        println!(
            "  n={:<8} build={:>10.3?}  10K searches={:>10.3?}",
            n, build_time, search_time
        );
    }
    println!("Notice: search time stays constant as text grows!");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count() {
        let fm = FMIndex::new(b"ABRACADABRA$");
        assert_eq!(fm.count(b"ABRA"), 2);
        assert_eq!(fm.count(b"A"), 5);
        assert_eq!(fm.count(b"Z"), 0);
    }

    #[test]
    fn test_locate() {
        let fm = FMIndex::new(b"ABRACADABRA$");
        let mut positions = fm.locate(b"ABRA");
        positions.sort();
        assert_eq!(positions, vec![0, 7]);
    }

    #[test]
    fn test_dna_pattern() {
        let fm = FMIndex::new(b"ATGATGATG$");
        assert_eq!(fm.count(b"ATG"), 3);
        assert_eq!(fm.count(b"GAT"), 2);
    }
}
