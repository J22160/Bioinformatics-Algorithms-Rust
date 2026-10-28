//! Suffix Arrays — Sorted Index of All Suffixes

use std::time::Instant;

/// Naive suffix array construction: sort indices by their suffixes.
pub fn suffix_array_naive(text: &[u8]) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..text.len()).collect();
    sa.sort_unstable_by(|&a, &b| text[a..].cmp(&text[b..]));
    sa
}

/// Binary search on suffix array to find all occurrences of pattern.
/// Returns (first, last+1) indices in the suffix array.
pub fn sa_search(text: &[u8], sa: &[usize], pattern: &[u8]) -> (usize, usize) {
    let n = sa.len();
    let m = pattern.len();

    // Find lower bound
    let lo = {
        let (mut left, mut right) = (0, n);
        while left < right {
            let mid = left + (right - left) / 2;
            let suffix = &text[sa[mid]..];
            let cmp_len = m.min(suffix.len());
            if &suffix[..cmp_len] < pattern {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        left
    };

    // Find upper bound
    let hi = {
        let (mut left, mut right) = (lo, n);
        while left < right {
            let mid = left + (right - left) / 2;
            let suffix = &text[sa[mid]..];
            let cmp_len = m.min(suffix.len());
            if &suffix[..cmp_len] <= pattern {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        left
    };

    (lo, hi)
}

/// Compute LCP (Longest Common Prefix) array using Kasai's algorithm.
/// LCP[i] = length of longest common prefix between SA[i] and SA[i-1].
pub fn build_lcp(text: &[u8], sa: &[usize]) -> Vec<usize> {
    let n = sa.len();
    let mut rank = vec![0usize; n];
    for i in 0..n {
        rank[sa[i]] = i;
    }

    let mut lcp = vec![0usize; n];
    let mut h = 0usize;
    for i in 0..n {
        if rank[i] > 0 {
            let j = sa[rank[i] - 1];
            while i + h < n && j + h < n && text[i + h] == text[j + h] {
                h += 1;
            }
            lcp[rank[i]] = h;
            if h > 0 { h -= 1; }
        }
    }
    lcp
}

fn main() {
    println!("=== Suffix Arrays: Genomic Motif Search ===\n");

    let text = b"ATCGATCGATCG$";
    let sa = suffix_array_naive(text);

    println!("Text: {}", std::str::from_utf8(text).unwrap());
    println!("Suffix Array: {:?}\n", sa);

    // Display sorted suffixes
    println!("Sorted suffixes:");
    for (rank, &pos) in sa.iter().enumerate() {
        let suffix = std::str::from_utf8(&text[pos..]).unwrap();
        println!("  SA[{}] = {} → {}", rank, pos, suffix);
    }

    // Search for patterns
    println!("\n--- Pattern Search ---");
    for pattern in &[b"ATCG" as &[u8], b"GAT", b"CGA", b"TTT"] {
        let (lo, hi) = sa_search(text, &sa, pattern);
        let positions: Vec<usize> = (lo..hi).map(|i| sa[i]).collect();
        println!(
            "  '{}': {} occurrences at {:?}",
            std::str::from_utf8(pattern).unwrap(),
            hi - lo,
            positions
        );
    }

    // LCP array
    let lcp = build_lcp(text, &sa);
    println!("\nLCP Array: {:?}", lcp);
    println!("(LCP enables O(m + log n) search instead of O(m log n))");

    // Benchmark
    println!("\n=== Benchmark ===");
    for &n in &[1_000, 10_000, 100_000] {
        let mut seq: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        seq.push(b'$');
        let start = Instant::now();
        let sa = suffix_array_naive(&seq);
        let elapsed = start.elapsed();
        let start2 = Instant::now();
        let (lo, hi) = sa_search(&seq, &sa, b"ACGTACGT");
        let search_time = start2.elapsed();
        println!(
            "  n={:<8} build={:>10.3?}  search={:>10.3?}  found={}",
            n, elapsed, search_time, hi - lo
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suffix_array() {
        let sa = suffix_array_naive(b"BANANA$");
        assert_eq!(sa, vec![6, 5, 3, 1, 0, 4, 2]);
    }

    #[test]
    fn test_search() {
        let text = b"ABRACADABRA$";
        let sa = suffix_array_naive(text);
        let (lo, hi) = sa_search(text, &sa, b"ABRA");
        assert_eq!(hi - lo, 2);
    }
}