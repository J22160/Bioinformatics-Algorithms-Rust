//! Selection Sort — Sorting Phred Quality Scores
//!
//! A simple O(n²) sorting algorithm that demonstrates the concept of
//! repeatedly finding the minimum element. In bioinformatics, sorting
//! quality scores helps identify low-quality reads for filtering.

use std::time::Instant;

/// Selection sort: repeatedly find the minimum of the unsorted portion
/// and swap it into the correct position.
///
/// # Performance Notes (Rust vs. Python)
///
/// Rust advantages here:
/// - In-place sorting on a contiguous `Vec<u8>` means excellent cache locality.
///   Each Phred score is a single byte, so 64 scores fit in one cache line.
/// - No bounds checking overhead in release mode (the compiler proves safety
///   at compile time, then removes the checks).
/// - The inner loop compiles to tight x86 instructions with no interpreter dispatch.
///
/// Even though selection sort is O(n²) in both languages, Rust is typically
/// 50-100x faster than Python for this kind of tight numerical loop because
/// Python's `for` loop has per-iteration interpreter overhead (~100ns/iteration)
/// while Rust's compiled loop runs at ~1ns/iteration.
pub fn selection_sort<T: PartialOrd>(arr: &mut [T]) {
    let n = arr.len();
    for i in 0..n {
        // Find index of minimum element in arr[i..n]
        let mut min_idx = i;
        for j in (i + 1)..n {
            if arr[j] < arr[min_idx] {
                min_idx = j;
            }
        }
        // Swap the found minimum into position i
        // Rust's slice::swap is safe and compiles to a simple register swap
        if min_idx != i {
            arr.swap(i, min_idx);
        }
    }
}

/// Filter reads below a quality threshold after sorting.
/// Returns indices of reads that pass the quality filter.
fn filter_low_quality(sorted_scores: &[u8], threshold: u8) -> usize {
    // After sorting, all scores below threshold are at the front.
    // Binary search finds the partition point in O(log n).
    sorted_scores.partition_point(|&score| score < threshold)
}

fn main() {
    println!("=== Selection Sort: Phred Quality Score Filtering ===\n");

    // Simulated Phred quality scores for sequencing reads
    // Phred scores typically range from 0-41 in Illumina data
    let mut scores: Vec<u8> = vec![
        38, 12, 27, 35, 8, 41, 19, 33, 5, 40, 28, 15, 37, 22, 9, 31,
        2, 39, 17, 36, 11, 30, 24, 7, 41, 13, 26, 34, 20, 3, 38, 16,
    ];

    println!("Original scores: {:?}", scores);
    println!("Number of reads: {}\n", scores.len());

    let start = Instant::now();
    selection_sort(&mut scores);
    let elapsed = start.elapsed();

    println!("Sorted scores:   {:?}", scores);
    println!("Sort time:       {:?}\n", elapsed);

    // Apply quality filter
    let threshold: u8 = 20;
    let cutoff = filter_low_quality(&scores, threshold);
    println!("Quality threshold: Q{}", threshold);
    println!("Reads FAILING filter (Q<{}): {} reads — {:?}", threshold, cutoff, &scores[..cutoff]);
    println!("Reads PASSING filter (Q≥{}): {} reads", threshold, scores.len() - cutoff);

    // Demonstrate scaling behavior
    println!("\n=== Scaling Benchmark ===");
    for &n in &[100, 1_000, 5_000, 10_000] {
        let mut data: Vec<u8> = (0..n).map(|i| ((i * 7 + 13) % 42) as u8).collect();
        let start = Instant::now();
        selection_sort(&mut data);
        let elapsed = start.elapsed();
        println!("  n={:<6}  time={:>10.3?}  (verify sorted: {})", n, elapsed, data.windows(2).all(|w| w[0] <= w[1]));
    }
    println!("\nNotice: doubling n roughly quadruples the time → O(n²) confirmed.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let mut arr: Vec<u8> = vec![];
        selection_sort(&mut arr);
        assert_eq!(arr, vec![]);
    }

    #[test]
    fn test_single() {
        let mut arr = vec![42u8];
        selection_sort(&mut arr);
        assert_eq!(arr, vec![42]);
    }

    #[test]
    fn test_already_sorted() {
        let mut arr = vec![1u8, 2, 3, 4, 5];
        selection_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_reverse_sorted() {
        let mut arr = vec![5u8, 4, 3, 2, 1];
        selection_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_phred_scores() {
        let mut scores = vec![38u8, 12, 27, 35, 8, 41, 19];
        selection_sort(&mut scores);
        assert_eq!(scores, vec![8, 12, 19, 27, 35, 38, 41]);
    }

    #[test]
    fn test_duplicates() {
        let mut arr = vec![3u8, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        selection_sort(&mut arr);
        assert_eq!(arr, vec![1, 1, 2, 3, 3, 4, 5, 5, 5, 6, 9]);
    }
}
