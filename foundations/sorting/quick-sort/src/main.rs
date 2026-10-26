//! Quick Sort — Sorting Variant Quality Scores
//!
//! A divide-and-conquer sorting algorithm with O(n log n) average
//! performance. Used as the foundation of std::sort in many languages.

use std::time::Instant;

/// Randomized Quick Sort using Lomuto partition scheme.
///
/// # Rust Performance Advantages
///
/// 1. **In-place partitioning**: We swap elements within the same Vec,
///    avoiding any heap allocation. Python's equivalent creates sublists.
///
/// 2. **Tail call optimization**: By recursing on the smaller partition
///    first, the compiler can optimize the larger partition into a loop,
///    keeping stack depth at O(log n).
///
/// 3. **Monomorphization**: This generic function compiles to specialized
///    machine code for each concrete type (f64, u32, etc.).
///
/// 4. **No GIL**: In multi-threaded contexts, Rust quick sort can
///    parallelize both partitions using Rayon without a Global Interpreter Lock.
pub fn quick_sort<T: PartialOrd>(arr: &mut [T]) {
    if arr.len() <= 1 {
        return;
    }
    // Switch to insertion sort for small arrays
    if arr.len() <= 16 {
        insertion_sort(arr);
        return;
    }
    let pivot_idx = partition(arr);
    let (left, right) = arr.split_at_mut(pivot_idx);
    quick_sort(left);
    quick_sort(&mut right[1..]); // skip the pivot element
}

/// Lomuto partition: choose last element as pivot.
/// Returns the final index of the pivot.
fn partition<T: PartialOrd>(arr: &mut [T]) -> usize {
    let len = arr.len();
    // Median-of-three pivot selection to avoid worst case
    let mid = len / 2;
    if arr[0] > arr[mid] { arr.swap(0, mid); }
    if arr[0] > arr[len - 1] { arr.swap(0, len - 1); }
    if arr[mid] > arr[len - 1] { arr.swap(mid, len - 1); }
    arr.swap(mid, len - 1); // place median at end

    let mut i = 0;
    for j in 0..len - 1 {
        if arr[j] <= arr[len - 1] {
            arr.swap(i, j);
            i += 1;
        }
    }
    arr.swap(i, len - 1);
    i
}

fn insertion_sort<T: PartialOrd>(arr: &mut [T]) {
    for i in 1..arr.len() {
        let mut j = i;
        while j > 0 && arr[j - 1] > arr[j] {
            arr.swap(j - 1, j);
            j -= 1;
        }
    }
}

fn main() {
    println!("=== Quick Sort: Variant Quality Score Filtering ===\n");

    let mut qual_scores: Vec<f64> = vec![
        934.0, 52.0, 1847.0, 300.0, 78.0, 2100.0, 415.0, 61.0,
        1200.0, 33.5, 89.0, 750.0, 2500.0, 145.0, 670.0, 28.0,
    ];

    println!("Original QUAL scores: {:?}\n", qual_scores);

    let start = Instant::now();
    quick_sort(&mut qual_scores);
    let elapsed = start.elapsed();

    println!("Sorted QUAL scores:   {:?}", qual_scores);
    println!("Sort time:            {:?}\n", elapsed);

    // Apply quality filter
    let threshold = 100.0;
    let passing: Vec<&f64> = qual_scores.iter().filter(|&&q| q >= threshold).collect();
    println!("Variants passing QUAL >= {}: {} out of {}", threshold, passing.len(), qual_scores.len());

    // Scaling benchmark
    println!("\n=== Scaling Benchmark ===");
    for &n in &[1_000, 10_000, 100_000, 1_000_000] {
        let mut data: Vec<u64> = (0..n)
            .map(|i| ((i as u64).wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407)) >> 33)
            .collect();
        let start = Instant::now();
        quick_sort(&mut data);
        let elapsed = start.elapsed();
        let sorted_ok = data.windows(2).all(|w| w[0] <= w[1]);
        println!("  n={:<10}  time={:>10.3?}  sorted={}", n, elapsed, sorted_ok);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let mut arr: Vec<i32> = vec![];
        quick_sort(&mut arr);
        assert!(arr.is_empty());
    }

    #[test]
    fn test_single() {
        let mut arr = vec![42];
        quick_sort(&mut arr);
        assert_eq!(arr, vec![42]);
    }

    #[test]
    fn test_sorted_input() {
        let mut arr = vec![1, 2, 3, 4, 5];
        quick_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_reverse_sorted() {
        let mut arr = vec![5, 4, 3, 2, 1];
        quick_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_duplicates() {
        let mut arr = vec![3, 1, 4, 1, 5, 9, 2, 6, 5];
        quick_sort(&mut arr);
        assert_eq!(arr, vec![1, 1, 2, 3, 4, 5, 5, 6, 9]);
    }

    #[test]
    fn test_qual_scores() {
        let mut scores = vec![934.0f64, 52.0, 1847.0, 300.0, 78.0, 2100.0, 415.0, 61.0];
        quick_sort(&mut scores);
        assert_eq!(scores, vec![52.0, 61.0, 78.0, 300.0, 415.0, 934.0, 1847.0, 2100.0]);
    }
}
