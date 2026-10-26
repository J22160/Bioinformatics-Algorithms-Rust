//! Merge Sort — Sorting Genomic Coordinates
//!
//! A divide-and-conquer sorting algorithm with guaranteed O(n log n)
//! performance. This is the basis of `samtools sort` and similar tools.

use std::time::Instant;

/// Merge sort implementation.
///
/// # Rust Performance Advantages
///
/// 1. **Allocation reuse**: We allocate one temporary buffer up front and pass
///    it through all recursive calls, avoiding repeated heap allocation.
///    Python's `sorted()` creates a new list every time.
///
/// 2. **Monomorphization**: Rust generates specialized machine code for each
///    concrete type (u64, f64, struct). Python uses generic object comparisons
///    with virtual dispatch overhead.
///
/// 3. **No GC pressure**: The temporary buffer is deterministically freed.
///    Python's garbage collector must track every intermediate list.
///
/// 4. **Copy semantics**: For primitive types like u64, copying is a single
///    MOV instruction. Python copies reference-counted pointers (slower).

pub fn merge_sort<T: Copy + PartialOrd>(arr: &mut [T]) {
    let n = arr.len();
    if n <= 1 {
        return;
    }
    // Allocate workspace once — this is the key optimization over naive versions
    let mut buffer = arr.to_vec();
    merge_sort_recursive(arr, &mut buffer, 0, n);
}

fn merge_sort_recursive<T: Copy + PartialOrd>(
    arr: &mut [T],
    buffer: &mut [T],
    left: usize,
    right: usize,
) {
    if right - left <= 1 {
        return;
    }

    // For small subarrays, insertion sort is faster due to lower overhead
    if right - left <= 32 {
        insertion_sort(&mut arr[left..right]);
        return;
    }

    let mid = left + (right - left) / 2;

    // Recursively sort both halves
    merge_sort_recursive(arr, buffer, left, mid);
    merge_sort_recursive(arr, buffer, mid, right);

    // Merge the sorted halves
    merge(arr, buffer, left, mid, right);
}

fn merge<T: Copy + PartialOrd>(
    arr: &mut [T],
    buffer: &mut [T],
    left: usize,
    mid: usize,
    right: usize,
) {
    // Copy the region to merge into the buffer
    buffer[left..right].copy_from_slice(&arr[left..right]);

    let mut i = left;    // Pointer into left half of buffer
    let mut j = mid;     // Pointer into right half of buffer
    let mut k = left;    // Write pointer into arr

    while i < mid && j < right {
        if buffer[i] <= buffer[j] {
            arr[k] = buffer[i];
            i += 1;
        } else {
            arr[k] = buffer[j];
            j += 1;
        }
        k += 1;
    }

    // Copy remaining elements from whichever half isn't exhausted
    while i < mid {
        arr[k] = buffer[i];
        i += 1;
        k += 1;
    }
    while j < right {
        arr[k] = buffer[j];
        j += 1;
        k += 1;
    }
}

/// Insertion sort for small subarrays — lower constant factor than merge sort
fn insertion_sort<T: PartialOrd>(arr: &mut [T]) {
    for i in 1..arr.len() {
        let mut j = i;
        while j > 0 && arr[j - 1] > arr[j] {
            arr.swap(j - 1, j);
            j -= 1;
        }
    }
}

/// Represents a genomic interval (e.g., from a BED file)
#[derive(Debug, Clone, Copy)]
struct GenomicInterval {
    chrom_id: u8,   // Chromosome encoded as number (saves memory vs String)
    start: u64,
    end: u64,
}

impl PartialEq for GenomicInterval {
    fn eq(&self, other: &Self) -> bool {
        self.chrom_id == other.chrom_id && self.start == other.start
    }
}

impl PartialOrd for GenomicInterval {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.chrom_id.cmp(&other.chrom_id) {
            std::cmp::Ordering::Equal => Some(self.start.cmp(&other.start)),
            ord => Some(ord),
        }
    }
}

fn main() {
    println!("=== Merge Sort: Sorting Genomic Coordinates ===\n");

    // ChIP-seq binding site positions on chromosome 7
    let mut positions: Vec<u64> = vec![
        155_799, 5_529, 87_160, 24_601, 117_482, 48_930, 72_215,
        3_100, 200_450, 99_999, 145_000, 67_800, 31_200, 180_300,
    ];

    println!("Unsorted binding sites (chr7): {:?}\n", positions);

    let start = Instant::now();
    merge_sort(&mut positions);
    let elapsed = start.elapsed();

    println!("Sorted binding sites (chr7):   {:?}", positions);
    println!("Time: {:?}\n", elapsed);

    // Scaling benchmark
    println!("=== Scaling Benchmark ===");
    for &n in &[1_000, 10_000, 100_000, 1_000_000, 10_000_000] {
        let mut data: Vec<u64> = (0..n)
            .map(|i| ((i as u64).wrapping_mul(2654435761) % 3_000_000_000))
            .collect();
        let start = Instant::now();
        merge_sort(&mut data);
        let elapsed = start.elapsed();
        let sorted_ok = data.windows(2).all(|w| w[0] <= w[1]);
        println!("  n={:<10}  time={:>12.3?}  sorted={}", n, elapsed, sorted_ok);
    }

    println!("\nNotice: 10x more data ≈ 10-12x more time → O(n log n) confirmed.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let mut arr: Vec<u64> = vec![];
        merge_sort(&mut arr);
        assert!(arr.is_empty());
    }

    #[test]
    fn test_sorted() {
        let mut arr = vec![1u64, 2, 3, 4, 5];
        merge_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_reverse() {
        let mut arr = vec![5u64, 4, 3, 2, 1];
        merge_sort(&mut arr);
        assert_eq!(arr, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_genomic_positions() {
        let mut positions = vec![155_799u64, 5_529, 87_160, 24_601, 117_482, 48_930, 72_215];
        merge_sort(&mut positions);
        assert_eq!(positions, vec![5_529, 24_601, 48_930, 72_215, 87_160, 117_482, 155_799]);
    }

    #[test]
    fn test_large_random() {
        let mut data: Vec<u64> = (0..10_000)
            .map(|i| (i as u64).wrapping_mul(6364136223846793005) % 1_000_000)
            .collect();
        merge_sort(&mut data);
        assert!(data.windows(2).all(|w| w[0] <= w[1]));
    }
}