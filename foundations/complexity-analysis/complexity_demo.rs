//! Complexity Analysis Demo
//!
//! Demonstrates O(n), O(n log n), and O(n^2) approaches to finding
//! the maximum gene expression value from an RNA-seq dataset.
//!
//! Run: rustc -O complexity_demo.rs && ./complexity_demo

use std::time::Instant;

/// O(n) — Single linear scan.
/// We visit each element exactly once, comparing it to our current best.
/// This mirrors how you might scan through microscopy images looking for
/// the brightest fluorescent signal.
fn find_max_linear(values: &[f64]) -> f64 {
    // `fold` iterates through the slice, carrying forward the running max.
    // f64::NEG_INFINITY ensures any real value will be larger.
    values.iter().fold(f64::NEG_INFINITY, |max, &v| {
        if v > max { v } else { max }
    })
}

/// O(n log n) — Sort and take the last element.
/// Sorting is more work than necessary here, but it demonstrates
/// the n·log(n) complexity class used by samtools sort and most
/// BAM processing pipelines.
fn find_max_sort(values: &mut [f64]) -> f64 {
    // Rust's sort is an adaptive merge sort (TimSort variant) — O(n log n)
    // guaranteed worst case, often faster on partially sorted data.
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    *values.last().unwrap()
}

/// O(n²) — Compare every pair.
/// This is deliberately wasteful. We check if each element is greater than
/// ALL others. This mirrors the brute-force motif search where you compare
/// every possible starting position against every sequence.
fn find_max_quadratic(values: &[f64]) -> f64 {
    let n = values.len();
    for i in 0..n {
        let mut is_maximum = true;
        for j in 0..n {
            if values[j] > values[i] {
                is_maximum = false;
                break; // Small optimization: still O(n²) worst case
            }
        }
        if is_maximum {
            return values[i];
        }
    }
    values[0] // Fallback (unreachable for non-empty input)
}

/// Generate synthetic gene expression values using a simple
/// linear congruential generator (avoids external dependencies).
fn generate_expression_data(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            // LCG parameters from Numerical Recipes
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            // Map to a TPM-like range [0, 100000)
            (state >> 33) as f64 / (1u64 << 31) as f64 * 100_000.0
        })
        .collect()
}

fn benchmark<F: Fn() -> f64>(name: &str, f: F) -> f64 {
    let start = Instant::now();
    let result = f();
    let elapsed = start.elapsed();
    println!("  {:<25} result={:.2}  time={:>10.3?}", name, result, elapsed);
    result
}

fn main() {
    println!("=== Complexity Analysis: Finding Maximum Gene Expression ===\n");

    // Test with increasing sizes to observe scaling behavior
    for &n in &[1_000, 10_000, 100_000, 1_000_000] {
        println!("--- n = {} ---", n);
        let data = generate_expression_data(n, 42);

        // O(n) linear scan
        let r1 = benchmark("O(n) linear scan", || find_max_linear(&data));

        // O(n log n) sort-based
        let mut data_copy = data.clone();
        let r2 = benchmark("O(n log n) sort", || find_max_sort(&mut data_copy));

        // O(n²) only for small sizes — would be too slow for large n
        if n <= 100_000 {
            let r3 = benchmark("O(n²) pairwise", || find_max_quadratic(&data));
            assert!((r1 - r3).abs() < 1e-10);
        } else {
            println!("  {:<25} SKIPPED (too slow)", "O(n²) pairwise");
        }

        assert!((r1 - r2).abs() < 1e-10);
        println!();
    }

    println!("=== Key Observations ===");
    println!("1. O(n) time roughly doubles when n doubles");
    println!("2. O(n log n) grows slightly faster than doubling");
    println!("3. O(n²) grows ~4x when n doubles — impractical for large datasets");
    println!("\nIn bioinformatics: always prefer O(n) or O(n log n) algorithms.");
    println!("This is why BWA uses FM-index search (O(m)) instead of naive O(nm) matching.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_methods_agree() {
        let data = generate_expression_data(500, 123);
        let mut sorted = data.clone();
        let r1 = find_max_linear(&data);
        let r2 = find_max_sort(&mut sorted);
        let r3 = find_max_quadratic(&data);
        assert!((r1 - r2).abs() < 1e-10);
        assert!((r1 - r3).abs() < 1e-10);
    }

    #[test]
    fn test_single_element() {
        assert_eq!(find_max_linear(&[42.0]), 42.0);
    }

    #[test]
    fn test_negative_values() {
        let data = vec![-5.0, -3.0, -10.0, -1.0];
        assert_eq!(find_max_linear(&data), -1.0);
    }
}
