//! The Change Problem — Minimum Fragments for Oligo Assembly
//!
//! Demonstrates dynamic programming by finding the minimum number
//! of pre-made DNA fragments needed to build an oligo of a target length.

use std::time::Instant;

/// Compute the minimum number of coins (fragments) to make the target amount.
/// Returns None if the amount cannot be made exactly.
///
/// # Rust Performance Advantages
///
/// 1. **Flat Vec<u32>**: The DP table is a contiguous block of 4-byte values.
///    Python uses a list of 28+ byte int objects scattered across the heap.
///
/// 2. **Bounds-checked at compile time**: In release mode, the compiler
///    proves that all array accesses are within bounds and removes checks.
///
/// 3. **No allocator pressure**: We allocate the DP table once. Python's
///    list comprehension allocates each int object individually.
///
/// 4. **SIMD-friendly inner loop**: The min comparison over d denominations
///    compiles to tight branch-free code.
pub fn min_coins(denominations: &[u32], target: u32) -> Option<u32> {
    let target = target as usize;
    let inf = u32::MAX;
    let mut dp = vec![inf; target + 1];
    dp[0] = 0;

    for amount in 1..=target {
        for &coin in denominations {
            let coin = coin as usize;
            if coin <= amount && dp[amount - coin] != inf {
                dp[amount] = dp[amount].min(dp[amount - coin] + 1);
            }
        }
    }

    if dp[target] == inf { None } else { Some(dp[target]) }
}

/// Also recover which coins were used (traceback).
pub fn min_coins_with_traceback(denominations: &[u32], target: u32) -> Option<(u32, Vec<u32>)> {
    let target = target as usize;
    let inf = u32::MAX;
    let mut dp = vec![inf; target + 1];
    let mut from = vec![0u32; target + 1]; // which coin was used
    dp[0] = 0;

    for amount in 1..=target {
        for &coin in denominations {
            let c = coin as usize;
            if c <= amount && dp[amount - c] != inf && dp[amount - c] + 1 < dp[amount] {
                dp[amount] = dp[amount - c] + 1;
                from[amount] = coin;
            }
        }
    }

    if dp[target] == inf {
        return None;
    }

    // Traceback
    let mut coins_used = Vec::new();
    let mut remaining = target;
    while remaining > 0 {
        coins_used.push(from[remaining]);
        remaining -= from[remaining] as usize;
    }

    Some((dp[target], coins_used))
}

fn main() {
    println!("=== Change Problem: Minimum BAC Clones for Tiling Path ===\n");

    // BAC insert sizes in kilobases
    let bac_sizes: Vec<u32> = vec![50, 120, 200];
    let target_kb: u32 = 1000;

    println!("BAC insert sizes: {:?} kb", bac_sizes);
    println!("Target coverage:  {} kb\n", target_kb);

    match min_coins_with_traceback(&bac_sizes, target_kb) {
        Some((count, clones)) => {
            println!("Minimum BAC clones needed: {}", count);
            println!("Clone sizes used: {:?}", clones);
            println!("Verification: sum = {} kb", clones.iter().sum::<u32>());
        }
        None => println!("Cannot achieve exactly {} kb with these BAC sizes.", target_kb),
    }

    // Second example: oligo synthesis
    println!("\n=== Oligo Fragment Assembly ===\n");
    let fragment_sizes: Vec<u32> = vec![1, 3, 7];
    let target_bases: u32 = 77;

    match min_coins_with_traceback(&fragment_sizes, target_bases) {
        Some((count, fragments)) => {
            println!("Fragment sizes: {:?}", fragment_sizes);
            println!("Target length:  {} bases", target_bases);
            println!("Min fragments:  {}", count);
            println!("Fragments used: {:?}", fragments);
        }
        None => println!("Cannot build {} bases with these fragments.", target_bases),
    }

    // Benchmark
    println!("\n=== Benchmark ===");
    let denominations = vec![1u32, 5, 10, 25, 50, 100];
    for &target in &[100, 1_000, 10_000, 100_000] {
        let start = Instant::now();
        let result = min_coins(&denominations, target);
        let elapsed = start.elapsed();
        println!("  target={:<8}  min_coins={:<6}  time={:>10.3?}",
                 target, result.unwrap_or(0), elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        assert_eq!(min_coins(&[1, 3, 7], 7), Some(1));
    }

    #[test]
    fn test_zero() {
        assert_eq!(min_coins(&[1, 5, 10], 0), Some(0));
    }

    #[test]
    fn test_impossible() {
        assert_eq!(min_coins(&[3, 7], 5), None);
    }

    #[test]
    fn test_standard() {
        assert_eq!(min_coins(&[1, 3, 7], 77), Some(11)); // 11 sevens
    }

    #[test]
    fn test_bac_tiling() {
        assert_eq!(min_coins(&[50, 120, 200], 1000), Some(5)); // 5 × 200
    }
}
