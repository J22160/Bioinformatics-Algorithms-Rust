//! Greedy Interval Scheduling — Sequencing Machine Optimization
//!
//! Given a set of experiments with start and end times, find the maximum
//! number of non-overlapping experiments. This mirrors the exon chaining
//! problem in gene prediction.

use std::time::Instant;

/// A genomic interval representing either an experiment time slot
/// or a candidate exon.
#[derive(Debug, Clone, Copy)]
struct Interval {
    id: usize,
    start: u64,
    end: u64,
}

/// Greedy interval scheduling: select maximum non-overlapping intervals.
///
/// # Algorithm
/// 1. Sort intervals by end time
/// 2. Greedily select the interval with earliest end time that doesn't overlap
///
/// # Proof of Optimality
/// This greedy choice is optimal because choosing the earliest-finishing
/// interval can never reduce the number of subsequent selections compared
/// to any other choice. (Exchange argument proof.)
///
/// # Rust Performance Advantages
///
/// 1. **In-place sorting**: `sort_unstable_by_key` sorts without allocating
///    extra memory and is typically faster than stable sort.
///    Python's `sorted()` always creates a new list.
///
/// 2. **Stack-allocated structs**: Each `Interval` is 24 bytes stored inline
///    in the Vec. Python stores references to heap-allocated objects.
///
/// 3. **Iterator chains**: `filter` and `collect` compile to a single pass
///    with no intermediate allocations. Python list comprehensions create
///    intermediate list objects.
pub fn max_non_overlapping(intervals: &mut [Interval]) -> Vec<Interval> {
    if intervals.is_empty() {
        return vec![];
    }

    // Sort by end time — the key greedy insight
    // sort_unstable is faster than sort because it doesn't preserve
    // order of equal elements (which we don't need)
    intervals.sort_unstable_by_key(|iv| iv.end);

    let mut selected = vec![intervals[0]];
    let mut last_end = intervals[0].end;

    for &iv in &intervals[1..] {
        if iv.start >= last_end {
            selected.push(iv);
            last_end = iv.end;
        }
    }

    selected
}

/// Weighted interval scheduling using dynamic programming.
/// This is more general: each interval has a weight (e.g., alignment score),
/// and we maximize total weight. This corresponds to the Exon Chaining problem.
///
/// Time: O(n log n), Space: O(n)
pub fn weighted_interval_scheduling(intervals: &mut [(Interval, f64)]) -> (f64, Vec<usize>) {
    if intervals.is_empty() {
        return (0.0, vec![]);
    }

    let n = intervals.len();
    intervals.sort_unstable_by(|a, b| a.0.end.cmp(&b.0.end));

    // For each interval, find the latest non-overlapping predecessor
    // using binary search (O(log n) per interval)
    let predecessors: Vec<Option<usize>> = (0..n)
        .map(|i| {
            let start = intervals[i].0.start;
            // Binary search for the rightmost interval ending <= start
            let pos = intervals[..i]
                .partition_point(|iv| iv.0.end <= start);
            if pos > 0 { Some(pos - 1) } else { None }
        })
        .collect();

    // DP: dp[i] = max weight using intervals 0..=i
    let mut dp = vec![0.0f64; n];
    dp[0] = intervals[0].1;

    for i in 1..n {
        let include = intervals[i].1
            + predecessors[i].map_or(0.0, |p| dp[p]);
        let exclude = dp[i - 1];
        dp[i] = include.max(exclude);
    }

    // Traceback to find which intervals were selected
    let mut selected = vec![];
    let mut i = n - 1;
    loop {
        let include = intervals[i].1
            + predecessors[i].map_or(0.0, |p| dp[p]);
        if include >= dp[i.saturating_sub(1)] || i == 0 {
            selected.push(intervals[i].0.id);
            match predecessors[i] {
                Some(p) => i = p,
                None => break,
            }
        } else {
            if i == 0 { break; }
            i -= 1;
        }
    }

    selected.reverse();
    (dp[n - 1], selected)
}

fn main() {
    println!("=== Greedy Interval Scheduling: Sequencing Machine Optimization ===\n");

    let mut experiments = vec![
        Interval { id: 0, start: 0, end: 6 },    // Experiment A
        Interval { id: 1, start: 1, end: 4 },    // Experiment B
        Interval { id: 2, start: 3, end: 5 },    // Experiment C
        Interval { id: 3, start: 5, end: 7 },    // Experiment D
        Interval { id: 4, start: 3, end: 9 },    // Experiment E
        Interval { id: 5, start: 5, end: 9 },    // Experiment F
        Interval { id: 6, start: 6, end: 10 },   // Experiment G
        Interval { id: 7, start: 8, end: 11 },   // Experiment H
    ];

    let labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];

    println!("All experiments:");
    for iv in &experiments {
        println!("  {} [{}, {})", labels[iv.id], iv.start, iv.end);
    }

    let selected = max_non_overlapping(&mut experiments);

    println!("\nGreedy selection (max non-overlapping):");
    for iv in &selected {
        println!("  {} [{}, {})", labels[iv.id], iv.start, iv.end);
    }
    println!("Total: {} experiments\n", selected.len());

    // Weighted version (exon chaining analogy)
    println!("=== Weighted Interval Scheduling (Exon Chaining) ===\n");

    let mut weighted_exons: Vec<(Interval, f64)> = vec![
        (Interval { id: 0, start: 100, end: 250 }, 15.0),   // Exon 1: score 15
        (Interval { id: 1, start: 200, end: 400 }, 25.0),   // Exon 2: score 25
        (Interval { id: 2, start: 350, end: 500 }, 20.0),   // Exon 3: score 20
        (Interval { id: 3, start: 450, end: 700 }, 30.0),   // Exon 4: score 30
        (Interval { id: 4, start: 600, end: 800 }, 18.0),   // Exon 5: score 18
    ];

    println!("Candidate exons (genomic position, alignment score):");
    for (iv, w) in &weighted_exons {
        println!("  Exon {} [{}, {})  score={:.1}", iv.id + 1, iv.start, iv.end, w);
    }

    let (total_weight, selected_ids) = weighted_interval_scheduling(&mut weighted_exons);

    println!("\nOptimal exon chain:");
    for &id in &selected_ids {
        let (iv, w) = weighted_exons.iter().find(|(iv, _)| iv.id == id).unwrap();
        println!("  Exon {} [{}, {})  score={:.1}", iv.id + 1, iv.start, iv.end, w);
    }
    println!("Total score: {:.1}", total_weight);

    // Benchmark
    println!("\n=== Scaling Benchmark ===");
    for &n in &[1_000, 10_000, 100_000, 1_000_000] {
        let mut intervals: Vec<Interval> = (0..n)
            .map(|i| {
                let start = (i as u64 * 7 + 3) % (n as u64 * 10);
                Interval { id: i, start, end: start + (i as u64 % 50) + 1 }
            })
            .collect();

        let start = Instant::now();
        let result = max_non_overlapping(&mut intervals);
        let elapsed = start.elapsed();
        println!("  n={:<10}  selected={}  time={:>10.3?}", n, result.len(), elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let mut intervals: Vec<Interval> = vec![];
        assert!(max_non_overlapping(&mut intervals).is_empty());
    }

    #[test]
    fn test_no_overlap() {
        let mut intervals = vec![
            Interval { id: 0, start: 0, end: 1 },
            Interval { id: 1, start: 1, end: 2 },
            Interval { id: 2, start: 2, end: 3 },
        ];
        assert_eq!(max_non_overlapping(&mut intervals).len(), 3);
    }

    #[test]
    fn test_all_overlap() {
        let mut intervals = vec![
            Interval { id: 0, start: 0, end: 10 },
            Interval { id: 1, start: 1, end: 9 },
            Interval { id: 2, start: 2, end: 8 },
        ];
        assert_eq!(max_non_overlapping(&mut intervals).len(), 1);
    }

    #[test]
    fn test_classic_example() {
        let mut intervals = vec![
            Interval { id: 0, start: 0, end: 6 },
            Interval { id: 1, start: 1, end: 4 },
            Interval { id: 2, start: 3, end: 5 },
            Interval { id: 3, start: 5, end: 7 },
            Interval { id: 4, start: 3, end: 9 },
            Interval { id: 5, start: 5, end: 9 },
            Interval { id: 6, start: 6, end: 10 },
            Interval { id: 7, start: 8, end: 11 },
        ];
        let result = max_non_overlapping(&mut intervals);
        assert_eq!(result.len(), 3);
    }
}