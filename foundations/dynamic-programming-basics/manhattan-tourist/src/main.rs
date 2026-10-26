//! Manhattan Tourist Problem — Foundation for Sequence Alignment

use std::time::Instant;

/// Solve the Manhattan Tourist Problem.
/// down: n rows × m cols of weights for vertical edges
/// right: (n+1) rows × (m-1) cols of weights for horizontal edges
/// Returns the maximum score path value.
pub fn manhattan_tourist(
    n: usize,
    m: usize,
    down: &[Vec<i32>],
    right: &[Vec<i32>],
) -> i32 {
    let mut s = vec![vec![0i32; m + 1]; n + 1];

    // Initialize first column (can only come from above)
    for i in 1..=n {
        s[i][0] = s[i - 1][0] + down[i - 1][0];
    }

    // Initialize first row (can only come from left)
    for j in 1..=m {
        s[0][j] = s[0][j - 1] + right[0][j - 1];
    }

    // Fill the rest
    for i in 1..=n {
        for j in 1..=m {
            let from_above = s[i - 1][j] + down[i - 1][j];
            let from_left = s[i][j - 1] + right[i][j - 1];
            s[i][j] = from_above.max(from_left);
        }
    }

    s[n][m]
}

fn main() {
    println!("=== Manhattan Tourist Problem ===\n");

    // Example grid (4 rows, 4 cols of intersections = 5×5 vertices)
    let down = vec![
        vec![3, 2, 4, 0, 1],
        vec![3, 2, 4, 2, 0],
        vec![0, 7, 3, 4, 4],
        vec![3, 3, 0, 2, 1],
    ];
    let right = vec![
        vec![1, 4, 5, 2],
        vec![0, 6, 5, 8],
        vec![4, 2, 5, 3],
        vec![1, 3, 2, 1],
        vec![3, 1, 1, 3],
    ];

    let n = 4; // number of rows of down edges
    let m = 4; // number of cols of right edges

    let start = Instant::now();
    let result = manhattan_tourist(n, m, &down, &right);
    let elapsed = start.elapsed();

    println!("Grid size: {}×{} vertices", n + 1, m + 1);
    println!("Maximum score path: {}", result);
    println!("Time: {:?}\n", elapsed);

    println!("This same algorithm, with different edge weights,");
    println!("becomes Needleman-Wunsch sequence alignment.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_2x2() {
        let down = vec![
            vec![1, 2],
            vec![3, 4],
        ];
        let right = vec![
            vec![5],
            vec![6],
            vec![7],
        ];
        let result = manhattan_tourist(2, 1, &down, &right);
        assert!(result > 0);
    }
}