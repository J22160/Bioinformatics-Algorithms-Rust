//! Shortest Superstring — Greedy Genome Assembly Heuristic

pub fn overlap_len(a: &str, b: &str) -> usize {
    let min_len = a.len().min(b.len());
    for k in (1..=min_len).rev() {
        if a.ends_with(&b[..k]) {
            return k;
        }
    }
    0
}

pub fn greedy_superstring(mut strings: Vec<String>) -> String {
    while strings.len() > 1 {
        let n = strings.len();
        let mut best_i = 0;
        let mut best_j = 1;
        let mut best_overlap = 0;

        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let ol = overlap_len(&strings[i], &strings[j]);
                    if ol > best_overlap {
                        best_overlap = ol;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        // Merge best pair
        let merged = format!("{}{}", &strings[best_i], &strings[best_j][best_overlap..]);

        // Remove in correct order (larger index first)
        let (first, second) = if best_i > best_j {
            (best_i, best_j)
        } else {
            (best_j, best_i)
        };
        strings.remove(first);
        strings.remove(second);
        strings.push(merged);
    }

    strings.pop().unwrap_or_default()
}

fn main() {
    println!("=== Shortest Superstring: Greedy Assembly ===\n");

    let reads = vec![
        "ATGCCA".to_string(), "GCCATG".to_string(),
        "CATGCC".to_string(), "TGCCAT".to_string(),
    ];

    println!("Reads: {:?}", reads);
    let concat_len: usize = reads.iter().map(|r| r.len()).sum();
    println!("Concatenation length: {}", concat_len);

    let superstring = greedy_superstring(reads.clone());
    println!("Superstring: {}", superstring);
    println!("Length: {} (saved {} chars)", superstring.len(),
             concat_len - superstring.len());

    // Verify
    for read in &reads {
        assert!(superstring.contains(read.as_str()),
                "Missing read {}", read);
    }
    println!("✓ All reads present in superstring");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_greedy() {
        let reads = vec!["ABC".to_string(), "BCD".to_string(), "CDE".to_string()];
        let result = greedy_superstring(reads);
        assert!(result.contains("ABC"));
        assert!(result.contains("BCD"));
        assert!(result.contains("CDE"));
        assert_eq!(result.len(), 5); // ABCDE
    }
}
