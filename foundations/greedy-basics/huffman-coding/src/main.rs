//! Huffman Coding — FASTQ Quality String Compression

use std::collections::{BinaryHeap, HashMap};
use std::cmp::Reverse;

#[derive(Debug)]
enum HuffTree {
    Leaf { ch: u8, freq: u64 },
    Internal { freq: u64, left: Box<HuffTree>, right: Box<HuffTree> },
}

impl HuffTree {
    fn freq(&self) -> u64 {
        match self {
            HuffTree::Leaf { freq, .. } => *freq,
            HuffTree::Internal { freq, .. } => *freq,
        }
    }
}

impl PartialEq for HuffTree {
    fn eq(&self, other: &Self) -> bool { self.freq() == other.freq() }
}
impl Eq for HuffTree {}
impl PartialOrd for HuffTree {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HuffTree {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.freq().cmp(&other.freq())
    }
}

/// Build a Huffman tree from byte frequencies.
pub fn build_huffman_tree(data: &[u8]) -> Option<HuffTree> {
    let mut freq_map: HashMap<u8, u64> = HashMap::new();
    for &b in data {
        *freq_map.entry(b).or_insert(0) += 1;
    }

    if freq_map.is_empty() {
        return None;
    }

    let mut heap: BinaryHeap<Reverse<HuffTree>> = freq_map
        .into_iter()
        .map(|(ch, freq)| Reverse(HuffTree::Leaf { ch, freq }))
        .collect();

    while heap.len() > 1 {
        let Reverse(left) = heap.pop().unwrap();
        let Reverse(right) = heap.pop().unwrap();
        let combined_freq = left.freq() + right.freq();
        heap.push(Reverse(HuffTree::Internal {
            freq: combined_freq,
            left: Box::new(left),
            right: Box::new(right),
        }));
    }

    heap.pop().map(|Reverse(tree)| tree)
}

/// Extract codes from the Huffman tree.
pub fn extract_codes(tree: &HuffTree) -> HashMap<u8, String> {
    let mut codes = HashMap::new();
    fn walk(node: &HuffTree, prefix: String, codes: &mut HashMap<u8, String>) {
        match node {
            HuffTree::Leaf { ch, .. } => {
                codes.insert(*ch, if prefix.is_empty() { "0".into() } else { prefix });
            }
            HuffTree::Internal { left, right, .. } => {
                walk(left, format!("{}0", prefix), codes);
                walk(right, format!("{}1", prefix), codes);
            }
        }
    }
    walk(tree, String::new(), &mut codes);
    codes
}

fn main() {
    println!("=== Huffman Coding: FASTQ Quality String Compression ===\n");

    // Simulated quality string from a FASTQ file
    let quality_string = b"IIIIIIIIIIFFFFFFFFFFFFFFFFFFF!!!!!!###########JJJJJJJJ";

    println!("Quality string length: {} bytes", quality_string.len());
    println!("Original size: {} bits\n", quality_string.len() * 8);

    if let Some(tree) = build_huffman_tree(quality_string) {
        let codes = extract_codes(&tree);

        println!("Huffman codes:");
        let mut sorted_codes: Vec<_> = codes.iter().collect();
        sorted_codes.sort_by_key(|(_, code)| code.len());
        for (ch, code) in &sorted_codes {
            println!("  '{}' → {} ({} bits)", **ch as char, code, code.len());
        }

        let compressed_bits: usize = quality_string
            .iter()
            .map(|b| codes[b].len())
            .sum();

        println!("\nCompressed size: {} bits", compressed_bits);
        println!("Compression ratio: {:.1}%",
                 (1.0 - compressed_bits as f64 / (quality_string.len() * 8) as f64) * 100.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_char() {
        let tree = build_huffman_tree(b"AAAA").unwrap();
        let codes = extract_codes(&tree);
        assert_eq!(codes.len(), 1);
        assert_eq!(codes[&b'A'], "0");
    }

    #[test]
    fn test_two_chars() {
        let tree = build_huffman_tree(b"AABB").unwrap();
        let codes = extract_codes(&tree);
        assert_eq!(codes.len(), 2);
        for (_, code) in &codes {
            assert_eq!(code.len(), 1);
        }
    }

    #[test]
    fn test_prefix_free() {
        let tree = build_huffman_tree(b"AAABBC").unwrap();
        let codes = extract_codes(&tree);
        let code_list: Vec<&String> = codes.values().collect();
        for i in 0..code_list.len() {
            for j in 0..code_list.len() {
                if i != j {
                    assert!(!code_list[j].starts_with(code_list[i].as_str()),
                            "Prefix-free property violated");
                }
            }
        }
    }
}
