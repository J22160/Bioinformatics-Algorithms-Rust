//! Pseudoalignment — Fast Read-Transcript Compatibility
//!
//! Determines which transcripts a read is compatible with
//! using k-mer matching, without full base-level alignment.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

/// K-mer index mapping each k-mer to its set of compatible transcripts.
pub struct KmerIndex {
    index: HashMap<Vec<u8>, Vec<usize>>,
    k: usize,
}

impl KmerIndex {
    /// Build a k-mer index from a set of transcript sequences.
    pub fn build(transcripts: &[Vec<u8>], k: usize) -> Self {
        let mut index: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();

        for (tid, seq) in transcripts.iter().enumerate() {
            if seq.len() < k { continue; }
            for i in 0..=(seq.len() - k) {
                let kmer = seq[i..i + k].to_vec();
                index.entry(kmer).or_default().push(tid);
            }
        }

        // Deduplicate transcript lists
        for tids in index.values_mut() {
            tids.sort_unstable();
            tids.dedup();
        }

        KmerIndex { index, k }
    }

    /// Pseudoalign a read: find the set of compatible transcripts.
    /// Returns the equivalence class (sorted transcript IDs).
    pub fn pseudoalign(&self, read: &[u8]) -> Option<Vec<usize>> {
        if read.len() < self.k {
            return None;
        }

        let mut compat: Option<HashSet<usize>> = None;

        for i in 0..=(read.len() - self.k) {
            let kmer = &read[i..i + self.k];
            if let Some(tids) = self.index.get(kmer) {
                let tid_set: HashSet<usize> = tids.iter().copied().collect();
                compat = Some(match compat {
                    Some(existing) => existing.intersection(&tid_set).copied().collect(),
                    None => tid_set,
                });
                if let Some(ref c) = compat {
                    if c.is_empty() { return None; }
                }
            }
        }

        compat.map(|set| {
            let mut v: Vec<usize> = set.into_iter().collect();
            v.sort_unstable();
            v
        })
    }
}

/// An equivalence class: a unique set of compatible transcripts + count.
#[derive(Debug, Clone)]
pub struct EquivalenceClass {
    pub transcript_ids: Vec<usize>,
    pub count: usize,
}

/// Collect reads into equivalence classes for efficient quantification.
pub fn collect_equivalence_classes(
    index: &KmerIndex,
    reads: &[Vec<u8>],
) -> Vec<EquivalenceClass> {
    let mut class_counts: HashMap<Vec<usize>, usize> = HashMap::new();

    for read in reads {
        if let Some(ec) = index.pseudoalign(read) {
            *class_counts.entry(ec).or_insert(0) += 1;
        }
    }

    class_counts.into_iter()
        .map(|(tids, count)| EquivalenceClass { transcript_ids: tids, count })
        .collect()
}

fn main() {
    println!("=== Pseudoalignment: Fast Read Classification ===\n");

    // Define a small transcriptome
    let transcript_names = vec!["GAPDH-001", "GAPDH-002", "ACTB-001", "TP53-001"];
    let transcripts: Vec<Vec<u8>> = vec![
        b"ATGCTAGCTAGCATCGATCGATCGATCGATCG".to_vec(),
        b"ATGCTAGCTAGCATCGATCG".to_vec(),              // prefix of first
        b"CCGGAATTCCGGAATTCCGG".to_vec(),
        b"TTGGCCAATTGGCCAATTGG".to_vec(),
    ];

    let k = 7;
    println!("Building k-mer index (k={})...", k);
    let start = Instant::now();
    let index = KmerIndex::build(&transcripts, k);
    println!("Index built in {:?} ({} unique k-mers)\n",
             start.elapsed(), index.index.len());

    // Simulate reads
    let reads: Vec<Vec<u8>> = vec![
        b"CTAGCTAGCATCGATCG".to_vec(),      // maps to GAPDH-001 and GAPDH-002
        b"GAATTCCGGAATTCCGG".to_vec(),       // maps to ACTB-001
        b"CCAATTGGCCAATTGG".to_vec(),        // maps to TP53-001
        b"ATCGATCGATCGATCG".to_vec(),        // unique to GAPDH-001
        b"CTAGCTAGCATCGATCG".to_vec(),       // duplicate of first
    ];

    println!("Pseudoaligning {} reads...\n", reads.len());
    for (i, read) in reads.iter().enumerate() {
        let ec = index.pseudoalign(read);
        let names: Vec<&str> = ec.as_ref()
            .map(|ids| ids.iter().map(|&id| transcript_names[id]).collect())
            .unwrap_or_default();
        println!("  Read {}: {} -> {:?}",
                 i, std::str::from_utf8(read).unwrap(), names);
    }

    // Equivalence classes
    let classes = collect_equivalence_classes(&index, &reads);
    println!("\nEquivalence classes ({} distinct):", classes.len());
    for ec in &classes {
        let names: Vec<&str> = ec.transcript_ids.iter()
            .map(|&id| transcript_names[id]).collect();
        println!("  {:?} × {}", names, ec.count);
    }

    println!("\n=== Key Insight ===");
    println!("{} reads compressed into {} equivalence classes", reads.len(), classes.len());
    println!("EM runs on classes (not reads) → massive speedup for millions of reads");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unique_mapping() {
        let transcripts = vec![b"AAACCCGGG".to_vec(), b"TTTGGGAAA".to_vec()];
        let index = KmerIndex::build(&transcripts, 5);
        let ec = index.pseudoalign(b"AAACCCGGG");
        assert_eq!(ec, Some(vec![0]));
    }

    #[test]
    fn test_multi_mapping() {
        let transcripts = vec![b"AAACCC".to_vec(), b"AAACCCGGG".to_vec()];
        let index = KmerIndex::build(&transcripts, 3);
        let ec = index.pseudoalign(b"AAACCC");
        assert!(ec.unwrap().len() == 2);
    }
}
