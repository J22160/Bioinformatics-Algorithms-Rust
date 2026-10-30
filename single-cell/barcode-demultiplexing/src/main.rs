//! Cell Barcode Demultiplexing for scRNA-seq
//!
//! Assigns sequencing reads to cells by matching observed barcodes
//! against a whitelist, with single-error correction.
//! This mirrors the approach used by Cell Ranger (10x Genomics).

use std::collections::{HashMap, HashSet};
use std::time::Instant;

/// Result of attempting to match a barcode against the whitelist.
#[derive(Debug, Clone, PartialEq)]
pub enum BarcodeMatch {
    /// Exact match found in whitelist
    Exact(String),
    /// Corrected to a whitelist barcode (original, corrected, quality)
    Corrected {
        observed: String,
        corrected: String,
        mismatch_pos: usize,
        posterior_prob: f64,
    },
    /// Ambiguous: multiple whitelist barcodes within distance 1
    Ambiguous(String),
    /// No match found within distance 1
    NoMatch(String),
}

/// A barcode demultiplexer using hash-based exact and 1-error matching.
///
/// # Rust Performance Advantages
///
/// 1. **HashSet for O(1) lookup**: The whitelist is stored in a
///    HashSet<Vec<u8>>, giving constant-time exact matching.
///    Python's set() has similar asymptotic complexity but higher
///    constant factors due to object overhead.
///
/// 2. **Byte-level barcode representation**: Barcodes are stored
///    as Vec<u8> (1 byte per base) rather than String objects.
///    This halves memory usage vs Python's str (2+ bytes per char
///    due to internal encoding).
///
/// 3. **In-place mutation for variant generation**: When checking
///    all 1-substitution variants, we mutate a single buffer in-place
///    rather than creating new strings. This avoids 48 heap allocations
///    per barcode (16 positions × 3 alternatives).
///
/// 4. **Cache-friendly sequential access**: Processing reads in order
///    through a contiguous Vec gives excellent CPU cache utilization.
pub struct BarcodeDemuxer {
    whitelist: HashSet<Vec<u8>>,
    barcode_counts: HashMap<Vec<u8>, u32>,
    barcode_len: usize,
}

impl BarcodeDemuxer {
    /// Build a demultiplexer from a whitelist of valid barcodes.
    pub fn new(whitelist_barcodes: &[&str]) -> Self {
        let barcode_len = whitelist_barcodes
            .first()
            .map(|b| b.len())
            .unwrap_or(16);

        let whitelist: HashSet<Vec<u8>> = whitelist_barcodes
            .iter()
            .map(|b| b.as_bytes().to_vec())
            .collect();

        BarcodeDemuxer {
            whitelist,
            barcode_counts: HashMap::new(),
            barcode_len,
        }
    }

    /// Attempt to match an observed barcode against the whitelist.
    pub fn match_barcode(&mut self, observed: &[u8], qualities: &[u8]) -> BarcodeMatch {
        // Step 1: Exact match (O(1) via HashSet)
        if self.whitelist.contains(observed) {
            *self.barcode_counts
                .entry(observed.to_vec())
                .or_insert(0) += 1;
            return BarcodeMatch::Exact(
                String::from_utf8_lossy(observed).to_string(),
            );
        }

        // Step 2: 1-error correction
        // Generate all possible 1-substitution variants and check whitelist
        let bases = [b'A', b'C', b'G', b'T'];
        let mut candidates: Vec<(Vec<u8>, usize)> = Vec::new();
        let mut buffer = observed.to_vec();

        for pos in 0..observed.len() {
            let original_base = buffer[pos];
            for &alt_base in &bases {
                if alt_base == original_base {
                    continue;
                }
                buffer[pos] = alt_base;
                if self.whitelist.contains(&buffer) {
                    candidates.push((buffer.clone(), pos));
                }
            }
            buffer[pos] = original_base; // restore
        }

        match candidates.len() {
            0 => BarcodeMatch::NoMatch(
                String::from_utf8_lossy(observed).to_string(),
            ),
            1 => {
                let (corrected, mismatch_pos) = &candidates[0];
                // Calculate posterior probability using base quality
                let qual = qualities.get(*mismatch_pos).copied().unwrap_or(30);
                let error_prob = 10.0f64.powf(-(qual as f64) / 10.0);
                // Simplified Bayesian: P(correct | error at this position)
                let posterior = 1.0 - error_prob;

                *self.barcode_counts
                    .entry(corrected.clone())
                    .or_insert(0) += 1;

                BarcodeMatch::Corrected {
                    observed: String::from_utf8_lossy(observed).to_string(),
                    corrected: String::from_utf8_lossy(corrected).to_string(),
                    mismatch_pos: *mismatch_pos,
                    posterior_prob: posterior,
                }
            }
            _ => BarcodeMatch::Ambiguous(
                String::from_utf8_lossy(observed).to_string(),
            ),
        }
    }

    /// Identify real cells using the "knee" in the UMI count distribution.
    /// Returns barcodes that pass the threshold.
    pub fn call_cells(&self, expected_cells: usize) -> Vec<(String, u32)> {
        let mut counts: Vec<(Vec<u8>, u32)> = self
            .barcode_counts
            .iter()
            .map(|(k, &v)| (k.clone(), v))
            .collect();

        // Sort by count descending
        counts.sort_by(|a, b| b.1.cmp(&a.1));

        // Simple knee detection: use expected_cells as initial estimate
        // In practice, Cell Ranger uses an adaptive algorithm
        let threshold = if counts.len() > expected_cells {
            // Use 1/10th of the count at the expected_cells rank
            counts[expected_cells.min(counts.len() - 1)].1 / 10
        } else {
            1
        };

        counts
            .into_iter()
            .filter(|(_, count)| *count > threshold)
            .map(|(bc, count)| (String::from_utf8_lossy(&bc).to_string(), count))
            .collect()
    }
}

fn main() {
    println!("=== Barcode Demultiplexing: scRNA-seq Cell Assignment ===\n");

    // Create a small whitelist for demonstration
    let whitelist = vec![
        "ACGTACGTACGTACGT",
        "TGCATGCATGCATGCA",
        "GGCCAATTGGCCAATT",
        "CCTTGGAACCTTGGAA",
        "AACCGGTTAACCGGTT",
    ];

    let mut demuxer = BarcodeDemuxer::new(&whitelist);

    // Simulate reads with various barcode qualities
    let test_reads: Vec<(&[u8], &[u8])> = vec![
        // Exact matches
        (b"ACGTACGTACGTACGT", &[30; 16]),
        (b"TGCATGCATGCATGCA", &[30; 16]),
        // 1-error (correctable)
        (b"ACGTACGTTCGTACGT", &[30, 30, 30, 30, 30, 30, 30, 30, 10, 30, 30, 30, 30, 30, 30, 30]),
        // No match
        (b"AAAAAAAAAAAAAAAA", &[30; 16]),
        // More exact matches to build counts
        (b"ACGTACGTACGTACGT", &[30; 16]),
        (b"ACGTACGTACGTACGT", &[30; 16]),
        (b"TGCATGCATGCATGCA", &[30; 16]),
    ];

    println!("Processing {} reads...\n", test_reads.len());

    let start = Instant::now();
    let mut exact = 0;
    let mut corrected = 0;
    let mut ambiguous = 0;
    let mut no_match = 0;

    for (barcode, quals) in &test_reads {
        let result = demuxer.match_barcode(barcode, quals);
        match &result {
            BarcodeMatch::Exact(bc) => {
                println!("  EXACT:     {} -> {}", std::str::from_utf8(barcode).unwrap(), bc);
                exact += 1;
            }
            BarcodeMatch::Corrected { observed, corrected: cor, mismatch_pos, posterior_prob } => {
                println!("  CORRECTED: {} -> {} (pos {}, P={:.4})",
                         observed, cor, mismatch_pos, posterior_prob);
                corrected += 1;
            }
            BarcodeMatch::Ambiguous(bc) => {
                println!("  AMBIGUOUS: {}", bc);
                ambiguous += 1;
            }
            BarcodeMatch::NoMatch(bc) => {
                println!("  NO MATCH:  {}", bc);
                no_match += 1;
            }
        }
    }
    let elapsed = start.elapsed();

    println!("\nSummary:");
    println!("  Exact:     {} ({:.1}%)", exact, exact as f64 / test_reads.len() as f64 * 100.0);
    println!("  Corrected: {} ({:.1}%)", corrected, corrected as f64 / test_reads.len() as f64 * 100.0);
    println!("  Ambiguous: {} ({:.1}%)", ambiguous, ambiguous as f64 / test_reads.len() as f64 * 100.0);
    println!("  No match:  {} ({:.1}%)", no_match, no_match as f64 / test_reads.len() as f64 * 100.0);
    println!("  Time: {:?}", elapsed);

    // Cell calling
    let cells = demuxer.call_cells(3);
    println!("\nCalled cells ({}):", cells.len());
    for (bc, count) in &cells {
        println!("  {} : {} UMIs", bc, count);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let mut d = BarcodeDemuxer::new(&["ACGT"]);
        let r = d.match_barcode(b"ACGT", &[30; 4]);
        assert!(matches!(r, BarcodeMatch::Exact(_)));
    }

    #[test]
    fn test_one_error_correction() {
        let mut d = BarcodeDemuxer::new(&["ACGT"]);
        let r = d.match_barcode(b"ACTT", &[30; 4]);
        assert!(matches!(r, BarcodeMatch::Corrected { .. }));
    }

    #[test]
    fn test_no_match() {
        let mut d = BarcodeDemuxer::new(&["ACGT"]);
        let r = d.match_barcode(b"TTTT", &[30; 4]);
        assert!(matches!(r, BarcodeMatch::NoMatch(_)));
    }
}
