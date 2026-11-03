//! Base Quality Score Recalibration (BQSR)
//!
//! Corrects systematic biases in base quality scores by comparing
//! reported vs. observed error rates at known non-variant sites.

use std::collections::HashMap;

/// A covariate key for recalibration.
#[derive(Debug, Hash, Eq, PartialEq, Clone)]
pub struct CovariateKey {
    pub reported_qual: u8,
    pub read_position: u8, // binned position in read
    pub dinucleotide: u16, // encoded 2-base context
}

/// Recalibration table entry.
#[derive(Debug, Default)]
pub struct RecalEntry {
    pub total_bases: u64,
    pub error_bases: u64,
}

impl RecalEntry {
    pub fn empirical_quality(&self) -> u8 {
        if self.total_bases == 0 { return 30; }
        let error_rate = (self.error_bases as f64 + 1.0) / (self.total_bases as f64 + 2.0);
        let q = -10.0 * error_rate.log10();
        q.round().min(60.0).max(2.0) as u8
    }
}

/// Build recalibration table from known non-variant sites.
pub fn build_recal_table(
    observations: &[(CovariateKey, bool)], // (key, is_error)
) -> HashMap<CovariateKey, RecalEntry> {
    let mut table: HashMap<CovariateKey, RecalEntry> = HashMap::new();
    for (key, is_error) in observations {
        let entry = table.entry(key.clone()).or_default();
        entry.total_bases += 1;
        if *is_error { entry.error_bases += 1; }
    }
    table
}

/// Apply recalibration to a quality score.
pub fn recalibrate(
    table: &HashMap<CovariateKey, RecalEntry>,
    key: &CovariateKey,
) -> u8 {
    table.get(key)
        .map(|entry| entry.empirical_quality())
        .unwrap_or(key.reported_qual)
}

fn main() {
    println!("=== BQSR: Base Quality Score Recalibration ===\n");

    // Simulate observations at known sites
    let mut observations = Vec::new();
    for q in [20u8, 25, 30, 35] {
        let expected_error = 10.0f64.powf(-(q as f64) / 10.0);
        // Simulate that Q20 bases actually have Q18 quality (slight miscalibration)
        let actual_error = expected_error * 1.5;
        let n = 10000;
        let n_errors = (n as f64 * actual_error) as u64;
        for i in 0..n {
            let key = CovariateKey { reported_qual: q, read_position: 50, dinucleotide: 0 };
            observations.push((key, (i as u64) < n_errors));
        }
    }

    let table = build_recal_table(&observations);

    println!("Recalibration results:");
    println!("{:<15} {:<15} {:<15}", "Reported Q", "Empirical Q", "Bases");
    for q in [20u8, 25, 30, 35] {
        let key = CovariateKey { reported_qual: q, read_position: 50, dinucleotide: 0 };
        let recal_q = recalibrate(&table, &key);
        let entry = table.get(&key).unwrap();
        println!("{:<15} {:<15} {:<15}", q, recal_q, entry.total_bases);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_quality() {
        let key = CovariateKey { reported_qual: 30, read_position: 0, dinucleotide: 0 };
        let obs: Vec<(CovariateKey, bool)> = (0..1000).map(|_| (key.clone(), false)).collect();
        let table = build_recal_table(&obs);
        let q = recalibrate(&table, &key);
        assert!(q >= 27); // should be very high
    }
}
