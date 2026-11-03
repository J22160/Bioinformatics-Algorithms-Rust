//! VCF Quality Filtering — Hard Filtering and FDR Correction

/// A variant call with associated quality metrics.
#[derive(Debug, Clone)]
pub struct VariantCall {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub qual: f64,
    pub depth: u32,
    pub quality_by_depth: f64,
    pub fisher_strand: f64,
    pub mapping_quality: f64,
}

/// Hard filtering thresholds (GATK best practices).
#[derive(Debug)]
pub struct HardFilterParams {
    pub min_qual: f64,
    pub min_depth: u32,
    pub min_qd: f64,
    pub max_fs: f64,
    pub min_mq: f64,
}

impl Default for HardFilterParams {
    fn default() -> Self {
        HardFilterParams {
            min_qual: 30.0,
            min_depth: 10,
            min_qd: 2.0,
            max_fs: 60.0,
            min_mq: 40.0,
        }
    }
}

/// Apply hard filters to a set of variant calls.
pub fn hard_filter(
    variants: &[VariantCall],
    params: &HardFilterParams,
) -> Vec<(usize, bool, Vec<String>)> {
    variants.iter().enumerate().map(|(idx, v)| {
        let mut reasons = Vec::new();
        if v.qual < params.min_qual { reasons.push("LowQual".into()); }
        if v.depth < params.min_depth { reasons.push("LowDepth".into()); }
        if v.quality_by_depth < params.min_qd { reasons.push("LowQD".into()); }
        if v.fisher_strand > params.max_fs { reasons.push("StrandBias".into()); }
        if v.mapping_quality < params.min_mq { reasons.push("LowMQ".into()); }
        let pass = reasons.is_empty();
        (idx, pass, reasons)
    }).collect()
}

/// Benjamini-Hochberg FDR correction on p-values.
pub fn bh_adjust(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    let mut indexed: Vec<(usize, f64)> = p_values.iter()
        .enumerate().map(|(i, &p)| (i, p)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut adjusted = vec![1.0; n];
    let mut min_so_far = 1.0f64;
    for rank in (0..n).rev() {
        let (orig_idx, p) = indexed[rank];
        let adj = (p * n as f64 / (rank + 1) as f64).min(1.0);
        min_so_far = min_so_far.min(adj);
        adjusted[orig_idx] = min_so_far;
    }
    adjusted
}

fn main() {
    println!("=== VCF Filtering: Hard Filters + FDR ===\n");

    let variants = vec![
        VariantCall { chrom: "chr1".into(), pos: 1000, ref_allele: "A".into(),
            alt_allele: "G".into(), qual: 500.0, depth: 45, quality_by_depth: 11.1,
            fisher_strand: 2.3, mapping_quality: 60.0 },
        VariantCall { chrom: "chr1".into(), pos: 2000, ref_allele: "C".into(),
            alt_allele: "T".into(), qual: 15.0, depth: 5, quality_by_depth: 3.0,
            fisher_strand: 1.0, mapping_quality: 55.0 },
        VariantCall { chrom: "chr1".into(), pos: 3000, ref_allele: "G".into(),
            alt_allele: "A".into(), qual: 200.0, depth: 30, quality_by_depth: 6.7,
            fisher_strand: 80.0, mapping_quality: 50.0 },
    ];

    let params = HardFilterParams::default();
    let results = hard_filter(&variants, &params);

    for (idx, pass, reasons) in &results {
        let v = &variants[*idx];
        let status = if *pass { "PASS" } else { &reasons.join(",") };
        println!("  {}:{} {}->{} QUAL={:.0} DP={} → {}",
                 v.chrom, v.pos, v.ref_allele, v.alt_allele,
                 v.qual, v.depth, status);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bh_correction() {
        let p = vec![0.01, 0.04, 0.03, 0.20];
        let adj = bh_adjust(&p);
        assert!(adj[0] <= adj[2]);
        for &a in &adj { assert!(a <= 1.0); }
    }
}
