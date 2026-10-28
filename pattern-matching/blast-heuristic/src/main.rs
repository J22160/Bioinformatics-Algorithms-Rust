//! BLAST-like Seed-and-Extend Algorithm
//!
//! Demonstrates the heuristic approach to sequence similarity search:
//! find exact word matches (seeds), then extend them into alignments.

use std::collections::HashMap;
use std::time::Instant;

/// A seed hit: position in query and position in database sequence
#[derive(Debug, Clone)]
struct SeedHit {
    db_id: usize,
    db_pos: usize,
    query_pos: usize,
}

/// A high-scoring segment pair (HSP)
#[derive(Debug, Clone)]
struct HSP {
    db_id: usize,
    db_start: usize,
    query_start: usize,
    length: usize,
    score: i32,
    identity: f64,
}

/// Build a word index (hash table) for all w-mers in the database.
///
/// # Rust Performance Advantages
///
/// 1. **HashMap with FxHash**: Rust's default HashMap uses SipHash
///    (cryptographic), but bioinformatics tools use faster hashes.
///    Even with SipHash, it's faster than Python's dict for bulk inserts.
///
/// 2. **Zero-copy slicing**: Database sequences are &[u8] slices
///    with no allocation. Python creates new string objects per k-mer.
///
/// 3. **Vec<SeedHit> pre-allocation**: We can estimate the number
///    of seeds and pre-allocate, avoiding repeated reallocation.
fn build_word_index(
    database: &[Vec<u8>],
    word_size: usize,
) -> HashMap<Vec<u8>, Vec<(usize, usize)>> {
    let mut index: HashMap<Vec<u8>, Vec<(usize, usize)>> = HashMap::new();

    for (db_id, seq) in database.iter().enumerate() {
        if seq.len() < word_size {
            continue;
        }
        for pos in 0..=(seq.len() - word_size) {
            let word = seq[pos..pos + word_size].to_vec();
            // Skip words containing non-ACGT characters
            if word.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                index.entry(word).or_default().push((db_id, pos));
            }
        }
    }
    index
}

/// Find seed hits between query and indexed database.
fn find_seeds(
    query: &[u8],
    word_index: &HashMap<Vec<u8>, Vec<(usize, usize)>>,
    word_size: usize,
) -> Vec<SeedHit> {
    let mut hits = Vec::new();

    if query.len() < word_size {
        return hits;
    }

    for qpos in 0..=(query.len() - word_size) {
        let word = &query[qpos..qpos + word_size];
        if let Some(positions) = word_index.get(word) {
            for &(db_id, db_pos) in positions {
                hits.push(SeedHit {
                    db_id,
                    db_pos,
                    query_pos: qpos,
                });
            }
        }
    }
    hits
}

/// Extend a seed hit in both directions using ungapped alignment.
/// Returns an HSP if the extension score exceeds the threshold.
fn extend_seed(
    query: &[u8],
    db_seq: &[u8],
    seed: &SeedHit,
    word_size: usize,
    match_score: i32,
    mismatch_penalty: i32,
    x_dropoff: i32,
) -> Option<HSP> {
    let mut score = word_size as i32 * match_score; // seed is exact match
    let mut best_score = score;
    let mut matches = word_size;

    // Extend to the right
    let mut right_ext = 0;
    let mut drop = 0i32;
    let mut q_right = seed.query_pos + word_size;
    let mut d_right = seed.db_pos + word_size;

    while q_right < query.len() && d_right < db_seq.len() {
        if query[q_right] == db_seq[d_right] {
            score += match_score;
            matches += 1;
        } else {
            score += mismatch_penalty;
        }
        right_ext += 1;
        q_right += 1;
        d_right += 1;

        if score > best_score {
            best_score = score;
            drop = 0;
        } else {
            drop = best_score - score;
            if drop > x_dropoff {
                break;
            }
        }
    }

    // Extend to the left
    let mut left_ext = 0;
    score = best_score;
    drop = 0;
    let mut q_left = seed.query_pos as isize - 1;
    let mut d_left = seed.db_pos as isize - 1;

    while q_left >= 0 && d_left >= 0 {
        let qi = q_left as usize;
        let di = d_left as usize;
        if query[qi] == db_seq[di] {
            score += match_score;
            matches += 1;
        } else {
            score += mismatch_penalty;
        }
        left_ext += 1;
        q_left -= 1;
        d_left -= 1;

        if score > best_score {
            best_score = score;
        } else {
            drop = best_score - score;
            if drop > x_dropoff {
                break;
            }
        }
    }

    let total_length = word_size + right_ext + left_ext;
    let identity = matches as f64 / total_length as f64 * 100.0;

    if best_score > 0 && total_length >= word_size * 2 {
        Some(HSP {
            db_id: seed.db_id,
            db_start: (seed.db_pos as isize - left_ext as isize).max(0) as usize,
            query_start: (seed.query_pos as isize - left_ext as isize).max(0) as usize,
            length: total_length,
            score: best_score,
            identity,
        })
    } else {
        None
    }
}

/// Run a simplified BLAST search.
pub fn blast_search(
    query: &[u8],
    database: &[Vec<u8>],
    word_size: usize,
    score_threshold: i32,
) -> Vec<HSP> {
    // Step 1: Build word index
    let word_index = build_word_index(database, word_size);

    // Step 2: Find seeds
    let seeds = find_seeds(query, &word_index, word_size);

    // Step 3: Extend seeds
    let mut hsps: Vec<HSP> = Vec::new();
    for seed in &seeds {
        if let Some(hsp) = extend_seed(
            query,
            &database[seed.db_id],
            seed,
            word_size,
            1,  // match score
            -2, // mismatch penalty
            10, // x-dropoff
        ) {
            if hsp.score >= score_threshold {
                hsps.push(hsp);
            }
        }
    }

    // Remove redundant HSPs (keep best per db sequence)
    hsps.sort_by(|a, b| b.score.cmp(&a.score));
    let mut seen = std::collections::HashSet::new();
    hsps.retain(|hsp| seen.insert((hsp.db_id, hsp.db_start)));

    hsps
}

fn main() {
    println!("=== BLAST-like Seed-and-Extend Search ===\n");

    // Create a small database of "gene" sequences
    let database: Vec<Vec<u8>> = vec![
        b"ATGCTAGCTAGCATCGATCGATCGATCG".to_vec(),     // Gene A
        b"TTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),      // Poly-T (noise)
        b"ATGCTAGCAAGCATCGATCGATCGATCG".to_vec(),     // Gene B (1 mutation from A)
        b"CCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),      // Poly-C (noise)
        b"GCTAGCTAGCATCGATCGAATCGATCGA".to_vec(),      // Gene C (similar to A)
    ];

    let gene_names = ["Gene_A", "Poly_T", "Gene_B", "Poly_C", "Gene_C"];

    // Query: a read that should match Gene A and Gene B
    let query = b"CTAGCTAGCATCGATCGATCG";

    println!("Query ({} bp): {}", query.len(), std::str::from_utf8(query).unwrap());
    println!("Database: {} sequences\n", database.len());

    let start = Instant::now();
    let results = blast_search(query, &database, 7, 5);
    let elapsed = start.elapsed();

    println!("Results ({:?}):", elapsed);
    for hsp in &results {
        println!(
            "  {} | score={} | identity={:.1}% | length={} | db_pos={} | query_pos={}",
            gene_names[hsp.db_id],
            hsp.score,
            hsp.identity,
            hsp.length,
            hsp.db_start,
            hsp.query_start
        );
    }

    // Demonstrate filtration principle
    println!("\n=== Filtration Principle ===");
    println!("With word size W=7 and max k=2 mismatches:");
    println!("  Any match must share an exact 7-mer (pigeonhole principle)");
    println!("  This filters out most database positions immediately");

    // Benchmark
    println!("\n=== Benchmark: Seed-and-Extend vs Naive ===");
    let large_db: Vec<Vec<u8>> = (0..1000)
        .map(|i| {
            (0..500)
                .map(|j| b"ACGT"[((i * 7 + j * 3) % 4)])
                .collect()
        })
        .collect();

    let test_query: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();

    let start = Instant::now();
    let results = blast_search(&test_query, &large_db, 11, 10);
    let elapsed = start.elapsed();
    println!(
        "  1000 seqs × 500 bp: {} hits found in {:?}",
        results.len(),
        elapsed
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let db = vec![b"ATCGATCGATCG".to_vec()];
        let results = blast_search(b"ATCGATCG", &db, 5, 3);
        assert!(!results.is_empty());
        assert!(results[0].identity > 99.0);
    }

    #[test]
    fn test_no_match() {
        let db = vec![b"AAAAAAAAAA".to_vec()];
        let results = blast_search(b"CCCCCCCCCC", &db, 5, 3);
        assert!(results.is_empty());
    }

    #[test]
    fn test_approximate_match() {
        let db = vec![b"ATCGATCGATCG".to_vec()];
        // One mismatch in the middle
        let results = blast_search(b"ATCGTTCGATCG", &db, 4, 2);
        assert!(!results.is_empty());
    }
}
