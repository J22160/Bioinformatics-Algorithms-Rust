//! Fractional Knapsack — Sequencing Budget Optimization

use std::time::Instant;

#[derive(Debug, Clone)]
struct Experiment {
    name: String,
    value: f64,
    cost: f64,
}

/// Solve the fractional knapsack problem.
/// Returns (total_value, Vec<(name, fraction_taken)>).
pub fn fractional_knapsack(
    experiments: &mut [Experiment],
    budget: f64,
) -> (f64, Vec<(String, f64)>) {
    // Sort by value-to-cost ratio (descending)
    experiments.sort_by(|a, b| {
        let ratio_a = a.value / a.cost;
        let ratio_b = b.value / b.cost;
        ratio_b.partial_cmp(&ratio_a).unwrap()
    });

    let mut remaining = budget;
    let mut total_value = 0.0;
    let mut selections = Vec::new();

    for exp in experiments.iter() {
        if remaining <= 0.0 {
            break;
        }
        if exp.cost <= remaining {
            // Take all of this experiment
            total_value += exp.value;
            remaining -= exp.cost;
            selections.push((exp.name.clone(), 1.0));
        } else {
            // Take a fraction
            let fraction = remaining / exp.cost;
            total_value += exp.value * fraction;
            remaining = 0.0;
            selections.push((exp.name.clone(), fraction));
        }
    }

    (total_value, selections)
}

fn main() {
    println!("=== Fractional Knapsack: Sequencing Budget Optimization ===\n");

    let mut experiments = vec![
        Experiment { name: "RNA-Seq".into(), value: 60.0, cost: 10.0 },
        Experiment { name: "WGS".into(), value: 100.0, cost: 20.0 },
        Experiment { name: "ChIP-Seq".into(), value: 120.0, cost: 30.0 },
        Experiment { name: "ATAC-Seq".into(), value: 80.0, cost: 25.0 },
        Experiment { name: "Hi-C".into(), value: 50.0, cost: 15.0 },
    ];

    let budget = 50.0;

    println!("Available experiments:");
    for exp in &experiments {
        println!("  {:<12} value={:<6.0} cost={:<6.0} ratio={:.2}",
                 exp.name, exp.value, exp.cost, exp.value / exp.cost);
    }
    println!("\nBudget: {}\n", budget);

    let start = Instant::now();
    let (total, selections) = fractional_knapsack(&mut experiments, budget);
    let elapsed = start.elapsed();

    println!("Optimal selection:");
    for (name, fraction) in &selections {
        println!("  {:<12} fraction={:.2}", name, fraction);
    }
    println!("\nTotal scientific value: {:.1}", total);
    println!("Time: {:?}", elapsed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let mut exps = vec![
            Experiment { name: "A".into(), value: 60.0, cost: 10.0 },
            Experiment { name: "B".into(), value: 100.0, cost: 20.0 },
            Experiment { name: "C".into(), value: 120.0, cost: 30.0 },
        ];
        let (value, _) = fractional_knapsack(&mut exps, 50.0);
        assert!((value - 240.0).abs() < 1e-6);
    }

    #[test]
    fn test_zero_budget() {
        let mut exps = vec![
            Experiment { name: "A".into(), value: 100.0, cost: 10.0 },
        ];
        let (value, _) = fractional_knapsack(&mut exps, 0.0);
        assert_eq!(value, 0.0);
    }
}
