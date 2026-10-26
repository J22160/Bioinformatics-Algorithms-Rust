#!/usr/bin/env python3
"""Fractional Knapsack — Sequencing Budget Optimization (Python)"""


def fractional_knapsack(experiments, budget):
    # Sort by value/cost ratio descending
    experiments.sort(key=lambda e: e["value"] / e["cost"], reverse=True)

    remaining = budget
    total_value = 0
    selections = []

    for exp in experiments:
        if remaining <= 0:
            break
        if exp["cost"] <= remaining:
            total_value += exp["value"]
            remaining -= exp["cost"]
            selections.append((exp["name"], 1.0))
        else:
            fraction = remaining / exp["cost"]
            total_value += exp["value"] * fraction
            remaining = 0
            selections.append((exp["name"], fraction))

    return total_value, selections


if __name__ == "__main__":
    experiments = [
        {"name": "RNA-Seq", "value": 60, "cost": 10},
        {"name": "WGS", "value": 100, "cost": 20},
        {"name": "ChIP-Seq", "value": 120, "cost": 30},
    ]
    value, sel = fractional_knapsack(experiments, 50)
    print(f"Total value: {value}")
    for name, frac in sel:
        print(f"  {name}: {frac:.2f}")
