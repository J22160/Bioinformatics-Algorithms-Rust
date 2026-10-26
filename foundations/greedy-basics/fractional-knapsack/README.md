# Fractional Knapsack

## What Is This?

You have a sequencing budget (capacity) and several sequencing
experiments to choose from. Each experiment has a cost and a scientific
value. Unlike the 0/1 knapsack, you can run a **fraction** of an
experiment (e.g., sequence only half the samples).

The **Fractional Knapsack** problem asks: how do you maximize total
scientific value while staying within budget?

## The Greedy Strategy

1. Compute the **value-to-cost ratio** for each experiment
2. Sort experiments by ratio (highest first)
3. Take as much as possible of the highest-ratio experiment
4. If budget remains, take the next highest-ratio experiment
5. Repeat until budget is exhausted

**This greedy approach is provably optimal** for the fractional case.

## Example

| Experiment | Value | Cost | Ratio |
|-----------|-------|------|-------|
| RNA-Seq   | 60    | 10   | 6.0   |
| WGS       | 100   | 20   | 5.0   |
| ChIP-Seq  | 120   | 30   | 4.0   |

Budget: 50

1. Take all of RNA-Seq (cost 10, value 60) → remaining budget: 40
2. Take all of WGS (cost 20, value 100) → remaining budget: 20
3. Take 20/30 of ChIP-Seq (value 120 × 20/30 = 80) → budget: 0

**Total value: 60 + 100 + 80 = 240**

## Complexity

- **Time**: O(n log n) for sorting
- **Space**: O(n)
