# Greedy Algorithms: Interval Scheduling

## What Is a Greedy Algorithm?

A **greedy algorithm** builds a solution step by step by choosing the option that looks best at the current stage. Instead of exploring every possible combination, it makes a locally optimal choice and proceeds with the remaining problem.

Greedy algorithms are efficient, but they do not always produce a globally optimal solution. Their correctness depends on whether the problem has the right mathematical properties.

For example:
- **Interval scheduling:** Selecting the maximum number of non-overlapping intervals using the earliest-finish-time strategy produces an optimal solution.
- **Change-making:** Greedily choosing the largest available coin denomination can fail to produce the minimum number of coins for some denomination systems.

This example explores interval scheduling, a classic problem that demonstrates when a greedy strategy works and why.

## The Interval Scheduling Problem

Imagine you have access to a single sequencing machine and several experiments that must run on it. Each experiment has a start time and an end time, and the machine can process only one experiment at a time.

Your objective is to **schedule the maximum possible number of experiments without any overlap**.

Formally, given a set of `n` intervals:

`I = {(s₁, f₁), (s₂, f₂), ..., (sₙ, fₙ)}`

where `sᵢ` is the start time and `fᵢ` is the finish time, select the largest subset of mutually compatible intervals.

Two intervals are compatible if one finishes at or before the other begins. We represent intervals as half-open ranges `[start, end)`, so an experiment ending at time `5` can be followed by another starting at time `5`.

**Important:** The objective is to maximize the *number* of experiments, not total experiment duration or scientific value.

## The Greedy Strategy: Earliest Finish Time

The key idea is to **always select the compatible experiment that finishes earliest**.

Why? An experiment that finishes earlier frees the machine sooner, leaving as much time as possible for the experiments that follow.

Follow these steps:

1. **Sort:** Arrange all experiments by increasing finish time.
2. **Initialize:** Start with an empty schedule and set the current finish time to negative infinity.
3. **Select:** Consider experiments in sorted order. Select an experiment if its start time is greater than or equal to the finish time of the last selected experiment.
4. **Update:** After selecting an experiment, update the current finish time to its finish time.
5. **Continue:** Skip experiments that overlap with the current schedule and repeat until every experiment has been considered.

The algorithm makes only one pass after sorting. It does not need to backtrack or reconsider earlier selections.

## Example: Scheduling Sequencing Experiments

Suppose the following experiments are available:

| Experiment | Start | Finish | Interval |
|---|---:|---:|---|
| A | 0 | 6 | `[0, 6)` |
| B | 1 | 4 | `[1, 4)` |
| C | 3 | 5 | `[3, 5)` |
| D | 5 | 7 | `[5, 7)` |
| E | 3 | 9 | `[3, 9)` |
| F | 5 | 9 | `[5, 9)` |
| G | 6 | 10 | `[6, 10)` |
| H | 8 | 11 | `[8, 11)` |

### Step 1: Sort by Finish Time

Sort the experiments by increasing finish time:

`B(1,4) → C(3,5) → A(0,6) → D(5,7) → E(3,9) → F(5,9) → G(6,10) → H(8,11)`

Experiments with the same finish time can appear in either order.

### Step 2: Select Compatible Experiments

**1. Select B `[1, 4)`**

B finishes at time `4`, earlier than every other experiment. It is the first greedy choice.

Current schedule: `B`

**2. Consider C `[3, 5)`**

C starts at time `3`, before B finishes at time `4`. The experiments overlap, so skip C.

**3. Consider A `[0, 6)`**

A also overlaps with B. Skip A.

**4. Consider D `[5, 7)`**

D starts at time `5`, which is after B finishes at time `4`. Select D.

Current schedule: `B → D`

**5. Consider E `[3, 9)`**

E starts before D finishes at time `7`. Skip E.

**6. Consider F `[5, 9)`**

F also overlaps with D. Skip F.

**7. Consider G `[6, 10)`**

G starts at time `6`, before D finishes at time `7`. Skip G.

**8. Consider H `[8, 11)`**

H starts at time `8`, after D finishes at time `7`. Select H.

Final schedule:

`B [1,4) → D [5,7) → H [8,11)`

**Maximum number of compatible experiments: 3**

### Step 3: Visualize the Schedule

```text
Time →  0  1  2  3  4  5  6  7  8  9  10 11

A       [-----------)
B          [------)
C                [----)
D                      [----)
E                [-----------)
F                      [--------)
G                         [--------)
H                               [------)

Selected:
B          [------)
D                      [----)
H                               [------)
```

The selected experiments do not overlap. Although other combinations are possible, none can contain more than three experiments for this input.

## Why Does This Greedy Strategy Work?

A greedy choice is not automatically optimal. For interval scheduling, however, the earliest-finish-time strategy can be proven optimal using an **exchange argument**.

1. **Choose the earliest-finishing interval.** Let this interval be `g`. By definition, no other interval finishes earlier than `g`.

2. **Compare it with an optimal solution.** Suppose an optimal schedule begins with interval `o`. Since `g` finishes no later than `o`, replacing `o` with `g` cannot prevent any interval that followed `o` from being scheduled.

3. **Preserve optimality.** The replacement produces another schedule with the same number of intervals, so there is an optimal schedule that begins with `g`.

4. **Repeat on the remaining intervals.** After selecting `g`, discard all intervals that overlap it and apply the same strategy to the compatible intervals that remain.

By repeating this reasoning, every greedy choice can be part of an optimal schedule. Therefore, the earliest-finish-time algorithm maximizes the number of non-overlapping intervals.

## Bioinformatics Connection

### Where Is Interval Scheduling Useful?

Many bioinformatics problems involve selecting compatible genomic intervals, including candidate exons, sequencing runs, and genomic features.

Consider a simplified **exon-selection problem**. A gene-prediction pipeline identifies candidate exons along a genomic region, and you want to select a set of non-overlapping candidates.

If every candidate exon has equal value and the only objective is to maximize the number of non-overlapping intervals, interval scheduling is a suitable model.

However, real exon chaining usually involves additional biological constraints, such as:
- Reading-frame consistency and splice-site compatibility.
- Strand orientation and exon order.
- Alignment scores or other evidence supporting each exon.
- Maximizing a total score rather than the number of exons.

In those settings, **weighted interval scheduling or more specialized dynamic programming algorithms** may be more appropriate. The basic earliest-finish-time greedy algorithm is not guaranteed to find the optimal biological solution when candidates have different scores or additional compatibility constraints.

## Complexity Analysis

Let `n` be the number of intervals.

- **Sorting:** `O(n log n)` to order intervals by finish time.
- **Selection:** `O(n)` to scan the sorted intervals once.
- **Total time:** `O(n log n)`, dominated by sorting.
- **Auxiliary space:** `O(1)` for the selection scan if the intervals are already sorted and the schedule is not copied; typically `O(n)` overall if sorting requires additional storage or the selected intervals are stored separately.

## Key Takeaways

- Greedy algorithms make locally optimal choices without exploring every possible solution.
- For interval scheduling, choosing the **earliest-finishing compatible interval** guarantees the maximum number of non-overlapping intervals.
- The algorithm is efficient, with `O(n log n)` time complexity.
- The same strategy should not be applied blindly to weighted intervals or biological problems with additional constraints.
- Understanding *why* a greedy strategy is optimal is just as important as understanding how to implement it.