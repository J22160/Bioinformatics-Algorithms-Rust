# Edit Distance (Levenshtein Distance)

## What Is This?

Suppose you have two DNA sequences and want to measure how different they are.

**Edit distance** is the minimum number of single-character operations required to transform one sequence into another. The three allowed operations are:

1. **Insertion:** Add a character.
2. **Deletion:** Remove a character.
3. **Substitution:** Replace one character with another.

For example, transforming `ACGT` into `GCAT` requires two substitutions:

- `A → G`
- `G → A`

Therefore, the Levenshtein distance is **2**.

## How It Works

Levenshtein distance uses **dynamic programming** to solve the problem by breaking it into smaller subproblems.

Imagine a matrix where one sequence runs across the top and the other runs down the left side. Each cell `D[i][j]` stores the minimum number of edits required to transform the first `i` characters of sequence A into the first `j` characters of sequence B.

### 1. Initialize the Matrix

Let sequence A have length `n` and sequence B have length `m`.

- Set `D[i][0] = i`, because transforming a sequence of length `i` into an empty sequence requires `i` deletions.
- Set `D[0][j] = j`, because transforming an empty sequence into a sequence of length `j` requires `j` insertions.
- Set `D[0][0] = 0`, because no edits are needed to transform an empty sequence into itself.

### 2. Fill the Matrix

For each pair of characters `A[i-1]` and `B[j-1]`, calculate the minimum cost of three possible operations:

- **Deletion:** `D[i-1][j] + 1`
- **Insertion:** `D[i][j-1] + 1`
- **Match or substitution:** `D[i-1][j-1] + cost`, where `cost = 0` if the characters match and `1` otherwise.

The recurrence is:

`D[i][j] = min(D[i-1][j] + 1, D[i][j-1] + 1, D[i-1][j-1] + cost)`

Here, `min` selects the lowest-cost operation. By filling the matrix row by row, each cell reuses results from previously solved subproblems.

### 3. Read the Result

The bottom-right cell, `D[n][m]`, contains the minimum edit distance between the two complete sequences.

**Example:** Calculate the edit distance between `ACGT` and `GCAT`.

|   | ∅ | G | C | A | T |
|---|---:|---:|---:|---:|---:|
| **∅** | 0 | 1 | 2 | 3 | 4 |
| **A** | 1 | 1 | 2 | 2 | 3 |
| **C** | 2 | 2 | 1 | 2 | 3 |
| **G** | 3 | 2 | 2 | 2 | 3 |
| **T** | 4 | 3 | 3 | 3 | 2 |

The bottom-right cell is `2`, so the minimum edit distance is **2**.

Because both sequences have the same length, two substitutions are sufficient. For sequences of different lengths, the optimal transformation may also involve insertions or deletions.

## Complexity

Let `n` and `m` be the lengths of the two sequences.

- **Time:** `O(n × m)` — each matrix cell is calculated in constant time.
- **Space:** `O(n × m)` — storing the full matrix. If only the distance is needed, rolling rows reduce auxiliary space to `O(min(n, m))`.

## Bioinformatics Problem

**Scenario:** You are analyzing CRISPR guide RNA sequences and want to identify candidate off-target sites that resemble a guide sequence.

Given:

- **Guide RNA:** `GAGTCCGAGCAGAAGAAGAA` (20 nt)
- **Candidate off-target:** `GAGTCCAAGCAGAAGATGAA` (20 nt)

The two sequences differ at positions 8 and 16 (using 1-based indexing). Therefore, their Levenshtein distance is **2**.

This candidate has two substitutions and falls within a hypothetical edit-distance threshold of 5.

**Important:** Edit distance alone does not establish whether a sequence is a biologically meaningful CRISPR off-target. Actual off-target assessment depends on factors such as the PAM sequence, mismatch positions, guide-target pairing, and the Cas system used. A mismatch threshold is a screening criterion, not a universal biological cutoff.

## Why This Matters

Edit distance is useful in several bioinformatics applications:

- **Sequence comparison:** Quantifies the minimum number of edits separating two sequences.
- **Read error correction:** Identifies reads that are similar to a consensus or expected sequence.
- **Approximate pattern matching:** Finds sequence regions that differ from a query by a limited number of edits.
- **CRISPR off-target screening:** Helps identify candidate sequences with a small edit distance, alongside PAM and mismatch-position analysis.

**Related algorithm:** Needleman–Wunsch also uses dynamic programming for global sequence alignment, but it generally optimizes an alignment score based on match, mismatch, and gap scoring parameters. Levenshtein distance is a specific minimum-edit-cost problem; the two methods are related, but they are not interchangeable in every scoring configuration.