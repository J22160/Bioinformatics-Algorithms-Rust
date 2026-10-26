#!/usr/bin/env python3
"""
Selection Sort — Sorting Phred Quality Scores (Python)

This is the Python equivalent of the Rust implementation.
Compare execution times to see why Rust is preferred for large-scale
bioinformatics sorting operations.

Usage: python selection_sort.py
"""

import time


def selection_sort(arr: list) -> None:
    """In-place selection sort. O(n²) time, O(1) extra space."""
    n = len(arr)
    for i in range(n):
        min_idx = i
        for j in range(i + 1, n):
            if arr[j] < arr[min_idx]:
                min_idx = j
        if min_idx != i:
            arr[i], arr[min_idx] = arr[min_idx], arr[i]


def main():
    print("=== Selection Sort: Phred Quality Score Filtering (Python) ===\n")

    scores = [38, 12, 27, 35, 8, 41, 19, 33, 5, 40, 28, 15, 37, 22, 9, 31,
              2, 39, 17, 36, 11, 30, 24, 7, 41, 13, 26, 34, 20, 3, 38, 16]

    print(f"Original scores: {scores}")

    start = time.perf_counter()
    selection_sort(scores)
    elapsed = time.perf_counter() - start

    print(f"Sorted scores:   {scores}")
    print(f"Sort time:       {elapsed:.6f}s\n")

    # Scaling benchmark
    print("=== Scaling Benchmark ===")
    for n in [100, 1_000, 5_000, 10_000]:
        data = [(i * 7 + 13) % 42 for i in range(n)]
        start = time.perf_counter()
        selection_sort(data)
        elapsed = time.perf_counter() - start
        is_sorted = all(data[i] <= data[i + 1] for i in range(len(data) - 1))
        print(f"  n={n:<6}  time={elapsed:>10.6f}s  (sorted: {is_sorted})")

    print()
    print("=== Why Rust Is Faster ===")
    print("For n=10,000: Python takes ~5-10 seconds, Rust takes ~50ms")
    print("Reason: Python's for-loop has ~100ns overhead per iteration")
    print("(bytecode dispatch, reference counting, type checking)")
    print("Rust compiles to native instructions: ~1ns per comparison")


if __name__ == "__main__":
    main()
