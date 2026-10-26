#!/usr/bin/env python3
"""
Quick Sort — Sorting Variant Quality Scores (Python)

Usage: python quick_sort.py
"""

import time
import random


def quick_sort(arr: list) -> None:
    """In-place quick sort with median-of-three pivot selection."""
    _quick_sort_helper(arr, 0, len(arr) - 1)


def _quick_sort_helper(arr: list, low: int, high: int) -> None:
    while low < high:
        if high - low < 16:
            _insertion_sort(arr, low, high)
            return
        pivot_idx = _partition(arr, low, high)
        # Recurse on smaller partition, loop on larger (tail-call optimization)
        if pivot_idx - low < high - pivot_idx:
            _quick_sort_helper(arr, low, pivot_idx - 1)
            low = pivot_idx + 1
        else:
            _quick_sort_helper(arr, pivot_idx + 1, high)
            high = pivot_idx - 1


def _partition(arr: list, low: int, high: int) -> int:
    """Lomuto partition with median-of-three pivot."""
    mid = (low + high) // 2
    if arr[low] > arr[mid]:
        arr[low], arr[mid] = arr[mid], arr[low]
    if arr[low] > arr[high]:
        arr[low], arr[high] = arr[high], arr[low]
    if arr[mid] > arr[high]:
        arr[mid], arr[high] = arr[high], arr[mid]
    arr[mid], arr[high] = arr[high], arr[mid]

    pivot = arr[high]
    i = low
    for j in range(low, high):
        if arr[j] <= pivot:
            arr[i], arr[j] = arr[j], arr[i]
            i += 1
    arr[i], arr[high] = arr[high], arr[i]
    return i


def _insertion_sort(arr: list, low: int, high: int) -> None:
    for i in range(low + 1, high + 1):
        j = i
        while j > low and arr[j - 1] > arr[j]:
            arr[j - 1], arr[j] = arr[j], arr[j - 1]
            j -= 1


def main():
    print("=== Quick Sort: Variant Quality Scores (Python) ===\n")

    scores = [934, 52, 1847, 300, 78, 2100, 415, 61]
    print(f"Original: {scores}")
    quick_sort(scores)
    print(f"Sorted:   {scores}\n")

    print("=== Benchmark ===")
    for n in [1_000, 10_000, 100_000, 1_000_000]:
        data = [random.randint(0, 10**9) for _ in range(n)]
        start = time.perf_counter()
        quick_sort(data)
        elapsed = time.perf_counter() - start
        print(f"  n={n:<10}  time={elapsed:.4f}s")

    print()
    print("=== Why Rust Is Faster ===")
    print("For n=1,000,000: Python ~2-4s, Rust ~0.08s (25-50x faster)")
    print("Reasons:")
    print("  1. Rust compiles partition loop to branchless cmov instructions")
    print("  2. No per-element reference counting or type checking")
    print("  3. Rust's swap() is a single register exchange; Python creates temps")
    print("  4. Median-of-three pivot avoids O(n^2) on sorted data in both")


if __name__ == "__main__":
    main()