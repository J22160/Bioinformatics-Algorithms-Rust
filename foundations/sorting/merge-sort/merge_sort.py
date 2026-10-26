#!/usr/bin/env python3
"""
Merge Sort — Sorting Genomic Coordinates (Python)

Usage: python merge_sort.py
"""

import time


def merge_sort(arr: list) -> list:
    """Recursive merge sort returning a new sorted list. O(n log n)."""
    if len(arr) <= 1:
        return arr

    mid = len(arr) // 2
    left = merge_sort(arr[:mid])
    right = merge_sort(arr[mid:])
    return merge(left, right)


def merge(left: list, right: list) -> list:
    """Merge two sorted lists into one sorted list."""
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result


def main():
    print("=== Merge Sort: Genomic Coordinates (Python) ===\n")

    positions = [155_799, 5_529, 87_160, 24_601, 117_482, 48_930, 72_215]
    print(f"Unsorted: {positions}")
    sorted_pos = merge_sort(positions)
    print(f"Sorted:   {sorted_pos}\n")

    print("=== Scaling Benchmark ===")
    for n in [1_000, 10_000, 100_000, 1_000_000]:
        data = [(i * 2654435761) % 3_000_000_000 for i in range(n)]
        start = time.perf_counter()
        sorted_data = merge_sort(data)
        elapsed = time.perf_counter() - start
        print(f"  n={n:<10}  time={elapsed:>10.4f}s")

    print()
    print("=== Why Rust Is Faster ===")
    print("For n=1,000,000: Python ~3-5s, Rust ~0.15s (20-30x faster)")
    print("Reasons:")
    print("  1. Python creates many temporary list objects (GC overhead)")
    print("  2. Each append() call may trigger reallocation")
    print("  3. Rust reuses a single pre-allocated buffer")
    print("  4. Rust's copy_from_slice compiles to memcpy (SIMD-accelerated)")
    print("  5. Rust switches to insertion sort for small subarrays")


if __name__ == "__main__":
    main()