#!/usr/bin/env python3
"""
Complexity Analysis Demo (Python)

Demonstrates O(n), O(n log n), and O(n^2) approaches to finding
the maximum gene expression value.

Usage: python complexity_demo.py
"""

import time
import random


def find_max_linear(values: list[float]) -> float:
    """O(n) — Single scan through the list."""
    current_max = float('-inf')
    for v in values:
        if v > current_max:
            current_max = v
    return current_max


def find_max_sort(values: list[float]) -> float:
    """O(n log n) — Sort and take the last element."""
    sorted_values = sorted(values)  # Python uses Timsort: O(n log n)
    return sorted_values[-1]


def find_max_quadratic(values: list[float]) -> float:
    """O(n^2) — Compare every element against all others."""
    n = len(values)
    for i in range(n):
        is_maximum = True
        for j in range(n):
            if values[j] > values[i]:
                is_maximum = False
                break
        if is_maximum:
            return values[i]
    return values[0]


def benchmark(name: str, func, *args) -> float:
    start = time.perf_counter()
    result = func(*args)
    elapsed = time.perf_counter() - start
    print(f"  {name:<25} result={result:.2f}  time={elapsed:.6f}s")
    return result


def generate_expression_data(n: int, seed: int = 42) -> list[float]:
    """Generate synthetic TPM-like expression values."""
    rng = random.Random(seed)
    return [rng.uniform(0, 100_000) for _ in range(n)]


def main():
    print("=== Complexity Analysis (Python) ===\n")

    for n in [1_000, 10_000, 100_000]:
        print(f"--- n = {n} ---")
        data = generate_expression_data(n)

        r1 = benchmark("O(n) linear scan", find_max_linear, data)
        r2 = benchmark("O(n log n) sort", find_max_sort, data)

        if n <= 50_000:
            r3 = benchmark("O(n^2) pairwise", find_max_quadratic, data)
        else:
            print(f"  {'O(n^2) pairwise':<25} SKIPPED (too slow in Python)")

        print()

    print("=== Rust vs. Python Performance Note ===")
    print("For n=1,000,000 the Rust O(n) scan completes in ~0.5ms.")
    print("Python's equivalent takes ~50ms — about 100x slower.")
    print()
    print("Why Rust is faster:")
    print("  1. Compiled to native machine code (no interpreter overhead)")
    print("  2. No garbage collector pauses")
    print("  3. Data stored contiguously in memory (cache-friendly Vec<f64>)")
    print("  4. LLVM auto-vectorizes the linear scan using SIMD instructions")
    print("  5. Zero-cost iterator abstractions compile to optimal loops")


if __name__ == "__main__":
    main()
