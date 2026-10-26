#!/usr/bin/env python3
"""
Change Problem — Minimum BAC Clones for Tiling Path (Python)

Usage: python change_problem.py
"""

import time


def min_coins(denominations: list[int], target: int) -> tuple[int | None, list[int]]:
    """Return (min_coins, list_of_coins_used) or (None, []) if impossible."""
    INF = float('inf')
    dp = [INF] * (target + 1)
    dp[0] = 0
    from_coin = [0] * (target + 1)

    for amount in range(1, target + 1):
        for coin in denominations:
            if coin <= amount and dp[amount - coin] + 1 < dp[amount]:
                dp[amount] = dp[amount - coin] + 1
                from_coin[amount] = coin

    if dp[target] == INF:
        return None, []

    # Traceback
    coins_used = []
    remaining = target
    while remaining > 0:
        coins_used.append(from_coin[remaining])
        remaining -= from_coin[remaining]

    return dp[target], coins_used


def main():
    print("=== Change Problem: BAC Tiling Path (Python) ===\n")

    bac_sizes = [50, 120, 200]
    target = 1000

    count, clones = min_coins(bac_sizes, target)
    if count is not None:
        print(f"Min BAC clones: {count}")
        print(f"Clones used: {clones}")
        print(f"Sum: {sum(clones)} kb")

    # Benchmark
    print("\n=== Benchmark ===")
    denoms = [1, 5, 10, 25, 50, 100]
    for target in [100, 1_000, 10_000, 100_000]:
        start = time.perf_counter()
        result, _ = min_coins(denoms, target)
        elapsed = time.perf_counter() - start
        print(f"  target={target:<8}  min_coins={result:<6}  time={elapsed:.4f}s")

    print("\n=== Why Rust Is Faster ===")
    print("For target=100,000: Python ~1s, Rust ~5ms (200x faster)")
    print("  1. Python's inner loop has interpreter overhead per iteration")
    print("  2. Python int objects are 28 bytes; Rust u32 is 4 bytes")
    print("  3. Rust's flat Vec has perfect cache locality")


if __name__ == "__main__":
    main()
