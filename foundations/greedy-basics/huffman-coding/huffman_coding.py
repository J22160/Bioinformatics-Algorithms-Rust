#!/usr/bin/env python3
"""Huffman Coding — FASTQ Quality String Compression (Python)"""

import heapq
from collections import Counter


class HuffNode:
    def __init__(self, char=None, freq=0, left=None, right=None):
        self.char = char
        self.freq = freq
        self.left = left
        self.right = right

    def __lt__(self, other):
        return self.freq < other.freq


def build_huffman_tree(data: bytes) -> HuffNode | None:
    freq = Counter(data)
    if not freq:
        return None

    heap = [HuffNode(char=ch, freq=f) for ch, f in freq.items()]
    heapq.heapify(heap)

    while len(heap) > 1:
        left = heapq.heappop(heap)
        right = heapq.heappop(heap)
        merged = HuffNode(freq=left.freq + right.freq, left=left, right=right)
        heapq.heappush(heap, merged)

    return heap[0]


def extract_codes(root: HuffNode) -> dict:
    codes = {}

    def walk(node, prefix=""):
        if node.char is not None:
            codes[node.char] = prefix if prefix else "0"
            return
        if node.left:
            walk(node.left, prefix + "0")
        if node.right:
            walk(node.right, prefix + "1")

    walk(root)
    return codes


if __name__ == "__main__":
    quality = b"IIIIIIIIIIFFFFFFFFFFFFFFFFFFF!!!!!!###########JJJJJJJJ"
    print(f"Original: {len(quality)} bytes = {len(quality) * 8} bits")

    tree = build_huffman_tree(quality)
    codes = extract_codes(tree)

    for ch, code in sorted(codes.items(), key=lambda x: len(x[1])):
        print(f"  '{chr(ch)}' -> {code} ({len(code)} bits)")

    compressed = sum(len(codes[b]) for b in quality)
    print(f"\nCompressed: {compressed} bits")
    print(f"Ratio: {(1 - compressed / (len(quality) * 8)) * 100:.1f}% reduction")
