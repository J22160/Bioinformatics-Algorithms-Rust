#!/usr/bin/env python3
"""Fitch's Algorithm for Small Parsimony (Python)"""

def fitch_bottom_up(tree):
    """Returns (score, set_of_chars) for a nested tuple tree."""
    if isinstance(tree, str):  # leaf
        return 0, {tree}
    left, right = tree
    ls, lset = fitch_bottom_up(left)
    rs, rset = fitch_bottom_up(right)
    inter = lset & rset
    if inter:
        return ls + rs, inter
    else:
        return ls + rs + 1, lset | rset

if __name__ == "__main__":
    # Tree: ((Human:A, Chimp:A), (Macaque:G, Baboon:G))
    tree = (('A', 'A'), ('G', 'G'))
    score, root_set = fitch_bottom_up(tree)
    print(f"Parsimony score: {score}")
    print(f"Root set: {root_set}")
