#!/usr/bin/env python3
"""BQSR (Python)"""
import math
from collections import defaultdict

def build_recal_table(observations):
    table = defaultdict(lambda: [0, 0])
    for key, is_err in observations:
        table[key][0] += 1
        if is_err: table[key][1] += 1
    return table

if __name__ == "__main__":
    table = defaultdict(lambda: [0, 0])
    for q in [20, 25, 30]:
        e = 10**(-q/10) * 1.5
        n = 10000
        for i in range(n):
            table[q][0] += 1
            if i < int(n * e): table[q][1] += 1
    for q in [20, 25, 30]:
        total, errs = table[q]
        emp_q = -10 * math.log10((errs+1)/(total+2))
        print(f"Q{q}: empirical Q{emp_q:.1f} ({total} bases, {errs} errors)")
