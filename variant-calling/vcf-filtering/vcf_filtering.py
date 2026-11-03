#!/usr/bin/env python3
"""VCF Hard Filtering (Python)"""

def hard_filter(variants, min_qual=30, min_dp=10, min_qd=2.0, max_fs=60):
    for v in variants:
        reasons = []
        if v['qual'] < min_qual: reasons.append("LowQual")
        if v['dp'] < min_dp: reasons.append("LowDP")
        if v['qd'] < min_qd: reasons.append("LowQD")
        if v['fs'] > max_fs: reasons.append("StrandBias")
        status = "PASS" if not reasons else ",".join(reasons)
        print(f"  {v['chrom']}:{v['pos']} {v['ref']}->{v['alt']} -> {status}")

if __name__ == "__main__":
    variants = [
        {"chrom":"chr1","pos":1000,"ref":"A","alt":"G","qual":500,"dp":45,"qd":11.1,"fs":2.3},
        {"chrom":"chr1","pos":2000,"ref":"C","alt":"T","qual":15,"dp":5,"qd":3.0,"fs":1.0},
    ]
    hard_filter(variants)
