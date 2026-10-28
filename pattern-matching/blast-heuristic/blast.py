#!/usr/bin/env python3
"""BLAST-like Seed-and-Extend (Python)"""

from collections import defaultdict


def build_word_index(database, w):
    index = defaultdict(list)
    for db_id, seq in enumerate(database):
        for pos in range(len(seq) - w + 1):
            word = seq[pos:pos+w]
            index[word].append((db_id, pos))
    return index


def blast_search(query, database, w=7, threshold=5):
    index = build_word_index(database, w)
    results = []

    for qpos in range(len(query) - w + 1):
        word = query[qpos:qpos+w]
        for db_id, db_pos in index.get(word, []):
            # Simple ungapped extension
            seq = database[db_id]
            score, matches, length = w, w, w

            # Extend right
            qi, di = qpos + w, db_pos + w
            while qi < len(query) and di < len(seq):
                if query[qi] == seq[di]:
                    score += 1; matches += 1
                else:
                    score -= 2
                length += 1; qi += 1; di += 1
                if score < 0:
                    break

            if score >= threshold:
                results.append({
                    "db_id": db_id, "score": score,
                    "identity": matches/length*100
                })

    # Deduplicate
    seen = set()
    unique = []
    for r in sorted(results, key=lambda x: -x["score"]):
        key = r["db_id"]
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


if __name__ == "__main__":
    db = [b"ATCGATCGATCGATCG", b"TTTTTTTTTTTTTTTT",
          b"ATCGATCAATCGATCG"]
    query = b"ATCGATCGATCG"
    results = blast_search(query, db)
    for r in results:
        print(f"  DB[{r['db_id']}]: score={r['score']} identity={r['identity']:.1f}%")
