#!/usr/bin/env python3
"""MetaBAT2's native pairwise similarity (composite_raw.lp) -> ASP affinity/3 facts.

composite_raw.lp holds, for every pair of contigs in MetaBAT2's graph, its composite score
(tetranucleotide frequency and coverage combined), a SIMILARITY in [0,1]. This writes
affinity(C1,C2,W) with W = round(sim * SCALE), integer ids from the dictionary, C1 < C2.

TOPK keeps each contig's TOPK strongest edges (an edge stays if it is in the top TOPK of either
endpoint). It keeps the program small enough to ground: the full graph has on the order of a
hundred edges per contig. refine.sh uses TOPK=5. TOPK=0 keeps all edges.

MODE selects the weight (refine.sh uses raw):
  raw  : weight = round(sim * SCALE)
  rank : weight = the edge's rank among the kept edges, spread evenly over 1..SCALE

Usage: build_affinity.py <master_dict.tsv> <composite_raw.lp> <out_affinity.lp> [SCALE] [TOPK] [MODE]
"""
import re, sys
from collections import defaultdict
dict_file, comp_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]
SCALE = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
TOPK  = int(sys.argv[5]) if len(sys.argv) > 5 else 0
MODE  = sys.argv[6] if len(sys.argv) > 6 else "raw"

name2id = {}
for line in open(dict_file):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2 and p[0] != "Clingo_ID":
        name2id[p[1]] = int(p[0])

rx = re.compile(r'composite\("([^"]+)","([^"]+)",([0-9.]+)\)')
raw = []
skipped = 0
for line in open(comp_file):
    m = rx.match(line)
    if not m: continue
    a, b, s = m.group(1), m.group(2), float(m.group(3))
    if a not in name2id or b not in name2id:
        skipped += 1; continue
    raw.append((a, b, s))

# top-k per contig (union of each endpoint's strongest TOPK edges)
if TOPK > 0:
    nb = defaultdict(list)
    for a, b, s in raw:
        nb[a].append((s, a, b)); nb[b].append((s, a, b))
    keep = set()
    for c, lst in nb.items():
        lst.sort(reverse=True)
        for s, a, b in lst[:TOPK]:
            keep.add((a, b))
    raw = [(a, b, s) for (a, b, s) in raw if (a, b) in keep]

# weight transform
if MODE == "rank":
    # uniformise: weight = global rank of sim among kept edges, scaled to 1..SCALE
    raw.sort(key=lambda t: t[2])           # ascending similarity
    n = len(raw)
    weighted = []
    for i, (a, b, s) in enumerate(raw):
        w = 1 + round(i / max(n - 1, 1) * (SCALE - 1))   # 1..SCALE, evenly spread
        weighted.append((a, b, w))
else:  # raw
    weighted = [(a, b, round(s * SCALE)) for (a, b, s) in raw]

out, zero = [], 0
for a, b, w in weighted:
    ca, cb = name2id[a], name2id[b]
    if w < 1:
        zero += 1; continue          # drop ~zero-weight edges (no signal)
    lo, hi = (ca, cb) if ca < cb else (cb, ca)
    out.append((lo, hi, w))

out.sort()
with open(out_file, "w") as f:
    f.write(f"% MetaBAT2-native composite affinity ({MODE}, scale {SCALE}, top-{TOPK or 'all'}); "
            f"affinity/3 fed directly, NO inversion\n")
    for a, b, w in out:
        f.write(f"affinity({a},{b},{w}).\n")
print(f"[+] {out_file}: {len(out)} affinity facts  top-{TOPK or 'all'}  "
      f"(dropped {zero} ~zero-sim, {skipped} unmapped)")
