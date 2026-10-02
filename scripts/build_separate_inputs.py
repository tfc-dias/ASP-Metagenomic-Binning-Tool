#!/usr/bin/env python3
"""Inputs for the separation operator (binning_separate.lp).

The builder only prepares data; every split decision is made by the ASP model. A contig's
affinity to the marker-bearing contigs (anchors) of its own bin is read from MetaBAT2's native
composite similarity (composite_raw.lp, sim in [0,1], higher = closer).

Emits (consumed by binning_separate.lp):
  seed(B,C).        the seed assignment (B = 1-based int id by sorted bin label).
  pin1(C).          the lowest-id anchor of each bin, fixed to sub-bin 1, so that relabelings
                    of the same split are not counted as different solutions.
  affinity(C,A,W).  C's K most similar anchors in its own bin, W = round(sim*SCALE).
                    Restricting to the nearest anchors keeps a large genome in the bin from
                    winning by sheer number of edges.
Markers (mg/2) and the contamination signal (comarker_sim/3) come from build_base.py and
build_comarker.py.

Usage: build_separate_inputs.py <dict> <mg> <composite_raw> <seed.tsv> <out.lp> [K=5] [SCALE=1000]
"""
import re, sys
from collections import defaultdict

DICT, MG, COMP, INB, OUT = sys.argv[1:6]
K     = int(sys.argv[6]) if len(sys.argv) > 6 else 5
SCALE = int(sys.argv[7]) if len(sys.argv) > 7 else 1000

# --- contig dictionary: name <-> clingo id ---
c2n, n2c = {}, {}
for line in open(DICT):
    p = line.rstrip("\n").split("\t")
    if p[0] == "Clingo_ID":
        continue
    cid = int(p[0]); c2n[cid] = p[1]; n2c[p[1]] = cid

# --- MetaBAT2 assignment: name -> label -> stable 1-based int id ---
mb = {}
for line in open(INB):
    if line.startswith("@") or not line.strip():
        continue
    c, b = line.rstrip("\n").split("\t")[:2]
    mb[c] = b
labels = sorted(set(mb.values()))
lab2id = {b: i + 1 for i, b in enumerate(labels)}
cid_bin = {n2c[c]: lab2id[b] for c, b in mb.items() if c in n2c}   # binned cid -> bin int
bin_cids = defaultdict(list)
for cid, b in cid_bin.items():
    bin_cids[b].append(cid)

# --- anchors: marker-bearing contigs ---
anchor = set()
for line in open(MG):
    m = re.match(r"mg\((\d+),", line)
    if m and int(m.group(1)) in cid_bin:
        anchor.add(int(m.group(1)))

# --- contig -> nearest anchors (same bin) from the native composite SIMILARITY ---
rx = re.compile(r'composite\("([^"]+)","([^"]+)",([0-9.]+)\)')
to_anchor = defaultdict(list)                  # contig -> [(sim, anchor)]
for line in open(COMP):
    m = rx.match(line)
    if not m:
        continue
    a, b = n2c.get(m.group(1)), n2c.get(m.group(2))
    if a is None or b is None:
        continue
    ba, bb = cid_bin.get(a), cid_bin.get(b)
    if ba is None or ba != bb:                 # same bin only
        continue
    s = float(m.group(3))
    if b in anchor and a != b:
        to_anchor[a].append((s, b))
    if a in anchor and b != a:
        to_anchor[b].append((s, a))

with open(OUT, "w") as f:
    f.write("% ASP separation inputs (build_separate_inputs.py)\n")
    f.write(f"% K={K} SCALE={SCALE}; {len(anchor)} anchors; from native composite similarity\n")
    for cid, b in sorted(cid_bin.items()):
        f.write(f"seed({b},{cid}).\n")
    for b, cids in sorted(bin_cids.items()):
        bin_anchors = [c for c in cids if c in anchor]
        f.write(f"pin1({min(bin_anchors) if bin_anchors else min(cids)}).\n")
    nedge = 0
    for c in sorted(to_anchor):
        # top-K most SIMILAR anchors (highest sim)
        for s, a in sorted(to_anchor[c], reverse=True)[:K]:
            w = round(s * SCALE)
            if w <= 0:
                continue
            f.write(f"affinity({c},{a},{w}).\n")
            nedge += 1
print(f"[+] {OUT}: {len(cid_bin)} seeded contigs in {len(labels)} bins, "
      f"{len(anchor)} anchors, {nedge} contig->anchor affinity edges (K={K}, SCALE={SCALE})")
