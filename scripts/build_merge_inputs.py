#!/usr/bin/env python3
"""Inputs for the ASP merge model (binning_merge.lp).

The third refinement operator, completing the set: where binning_recruit.lp ADDS
unbinned contigs and binning_separate.lp SPLITS chimeric bins, this MERGES bins a
fast binner OVER-SPLIT (one genome scattered across two/three bins). This builder
ONLY marshals data; ALL merge decisions are made by the ASP model.

Candidate driver = native top-5 cross-edge COUNT between bins (NOT total distance
weight, which is size-biased toward big bins; the per-contig top-5 cap makes the
count size-robust). The model then keeps only MUTUAL-nearest-neighbour pairs and
gates every merge through the affinity-aware contamination guard.

Emits:
  seed(B,C).        the seed binning (B = integer bin id, 1-based by sorted label).
  binaff(B1,B2,N).  B1<B2; N = number of native top-5 affinity edges that cross
                    between bin B1 and bin B2 (the cross-bin connectivity). The
                    model symmetrises, takes each bin's argmax (best partner),
                    keeps mutual pairs, and merges the guard-safe ones.
Markers (mg/2) and the contamination signal (comarker_sim/3) come from the existing
fact files, reused as-is (the model computes the guard itself).

Usage: build_merge_inputs.py <dict> <native_affinity.lp> <seed.tsv> <out.lp>
"""
import re, sys
from collections import defaultdict

DICT, AFF, INB, OUT = sys.argv[1:5]

# --- contig dictionary: original name <-> clingo id ---
c2n, n2c = {}, {}
for line in open(DICT):
    p = line.rstrip("\n").split("\t")
    if p[0] == "Clingo_ID":
        continue
    cid = int(p[0]); c2n[cid] = p[1]; n2c[p[1]] = cid

# --- seed binning: name -> label, label -> stable 1-based int id ---
mb = {}
for line in open(INB):
    if line.startswith("@") or not line.strip():
        continue
    c, b = line.rstrip("\n").split("\t")[:2]
    mb[c] = b
labels = sorted(set(mb.values()))
lab2id = {b: i + 1 for i, b in enumerate(labels)}
cid_bin = {n2c[c]: lab2id[b] for c, b in mb.items() if c in n2c}   # binned cid -> bin int

# --- native top-5 affinity edges -> cross-bin edge COUNT ---
xcount = defaultdict(int)
rx = re.compile(r"affinity\((\d+),\s*(\d+),")
for line in open(AFF):
    m = rx.match(line)
    if not m:
        continue
    a, b = int(m.group(1)), int(m.group(2))
    ba, bb = cid_bin.get(a), cid_bin.get(b)
    if ba is None or bb is None or ba == bb:
        continue
    key = (ba, bb) if ba < bb else (bb, ba)
    xcount[key] += 1

with open(OUT, "w") as f:
    f.write("% ASP merge inputs (build_merge_inputs.py)\n")
    f.write(f"% native top-5 cross-edge counts; bin int ids: {lab2id}\n")
    for cid, b in sorted(cid_bin.items()):
        f.write(f"seed({b},{cid}).\n")
    for (b1, b2), n in sorted(xcount.items()):
        f.write(f"binaff({b1},{b2},{n}).\n")
print(f"[+] {OUT}: {len(cid_bin)} seeded contigs in {len(labels)} bins, "
      f"{len(xcount)} cross-bin candidate pairs")
