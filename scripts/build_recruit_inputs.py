#!/usr/bin/env python3
"""Recruitment inputs from the NATIVE composite graph (generic, dump-free).

Pins an ARBITRARY current binning (MetaBAT2, or a refined sep/merge prediction) and
lets binning_recruit.lp place the still-UNBINNED (recruitable) contigs. Recruit edges
come from the native composite similarity (composite_raw.lp), not from the intractable
O(n^2) abundance and tetranucleotide dumps.

Emits:
  <out_recruitable.lp>  recruitable(C): every contig of the dictionary the pinned binning
                        leaves unbinned (the contamination guard is built for these).
  <out.lp>              pinned(B,C).        the pinned binning (B = int id, first-seen order).
                        recruitable(C).     the contigs to (optionally) place.
                        affinity(Cu,Cb,W).  recruitable Cu -> pinned Cb, W = round(sim*1000).

comarker_sim (the contamination guard) is supplied SEPARATELY to the model (built by
build_comarker.py --mode guard), reused as-is.

Usage:
  build_recruit_inputs.py <dict.tsv> <composite_raw.lp> <pinned_binning.tsv> <out_recruitable.lp> <out.lp>
"""
import re
import sys

if len(sys.argv) != 6:
    sys.exit(__doc__)
DICT, COMP, PINNED, RECR, OUT = sys.argv[1:6]

name2id = {}
with open(DICT) as f:
    next(f, None)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            name2id[p[1]] = int(p[0])

# pinned binning: contig name -> label -> stable 1-based int id
lab2int, nextb, binned_bin = {}, 1, {}
for line in open(PINNED):
    if line.startswith("@") or not line.strip():
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 2 or p[0] not in name2id:
        continue
    lab = p[1]
    if lab not in lab2int:
        lab2int[lab] = nextb
        nextb += 1
    binned_bin[name2id[p[0]]] = lab2int[lab]

recruitable = set(name2id.values()) - set(binned_bin)
with open(RECR, "w") as f:
    f.write("% recruitable(Contig): contigs the pinned binning leaves unbinned\n")
    for c in sorted(recruitable):
        f.write(f"recruitable({c}).\n")

crx = re.compile(r'composite\("([^"]+)","([^"]+)",([0-9.eE+-]+)\)')
edges = []
for line in open(COMP):
    m = crx.match(line)
    if not m:
        continue
    a, b, s = m.group(1), m.group(2), float(m.group(3))
    if a not in name2id or b not in name2id:
        continue
    ia, ib = name2id[a], name2id[b]
    w = round(s * 1000)
    if w <= 0:
        continue
    if ia in recruitable and ib in binned_bin:
        edges.append((ia, ib, w))
    elif ib in recruitable and ia in binned_bin:
        edges.append((ib, ia, w))

with open(OUT, "w") as f:
    f.write("% recruit inputs (native-composite recruit edges); pinned = current binning\n")
    for c, b in sorted(binned_bin.items()):
        f.write(f"pinned({b},{c}).\n")
    for c in sorted(recruitable):
        f.write(f"recruitable({c}).\n")
    for c, cb, w in edges:
        f.write(f"affinity({c},{cb},{w}).\n")

with_edge = len({c for c, _, _ in edges})
print(f"[+] {OUT}: {len(binned_bin)} pinned in {len(lab2int)} bins, "
      f"{len(recruitable)} recruitable, {len(edges)} recruit edges "
      f"({with_edge} recruitable have >=1 edge)")
