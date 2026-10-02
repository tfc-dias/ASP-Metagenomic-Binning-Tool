#!/usr/bin/env python3
"""Restrict mg/2 facts to the contigs mentioned in a seed/pinned fact file.

The separation model grounds sharemarker(C1,C2) :- mg(C1,G), mg(C2,G) over ALL
marker-bearing contigs. On a real metagenome that is millions of pairs, most of them
irrelevant: separation only reasons about SEEDED (binned) contigs (sub/2 exists only
for them), so markers on unbinned contigs cannot affect any split. Feeding a
seed-restricted mg is therefore SEMANTICALLY IDENTICAL and keeps grounding small.

Usage: mg_subset.py <mg.lp> <seed_or_pinned.lp> <out_mg.lp>
"""
import re
import sys

if len(sys.argv) != 4:
    sys.exit(__doc__)
MG, REF, OUT = sys.argv[1:4]

keep = set()
for line in open(REF):
    m = re.match(r"(?:seed|pinned)\(\d+,(\d+)\)", line.strip())
    if m:
        keep.add(int(m.group(1)))

n = 0
with open(OUT, "w") as o:
    o.write(f"% mg restricted to {len(keep)} seeded contigs (grounding safety; semantics unchanged)\n")
    for line in open(MG):
        m = re.match(r"mg\((\d+),\d+\)", line)
        if m and int(m.group(1)) in keep:
            o.write(line)
            n += 1
print(f"[+] {OUT}: {n} mg facts over {len(keep)} seeded contigs")
