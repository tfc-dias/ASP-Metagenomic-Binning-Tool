#!/usr/bin/env python3
"""Per-bin sub-bin budget for the separation operator (binning_separate.lp).

Separation may only split a seed bin into as many sub-bins as its marker conflicts require.
For each bin, the cannot-link graph has the bin's marker-bearing contigs as nodes and an edge
between every two that share a marker and are dissimilar (D >= simthresh). The number of
sub-bins needed is the colouring number of that graph; a greedy colouring gives a count that
always admits a valid split, so no bin can make the whole problem unsatisfiable. Bins with no
such edge get a budget of 1 (the model's default) and cannot split.

The budget uses only the marker and distance facts and the same simthresh the model uses.

Emits binmax(B,K) for K >= 2. No cap by default; the model labels sub-bins B*100+S, so K up to
99 is representable.

Usage: build_binmax.py <seed(separate_inputs).lp> <comarker_within.lp> <out.lp> [simthresh=50] [slotcap=0 (none)]
"""
import re, sys
from collections import defaultdict

SEED, WITHIN, OUT = sys.argv[1:4]
ST      = int(sys.argv[4]) if len(sys.argv) > 4 else 50
SLOTCAP = int(sys.argv[5]) if len(sys.argv) > 5 else 0   # 0 = no cap

# seed bin per contig
binof = {}
for l in open(SEED):
    m = re.match(r"seed\((\d+),(\d+)\)", l)
    if m:
        binof[int(m.group(2))] = int(m.group(1))

# cannot-link edges per bin (within-bin same-marker pairs with D >= simthresh)
edges = defaultdict(list)   # bin -> [(c1,c2)]
for l in open(WITHIN):
    m = re.match(r"comarker_sim\((\d+),(\d+),(\d+)\)", l)
    if not m:
        continue
    c1, c2, D = int(m.group(1)), int(m.group(2)), int(m.group(3))
    b = binof.get(c1)
    if b is not None and binof.get(c2) == b and D >= ST:
        edges[b].append((c1, c2))

def greedy_colours(elist):
    adj = defaultdict(set)
    for a, b in elist:
        adj[a].add(b); adj[b].add(a)
    colour = {}
    # colour high-degree nodes first (fewer colours)
    for n in sorted(adj, key=lambda x: -len(adj[x])):
        used = {colour[m] for m in adj[n] if m in colour}
        c = 0
        while c in used:
            c += 1
        colour[n] = c
    return (max(colour.values()) + 1) if colour else 1

rows = []
capped = 0
for b, elist in edges.items():
    k = greedy_colours(elist)
    if SLOTCAP and k > SLOTCAP:
        k = SLOTCAP; capped += 1
    if k >= 2:
        rows.append((b, k))

with open(OUT, "w") as f:
    f.write(f"% per-bin sub-slot budget (greedy colouring of cannot-link graph, simthresh={ST}, cap={SLOTCAP or 'none'})\n")
    for b, k in sorted(rows):
        f.write(f"binmax({b},{k}).\n")

hist = defaultdict(int)
for _, k in rows:
    hist[k] += 1
print(f"[+] {OUT}: {len(rows)} bins need >=2 slots"
      + (f" ({capped} capped at {SLOTCAP})" if capped else ""))
print("  budget histogram:", dict(sorted(hist.items())))
