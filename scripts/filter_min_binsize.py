#!/usr/bin/env python3
"""Post-hoc minimum-bin-size filter (reporting parity with MetaBAT2).

MetaBAT2 emits no bin below its default `-s 200000` bp; our ASP refinement does
not filter, so separation can leave sub-threshold single-contig scraps. This
drops any bin whose total length < MINBP; its contigs become UNBINNED (simply
omitted), exactly as MetaBAT2 treats sub-threshold material -- keeping the
comparison apples-to-apples. The threshold is MetaBAT2's own documented default,
NOT tuned against AMBER.

Usage: filter_min_binsize.py <pred.tsv> <dict.tsv> <out.tsv> [minbp=200000]
  pred.tsv : AMBER-format prediction (@ headers + SEQUENCEID<TAB>BINID)
  dict.tsv : contig dictionary from build_base.py (id, contig name, length; 1 header row)
"""
import sys, csv
from collections import defaultdict

pred, dictf, out = sys.argv[1], sys.argv[2], sys.argv[3]
MINBP = int(sys.argv[4]) if len(sys.argv) > 4 else 200000

# orig contig -> length
L = {}
with open(dictf) as fh:
    r = csv.reader(fh, delimiter="\t"); next(r)
    for row in r:
        L[row[1]] = int(row[2])

# read prediction: keep headers, collect (contig, bin) data rows
headers, rows = [], []
for line in open(pred):
    if line.startswith("@"):
        headers.append(line.rstrip("\n"))
    elif line.strip():
        c, b = line.rstrip("\n").split("\t")[:2]
        rows.append((c, b))

# bin -> total bp
binbp = defaultdict(int)
for c, b in rows:
    binbp[b] += L.get(c, 0)

dropped_bins = {b for b, bp in binbp.items() if bp < MINBP}
kept = [(c, b) for c, b in rows if b not in dropped_bins]

with open(out, "w") as fh:
    for h in headers:
        fh.write(h + "\n")
    for c, b in kept:
        fh.write(f"{c}\t{b}\n")

n_bins_before = len(binbp)
n_bins_after = n_bins_before - len(dropped_bins)
dropped_bp = sum(binbp[b] for b in dropped_bins)
dropped_contigs = sum(1 for c, b in rows if b in dropped_bins)
print(f"{pred}")
print(f"  min bin size      : {MINBP} bp (MetaBAT2 default)")
print(f"  bins              : {n_bins_before} -> {n_bins_after}  (dropped {len(dropped_bins)})")
print(f"  contigs unbinned  : {dropped_contigs}  ({dropped_bp} bp)")
print(f"  dropped bins (bp) : {sorted(round(binbp[b]/1000,1) for b in dropped_bins)} kb")
print(f"  wrote             : {out}")
