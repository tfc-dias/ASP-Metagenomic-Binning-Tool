#!/usr/bin/env python3
"""One summary row per refinement stage, written to summary.tsv by refine.sh.

For a stage that turned <in.tsv> into <out.tsv> (and <filtered.tsv> after the minimum bin size),
reports what changed and how the solver finished. Everything is read from files the stage left
behind, so a stage reused from an earlier run is summarised the same way.

Columns:
  stage              the chain so far (e.g. separate_merge)
  operator           separate | merge | recruit
  bins_in/out        bins before and after the operator
  contigs_in/out     binned contigs before and after
  bins_split         input bins whose contigs went to more than one output bin
  bins_joined        output bins made from more than one input bin
  contigs_recruited  contigs binned by this stage that were unbinned before
  bins_kept, contigs_kept   after the minimum bin size filter (final stage only; blank otherwise)
  unmeasured_as_dissimilar  marker pairs MetaBAT2 could not measure, counted as dissimilar
  solver             optimal (proved best) | solved (no objective; complete) |
                     time limit (best found, not proved) | unknown
  seconds            wall time of the stage (blank if not recorded)

Usage:
  stage_summary.py --header
  stage_summary.py <stage> <operator> <in.tsv> <out.tsv> <filtered.tsv or -> <stage_dir>
"""
import glob
import os
import sys
from collections import defaultdict

COLS = ["stage", "operator", "bins_in", "bins_out", "contigs_in", "contigs_out", "bins_split",
        "bins_joined", "contigs_recruited", "bins_kept", "contigs_kept",
        "unmeasured_as_dissimilar", "solver", "seconds"]

if sys.argv[1:] == ["--header"]:
    print("\t".join(COLS))
    sys.exit(0)
if len(sys.argv) != 7:
    sys.exit(__doc__)
stage, op, fin, fout, ffilt, d = sys.argv[1:7]


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@") or not line.strip():
                continue
            c, b = line.rstrip("\n").split("\t")[:2]
            out[c] = b
    return out


a, b = load(fin), load(fout)
k = load(ffilt) if ffilt != "-" else None     # only the final stage is filtered

to_out = defaultdict(set)      # input bin -> output bins
from_in = defaultdict(set)     # output bin -> input bins
for c, bo in b.items():
    if c in a:
        to_out[a[c]].add(bo)
        from_in[bo].add(a[c])

filled = 0
for stats in glob.glob(os.path.join(d, "*.stats")):
    for line in open(stats):
        key, val = line.split()
        if key == "filled":
            filled += int(val)

solver = "unknown"
clingo = os.path.join(d, "clingo.txt")
if os.path.exists(clingo):
    text = open(clingo).read()
    if "OPTIMUM FOUND" in text:
        solver = "optimal"
    elif "TIME LIMIT" in text:
        solver = "time limit"
    elif "SATISFIABLE" in text:
        solver = "solved"

seconds = ""
sec_file = os.path.join(d, "seconds.txt")
if os.path.exists(sec_file):
    seconds = open(sec_file).read().strip()

row = [stage, op, len(set(a.values())), len(set(b.values())), len(a), len(b),
       sum(1 for s in to_out.values() if len(s) > 1),
       sum(1 for s in from_in.values() if len(s) > 1),
       len(set(b) - set(a)), len(set(k.values())) if k is not None else "",
       len(k) if k is not None else "", filled, solver, seconds]
print("\t".join(str(x) for x in row))
