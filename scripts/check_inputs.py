#!/usr/bin/env python3
"""Check that the inputs of a run belong together, before anything long starts.

Problems this catches, each of which would otherwise surface only partway through a run:
  * the depth file does not list the assembly's contigs, or lists them in another order
    (MetaBAT2 stops with an error in that case);
  * the marker file or the seed uses contig names that do not match the FASTA headers;
  * duplicate contig names in the assembly.

Contig names are the first word of each FASTA header. A marker or seed name also matches if the
part before '_LN:' does (some assemblies store the full header with spaces turned into '_').

Prints one line per check and exits with status 1 if any check fails. Warnings do not fail.

Usage:
  check_inputs.py --fasta F --depth D --minlen 1500 --markers pfam|busco|table --marker-file M
                  [--marker-file M2 ...] [--seed DIR_or_TABLE]
"""
import argparse
import gzip
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--depth", required=True)
ap.add_argument("--minlen", type=int, default=1500)
ap.add_argument("--markers", required=True, choices=["pfam", "busco", "table"])
ap.add_argument("--marker-file", dest="marker_files", action="append", required=True)
ap.add_argument("--seed")
a = ap.parse_args()

FASTA_EXT = (".fa", ".fasta", ".fna", ".fas")
FASTA_EXT = FASTA_EXT + tuple(e + ".gz" for e in FASTA_EXT)
failed = False


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def ok(msg):
    print(f"  ok       {msg}")


def warn(msg):
    print(f"  WARNING  {msg}")


def fail(msg):
    global failed
    failed = True
    print(f"  ERROR    {msg}")


# ---- assembly ----------------------------------------------------------------------------
order, lengths = [], {}
name, ln = None, 0
with open_text(a.fasta) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                order.append(name); lengths[name] = ln
            name, ln = line[1:].split()[0] if line[1:].split() else "", 0
        else:
            ln += len(line.strip())
    if name is not None:
        order.append(name); lengths[name] = ln
kept = {c for c in order if lengths[c] >= a.minlen}
if len(lengths) != len(order):
    fail(f"assembly: {len(order) - len(lengths)} duplicate contig names")
elif not kept:
    fail(f"assembly: no contig of at least {a.minlen} bp")
else:
    ok(f"assembly: {len(order):,} contigs, {len(kept):,} of at least {a.minlen} bp")


def resolve(n):
    if n in lengths:
        return n
    short = n.split("_LN:", 1)[0]
    return short if short in lengths else None


# ---- depth -------------------------------------------------------------------------------
dnames = []
with open_text(a.depth) as f:
    header = f.readline().rstrip("\n").split("\t")
    for line in f:
        if line.strip():
            dnames.append(line.split("\t", 1)[0])
if not header or header[0] != "contigName":
    warn("depth: first column header is not 'contigName' (is this a jgi_summarize_bam_contig_depths file?)")
known = [d for d in dnames if d in lengths]
missing = kept - set(dnames)
pos = {c: i for i, c in enumerate(order)}
out_of_order = any(pos[x] > pos[y] for x, y in zip(known, known[1:]))
if not known:
    fail("depth: none of its contig names occur in the assembly")
elif out_of_order:
    fail("depth: contigs are not in the same order as in the assembly (MetaBAT2 requires that)")
elif missing:
    warn(f"depth: {len(missing):,} contigs of at least {a.minlen} bp are missing from it "
         f"(e.g. {sorted(missing)[0]}); MetaBAT2 will have no coverage for them")
else:
    ok(f"depth: covers all {len(kept):,} contigs, same order as the assembly, "
       f"{max(len(header) - 3, 0) // 2 or 1} sample(s)")

# ---- markers -----------------------------------------------------------------------------
mnames, total = set(), 0
for path in a.marker_files:
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if a.markers == "pfam" and len(p) >= 2:
                total += 1
                if p[1].strip() not in ("{}", ""):
                    mnames.add(p[0])
            elif a.markers == "busco" and len(p) >= 3 and p[1] in ("Complete", "Duplicated"):
                total += 1
                mnames.add(p[2])
            elif a.markers == "table" and not line.startswith("@") and len(p) >= 2 and p[1].strip():
                total += 1
                mnames.add(p[0])
matched = {resolve(n) for n in mnames} - {None}
if not mnames:
    fail(f"markers: no marker hits read from {', '.join(a.marker_files)} "
         f"(is the format right: --markers-{a.markers}?)")
elif not matched:
    fail("markers: none of the contigs carrying markers occur in the assembly (contig names differ)")
else:
    unmatched = len(mnames) - len(matched)
    msg = f"markers: {len(matched):,} contigs carry markers"
    if unmatched > len(mnames) / 2:
        warn(msg + f"; {unmatched:,} more name contigs not in the assembly. Check the names.")
    else:
        ok(msg + (f" ({unmatched:,} names not in the assembly, ignored)" if unmatched else ""))
    if len(matched & kept) < 2:
        warn("markers: fewer than 2 marker-bearing contigs pass the length filter; the operators "
             "will have nothing to work with")

# ---- seed --------------------------------------------------------------------------------
if a.seed:
    snames, nbins = [], 0
    if os.path.isdir(a.seed):
        for fn in sorted(os.listdir(a.seed)):
            if not fn.lower().endswith(FASTA_EXT):
                continue
            if any(w in fn.lower() for w in ("unbinned", "lowdepth", "tooshort")):
                continue
            nbins += 1
            with open_text(os.path.join(a.seed, fn)) as f:
                snames += [l[1:].split()[0] for l in f if l.startswith(">") and l[1:].split()]
    else:
        bins = set()
        with open_text(a.seed) as f:
            for line in f:
                if line.startswith(("@", "#")) or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    snames.append(p[0]); bins.add(p[1])
        nbins = len(bins)
    res = [resolve(n) for n in snames]
    matched = [r for r in res if r is not None]
    in_kept = [r for r in matched if r in kept]
    if not snames:
        fail("seed: no bins or contigs found")
    elif not matched:
        fail("seed: none of its contig names occur in the assembly")
    else:
        msg = (f"seed: {nbins:,} bins, {len(snames):,} contigs; {len(in_kept):,} pass the "
               f"length filter")
        if len(matched) < len(snames) / 2:
            warn(msg + f"; {len(snames) - len(matched):,} names not in the assembly. Check the names.")
        else:
            ok(msg)

sys.exit(1 if failed else 0)
