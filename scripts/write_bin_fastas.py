#!/usr/bin/env python3
"""Write binnings as folders of FASTA files, one file per bin.

Most downstream tools take bins as separate FASTA files. This reads the assembly once (plain or
.gz) and, for each binning given, writes <out_dir>/<bin>.fa holding that bin's contigs. An
existing <out_dir> is replaced.

Usage:
  write_bin_fastas.py <assembly.fasta[.gz]> <binning.tsv> <out_dir> [<binning.tsv> <out_dir> ...]
"""
import gzip
import os
import re
import shutil
import sys

if len(sys.argv) < 4 or len(sys.argv) % 2 != 0:
    sys.exit(__doc__)
FASTA = sys.argv[1]
jobs = list(zip(sys.argv[2::2], sys.argv[3::2]))


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def safe(label):
    return re.sub(r"[^A-Za-z0-9._-]", "_", label)


# contig -> list of (job index, bin)
where = {}
for j, (pred, out_dir) in enumerate(jobs):
    with open_text(pred) as f:
        for line in f:
            if line.startswith("@") or not line.strip():
                continue
            c, b = line.rstrip("\n").split("\t")[:2]
            where.setdefault(c, []).append((j, b))
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

written = [0] * len(jobs)
bins_seen = [set() for _ in jobs]


def flush(name, lines):
    """Append one contig's record to every bin file it belongs to (one file open at a time,
    so the number of bins is not limited by the open-file limit)."""
    for j, b in where.get(name, []):
        with open(os.path.join(jobs[j][1], safe(b) + ".fa"), "a") as h:
            h.writelines(lines)
        written[j] += 1
        bins_seen[j].add(b)


name, lines = None, []
with open_text(FASTA) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                flush(name, lines)
            name = line[1:].split()[0] if line[1:].split() else ""
            lines = [line] if name in where else []
        elif lines:
            lines.append(line)
if name is not None:
    flush(name, lines)

for j, (pred, out_dir) in enumerate(jobs):
    print(f"[+] {out_dir}: {len(bins_seen[j])} bins, {written[j]} contigs")
