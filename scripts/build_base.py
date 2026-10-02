#!/usr/bin/env python3
"""Base facts for any assembly: the contig dictionary, the lengths and the marker facts.

One builder for every dataset. The only thing that differs between dataset families is the
FORMAT the markers arrive in, which is an adapter here, not a separate pipeline:

  --markers pfam    marker_gene_stats.tsv, the column format written by CheckM (marker
                    accessions per gene, per contig).
  --markers busco   one or more BUSCO full_table.tsv (one per lineage dataset, e.g. bacteria
                    and archaea); hits are merged and deduplicated.
  --markers table   one or more plain tables, one <contig> <TAB> <marker> pair per line (lines
                    starting with '#' or '@' are skipped; extra columns are ignored), for
                    annotations from any other tool.

Everything downstream reads the integer fact space written here, so the two families go
through exactly the same builders from this point on.

Emits (into OUTDIR):
  dict.tsv        Clingo_ID<TAB>Original_ID<TAB>Length   (integer id <-> contig name; contigs >= MINLEN)
  mg.lp           mg(ContigInt, MarkerInt).              (marker presence, deduplicated per contig)
  markers_map.tsv MarkerInt<TAB>Marker                   (marker id legend)

Contig ids are assigned in FASTA order, 0-based, to contigs >= MINLEN (matching MetaBAT2's
own -m filter, so the seed, the native composite and the dictionary share one contig set).

Marker ids are assigned 0-based in first-seen order and carry no meaning; the encodings only
ever ask whether two contigs share one.

Usage:
Any input file may be gzip-compressed (.gz).

  build_base.py --fasta assembly.fasta --markers pfam  --marker-file marker_gene_stats.tsv \
                --out FACTS_DIR [--minlen 1500]
  build_base.py --fasta assembly.fasta --markers busco --marker-file bac/full_table.tsv \
                --marker-file arc/full_table.tsv --out FACTS_DIR
"""
import argparse
import ast
import gzip
import os

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--markers", required=True, choices=["pfam", "busco", "table"],
                help="format of the marker files (pfam = marker_gene_stats.tsv, "
                     "busco = full_table.tsv, table = contig<TAB>marker)")
ap.add_argument("--marker-file", dest="marker_files", action="append", required=True,
                help="repeat once per file; BUSCO needs one per lineage")
ap.add_argument("--out", required=True, help="output directory for the fact files")
ap.add_argument("--minlen", type=int, default=1500)
a = ap.parse_args()

os.makedirs(a.out, exist_ok=True)

# BUSCO statuses that count as marker PRESENCE. Duplicated is kept on purpose: it is the
# contamination signal the guard reads.
BUSCO_KEEP = {"Complete", "Duplicated"}


def open_text(path):
    """Open a text file for reading, gzip-compressed (.gz) or not."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fasta_lengths(path):
    """Stream (name, length); name = first whitespace token of the header."""
    header, ln = None, 0
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, ln
                header = line[1:].split()[0]
                ln = 0
            else:
                ln += len(line.strip())
        if header is not None:
            yield header, ln


# ---- 1) contig -> integer id (FASTA order, >= MINLEN), dict.tsv -----------------------
name2id = {}
dict_path = os.path.join(a.out, "dict.tsv")
kept = total = 0
with open(dict_path, "w") as fd:
    fd.write("Clingo_ID\tOriginal_ID\tLength\n")
    for name, length in fasta_lengths(a.fasta):
        total += 1
        if length < a.minlen:
            continue
        cid = kept
        name2id[name] = cid
        fd.write(f"{cid}\t{name}\t{length}\n")
        kept += 1
print(f"[+] {dict_path}: {kept}/{total} contigs kept (>= {a.minlen} bp)")


def resolve(name):
    """Contig name -> clingo id, or None if the contig was filtered out or is absent.

    Some assemblies store the FULL FASTA header with spaces turned into underscores, e.g.
    'edge_10028_LN:i:1626_RC:i:28_XC:f:...'; the contig id is the part before the first
    '_LN:' annotation. Tried only when the exact name misses.
    """
    cid = name2id.get(name)
    if cid is None:
        cid = name2id.get(name.split("_LN:", 1)[0])
    return cid


# ---- 2) markers -> mg/2 ----------------------------------------------------------------
marker2id = {}
pairs = set()                       # (contig_int, marker_int), deduplicated
parse_fail = 0

if a.markers == "pfam":
    # One contig per line: <contig>\t{'<orf>': {'<PFAM>': [[s,e], ...], ...}, ...}  ('{}' = none)
    for path in a.marker_files:
        for line in open_text(path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            cid = resolve(parts[0])
            if cid is None:         # contig below MINLEN or not in the assembly
                continue
            try:
                d = ast.literal_eval(parts[1])
            except (ValueError, SyntaxError):
                parse_fail += 1
                continue
            if not d:
                continue
            for orf, mm in d.items():           # orf -> {PFAM: [[ranges]]}
                if not isinstance(mm, dict):
                    continue
                for pf in mm:
                    mid = marker2id.setdefault(pf, len(marker2id))
                    pairs.add((cid, mid))
elif a.markers == "table":
    # plain table: <contig> <TAB> <marker> [extra columns ignored]
    for path in a.marker_files:
        for line in open_text(path):
            if line.startswith(("#", "@")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1].strip():
                parse_fail += 1
                continue
            cid = resolve(p[0])
            if cid is None:
                continue
            mid = marker2id.setdefault(p[1].strip(), len(marker2id))
            pairs.add((cid, mid))
else:
    # BUSCO full_table.tsv: Busco_id, Status, Sequence, ...
    for path in a.marker_files:
        if not os.path.exists(path):
            print(f"[!] missing {path}, skipped")
            continue
        for line in open_text(path):
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[1] not in BUSCO_KEEP:
                continue
            cid = resolve(p[2])
            if cid is None:
                continue
            mid = marker2id.setdefault(p[0], len(marker2id))
            pairs.add((cid, mid))

mg_path = os.path.join(a.out, "mg.lp")
with open(mg_path, "w") as fm:
    fm.write(f"% marker presence: mg(ContigInt, MarkerInt). Source format: {a.markers}.\n")
    for c, g in sorted(pairs):
        fm.write(f"mg({c},{g}).\n")

map_path = os.path.join(a.out, "markers_map.tsv")
with open(map_path, "w") as fx:
    fx.write("MarkerInt\tMarker\n")
    for name, mid in sorted(marker2id.items(), key=lambda kv: kv[1]):
        fx.write(f"{mid}\t{name}\n")

marker_contigs = len({c for c, _ in pairs})
print(f"[+] {mg_path}: {len(pairs)} mg facts over {marker_contigs} marker-bearing contigs, "
      f"{len(marker2id)} distinct markers"
      + (f"  ({parse_fail} unparseable marker rows skipped)" if parse_fail else ""))
print(f"[+] {map_path}")
