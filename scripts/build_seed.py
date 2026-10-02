#!/usr/bin/env python3
"""Read a seed binning, from ANY binner, into the facts the operators refine.

The seed can be given in either of two forms:

  a folder of bin FASTA files   one file per bin (.fa .fasta .fna .fas, optionally .gz), any names; the bin's
                                name is the file name without its extension. Files whose name
                                contains 'unbinned', 'lowDepth' or 'tooShort' (MetaBAT2's leftover
                                files) are skipped: those contigs are NOT binned.
  a table                       one contig per line, tab-separated: <contig> <TAB> <bin> [extra
                                columns ignored]. Lines starting with '@' or '#' are skipped, so
                                a bioboxes binning file works as it is. May be gzip-compressed (.gz).

Contig names must match the FASTA headers of the assembly (first word of the header). Contigs
the length filter removed are ignored; contigs of the assembly that the seed does not bin become
'recruitable', the pool the recruitment operator draws from.

Emits:
  <out_pred>         the seed as an AMBER/bioboxes binning (binned contigs only).
  <out_seed>         seed(BinInt, ContigInt).   (bins numbered 1.. in sorted name order)
  <out_recruitable>  recruitable(ContigInt).

Usage:
  build_seed.py <dict.tsv> <seed_dir_or_table> <out_pred.tsv> <out_seed.lp> <out_recruitable.lp> <sample_id>
"""
import gzip
import os
import sys

if len(sys.argv) != 7:
    sys.exit(__doc__)
DICT, SEED, OUT_PRED, OUT_SEED, OUT_RECR, SAMPLE = sys.argv[1:7]

FASTA_EXT = (".fa", ".fasta", ".fna", ".fas")
FASTA_EXT = FASTA_EXT + tuple(e + ".gz" for e in FASTA_EXT)
SKIP_WORDS = ("unbinned", "lowdepth", "tooshort")


def open_text(path):
    """Open a text file for reading, gzip-compressed (.gz) or not."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

name2id, name2len = {}, {}
with open(DICT) as f:
    next(f, None)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            name2id[p[1]] = int(p[0])
            name2len[p[1]] = int(p[2])


def resolve(name):
    """Seed contig name -> dictionary name, or None. Same rule as build_base.py."""
    if name in name2id:
        return name
    short = name.split("_LN:", 1)[0]
    return short if short in name2id else None


# ---- read the seed: list of (contig name as given, bin label) ------------------------------
pairs = []
if os.path.isdir(SEED):
    files = sorted(fn for fn in os.listdir(SEED) if fn.lower().endswith(FASTA_EXT))
    skipped_files = [fn for fn in files if any(w in fn.lower() for w in SKIP_WORDS)]
    files = [fn for fn in files if fn not in skipped_files]
    if not files:
        sys.exit(f"ERROR: no bin FASTA files ({', '.join(FASTA_EXT)}) in {SEED}")
    for fn in files:
        label = os.path.splitext(fn[:-3] if fn.lower().endswith(".gz") else fn)[0]
        with open_text(os.path.join(SEED, fn)) as fh:
            for line in fh:
                if line.startswith(">"):
                    pairs.append((line[1:].split()[0], label))
    if skipped_files:
        print(f"[i] skipped {len(skipped_files)} leftover file(s): {', '.join(skipped_files)}")
elif os.path.isfile(SEED):
    for line in open_text(SEED):
        if line.startswith(("@", "#")) or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            sys.exit(f"ERROR: {SEED}: expected <contig><TAB><bin>, got: {line.rstrip()}")
        pairs.append((p[0], p[1]))
else:
    sys.exit(f"ERROR: seed not found: {SEED}")

# ---- map onto the dictionary ---------------------------------------------------------------
binned = {}                       # dictionary contig name -> bin label
unknown, conflicts = 0, []
for c, label in pairs:
    name = resolve(c)
    if name is None:              # filtered by length, or a name the assembly does not have
        unknown += 1
        continue
    if name in binned and binned[name] != label:
        conflicts.append((name, binned[name], label))
    binned[name] = label

if conflicts:
    c, a, b = conflicts[0]
    sys.exit(f"ERROR: {len(conflicts)} contig(s) are in more than one bin (e.g. {c} in {a} and {b}). "
             f"The operators need each contig in at most one bin.")
if not binned:
    sys.exit(f"ERROR: none of the {len(pairs)} seed contigs matches a contig of the assembly. "
             f"Check that the seed uses the same contig names as the FASTA headers.")
if unknown:
    print(f"[i] {unknown} of {len(pairs)} seed contigs are not in the filtered assembly "
          f"(shorter than the length filter, or unknown names); they are ignored")
    if unknown > len(pairs) / 2:
        print(f"[!] WARNING: more than half of the seed's contigs did not match. If the seed was "
              f"made on this assembly, its contig names probably differ from the FASTA headers.")

labels = sorted(set(binned.values()))
lab2int = {b: i + 1 for i, b in enumerate(labels)}

os.makedirs(os.path.dirname(OUT_PRED) or ".", exist_ok=True)
with open(OUT_PRED, "w") as fp, open(OUT_SEED, "w") as fs:
    fp.write(f"@Version:0.9.0\n@SampleID:{SAMPLE}\n@@SEQUENCEID\tBINID\tLENGTH\n")
    fs.write(f"% seed(Bin,Contig) for {SAMPLE}\n")
    for c, b in binned.items():
        fp.write(f"{c}\t{b}\t{name2len[c]}\n")
        fs.write(f"seed({lab2int[b]},{name2id[c]}).\n")

recruitable = [c for c in name2id if c not in binned]
with open(OUT_RECR, "w") as fr:
    fr.write("% recruitable(Contig): filtered contigs the seed left unbinned\n")
    for c in recruitable:
        fr.write(f"recruitable({name2id[c]}).\n")

print(f"[+] seed: {len(binned)} binned contigs in {len(labels)} bins; "
      f"{len(recruitable)} unbinned (recruitable); {len(name2id)} contigs after the length filter")
