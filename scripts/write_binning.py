#!/usr/bin/env python3
"""Write an operator's result as a binning file.

Reads clingo's output, keeps the last (best) answer, and turns its att(Bin,Contig) atoms
(integer ids) back into the original contig names, using the dictionary from build_base.py.
The file is in the bioboxes binning format (@ headers, then contig <TAB> bin); --sample sets its
@SampleID header (default SAMPLE; refine.sh passes --name).

Usage:
  write_binning.py <dict.tsv> --pred <clingo_output.txt> <out.tsv> [--sample ID]
"""
import re
import sys


def load_dict(path):
    """Clingo_ID -> (Original_ID, Length). Skips the header row."""
    mapping = {}
    with open(path) as f:
        next(f, None)  # header: Clingo_ID  Original_ID  Length
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                mapping[parts[0]] = (parts[1], parts[2])
    return mapping


def last_answer_block(text):
    # clingo prints one "Answer:" block per improved model; keep only the final.
    idx = text.rfind("Answer:")
    return text[idx:] if idx != -1 else text


def att_pairs(text):
    return re.findall(r"att\(\s*(\d+)\s*,\s*(\d+)\s*\)", text)


def main():
    argv = sys.argv[1:]
    sample_id = None
    if "--sample" in argv:
        i = argv.index("--sample")
        sample_id = argv[i + 1]
        del argv[i:i + 2]
    dict_file = argv[0]
    mode = argv[1]
    att_file, out_tsv = argv[2], argv[3]
    if sample_id is None:
        sample_id = "SAMPLE"
    mapping = load_dict(dict_file)

    with open(att_file) as f:
        text = f.read()
    if mode != "--pred":
        sys.exit(f"unknown mode {mode!r} (use --pred)")
    text = last_answer_block(text)

    pairs = sorted(set(att_pairs(text)), key=lambda p: (int(p[0]), int(p[1])))
    missing = 0
    with open(out_tsv, "w") as out:
        out.write(f"@Version:0.9.0\n@SampleID:{sample_id}\n")
        out.write("@@SEQUENCEID\tBINID\n")
        written = 0
        for bin_id, contig in pairs:
            if contig not in mapping:
                missing += 1
                continue
            name, _ = mapping[contig]
            out.write(f"{name}\tBin_{bin_id}\n")
            written += 1
    note = f" ({missing} contigs not in dictionary, skipped)" if missing else ""
    print(f"[+] {out_tsv}: {written} rows{note}")


if __name__ == "__main__":
    main()
