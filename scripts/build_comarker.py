#!/usr/bin/env python3
"""comarker_sim(C1,C2,D) from MetaBAT2's own distance functions.

For every pair of contigs that share a marker gene, the distance MetaBAT2 itself puts between
them. This is the contamination signal all three operators read: a marker shared by SIMILAR
contigs is a duplication inside one genome; shared by DISSIMILAR contigs, it is two genomes in
one bin.

The distance is MetaBAT2's cal_tnf_dist + cal_abd_dist on exactly the pairs asked for, computed
either by mb2_features.py (a Python port of those functions; --features) or by the modified
MetaBAT2 binary in its ASP_DUMP_PAIRS mode (--metabat2). D = tnf_prob*100 + abd_prob*100, so it
lies in 0..200 on every dataset, the scale the threshold is derived on.

Modes (which pairs): within = same-bin pairs of --pred (separation); merge = different-bin pairs
of --pred (merge guard, and the threshold derivation); guard = each --recruitable contig against
every contig sharing a marker with it (recruitment); all = every pair (not used by refine.sh).

--scale is accepted for compatibility and ignored: the scale is fixed at 0..200.

Usage:
  build_comarker.py --mode {within,merge,guard,all} --dict D --mg MG
                    (--features FEATURES.npz | --metabat2 BIN --fasta F --depth DP [--minlen N])
                    (--pred PRED.tsv | --recruitable RECR.lp) --out OUT.lp
  (FEATURES.npz is written by `mb2_features.py features`.)
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--mode", required=True, choices=["within", "merge", "guard", "all"])
ap.add_argument("--dict", required=True)
ap.add_argument("--mg", required=True)
ap.add_argument("--features", help="features.npz from mb2_features.py (Python backend)")
ap.add_argument("--metabat2", help="the modified MetaBAT2 binary (MetaBAT2 backend)")
ap.add_argument("--fasta", help="assembly (MetaBAT2 backend)")
ap.add_argument("--depth", help="depth table (MetaBAT2 backend)")
ap.add_argument("--minlen", type=int, default=1500, help="minimum contig length (MetaBAT2 backend)")
ap.add_argument("--pred")
ap.add_argument("--recruitable")
ap.add_argument("--out", required=True)
ap.add_argument("--scale", type=int, default=1000, help="accepted for compatibility, IGNORED")
ap.add_argument("--miss-fill", dest="miss_fill", type=int, default=None,
                help="D for pairs that remain unmeasurable even with the coverage filter "
                     "relaxed. Default: none, such pairs are skipped in every mode. See the "
                     "note by MISS_FILL below.")
ap.add_argument("--relax-cv", action="store_true",
                help="Measure at minCV=minCVSum=0 so contigs MetaBAT2 excluded for low "
                     "coverage still get a distance. OFF by default: minCV is also the "
                     "signal gate inside cal_abd_dist, so relaxing it changes the distance "
                     "of pairs MetaBAT2 kept anyway (on a multi-sample dataset it moved 86%% "
                     "of them). The run refuses to write if any measured pair changes.")
a = ap.parse_args()
if bool(a.features) == bool(a.metabat2):
    ap.error("give either --features (Python backend) or --metabat2 (MetaBAT2 backend)")
if a.metabat2 and not (a.fasta and a.depth):
    ap.error("--metabat2 needs --fasta and --depth")

if a.scale != 1000:
    print(f"[{a.mode}] NOTE: --scale {a.scale} ignored; native D is fixed on 0..200.",
          file=sys.stderr)

# ---- dict: name <-> int ----
id2name, name2id = {}, {}
with open(a.dict) as f:
    next(f, None)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            id2name[int(p[0])] = p[1]
            name2id[p[1]] = int(p[0])

# ---- marker -> contigs carrying it ----
contigs_of_marker = defaultdict(list)
n_mg = 0
for line in open(a.mg):
    # tolerant of whitespace and CRLF: the 30-genome facts are "mg(4642, 1).\r\n" while the
    # others are "mg(12478,1).\n". A strict re.match matches NOTHING on the former and
    # silently yields zero facts.
    m = re.search(r"mg\(\s*(\d+)\s*,\s*(\d+)\s*\)", line)
    if m:
        contigs_of_marker[int(m.group(2))].append(int(m.group(1)))
        n_mg += 1
if n_mg == 0:
    sys.exit(f"[{a.mode}] parsed 0 mg/2 facts from {a.mg} — refusing to write an empty guard. "
             f"An empty file would make the contamination guard silently allow everything, "
             f"so this is an error, not an empty result.")

# ---- select the pairs for this mode ----
pairs = []                        # list of (c1, c2); guard keeps c1 = recruitable
if a.mode == "all":
    # every co-marker pair, no bin structure assumed (not used by the refinement pipeline;
    # kept for whole-partition models, which start with no seed partition to filter against).
    seen = set()
    for g, cs in contigs_of_marker.items():
        cs = sorted(cs)
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                key = (cs[i], cs[j])
                if key not in seen:
                    seen.add(key)
                    pairs.append(key)
elif a.mode in ("within", "merge"):
    if not a.pred:
        raise SystemExit(f"--mode {a.mode} needs --pred <binning.tsv>")
    cid_bin = {}
    for line in open(a.pred):
        if line.startswith("@") or not line.strip():
            continue
        c, b = line.rstrip("\n").split("\t")[:2]
        if c in name2id:
            cid_bin[name2id[c]] = b
    same = (a.mode == "within")
    seen = set()
    for g, cs in contigs_of_marker.items():
        cs = [c for c in cs if c in cid_bin]
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                c1, c2 = cs[i], cs[j]
                if (cid_bin[c1] == cid_bin[c2]) == same:
                    key = (min(c1, c2), max(c1, c2))
                    if key not in seen:
                        seen.add(key)
                        pairs.append(key)
else:  # guard: recruitable-anchored, DIRECTIONAL (c1 must stay the recruitable endpoint)
    if not a.recruitable:
        raise SystemExit("--mode guard needs --recruitable <recruitable.lp>")
    recruitable = set()
    for line in open(a.recruitable):
        m = re.match(r"recruitable\((\d+)\)", line)
        if m:
            recruitable.add(int(m.group(1)))
    seen = set()
    for g, cs in contigs_of_marker.items():
        rec_here = [c for c in cs if c in recruitable]
        for c1 in rec_here:
            for c2 in cs:
                if c1 != c2 and (c1, c2) not in seen:
                    seen.add((c1, c2))
                    pairs.append((c1, c2))

print(f"[{a.mode}] co-marker pairs: {len(pairs)} (from {n_mg} mg facts)")
if not pairs:
    if a.mode in ("merge", "guard"):
        print(f"[{a.mode}] WARNING: 0 pairs from {n_mg} mg facts. An empty file makes the "
              f"contamination guard VACUOUS — it will silently permit everything. Check that "
              f"--pred/--recruitable name the same contigs as --dict.", file=sys.stderr)
    open(a.out, "w").write(
        f"% {a.mode} comarker_sim: no pairs\n")
    print(f"[+] {a.out}: 0 facts")
    raise SystemExit(0)

# ---- MetaBAT2's distance for exactly those pairs ----
name_pairs = [(id2name[c1], id2name[c2]) for c1, c2 in pairs]

# MetaBAT2 excludes contigs whose total effective mean coverage is below minCVSum
# (metabat2.cpp:599, default 1.0). Those contigs are long enough for the ASP filter and
# carry markers, so the guard needs their distance, but none can be computed for
# a contig that was never loaded.
#
# The catch: minCV is ALSO the gate inside cal_abd_dist (metabat2.cpp:1652) deciding
# whether a sample carries signal for a pair, and minCVSum = max(minCV, minCVSum)
# (line 168), so the contigs cannot be kept without also lowering that gate. In
# principle that changes every pair's distance, not just the recovered ones.
#
# The effect is dataset-dependent. On a single-sample assembly relaxing recovered 230 of
# 1,117 pairs and changed none of the others; on a four-sample assembly it recovered none
# and changed 86 % of them, because low-coverage samples then count as informative and the
# averaged abundance distance moves for almost every pair. So relaxing is NOT safe in general. Default: run at MetaBAT2's own settings and let the
# MISS_FILL policy below cover whatever it could not measure. --relax-cv is opt-in and
# REFUSES to write if the invariance does not hold.
def measure(relaxed):
    """{(name1, name2): D} for the pairs MetaBAT2 can measure, at its own coverage settings or,
    if relaxed, at minCV = minCVSum = 0."""
    if a.features:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from mb2_features import pair_distances
        cv = 0.0 if relaxed else 1.0
        return pair_distances(a.features, name_pairs, min_cv=cv, min_cv_sum=cv)
    # the modified binary writes comarker_raw.lp into its working directory, so it runs in a
    # temporary one and every input path is absolute
    rx = re.compile(r'comarker_sim\("([^"]+)","([^"]+)",(\d+)\)')
    with tempfile.TemporaryDirectory() as td:
        plist = os.path.join(td, "pairs.tsv")
        with open(plist, "w") as f:
            for n1, n2 in name_pairs:
                f.write(f"{n1}\t{n2}\n")
        cmd = [os.path.abspath(a.metabat2),
               "-i", os.path.abspath(a.fasta), "-a", os.path.abspath(a.depth),
               "-o", os.path.join(td, "out", "bin"), "-m", str(a.minlen), "--seed", "1", "-q"]
        if relaxed:
            cmd += ["-x", "0", "--minCVSum", "0"]
        r = subprocess.run(cmd, cwd=td, env=dict(os.environ, ASP_DUMP_PAIRS=plist),
                           capture_output=True, text=True)
        raw = os.path.join(td, "comarker_raw.lp")
        if not os.path.exists(raw):
            sys.exit(f"[{a.mode}] ASP_DUMP_PAIRS produced no comarker_raw.lp\n"
                     f"  cmd: {' '.join(cmd)}\n  stderr tail:\n{r.stderr[-2000:]}")
        out = {}
        for line in open(raw):
            m = rx.search(line)
            if m:
                out[(m.group(1), m.group(2))] = int(m.group(3))
        return out


strict = measure(False)
if not a.relax_cv:
    dist = strict
else:
    dist = measure(True)
    changed = [k for k in set(strict) & set(dist) if strict[k] != dist[k]]
    recovered = len(set(dist) - set(strict))
    print(f"[{a.mode}] --relax-cv: {recovered} pairs recovered, {len(changed)} of the "
          f"{len(set(strict) & set(dist))} MetaBAT2 kept anyway CHANGED")
    if changed:
        sys.exit(f"[{a.mode}] REFUSING to write: relaxing minCV moved {len(changed)} pairs "
                 f"that MetaBAT2 measured at its own settings, so these are no longer "
                 f"MetaBAT2's numbers. Drop --relax-cv and let --miss-fill cover the gap.")

# ---- policy for pairs MetaBAT2 could not measure ----
# MetaBAT2 excludes contigs whose total effective mean coverage is below minCVSum
# (metabat2.cpp:599), so a few pairs the marker facts name have no distance. They are skipped.
#
# Such a contig is also absent from composite_raw.lp, so it has no affinity edge and
# binning_recruit.lp can never place it. When the seed comes from MetaBAT2 it is not binned
# either, so it is never in a bin and its pairs could not affect any guard, whatever D they got.
#
# With a seed from ANOTHER binner this no longer holds: that binner may have binned the contig.
# Its markers are then invisible to the guards (no distance), which can let a merge or a
# recruitment through that a measured distance would have blocked. --miss-fill 200 makes such
# pairs count as dissimilar instead (the cautious choice). refine.sh passes it for the merge and
# recruitment guards only; without it (the default here) such pairs are skipped.
MISS_FILL = a.miss_fill

out, skipped, filled = [], 0, 0
for c1, c2 in pairs:
    d = dist.get((id2name[c1], id2name[c2]))
    if d is None:                      # MetaBAT2 never loaded this contig (low coverage)
        if MISS_FILL is None:
            skipped += 1
            continue
        d = MISS_FILL
        filled += 1
    out.append((c1, c2, d))

with open(a.out, "w") as f:
    f.write(f"% comarker_sim(C1,C2,D) [{a.mode}]: MetaBAT2-NATIVE cal_tnf_dist + cal_abd_dist,\n"
            f"% D = tnf*100 + abd*100 on 0..200 (low = similar).\n")
    for c1, c2, d in out:
        f.write(f"comarker_sim({c1},{c2},{d}).\n")
with open(a.out + ".stats", "w") as f:       # read by stage_summary.py
    f.write(f"measured\t{len(out) - filled}\nfilled\t{filled}\nskipped\t{skipped}\n")
print(f"[+] {a.out}: {len(out)} comarker_sim facts "
      f"({len(out) - filled} measured, {filled} filled with {MISS_FILL}, {skipped} skipped)")
