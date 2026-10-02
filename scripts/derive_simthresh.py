#!/usr/bin/env python3
"""Derive the co-marker distance threshold from the seed binning, without ground truth.

The operators call two contigs that carry the same marker "dissimilar" when their distance
reaches simthresh. Telling same-genome from different-genome pairs directly would need a
ground truth. Instead, the seed binning stands in for one: pairs it places in the SAME bin are
probably same-genome, pairs it places in DIFFERENT bins probably not. The threshold is the cut
in 0..200 that disagrees least with that grouping (ties go to the lower cut).

Values derived on the dissertation's datasets: CAMI TOY LOW 74, CAMI TOY MEDIUM 36, strong100
25, and 17-18 on three real wastewater and soil assemblies. Where a fixed 50 was also tried, the
operators' output was unchanged or nearly so.

KNOWN LIMITATION. This always returns a number, including on datasets where the two
populations overlap entirely and no threshold can separate them (strong100). Two indicators
were tried to detect that case without a ground truth and both failed.
The consequence is benign in direction: where the signal is absent the guard becomes
over-strict, so the operators do LESS than they could, never something wrong.

Usage:
  derive_simthresh.py <comarker_within.lp> <comarker_merge.lp> [--verbose]

  within : co-marker pairs the seed puts in the SAME bin      (stand-in for same-genome)
  merge  : co-marker pairs the seed puts in DIFFERENT bins    (stand-in for cross-genome)
Prints the threshold on stdout, so it can be used as:  -c simthresh=$(derive_simthresh.py ...)
"""
import re
import sys

SCALE_MAX = 200          # native co-marker distance is tnf*100 + abd*100
FALLBACK = 50            # the historical carried value, used only if derivation is impossible


def load(path):
    out = []
    rx = re.compile(r"comarker_sim\(\s*\d+\s*,\s*\d+\s*,\s*(\d+)\s*\)")
    for line in open(path):
        m = rx.search(line)
        if m:
            out.append(int(m.group(1)))
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    verbose = "--verbose" in sys.argv
    if len(args) != 2:
        sys.exit(__doc__)
    same, cross = load(args[0]), load(args[1])

    # The co-marker distance is tnf*100 + abd*100, so it cannot exceed 200. Anything above that
    # means the file was not built by build_comarker.py, and a threshold derived on another
    # scale would be meaningless for the encodings, so refuse.
    over = [v for v in same + cross if v > SCALE_MAX]
    if over:
        sys.exit(f"[derive_simthresh] ERROR: {len(over)} distances exceed {SCALE_MAX} "
                 f"(max {max(over)}), so they are not on the 0..{SCALE_MAX} scale of MetaBAT2's "
                 f"distances. Rebuild them with build_comarker.py.")

    if not same or not cross:
        # One population empty: there is nothing to calibrate against, so fall back to the
        # historical value and say so loudly. This is a deliberate guard, not a silent default.
        #
        # What an empty same-bin population MEANS: the seed binner placed no two marker-sharing
        # contigs in the same bin at all. Separation consumes exactly that file, so it would
        # abstain regardless of any threshold. Merge and recruitment read different files and
        # would still run, but without local calibration. It did not occur on any dataset in the
        # dissertation.
        print(FALLBACK)
        print("=" * 78, file=sys.stderr)
        print(f"[derive_simthresh] FALLBACK IN USE: simthresh = {FALLBACK}, NOT derived.",
              file=sys.stderr)
        print(f"[derive_simthresh]   same-bin pairs: {len(same)}   cross-bin pairs: {len(cross)}",
              file=sys.stderr)
        print("[derive_simthresh]   One population is empty, so there is nothing to calibrate",
              file=sys.stderr)
        print("[derive_simthresh]   against. 50 is the value read on the CAMI TOY LOW benchmark,",
              file=sys.stderr)
        print("[derive_simthresh]   carried over. Report it as such, not as derived.",
              file=sys.stderr)
        print("=" * 78, file=sys.stderr)
        return

    # the cut minimising disagreement with the seed's own grouping
    best_err, best_t = min(
        (sum(1 for x in same if x >= t) + sum(1 for x in cross if x < t), t)
        for t in range(0, SCALE_MAX + 1))
    print(best_t)

    fp = sum(1 for x in same if x >= best_t)
    fn = sum(1 for x in cross if x < best_t)
    print(f"[derive_simthresh] same-bin pairs {len(same)}, cross-bin pairs {len(cross)} "
          f"-> simthresh={best_t}", file=sys.stderr)
    print(f"[derive_simthresh]   same-bin pairs above the cut: {fp} ({100*fp/len(same):.1f} %); "
          f"cross-bin below: {fn} ({100*fn/len(cross):.2f} %)", file=sys.stderr)
    if verbose:
        s = sorted(same); c = sorted(cross)
        print(f"[derive_simthresh]   median same-bin {s[len(s)//2]}, "
              f"median cross-bin {c[len(c)//2]}", file=sys.stderr)


if __name__ == "__main__":
    main()
