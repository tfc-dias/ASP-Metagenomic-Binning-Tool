# Modified MetaBAT2

`metabat2_asp.patch` modifies a single MetaBAT2 source file, `src/metabat2.cpp`. It adds operating
modes in which MetaBAT2 computes its pairwise contig features, writes them to a file and exits
before binning. No other part of MetaBAT2 is changed, and its binning algorithm is unaffected.

It is used only with `refine.sh --backend metabat2`. By default, `refine.sh` computes the same
features with `scripts/mb2_features.py`, a Python port of MetaBAT2's functions that was validated
against this patch. To use this backend, set `METABAT2_MODIFIED` in `config.sh` to the binary built
below.

## Installation

The patch applies to MetaBAT2 at commit `c869c524d0f131d60a03be64bd26b89738160652` (May 2026,
shortly after the 2.18 release). It does not compile against the 2.18 release itself, which lacks
a declaration used by the modified file.

```
git clone https://bitbucket.org/berkeleylab/metabat.git metabat-asp
cd metabat-asp
git checkout c869c524d0f131d60a03be64bd26b89738160652
git apply /path/to/metabat2_modified/metabat2_asp.patch
mkdir build && cd build && cmake .. && make
```

The resulting binary is `build/src/metabat2`.

MetaBAT2's own build requirements apply: Boost, cmake, a C++17 compiler, and autoconf and automake
(htslib is downloaded during the build). Further details are given in MetaBAT2's `README.md` and
`INSTALL.md`.

## Operating modes

The mode is selected
by an environment variable. In every mode, MetaBAT2 reads its standard inputs (`-i` assembly, `-a`
depth, `-m` minimum contig length, 1500 in the pipeline), computes its features, writes one fact
file to the working directory and exits before binning.

| Mode | Invoked by | Output |
|---|---|---|
| Default (no variable set) | `refine.sh --backend metabat2`, once per run | `composite_raw.lp`: MetaBAT2's native similarity score, combining tetranucleotide frequency and coverage, for every pair of contigs in its similarity graph. All affinity facts used by the operators are derived from this file. |
| `ASP_DUMP_PAIRS=<pair file>` | `scripts/build_comarker.py`, once per operator input | `comarker_raw.lp`: MetaBAT2's tetranucleotide and abundance distances for the supplied contig pairs only (the pairs that share a marker gene), as `D = tnf×100 + abd×100`, on the range 0 (identical) to 200 (maximally different). These distances are compared with the dissimilarity threshold. |
| `ASP_DUMP_RAW=1` | Not used by the pipeline | The tetranucleotide and abundance distances for all pairs of contigs. Its cost is quadratic in the number of contigs, which makes it impractical for large assemblies. |

The `-o` argument is supplied only because MetaBAT2 requires it; no output is written to that
location. Pairs involving a contig that MetaBAT2 did not load, because it is shorter than the
minimum length or has insufficient coverage, receive no distance.

