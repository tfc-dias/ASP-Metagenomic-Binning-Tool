# ASP refinement of metagenomic binnings

This repository provides a refinement method for metagenomic binnings based on Answer Set
Programming (ASP), implemented with clingo. The method takes an existing binning, referred to as
the seed (produced, for example, by MetaBAT2), and corrects it with
three declarative operators, each addressing one class of binning error:

| Operator | Error addressed | Principle |
|---|---|---|
| **Separation** | Chimeric bins (two genomes in one bin) | Two dissimilar contigs that carry the same single-copy marker gene incur a penalty for remaining in the same bin; a bin is split where this evidence outweighs the affinity that holds it together. |
| **Merge** | Over-split genomes (one genome in several bins) | Two bins that are mutual nearest neighbours are joined, unless a marker gene carried by dissimilar contigs in both bins indicates that they belong to different genomes. |
| **Recruitment** | Unbinned contigs | Each unbinned contig may be assigned to the bin it is most similar to, unless the assignment would place a marker gene on two dissimilar contigs. |

The operators to be applied, and their order, are selectable (see
[Operator selection](#operator-selection)); by default, all three are applied in the order above.
Dissimilarity is decided by a distance threshold derived from the seed binning itself, so the
method requires neither reference genomes nor a ground truth.

## Requirements

| Software | Notes |
|---|---|
| clingo 5.4 or later | Available from conda-forge (`conda install -c conda-forge clingo`). |
| Python 3.10 or later, with numpy and scipy | `pip install numpy scipy`, or from conda-forge. |

Nothing needs to be compiled. The modified MetaBAT2 is optional (see [Backends](#backends)).

## MetaBAT2's features

The operators rely on two quantities defined by MetaBAT2: its pairwise contig similarity, which
combines tetranucleotide frequency and coverage, and its distance between contigs that share a
marker gene. The seed binning determines the initial bins; these quantities determine how the
operators modify them, so they are used whichever binner produced the seed.

`scripts/mb2_features.py` computes both in Python. It is a port of MetaBAT2's own functions
(tetranucleotide frequencies, `cal_tnf_dist`, `cal_abd_dist`, the abundance correlation, the
construction of its similarity graph and its edge score), at the commit
`c869c524d0f131d60a03be64bd26b89738160652`, with MetaBAT2's default settings, and was validated
against MetaBAT2 compiled with the patch in `metabat2_modified/`. `refine.sh` runs it
automatically; it never produces a binning.

### Backends

The features can be computed in two ways, selected with `--backend`:

| `--backend` | Computes the features with | Requires |
|---|---|---|
| `python` (default) | `scripts/mb2_features.py` | numpy and scipy |
| `metabat2` | MetaBAT2 compiled with the patch in `metabat2_modified/`, which exports the same quantities | The modified MetaBAT2, built as described in `metabat2_modified/README.md`, and `METABAT2_MODIFIED` set in `config.sh` |

Both give the same features and, on the assemblies tested, identical binnings. The Python
backend needs no compilation; the MetaBAT2 backend is the one the dissertation's results were
produced with. The backend is recorded in `run_info.txt`, and an output directory cannot be
resumed with a different backend.

## Configuration

No configuration is needed when `clingo` and `python3` are on the system path. Otherwise, their
locations are set in a configuration file created from the template provided:

```
cp config.sh.example config.sh
```

## Usage

```
./refine.sh --name mysample \
            --fasta assembly.fasta \
            --depth depth.txt \
            --markers-pfam marker_gene_stats.tsv \
            --seed my_bins/ \
            --write-bins
```

The option `--check` validates the configuration and the inputs without running the refinement.
The same validation is performed at the start of every run.

### Inputs and options

| Option | Description | Notes |
|---|---|---|
| `--name` | Label for the run. | Also written as `@SampleID` in the output binnings. |
| `--fasta` | Assembly. | FASTA format, plain or gzip-compressed (`.gz`). |
| `--depth` | Coverage depth. | Format of MetaBAT2's `jgi_summarize_bam_contig_depths`, with contigs in the same order as in the assembly. Produced beforehand from the sequencing reads. |
| `--markers-pfam`, `--markers-busco` or `--markers-table` | Single-copy marker gene annotation. | Produced beforehand with the annotation tool of choice; see [Marker gene annotation](#marker-gene-annotation). |
| `--seed` | Seed binning. | Either a directory of bin FASTA files (`.fa`, `.fasta`, `.fna` or `.fas`, plain or gzip-compressed, one file per bin, arbitrary names) or a table with one `contig <TAB> bin` pair per line (lines beginning with `@` or `#` are ignored, so a bioboxes binning file is accepted; it may be gzip-compressed). Required: the seed must be produced beforehand, with any binner. |
| `--operators` | Operators to apply, and their order. | Default `separate,merge,recruit`; see [Operator selection](#operator-selection). |
| `--out` | Output directory. | Default `refine_out/<name>`. |
| `--write-bins` | Also write the final binning as a directory of FASTA files, one per bin. | |
| `--backend` | How MetaBAT2's features are computed: `python` or `metabat2`. | Default `python`; see [Backends](#backends). |
| `--check` | Validate the configuration and the inputs, then stop. | |
| `--version` | Print the version. | |

Contig names in the seed and in the marker file must match the FASTA headers (their first word).

The following environment variables are optional: `SEP_TL`, `MRG_TL` and `REC_TL`, the time limit
of each operator in seconds (default 600, 180 and 900); `THREADS`, the number of solver threads
(default 8); and `MINLEN`, the minimum contig length in base pairs (default 1500).

A run repeated with the same `--out` resumes from where it stopped: completed stages are skipped
and an interrupted stage is recomputed. A run whose seed or parameters differ from those recorded
in `--out` is rejected; a new output directory must then be used. Time limits and the number of
threads may differ between runs.

The method's parameters, and the options that modify them, are described under
[Parameters](#parameters).

### Marker gene annotation

The single-copy marker gene annotation must be produced beforehand, with the annotation tool of
choice, and supplied in one of three formats. The method uses only one piece of information from it:
which marker genes each contig carries. Contigs shorter than the minimum contig length are
ignored, and a marker gene is counted once per contig.

| Option | Format | Content used |
|---|---|---|
| `--markers-pfam` | `marker_gene_stats.tsv`, as written by CheckM: one line per contig, holding the contig name and, for each gene on the contig, the marker accessions it matched. | Every marker accession found on a contig. The marker set is the one CheckM was run with. |
| `--markers-busco` | `full_table.tsv`, as written by BUSCO: one line per BUSCO gene and hit, with the gene identifier, its status and the contig. The option may be repeated to combine several lineage datasets (for example, bacteria and archaea). | Genes with status `Complete` or `Duplicated`. Duplicated genes are retained because they are the occurrences the contamination criterion examines; `Fragmented` and `Missing` entries are ignored. |
| `--markers-table` | A plain table with one `contig <TAB> marker` pair per line, where `marker` is any identifier of a single-copy marker gene. Lines beginning with `#` or `@` are ignored, as are additional columns. The option may be repeated. This format accepts annotations from any tool. | Every pair listed. |

Any of these files may be gzip-compressed, and only one format may be used per run.

### Operator selection

`--operators` takes a comma-separated list of `separate`, `merge` and `recruit`, each appearing at
most once, and applies them in the order given. Each operator is applied to the output of the
previous one. Each stage's output is named after the sequence of operators applied up to it, and
the final stage's output is also given after the minimum bin size filter:

| `--operators` | Outputs |
|---|---|
| `separate,merge,recruit` (default) | Full refinement. Final binning: `pred_separate_merge_recruit_min200k.tsv`. |
| `separate,merge` | Structural refinement only, without recruitment. Final binning: `pred_separate_merge_min200k.tsv`. |
| `recruit` | Recruitment of unbinned contigs only. Final binning: `pred_recruit_min200k.tsv`. |
| `merge` | Merging of over-split genomes only. Final binning: `pred_merge_min200k.tsv`. |

Separation is applied first so that the fragments it
separates can subsequently be merged with the bins to which they belong; recruitment is applied
last so that contigs are assigned to bins whose boundaries are already established. Other orders
are permitted but produce different results.

Several sequences may be computed in the same output directory; stages shared between them are
computed only once (for example, `--operators separate,recruit` after the default sequence reuses
the separation stage). The distance threshold is always derived from the seed, and is therefore
the same for every sequence.

## Outputs

All outputs are written to `refine_out/<name>/`, or to the directory given by `--out`:

| File | Content |
|---|---|
| `pred_separate_merge_recruit_min200k.tsv` | **The final binning**: the output of the last stage after the minimum bin size filter (named after the sequence of operators applied; default sequence shown). |
| `bins_separate_merge_recruit_min200k/` | With `--write-bins`, the final binning as a directory of FASTA files, one `.fa` file per bin. |
| `pred_seed.tsv` | The seed binning, restricted to the contigs the method operates on. |
| `pred_separate.tsv`, `pred_separate_merge.tsv`, `pred_separate_merge_recruit.tsv` | The output of each stage, before the minimum bin size filter, recorded for inspection. |
| `summary.tsv` | The changes made by each stage, and whether the solver proved its result optimal (described below). |
| `run_info.txt` | The versions of the tool, clingo, Python, and numpy and scipy or the modified MetaBAT2 (depending on the backend), the command line and all parameters used. |
| `simthresh_used.txt` | The distance threshold used, and whether it was derived or set manually. |
| `run.log` | The complete log of the run. |
| `work/` | Intermediate facts and solver outputs, with one directory per stage under `work/stages/`. |

The binnings are tab-separated `contig <TAB> bin` files in the bioboxes binning format. Only the
final stage is filtered by bin size; a filtered result after an intermediate stage is obtained by
running the corresponding shorter sequence (for example, `--operators separate,merge` for the
structural refinement), which reuses the stages already computed.

The tool does not evaluate its results. They can be assessed with any suitable method, for
example AMBER against a ground truth, or CheckM or CheckM2 in the absence of one (applied to the
`bins_*` directory).

### Stage summary

`summary.tsv` contains one row per stage:

| Column | Description |
|---|---|
| `stage`, `operator` | The sequence of operators up to this stage (for example, `separate_merge`) and the operator applied at it. |
| `bins_in`, `bins_out`, `contigs_in`, `contigs_out` | Number of bins and of binned contigs before and after the stage. |
| `bins_split` | Input bins whose contigs were assigned to more than one output bin. |
| `bins_joined` | Output bins formed from more than one input bin. |
| `contigs_recruited` | Contigs binned at this stage that were previously unbinned. |
| `bins_kept`, `contigs_kept` | Number of bins and of binned contigs after the minimum bin size filter (final stage only). |
| `unmeasured_as_dissimilar` | Marker gene pairs that MetaBAT2 could not measure, treated as dissimilar (see [Notes and limitations](#notes-and-limitations)). |
| `solver` | `optimal`: the result was proved optimal. `solved`: merge has no objective to optimise, so its result is complete. `time limit`: the best result found before the time limit, not proved optimal. |
| `seconds` | Wall-clock time of the stage. |

## Example: strong100

`strong100` is a small synthetic community (1,150 contigs, 81 genomes) from the GraphMB datasets,
distributed with its own depth file and marker genes. The dataset can be obtained from
[Zenodo record 6122610](https://zenodo.org/records/6122610) (CC-BY-4.0) as `strong100.zip`.

A seed binning must first be produced with a binner of choice. With MetaBAT2, as in the
dissertation:

```
metabat2 -i strong100/assembly.fasta -a strong100/assembly_depth.txt \
         -o strong100_seed/bin -m 1500 --seed 1
```

The seed is then refined as follows:

```
./refine.sh --name strong100 \
            --fasta strong100/assembly.fasta \
            --depth strong100/assembly_depth.txt \
            --markers-pfam strong100/marker_gene_stats.tsv \
            --seed strong100_seed/
```

The same Zenodo record contains the wastewater and soil assemblies (`hjor`, `aale`, `soil` and
others), each with its depth file and marker genes, in the same layout.

## Parameters

| Parameter | Default | Rationale | Option |
|---|---|---|---|
| Distance threshold | Derived from the seed | The cut that best separates contig pairs placed in the same seed bin from pairs placed in different seed bins. If it cannot be derived, the value 50 is used, and this is stated in `simthresh_used.txt` and in the log. | `--simthresh N` (0 to 200) |
| Separation penalty weight | 1 | Penalty and affinity are compared on their own scales. With this setting, separation can fail to perform a split but cannot introduce a spurious one. Higher values make splits easier to obtain; lower values, harder. | `--sepmul N` |
| Merge minimum connectivity | 0 | The minimum number of crossing affinity edges required for two bins to be merged. At 0, the marker gene criterion alone decides merges. | `--minconn N` |
| Minimum bin size | 200,000 bp | MetaBAT2's own minimum. Bins smaller than this are removed from the final binning, whose file name states the value used: `_min200k` by default, `_min100k` for `--min-bin-size 100000`. | `--min-bin-size N` |
| Recruitment tolerance | 1 | A single marker gene conflict with a dissimilar contig prevents an assignment. | Not configurable |
| Unmeasured marker pairs | Distance 200 in merge and recruitment | Marker gene pairs that MetaBAT2 cannot measure are treated as dissimilar in these operators (see [Notes and limitations](#notes-and-limitations)). | Not configurable |

The options are intended for sensitivity analyses and benchmarking rather than routine use. Their
use is recorded in the run's output directory, in `run.log` and in `run_info.txt` (see
[Outputs](#outputs)).

## Large assemblies

On assemblies of tens of thousands of contigs, separation and recruitment may reach their time
limits before proving their results optimal; `summary.tsv` indicates which stages did so, and their
result is the best found. Recruitment reaches a good solution early, so its default time limit is
generous. Separation may return different, near-equivalent partitions under different time limits
on the most difficult assemblies; increasing `SEP_TL` does not necessarily remove this variation,
and the time limit used should be reported.

## Notes and limitations

* **Contigs shorter than 1,500 bp are excluded**, even when the seed assigns them to a bin, because
  MetaBAT2 computes no features for them. They are absent from every output, including
  `pred_seed.tsv`.
* **Recruitment can only assign contigs for which MetaBAT2 provides features.** Contigs excluded
  by MetaBAT2 for very low coverage, and contigs with no similar contig in any bin, remain unbinned.
* **Contigs that MetaBAT2 cannot measure are treated as conflicting in merge and recruitment.**
  MetaBAT2 excludes contigs with very low coverage, so no distance is available between such a
  contig and the other contigs that share its marker genes. A seed produced by another binner may
  nevertheless assign such a contig to a bin. To prevent its marker genes from being disregarded,
  **merge and recruitment treat every unmeasured pair as dissimilar** (distance 200, the maximum).
  This is the conservative choice: it may prevent a correct merge or assignment in the absence of
  evidence, but it never permits one because evidence is missing. Separation and the threshold
  derivation do not use these pairs. For a seed produced by MetaBAT2 the setting has no effect,
  since such contigs are then never binned. `summary.tsv` reports the number of such pairs at each
  stage.
* **Bins smaller than 200 kb are removed from the final binning**, and their contigs become
  unbinned. This is the minimum bin size applied by MetaBAT2, and is needed because separation can
  produce small fragments. The filter is applied once, after the last operator; the operators
  themselves work on unfiltered binnings.
* **An operator may reach its time limit** before proving its result optimal, in which case the
  best result found is returned (see [Large assemblies](#large-assemblies)).
* **Equally optimal solutions may differ between runs.** When two solutions have exactly the same
  objective value, the one returned depends on the scheduling of the solver's parallel threads, so
  repeated runs may assign a contig differently even when both prove their result optimal. On
  `strong100`, two runs of `--operators recruit` assigned 1 of 301 recruited contigs to different
  bins, with identical optimal objective values. Both results are valid; the method does not
  prefer either.

## Troubleshooting

| Message | Cause and resolution |
|---|---|
| `clingo not found`, `Python 3.10+ not found` | A program path is incorrect. Set `CLINGO` or `PYTHON` in `config.sh` accordingly. |
| `numpy and scipy are required` | Install them for the Python that runs the tool (`pip install numpy scipy`), or use `--backend metabat2`. |
| `METABAT2_MODIFIED is not set`, `modified MetaBAT2 not found` (with `--backend metabat2`) | Set `METABAT2_MODIFIED` in `config.sh` to the modified MetaBAT2 binary. |
| `... is not the MODIFIED MetaBAT2` | `METABAT2_MODIFIED` points to an unmodified MetaBAT2, or the patch was not applied when it was built. Rebuild it as described in `metabat2_modified/README.md`; check that `git apply` reported no error. |
| `error while loading shared libraries` when MetaBAT2 runs | The modified MetaBAT2 cannot locate a library it was built against. Set `METABAT2_MODIFIED_LIBS` in `config.sh` to the directory containing it. |
| `modified MetaBAT2 did not write composite_raw.lp` | MetaBAT2 itself failed. Its error message appears immediately above in `run.log`. |
| `depth: contigs are not in the same order as in the assembly` | The depth file was produced from a different version of the assembly, or reordered afterwards. Regenerate it from the same FASTA file. |
| `none of its contig names occur in the assembly` (seed or markers) | The contig names differ from the FASTA headers (for example, renamed contigs, or full headers instead of their first word). The same names must be used. |
| `the contig features could not be computed` | The assembly or the depth file could not be read. The error message appears immediately above in `run.log`. |
| `... was made with different settings ... Use a new --out` | The output directory belongs to a run with a different seed or different parameters. Use a new output directory. |
| `the prediction has N contigs, expected at least M` | clingo was stopped by its time limit while writing its result. Increase the time limit of that operator (`SEP_TL`, `MRG_TL` or `REC_TL`) and run again; completed stages are retained. |
| `simthresh COULD NOT BE DERIVED` | The seed contains no marker-sharing contigs within the same bin, or none across different bins, so the threshold cannot be derived. The value 50 is used, and should be reported as such. |
| `solver: time limit` in `summary.tsv` | The stage returned the best result found, not a proved optimum (see [Large assemblies](#large-assemblies)). |

## Repository structure

```
refine.sh             Main command.
config.sh.example     Template for the configuration file (config.sh).
asp_encodings/        The three operators, as clingo programs.
scripts/              Auxiliary scripts invoked by refine.sh; not intended to be run directly.
metabat2_modified/    Patch to MetaBAT2 that exports its features (for --backend metabat2).
```
