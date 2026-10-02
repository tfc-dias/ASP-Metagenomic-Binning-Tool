#!/usr/bin/env bash
# =============================================================================
# ASP refinement of a seed binning with any of three operators:
#   separate  splits bins that hold two genomes
#   merge     joins bins that hold parts of one genome
#   recruit   places contigs the seed left unbinned
#
# Writes, in the bioboxes binning format, the binning after each operator (named after the chain so
# far, e.g. pred_separate_merge.tsv) and, for the final stage, the result after the minimum bin size
# filter: with the default stack, pred_separate_merge_recruit_min200k.tsv.
#
# Usage:
#   refine.sh --name NAME --fasta ASSEMBLY.fasta --depth DEPTH.txt \
#             (--markers-pfam FILE | --markers-busco FILE [...] | --markers-table FILE [...]) \
#             --seed BINS_DIR_or_TABLE [--operators LIST] [--out DIR] [--write-bins] [--check]
#
#   --seed     the seed binning to refine (required): a folder of bin FASTA files, or a
#              <contig><TAB><bin> table.
#   --operators  which operators to run, comma-separated, in that order (each at most once).
#              Default separate,merge,recruit. E.g. --operators recruit, or --operators separate,merge.
#   --out      output folder (default ./refine_out/NAME). Re-running with the same --out resumes:
#              any stage whose output exists is skipped. Use a new --out for a different seed or
#              different settings.
#   --write-bins  also write the final binning as a folder of FASTA files, one per bin.
#   --check    only check the configuration and the inputs, then stop.
#
# Advanced (change the method's defaults; for sensitivity studies, not for normal use):
#   --simthresh N     use N as the distance threshold instead of deriving it from the seed
#   --sepmul N        separation penalty weight (default 1)
#   --minconn N       merge: minimum crossing edges between two bins (default 0)
#   --min-bin-size N  minimum bin size in bp for the filtered outputs (default 200000)
#
# Optional environment variables: SEP_TL MRG_TL REC_TL (clingo time limits, s; default 600/180/900),
# THREADS (8), MINLEN (1500, minimum contig length), REFINE_CONFIG (config file; default ./config.sh
# next to this script).
# =============================================================================
set -u
VERSION="1.0.0"
usage(){ sed -n '3,/^# =====/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 1; }

NAME="" FASTA="" DEPTH="" SEED="" OUT="" OPERATORS="separate,merge,recruit"
WRITE_BINS=0 CHECK_ONLY=0 ST_OVERRIDE="" SEPMUL=1 MINCONN=0 MINBIN=200000
CMDLINE="$0 $*"
MARKERS="" MARK_ARGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --name)          NAME="$2"; shift 2 ;;
    --fasta)         FASTA="$2"; shift 2 ;;
    --depth)         DEPTH="$2"; shift 2 ;;
    --seed)          SEED="$2"; shift 2 ;;
    --out)           OUT="$2"; shift 2 ;;
    --operators)     OPERATORS="$2"; shift 2 ;;
    --write-bins)    WRITE_BINS=1; shift ;;
    --check)         CHECK_ONLY=1; shift ;;
    --simthresh)     ST_OVERRIDE="$2"; shift 2 ;;
    --sepmul)        SEPMUL="$2"; shift 2 ;;
    --minconn)       MINCONN="$2"; shift 2 ;;
    --min-bin-size)  MINBIN="$2"; shift 2 ;;
    --version)       echo "refine.sh $VERSION"; exit 0 ;;
    --markers-pfam|--markers-busco|--markers-table)
                     fmt="${1#--markers-}"
                     [ -z "$MARKERS" ] || [ "$MARKERS" = "$fmt" ] || { echo "use one marker format per run (--markers-pfam, --markers-busco or --markers-table)"; exit 1; }
                     MARKERS="$fmt"; MARK_ARGS+=(--marker-file "$2"); shift 2 ;;
    -h|--help)       usage ;;
    *)               echo "unknown option: $1"; usage ;;
  esac
done
[ -n "$NAME" ] && [ -n "$FASTA" ] && [ -n "$DEPTH" ] && [ -n "$MARKERS" ] || usage
[ -n "$SEED" ] || { echo "--seed is required: the binning to refine (a folder of bin FASTA files, or a"
                    echo "contig<TAB>bin table), produced beforehand with any binner."; exit 1; }
IFS=',' read -r -a OPS <<< "$OPERATORS"
[ "${#OPS[@]}" -gt 0 ] || { echo "--operators is empty"; exit 1; }
for op in "${OPS[@]}"; do
  case "$op" in separate|merge|recruit) ;; *) echo "unknown operator '$op' (use separate, merge, recruit)"; exit 1 ;; esac
done
[ "$(printf '%s\n' "${OPS[@]}" | sort | uniq -d)" = "" ] || { echo "--operators: each operator at most once"; exit 1; }
isint(){ case "$1" in ''|*[!0-9]*) return 1 ;; esac; }
[ -z "$ST_OVERRIDE" ] || { isint "$ST_OVERRIDE" && [ "$ST_OVERRIDE" -le 200 ]; } || { echo "--simthresh must be an integer in 0..200"; exit 1; }
for v in "$SEPMUL" "$MINCONN" "$MINBIN"; do isint "$v" || { echo "--sepmul, --minconn and --min-bin-size take non-negative integers"; exit 1; }; done
# filtered outputs are named after the minimum bin size, e.g. _min200k
if [ $((MINBIN % 1000)) -eq 0 ]; then FSUF="_min$((MINBIN / 1000))k"; else FSUF="_min${MINBIN}bp"; fi

# Every path must be ABSOLUTE: the modified MetaBAT2 runs from inside the work folder.
for f in "$FASTA" "$DEPTH" "$SEED"; do [ -e "$f" ] || { echo "not found: $f"; exit 1; }; done
for i in "${!MARK_ARGS[@]}"; do
  [ "${MARK_ARGS[$i]}" = --marker-file ] || [ -e "${MARK_ARGS[$i]}" ] || { echo "not found: ${MARK_ARGS[$i]}"; exit 1; }
done
abspath(){ echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
FASTA="$(abspath "$FASTA")"; DEPTH="$(abspath "$DEPTH")"
SEED="$(abspath "$SEED")"
for i in "${!MARK_ARGS[@]}"; do
  [ "${MARK_ARGS[$i]}" = --marker-file ] || MARK_ARGS[$i]="$(abspath "${MARK_ARGS[$i]}")"
done

# ---- configuration -------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${REFINE_CONFIG:-$HERE/config.sh}"
[ -f "$CONFIG" ] || { echo "missing $CONFIG (copy config.sh.example to config.sh and edit it)"; exit 1; }
# shellcheck source=/dev/null
. "$CONFIG"
# Settings from config.sh. Only METABAT2_MODIFIED has no default; CLINGO and PYTHON fall back to
# the programs of that name on the system path.
CLINGO="${CLINGO:-clingo}"
PY="${PYTHON:-python3}"
MODMB2="${METABAT2_MODIFIED:-}"
MODMB2_LIBS="${METABAT2_MODIFIED_LIBS:-}"
SEP_TL="${SEP_TL:-600}"; MRG_TL="${MRG_TL:-180}"; REC_TL="${REC_TL:-900}"; THREADS="${THREADS:-8}"
MINLEN="${MINLEN:-1500}"

S="$HERE/scripts"
SEP_MODEL="$HERE/asp_encodings/binning_separate.lp"
MRG_MODEL="$HERE/asp_encodings/binning_merge.lp"
REC_MODEL="$HERE/asp_encodings/binning_recruit.lp"

# ---- pre-flight: tools and inputs ----------------------------------------------
echo "Checking the configuration and the inputs:"
preflight_fail=0
command -v "$CLINGO" >/dev/null 2>&1 && echo "  ok       clingo: $("$CLINGO" --version 2>/dev/null | head -1)" \
  || { echo "  ERROR    clingo not found ($CLINGO); set CLINGO in config.sh"; preflight_fail=1; }
"$PY" -c 'import sys; sys.exit(sys.version_info < (3, 10))' 2>/dev/null && echo "  ok       python: $("$PY" --version 2>&1)" \
  || { echo "  ERROR    Python 3.10+ not found ($PY); set PYTHON in config.sh"; preflight_fail=1; }
if [ -z "$MODMB2" ]; then
  echo "  ERROR    METABAT2_MODIFIED is not set in config.sh; build the modified MetaBAT2 (metabat2_modified/README.md) and set it"; preflight_fail=1
elif [ -x "$MODMB2" ]; then
  echo "  ok       modified MetaBAT2: $MODMB2"
else echo "  ERROR    modified MetaBAT2 not found at $MODMB2; check METABAT2_MODIFIED in config.sh"; preflight_fail=1; fi
if [ "$preflight_fail" = 0 ]; then
  "$PY" "$HERE/scripts/check_inputs.py" --fasta "$FASTA" --depth "$DEPTH" --minlen "$MINLEN" \
       --markers "$MARKERS" "${MARK_ARGS[@]}" --seed "$SEED" || preflight_fail=1
fi
[ "$preflight_fail" = 0 ] || { echo "Fix the errors above, then run again."; exit 1; }
if [ "$CHECK_ONLY" = 1 ]; then echo "All checks passed."; exit 0; fi

# ---- layout --------------------------------------------------------------------
OUT="${OUT:-$PWD/refine_out/$NAME}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
RES="$OUT"; WORK="$OUT/work"; F="$WORK/facts"
mkdir -p "$F" "$WORK/moddump"
LOG="$RES/run.log"

log(){ echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }
have(){ [ -s "$1" ]; }                 # output exists and non-empty
die(){ log "FATAL: $*"; exit 1; }
asp_ok(){ grep -qE "Answer: |SATISFIABLE|OPTIMUM FOUND" "$1"; }

# Number of contigs in a prediction (data lines, without the @ headers).
npred(){ [ -s "$1" ] && grep -vc '^@' "$1" || echo 0; }
# GUARD AGAINST TRUNCATION. If clingo is killed by --time-limit while printing a model, the
# last "Answer:" block is cut off and write_binning.py, which reads the LAST one, writes a
# prediction that is missing contigs without complaining. asp_ok does not catch this, because
# the "Answer:" line is still there. It has happened (one run lost 144 of 6,693 atoms), and the
# risk is real: on a large assembly separation can print hundreds of models before the limit.
# No stage may LOSE contigs relative to the previous one: separation and merge only relabel,
# recruitment only adds.
check_n(){  # $1=label  $2=new prediction  $3=minimum expected
  local got; got=$(npred "$2")
  if [ "$got" -lt "$3" ]; then
    die "$1: the prediction has $got contigs, expected at least $3. The model was probably
       truncated when --time-limit fired. Raise that stage's time limit and run again."
  fi
  log "   $1: $got contigs (expected >= $3) OK"
}

# A resumed run must use the same seed and the same settings, or the cached stages would belong
# to another run. (Time limits and threads are not part of this: they do not change the problem.)
SEED_DESC="$SEED"
SETTINGS="seed=$SEED_DESC fasta=$FASTA minlen=$MINLEN markers=$MARKERS simthresh=${ST_OVERRIDE:-derived} sepmul=$SEPMUL minconn=$MINCONN"
if have "$RES/settings.txt" && [ "$(cat "$RES/settings.txt")" != "$SETTINGS" ]; then
  echo "$OUT was made with different settings:"
  echo "  before: $(cat "$RES/settings.txt")"
  echo "  now:    $SETTINGS"
  echo "Use a new --out."; exit 1
fi
echo "$SETTINGS" > "$RES/settings.txt"

log "===== START $NAME ====="
log "fasta: $FASTA | depth: $DEPTH | markers: $MARKERS | seed: $SEED_DESC"
log "operators: ${OPERATORS//,/ -> } | minlen=$MINLEN | threads=$THREADS | TL separate=$SEP_TL merge=$MRG_TL recruit=$REC_TL"
log "out: $OUT"
[ "$SEPMUL$MINCONN$MINBIN${ST_OVERRIDE}" = "10200000" ] || \
  log "NOTE: non-default settings: simthresh=${ST_OVERRIDE:-derived} sepmul=$SEPMUL minconn=$MINCONN min-bin-size=$MINBIN"

# Everything needed to trace or reproduce this run.
{
  echo "refine.sh version   $VERSION"
  echo "date                $(date '+%F %T')"
  echo "command             $CMDLINE"
  echo "clingo              $("$CLINGO" --version 2>/dev/null | head -1)"
  echo "python              $("$PY" --version 2>&1)"
  echo "modified MetaBAT2   $("$MODMB2" -h 2>&1 | grep -m1 -o 'version [^)]*')  ($MODMB2)"
  echo "seed                $SEED_DESC"
  echo "assembly            $FASTA"
  echo "depth               $DEPTH"
  echo "markers             $MARKERS: $(for x in "${MARK_ARGS[@]}"; do [ "$x" = --marker-file ] || printf '%s ' "$x"; done)"
  echo "operators           $OPERATORS"
  echo "min contig length   $MINLEN"
  echo "sepmul              $SEPMUL"
  echo "minconn             $MINCONN"
  echo "unmeasured pairs    counted as dissimilar (200) in merge and recruitment"
  echo "min bin size        $MINBIN bp"
  echo "time limits         separate ${SEP_TL}s, merge ${MRG_TL}s, recruit ${REC_TL}s; threads $THREADS"
} > "$RES/run_info.txt"

# ---- 0) base facts: dict + mg --------------------------------------------------
if have "$F/dict.tsv" && have "$F/mg.lp"; then log "0) base facts — skip (exist)"
else log "0) base facts (dict + mg)"
  $PY "$S/build_base.py" --fasta "$FASTA" --markers "$MARKERS" "${MARK_ARGS[@]}" \
       --out "$F" --minlen "$MINLEN" 2>&1 | tee -a "$LOG" || die "base"; fi
DICT="$F/dict.tsv"; MG="$F/mg.lp"

# ---- 1) the seed ---------------------------------------------------------------
if have "$RES/pred_seed.tsv"; then log "1) seed — skip (exist)"
else
  log "1) seed from $SEED"
  SEED_IN="$SEED"
  $PY "$S/build_seed.py" "$DICT" "$SEED_IN" \
       "$RES/pred_seed.tsv" "$F/seed.lp" "$F/recruitable.lp" "$NAME" 2>&1 | tee -a "$LOG"
  have "$RES/pred_seed.tsv" || die "seed parse"
fi

# ---- 2) modified MetaBAT2 -> native composite affinity -------------------------
# The features come from MetaBAT2 whatever binner made the seed: the operators read its
# pairwise score and its distance functions.
COMP="$WORK/moddump/composite_raw.lp"
if have "$COMP"; then log "2) native composite — skip (exist)"
else
  log "2) modified MetaBAT2 (composite_raw.lp dump)"
  # the modified binary may need its build's shared libraries; set METABAT2_MODIFIED_LIBS in config.sh if so
  ( cd "$WORK/moddump" && LD_LIBRARY_PATH="${MODMB2_LIBS:-}:${LD_LIBRARY_PATH:-}" \
       "$MODMB2" -i "$FASTA" -a "$DEPTH" -o "$WORK/moddump/bin" -m "$MINLEN" --seed 1 ) \
       >>"$LOG" 2>&1 || true
  have "$COMP" || die "modified MetaBAT2 did not write composite_raw.lp"
fi
AFF="$F/affinity_native_top5.lp"
have "$AFF" || $PY "$S/build_affinity.py" "$DICT" "$COMP" "$AFF" 1000 5 2>&1 | tee -a "$LOG"

# ---- 3) the distance threshold, from the SEED ---------------------------------
# NB comarker_sim comes from MetaBAT2's OWN cal_tnf_dist + cal_abd_dist via the modified
# binary's pair-dump mode, so D is on the 0..200 scale the threshold is read on.
MODLD="${MODMB2_LIBS:-}:${LD_LIBRARY_PATH:-}"   # libraries for the modified binary, scoped to its calls
comarker(){  # $1 mode  $2 out  $3.. extra args
  local mode="$1" out="$2"; shift 2
  LD_LIBRARY_PATH="$MODLD" $PY "$S/build_comarker.py" --mode "$mode" \
       --dict "$DICT" --mg "$MG" --fasta "$FASTA" --depth "$DEPTH" \
       --metabat2 "$MODMB2" --minlen "$MINLEN" --out "$out" "$@" 2>&1 | tee -a "$LOG"
}
# SIMTHRESH IS DERIVED, not carried over. The threshold is the cut that best separates the
# co-marker pairs the SEED puts in the same bin from those it puts in different bins, with no
# ground truth; see scripts/derive_simthresh.py. It is computed once, on the seed, whichever
# operators run, so every operator in every stack uses the same value.
log "3) distance threshold from the seed"
have "$F/comarker_within_seed.lp" || comarker within "$F/comarker_within_seed.lp" --pred "$RES/pred_seed.tsv"
if [ -n "$ST_OVERRIDE" ]; then
  ST="$ST_OVERRIDE"
  log "simthresh set by --simthresh: $ST (not derived)"
  echo "$ST (set by --simthresh, not derived)" > "$RES/simthresh_used.txt"
else
  have "$F/comarker_merge_seed.lp" || comarker merge "$F/comarker_merge_seed.lp" --pred "$RES/pred_seed.tsv"
  ST=$($PY "$S/derive_simthresh.py" "$F/comarker_within_seed.lp" "$F/comarker_merge_seed.lp" 2>&1 | tail -1)
  case "$ST" in
    ''|*[!0-9]*) ST=50
       log "!!! simthresh COULD NOT BE DERIVED on this dataset: using the historical value 50 !!!"
       log "!!! Report it as carried over, not as derived.                                   !!!"
       echo "50 (FALLBACK: could not be derived on this dataset)" > "$RES/simthresh_used.txt" ;;
    *) log "simthresh derived for this dataset: $ST"
       echo "$ST (derived from the seed)" > "$RES/simthresh_used.txt" ;;
  esac
fi

echo "distance threshold  $(cat "$RES/simthresh_used.txt")" >> "$RES/run_info.txt"

# UNMEASURED PAIRS COUNT AS DISSIMILAR in the merge and recruitment guards (--miss-fill 200).
# MetaBAT2 gives no distance for a contig it does not load (very low coverage). A seed from
# another binner may still bin such a contig, and with no distance its markers could not block
# a merge or a placement. Filling those pairs with the maximum distance, 200, makes the guards
# treat them as a clash: the cautious choice. It can block a correct merge or placement on no
# real evidence. It is NOT used for separation (a made-up distance must not force a split) nor
# for the threshold derivation (a made-up distance must not move the threshold). With a seed
# made by MetaBAT2 it changes nothing, since such contigs are then never in a bin.
MISS_FILL=200

# ---- 4) the operators, in the order given by --operators -----------------------
# Each operator reads the previous stage's binning ($1) and writes the next ($2). Its facts and
# solver output go to its own folder ($3), named after the chain so far (e.g. separate_merge), so a
# later run with a different --operators reuses every stage it shares with an earlier one.

run_separate(){  # $1 in.tsv  $2 out.tsv  $3 stage dir
  local in="$1" out="$2" d="$3"
  log "   separation inputs"
  $PY "$S/build_separate_inputs.py" "$DICT" "$MG" "$COMP" "$in" "$d/separate_inputs.lp" 2>&1 | tee -a "$LOG"
  if [ "$in" = "$RES/pred_seed.tsv" ]; then cp "$F/comarker_within_seed.lp" "$d/comarker_within.lp"
  else comarker within "$d/comarker_within.lp" --pred "$in"; fi
  $PY "$S/mg_subset.py" "$MG" "$d/separate_inputs.lp" "$d/mg_binned.lp" 2>&1 | tee -a "$LOG"
  # sub-bin budget: colouring number of each bin's cannot-link graph at the derived threshold
  $PY "$S/build_binmax.py" "$d/separate_inputs.lp" "$d/comarker_within.lp" "$d/binmax.lp" $ST 2>&1 | tee -a "$LOG"
  log "   run separation (TL=${SEP_TL}s)"
  "$CLINGO" "$SEP_MODEL" "$d/separate_inputs.lp" "$d/mg_binned.lp" "$d/comarker_within.lp" "$d/binmax.lp" \
       -c simthresh=$ST -c sepmul=$SEPMUL -t "$THREADS" --quiet=1 --time-limit="$SEP_TL" > "$d/clingo.txt" 2>&1 || true
  asp_ok "$d/clingo.txt" || die "separation produced no model (see $d/clingo.txt)"
  $PY "$S/write_binning.py" "$DICT" --pred "$d/clingo.txt" "$out" --sample "$NAME" 2>&1 | tee -a "$LOG"
  check_n "separation" "$out" "$(npred "$in")"
}

run_merge(){  # $1 in.tsv  $2 out.tsv  $3 stage dir
  local in="$1" out="$2" d="$3"
  log "   merge inputs + guard"
  $PY "$S/build_merge_inputs.py" "$DICT" "$AFF" "$in" "$d/merge_inputs.lp" 2>&1 | tee -a "$LOG"
  comarker merge "$d/comarker_merge.lp" --pred "$in" --miss-fill "$MISS_FILL"
  log "   merge guard: unmeasured pairs counted as dissimilar (D=$MISS_FILL); see 'filled' above"
  log "   run merge (TL=${MRG_TL}s)"
  "$CLINGO" "$MRG_MODEL" "$d/merge_inputs.lp" "$MG" "$d/comarker_merge.lp" \
       -c simthresh=$ST -c minconn=$MINCONN -t "$THREADS" --quiet=1 --time-limit="$MRG_TL" > "$d/clingo.txt" 2>&1 || true
  asp_ok "$d/clingo.txt" || die "merge produced no model (see $d/clingo.txt)"
  $PY "$S/write_binning.py" "$DICT" --pred "$d/clingo.txt" "$out" --sample "$NAME" 2>&1 | tee -a "$LOG"
  check_n "merge" "$out" "$(npred "$in")"
}

run_recruit(){  # $1 in.tsv  $2 out.tsv  $3 stage dir
  local in="$1" out="$2" d="$3"
  log "   recruitment inputs + guard"
  # the contigs still unbinned in the input are the ones recruitment may place
  $PY "$S/build_recruit_inputs.py" "$DICT" "$COMP" "$in" "$d/recruitable.lp" "$d/recruit_inputs.lp" 2>&1 | tee -a "$LOG"
  comarker guard "$d/comarker_guard.lp" --recruitable "$d/recruitable.lp" --miss-fill "$MISS_FILL"
  log "   recruitment guard: unmeasured pairs counted as dissimilar (D=$MISS_FILL); see 'filled' above"
  log "   run recruitment (TL=${REC_TL}s)"
  "$CLINGO" "$REC_MODEL" "$d/recruit_inputs.lp" "$MG" "$d/comarker_guard.lp" \
       -c threshold_conta=1 -c conta_simthresh=$ST -c minaff=0 -t "$THREADS" --quiet=1 --time-limit="$REC_TL" \
       > "$d/clingo.txt" 2>&1 || true
  asp_ok "$d/clingo.txt" || die "recruitment produced no model (see $d/clingo.txt)"
  $PY "$S/write_binning.py" "$DICT" --pred "$d/clingo.txt" "$out" --sample "$NAME" 2>&1 | tee -a "$LOG"
  check_n "recruitment" "$out" "$(npred "$in")"
}

IN="$RES/pred_seed.tsv"; CHAIN=""; STAGES=()
step=4
for op in "${OPS[@]}"; do
  CHAIN="${CHAIN:+${CHAIN}_}$op"
  OUTP="$RES/pred_$CHAIN.tsv"
  if have "$OUTP"; then log "$step) $op (stage $CHAIN) — skip (exist)"
  else
    log "$step) $op (stage $CHAIN)"
    D="$WORK/stages/$CHAIN"; rm -rf "$D"; mkdir -p "$D"
    # written under a temporary name and renamed only once it has passed its checks, so an
    # interrupted stage is redone on the next run instead of being taken as finished
    t0=$SECONDS
    "run_$op" "$IN" "$OUTP.partial" "$D"
    mv "$OUTP.partial" "$OUTP"
    echo $((SECONDS - t0)) > "$D/seconds.txt"
  fi
  STAGES+=("$CHAIN:$op:$IN:$OUTP")
  IN="$OUTP"; step=$((step + 1))
done

# ---- minimum bin size filter, on the final stage only ---------------------------
# (default 200 kb, MetaBAT2's own minimum). Intermediate stages are kept unfiltered: each stage
# is applied to the unfiltered output of the previous one, and the filter belongs to the result.
FINAL="pred_${CHAIN}${FSUF}.tsv"
log "$step) minimum bin size filter ($MINBIN bp) on the final stage"
$PY "$S/filter_min_binsize.py" "$IN" "$DICT" "$RES/$FINAL" "$MINBIN" 2>&1 | tee -a "$LOG"

# ---- summary of what each stage did ---------------------------------------------
{
  echo "# $NAME | operators: ${OPERATORS//,/ -> } | distance threshold: $(cat "$RES/simthresh_used.txt")"
  $PY "$S/stage_summary.py" --header
  for st in "${STAGES[@]}"; do
    IFS=':' read -r chain op in outp <<< "$st"
    filt="-"; [ "$chain" = "$CHAIN" ] && filt="$RES/$FINAL"   # the filter applies to the final stage only
    $PY "$S/stage_summary.py" "$chain" "$op" "$in" "$outp" "$filt" "$WORK/stages/$chain"
  done
} > "$RES/summary.tsv"
head -1 "$RES/summary.tsv" | tee -a "$LOG"
tail -n +2 "$RES/summary.tsv" | { column -t -s $'\t' 2>/dev/null || cat; } | tee -a "$LOG"
grep -q "time limit" "$RES/summary.tsv" && \
  log "NOTE: at least one stage stopped at its time limit; its result is the best found, not proved best."

# ---- optional: final binning as FASTA files -------------------------------------
BINDIR="bins_${FINAL#pred_}"; BINDIR="${BINDIR%.tsv}"      # pred_X_min200k.tsv -> bins_X_min200k/
if [ "$WRITE_BINS" = 1 ]; then
  log "writing the final binning as FASTA files"
  $PY "$S/write_bin_fastas.py" "$FASTA" "$RES/$FINAL" "$RES/$BINDIR" 2>&1 | tee -a "$LOG"
fi

log "===== DONE $NAME =====  (in $RES)"
log "   final binning:  $FINAL"
[ "$WRITE_BINS" = 1 ] && log "   final bins:     $BINDIR/"
log "   seed and stage outputs (unfiltered): pred_seed.tsv, $(for st in "${STAGES[@]}"; do printf 'pred_%s.tsv ' "${st%%:*}"; done)"
log "   summary: summary.tsv; settings and versions: run_info.txt"
