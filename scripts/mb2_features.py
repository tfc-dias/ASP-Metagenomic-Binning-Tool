#!/usr/bin/env python3
"""MetaBAT2's pairwise contig features, reimplemented in Python (numpy + scipy).

The refinement operators read two things from MetaBAT2, and this module computes both without
the MetaBAT2 binary:

  composite   MetaBAT2's native edge score on its TNF-pruned similarity graph (the score it
              clusters on), for every edge of that graph -> composite_raw.lp, the same file the
              patched binary wrote in its default mode.
  pairs       cal_tnf_dist + cal_abd_dist for a supplied list of contig pairs, as
              D = int(tnf*100) + int(abd*100) on 0..200 (what the patched binary wrote in its
              ASP_DUMP_PAIRS mode). Used by build_comarker.py through pair_distances().

It is a port of metabat2.cpp at the commit the patch targets (c869c52), default settings, and
follows the original step by step, including where MetaBAT2 rounds to single precision (the TNF
matrix, the depth matrix, log-lengths, the normal distributions of cal_abd_dist, the stored
edge scores). The pTNF cutoff search draws random numbers; it uses an exact emulation of glibc's
rand(), seeded like MetaBAT2's --seed, so the cutoff is the one the binary picks. What is NOT
reproduced bit for bit: the order of floating-point additions inside a few sums (differences
around 1e-16), which can move a stored score by one unit in its last float digit, and the choice
between two neighbours with exactly the same similarity when only one fits in a contig's 200
edges.

Usage:
  mb2_features.py features  --fasta A.fasta --depth DEPTH.txt --minlen 1500 --out features.npz
  mb2_features.py composite --features features.npz --out composite_raw.lp [--threads N] [--seed 1]
  mb2_features.py pairs     --features features.npz --pairs PAIRS.tsv --out comarker_raw.lp
"""
import os
# one BLAS thread per worker process; parallelism comes from the process pool
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import gzip
import math
import sys

import numpy as np
from scipy.special import erfc

F32, F64 = np.float32, np.float64

# ---- MetaBAT2 defaults (metabat2.cpp / metabat2.h) ---------------------------------------
MAXP = 95 / 100.          # --maxP
MAX_EDGES = 200           # --maxEdges
MIN_CV = 1.0              # --minCV
MIN_CV_SUM = 1.0          # --minCVSum
MIN_SAMPLE = 3            # abundance correlation used from 3 samples on

# ---- tetranucleotides (metabat2.h) -------------------------------------------------------
TN = ["GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT",
      "CAAA", "CCAG", "GGAC", "ATTA", "GATC", "CCTC", "CTAA", "ACTA", "AGGC",
      "GCAA", "CCGC", "CGCC", "AAAC", "ACTC", "ATCC", "GACC", "GAGA", "ATAG",
      "ATCA", "CAGA", "AGTA", "ATGA", "AAAT", "TTAA", "TATA", "AGTG", "AGCT",
      "CCAC", "GGCC", "ACCC", "GGGA", "GCGC", "ATAC", "CTGA", "TAGA", "ATAT",
      "GTCA", "CTCC", "ACAA", "ACCT", "TAAA", "AACG", "CGAG", "AGGG", "ATCG",
      "ACGC", "TCAA", "CTAC", "CTCA", "GACA", "GGAA", "CTTC", "GCCC", "CTGC",
      "TGCA", "GGCA", "CACG", "GAGC", "AACT", "CATG", "AATT", "ACAG", "AGAT",
      "ATAA", "CATC", "GCCA", "TCGA", "CACA", "CAAC", "AAGG", "AGCA", "ATGG",
      "ATTC", "GTGA", "ACCG", "GATA", "GCTA", "CGTC", "CCCG", "AAGC", "CGTA",
      "GTAC", "AGGA", "AATG", "CACC", "CAGC", "CGGC", "ACAC", "CCGG", "CCGA",
      "CCCC", "TGAA", "AACA", "AGAG", "CCCA", "CGGA", "TACA", "ACCA", "ACGT",
      "GAAC", "GTAA", "ATGC", "GTTA", "TCCA", "CAGG", "ACTG", "AAAG", "AAGA",
      "CAAG", "GCGA", "AACC", "ACGG", "CCAA", "CTTA", "AGAC", "AGCG", "GAAA",
      "AATC", "ATTG", "GCAC", "CCTA", "CGAC", "CTAG", "AGAA", "CGCA", "CGCG",
      "AATA"]
TNP = {"ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC",
       "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC"}
NTNF = len(TN)            # 136


def _tn_lookup():
    """TNLookup of main(): raw 4-mer number (first base least significant) -> TNF column."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    tnmap = {t: i for i, t in enumerate(TN)}
    lut = np.full(257, NTNF, dtype=np.int64)
    bases = "ACGT"
    for i0 in range(4):
        for i1 in range(4):
            for i2 in range(4):
                for i3 in range(4):
                    s = bases[i0] + bases[i1] + bases[i2] + bases[i3]
                    num = i0 + 4 * i1 + 16 * i2 + 64 * i3
                    if s in tnmap:
                        lut[num] = tnmap[s]
                        continue
                    rc = "".join(comp[c] for c in reversed(s))
                    lut[num] = NTNF if rc in TNP else tnmap[rc]
    return lut


TN_LUT = _tn_lookup()
BASE_NUM = np.full(256, 4, dtype=np.int64)
for _c, _n in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    BASE_NUM[ord(_c)] = _n
    BASE_NUM[ord(_c.lower())] = _n


def tnf_vector(seq):
    """generateTNF: canonical tetranucleotide counts, L2-normalised, stored as float."""
    b = BASE_NUM[np.frombuffer(seq, dtype=np.uint8)]
    bad = b == 4
    b = np.where(bad, 0, b)
    code = b[:-3] + 4 * b[1:-2] + 16 * b[2:-1] + 64 * b[3:]
    code[bad[:-3] | bad[1:-2] | bad[2:-1] | bad[3:]] = 256
    counts = np.bincount(TN_LUT[code], minlength=NTNF + 1)[:NTNF].astype(F32)
    sq = (counts * counts).astype(F64)              # float product, summed in double
    rsum = math.sqrt(np.add.accumulate(sq)[-1])
    return (counts.astype(F64) / rsum).astype(F32)


# ---- input -------------------------------------------------------------------------------
def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if gz else open(path, "rb")


def read_fasta(path, minlen):
    """Contigs of at least minlen bp, in file order: names, lengths, TNF rows."""
    names, lens, tnfs, seen = [], [], [], set()

    def flush(name, chunks):
        if name is None:
            return
        seq = b"".join(chunks)
        if len(seq) >= minlen:
            if name in seen:
                sys.exit(f"[Error!] found duplicate contig: {name}")
            seen.add(name)
            names.append(name)
            lens.append(len(seq))
            tnfs.append(tnf_vector(seq))

    name, chunks = None, []
    with _open(path) as f:
        for line in f:
            if line.startswith(b">"):
                flush(name, chunks)
                name = line[1:].split(None, 1)[0].decode() if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
        flush(name, chunks)
    tnf = np.vstack(tnfs) if tnfs else np.zeros((0, NTNF), F32)
    return names, np.array(lens, dtype=np.int64), tnf


def read_depth(path, names):
    """jgi_summarize_bam_contig_depths table: per large contig, mean and variance per sample.

    Returns ABD and ABD_VAR as double (rounded to float where used), and a flag for rows that
    pass MetaBAT2's sanity checks. The minCV / minCVSum filter is applied later, by load().
    """
    index = {n: i for i, n in enumerate(names)}
    n = len(names)
    with open(path) as f:
        f.readline()
        first = f.readline().rstrip("\r\n")
    ncol = len(first.split("\t")) if first else 0
    nabd = ncol - 3
    if nabd <= 0 or nabd % 2:
        sys.exit("[Error!] Number of columns (excluding the first column) in abundance data "
                 "file is not even.")
    nabd //= 2
    abd = np.zeros((n, nabd), F64)
    var = np.zeros((n, nabd), F64)
    good = np.zeros(n, bool)
    got = np.zeros(n, bool)
    num = 0
    with open(path) as f:
        f.readline()
        for row in f:
            row = row.rstrip("\r\n")
            if not row:
                continue
            fields = row.split("\t")
            label = fields[0].split(" ")[0]
            i = index.get(label)
            if i is None:
                continue                      # short contig or one not in the assembly
            if i != num:
                sys.exit(f"[Error!] the order of contigs in abundance file is not the same as "
                         f"the assembly file: {label}")
            num += 1
            got[i] = True
            vals = fields[3:]
            ok = len(vals) == 2 * nabd
            if ok:
                m = np.array(vals[0::2], F64)
                v = np.array(vals[1::2], F64)
                ok = bool(np.all(m <= 1e7) and np.all(m >= 0)
                          and np.all(v <= 1e14) and np.all(v >= 0))
                abd[i], var[i] = m, v
            good[i] = ok
    if not got.all():
        miss = names[int(np.argmin(got))]
        sys.exit(f"[Error!] {int((~got).sum())} contigs of the assembly are missing from the "
                 f"depth file (first: {miss})")
    return abd, var, good


def build_features(fasta, depth, minlen, out):
    names, lens, tnf = read_fasta(fasta, minlen)
    if not names:
        sys.exit(f"[Error!] no contigs >= {minlen} bp in {fasta}")
    abd, var, good = read_depth(depth, names)
    np.savez(out, names=np.array(names), lens=lens, tnf=tnf, abd=abd, var=var, good=good,
             minlen=minlen)
    print(f"[features] {len(names)} contigs >= {minlen} bp, {abd.shape[1]} samples -> {out}")


class Features:
    """The contigs MetaBAT2 keeps at the given minCV / minCVSum, with its stored types."""

    def __init__(self, path, min_cv=MIN_CV, min_cv_sum=MIN_CV_SUM):
        z = np.load(path, allow_pickle=False)
        min_cv_sum = max(min_cv, min_cv_sum)
        abd = z["abd"]
        mean_sum = np.where(abd >= min_cv, abd, 0.0).sum(axis=1)
        keep = z["good"] & (mean_sum >= min_cv_sum)
        self.min_cv = min_cv
        self.names = [str(x) for x in z["names"][keep]]
        self.index = {n: i for i, n in enumerate(self.names)}
        self.n = len(self.names)
        self.tnf = z["tnf"][keep]                                   # float
        self.abd = abd[keep].astype(F32)                            # float
        self.var = z["var"][keep].astype(F32)
        self.nabd = self.abd.shape[1]
        lens = z["lens"][keep]
        logsize = np.log10(np.minimum(lens, 500000).astype(F64)).astype(F32)
        # powers of the log-length as _cal_tnf_dist forms them: lw?2 is a float product,
        # the higher powers are double products
        P = np.zeros((self.n, 8), F64)
        P[:, 1] = logsize
        P[:, 2] = (logsize * logsize).astype(F64)
        for k in range(3, 8):
            P[:, k] = P[:, k - 1] * P[:, 1]
        self.P = P
        self.logsize = logsize
        self.nz = (self.abd.astype(F64) > min_cv).any(axis=1)      # is_nz, per contig


# ---- distance functions ------------------------------------------------------------------
FLOOR_PROB = 0.1
FLOOR_PRE = math.log((1.0 / FLOOR_PROB) - 1.0)


# The two logistic models of _cal_tnf_dist, b + c*d each, as polynomials in the log-lengths of the
# shorter (x) and the longer (y) contig of a pair. Each is listed term by term in MetaBAT2's own
# order of summation, so that tnf_prob() reproduces its rounding: (kind, power, coefficient) with
# kind "x" or "y" (coef * x^k), "xy" (coef * x^k * y^k), the first entry being the constant.
TNF_MODEL = {
    "b1": [("1", 0, 46349.1624324381), ("x", 1, -76092.3748553155), ("y", 1, -639.918334183),
           ("x", 2, 53873.3933743949), ("y", 2, -156.6547554844), ("x", 3, -21263.6010657275),
           ("y", 3, 64.7719132839), ("x", 4, 5003.2646455284), ("y", 4, -8.5014386744),
           ("x", 5, -700.5825500292), ("y", 5, 0.3968284526), ("x", 6, 54.037542743),
           ("x", 7, -1.7713972342), ("xy", 1, 474.0850141891), ("xy", 2, -23.966597785),
           ("xy", 3, 0.7800219061), ("xy", 4, -0.0138723693), ("xy", 5, 0.0001027543)],
    "c1": [("1", 0, -443565.465710869), ("x", 1, 718862.10804858), ("y", 1, 5114.1630934534),
           ("x", 2, -501588.206183097), ("y", 2, 784.4442123743), ("x", 3, 194712.394138513),
           ("y", 3, -377.9645994741), ("x", 4, -45088.7863182741), ("y", 4, 50.5960513287),
           ("x", 5, 6220.3310639927), ("y", 5, -2.3670776453), ("x", 6, -473.269785487),
           ("x", 7, 15.3213264134), ("xy", 1, -3282.8510348085), ("xy", 2, 164.0438603974),
           ("xy", 3, -5.2778800755), ("xy", 4, 0.0929379305), ("xy", 5, -0.0006826817)],
    "b2": [("1", 0, 6770.9351457442), ("x", 1, -5933.7589419767), ("y", 1, -2976.2879986855),
           ("x", 2, 3279.7524685865), ("y", 2, 1602.7544794819), ("x", 3, -967.2906583423),
           ("y", 3, -462.0149190219), ("x", 4, 159.8317289682), ("y", 4, 74.4884405822),
           ("x", 5, -14.0267151808), ("y", 5, -6.3644917671), ("x", 6, 0.5108811613),
           ("y", 6, 0.2252455343), ("xy", 2, 0.965040193), ("xy", 3, -0.0546309127),
           ("xy", 4, 0.0012917084), ("xy", 5, -1.14383e-05)],
    "c2": [("1", 0, 39406.5712626297), ("x", 1, -77863.1741143294), ("y", 1, 9586.8761567725),
           ("x", 2, 55360.1701572325), ("y", 2, -5825.2491611377), ("x", 3, -21887.8400068324),
           ("y", 3, 1751.6803621934), ("x", 4, 5158.3764225203), ("y", 4, -290.1765894829),
           ("x", 5, -724.0348081819), ("y", 5, 25.364646181), ("x", 6, 56.0522105105),
           ("y", 6, -0.9172073892), ("x", 7, -1.8470088417), ("xy", 1, 449.4660736502),
           ("xy", 2, -24.4141920625), ("xy", 3, 0.8465834103), ("xy", 4, -0.0158943762),
           ("xy", 5, 0.0001235384)],
}


def _poly(name, X, Y):
    """One model coefficient, summed in MetaBAT2's order."""
    terms = TNF_MODEL[name]
    acc = terms[0][2]
    for kind, k, coef in terms[1:]:
        if kind == "x":
            acc = acc + coef * X[k]
        elif kind == "y":
            acc = acc + coef * Y[k]
        else:
            acc = acc + coef * X[k] * Y[k]
    return acc


def tnf_prob(d, X, Y):
    """_cal_tnf_dist's two logistic models. X / Y: power columns of the shorter / longer contig."""
    with np.errstate(over="ignore"):
        pre = -(_poly("b1", X, Y) + _poly("c1", X, Y) * d)
        prob = np.where(pre <= FLOOR_PRE, FLOOR_PROB, 1.0 / (1 + np.exp(pre)))
        pre2 = -(_poly("b2", X, Y) + _poly("c2", X, Y) * d)
        prob2 = np.where(pre2 <= FLOOR_PRE, 1.0 / (1 + np.exp(pre2)), FLOOR_PROB)
    return np.where(prob >= FLOOR_PROB, prob2, prob)


def _minmax_powers(ft, i, j):
    """Power columns of the shorter (X) and longer (Y) contig of each pair (i, j)."""
    Pi, Pj = ft.P[i], ft.P[j]
    shorter_i = (ft.logsize[i] <= ft.logsize[j])[..., None]
    X = np.where(shorter_i, Pi, Pj)
    Y = np.where(shorter_i, Pj, Pi)
    return np.moveaxis(X, -1, 0), np.moveaxis(Y, -1, 0)


def tnf_dist_pairs(ft, i, j):
    """cal_tnf_dist for pair arrays i, j, as the binary computes it (float differences)."""
    diff = (ft.tnf[i] - ft.tnf[j]).astype(F64)
    d = np.sqrt((diff * diff).sum(axis=1))
    X, Y = _minmax_powers(ft, i, j)
    return tnf_prob(d, X, Y)


ROOT2_F = F32(math.sqrt(2.0))


def _norm_cdf(mean, sd, x):
    """boost::math::cdf(normal_distribution<float>, x): float arithmetic, erfc in double."""
    xf = x.astype(F32)
    diff = (xf - mean) / (sd * ROOT2_F)
    return erfc(-diff.astype(F64)).astype(F32) / F32(2)


def abd_dist_sample(ft, i, j, s):
    """cal_abd_dist for one sample s over pair arrays: (distance, sample counts)."""
    m1, m2 = ft.abd[i, s], ft.abd[j, s]
    nz = (m1.astype(F64) > ft.min_cv) | (m2.astype(F64) > ft.min_cv)
    m1 = np.maximum(m1, F32(1e-6))
    m2 = np.maximum(m2, F32(1e-6))
    v1, v2 = ft.var[i, s], ft.var[j, s]
    sd1 = np.sqrt(np.where(v1 < 1, F32(1), v1))
    sd2 = np.sqrt(np.where(v2 < 1, F32(1), v2))
    M1, M2 = m1.astype(F64), m2.astype(F64)
    V1 = sd1.astype(F64) ** 2
    V2 = sd2.astype(F64) ** 2
    vd = V1 - V2
    eqv = np.abs(vd) < 1e-4
    with np.errstate(divide="ignore", invalid="ignore"):
        mm = M1 - M2
        tmp = np.sqrt(V1 * V2 * (mm * mm - 2 * vd * np.log(np.sqrt(V2 / V1))))
        m1v2, m2v1 = M1 * V2, M2 * V1
        k1 = (tmp - m1v2 + m2v1) / vd
        k2 = (tmp + m1v2 - m2v1) / (-vd)
    mid = (M1 + M2) / 2
    k1 = np.where(eqv, mid, k1)
    k2 = np.where(eqv, mid, k2)
    k1, k2 = np.minimum(k1, k2), np.maximum(k1, k2)
    sw = V1 > V2                                   # std::swap(p1, p2)
    a_m, a_s = np.where(sw, m2, m1), np.where(sw, sd2, sd1)
    b_m, b_s = np.where(sw, m1, m2), np.where(sw, sd1, sd2)
    with np.errstate(invalid="ignore"):
        one = np.abs(_norm_cdf(a_m, a_s, k1) - _norm_cdf(b_m, b_s, k1))
        two = np.abs(_norm_cdf(a_m, a_s, k2) - _norm_cdf(a_m, a_s, k1)
                     + _norm_cdf(b_m, b_s, k1) - _norm_cdf(b_m, b_s, k2))
    d = np.where(k1 == k2, one, two).astype(F64)
    d = np.where(m1 != m2, d, 0.0)
    return np.minimum(np.maximum(d, 1e-6), 1.0 - 1e-6), nz


def abd_corr(ft, i, j):
    """_cal_abd_corr: Pearson correlation of the sample means, in its streaming form."""
    sxx = np.zeros(len(i)); syy = np.zeros(len(i)); sxy = np.zeros(len(i))
    mx = ft.abd[i, 0].astype(F64)
    my = ft.abd[j, 0].astype(F64)
    for s in range(1, ft.nabd):
        ip1 = float(s + 1)
        ratio = s / ip1
        dx = ft.abd[i, s].astype(F64) - mx
        dy = ft.abd[j, s].astype(F64) - my
        sxx += dx * dx * ratio
        syy += dy * dy * ratio
        sxy += dx * dy * ratio
        mx += dx / ip1
        my += dy / ip1
    with np.errstate(divide="ignore", invalid="ignore"):
        return sxy / (np.sqrt(sxx) * np.sqrt(syy))


def pair_distances_ft(ft, pairs):
    """{(name1, name2): D} for the pairs both of whose contigs MetaBAT2 keeps."""
    ok = [(a, b) for a, b in pairs if a in ft.index and b in ft.index and a != b]
    if not ok:
        return {}
    i = np.array([ft.index[a] for a, _ in ok])
    j = np.array([ft.index[b] for _, b in ok])
    out = {}
    CH = 200000
    for s0 in range(0, len(ok), CH):
        ii, jj = i[s0:s0 + CH], j[s0:s0 + CH]
        tnf_int = (tnf_dist_pairs(ft, ii, jj) * 100).astype(np.int64)
        total = np.zeros(len(ii)); valid = np.zeros(len(ii), np.int64)
        for s in range(ft.nabd):
            d, nz = abd_dist_sample(ft, ii, jj, s)
            total = np.where(nz, total + d, total)
            valid += nz
        with np.errstate(divide="ignore", invalid="ignore"):
            abd_prob = np.where(valid > 0, total / valid, 1.0)
        abd_int = (abd_prob * 100).astype(np.int64)
        for k, D in enumerate(tnf_int + abd_int):
            out[ok[s0 + k]] = int(D)
    return out


def pair_distances(features, pairs, min_cv=MIN_CV, min_cv_sum=MIN_CV_SUM):
    return pair_distances_ft(Features(features, min_cv, min_cv_sum), pairs)


# ---- the composite graph -----------------------------------------------------------------
class GlibcRand:
    """glibc rand() after srand(seed) (TYPE_3 additive feedback generator)."""

    def __init__(self, seed):
        seed &= 0xFFFFFFFF
        if seed == 0:
            seed = 1
        r = [0] * 34
        r[0] = seed - 2 ** 32 if seed >= 2 ** 31 else seed
        for k in range(1, 31):
            hi = int(r[k - 1] / 127773)
            lo = r[k - 1] - hi * 127773
            w = 16807 * lo - 2836 * hi
            r[k] = w + 2147483647 if w < 0 else w
        for k in range(31, 34):
            r[k] = r[k - 31]
        self.r = [x & 0xFFFFFFFF for x in r]
        for _ in range(310):
            self.rand()

    def rand(self):
        r = self.r
        v = (r[-31] + r[-3]) & 0xFFFFFFFF
        r.append(v)
        if len(r) > 64:
            del r[:len(r) - 34]
        return v >> 1


_FT = None            # shared with the worker processes (fork)
_K = 0


def _screen_tables(ft):
    """Per-contig pieces of the TNF models for the screening pass. Each model coefficient of a
    pair is A(shorter) + B(longer) + sum_k w_k * x^k * y^k: A and B are computed once per
    contig, and the cross sum is a matrix product."""
    ft.A, ft.B, ft.W = {}, {}, {}
    for name, terms in TNF_MODEL.items():
        a = np.full(ft.n, terms[0][2])
        b = np.zeros(ft.n)
        w = np.zeros(5)
        for kind, k, coef in terms[1:]:
            if kind == "x":
                a += coef * ft.P[:, k]
            elif kind == "y":
                b += coef * ft.P[:, k]
            else:
                w[k - 1] = coef
        ft.A[name], ft.B[name], ft.W[name] = a, b, w
    ft.P15 = np.ascontiguousarray(ft.P[:, 1:6])


def _screen_key(ft, r0, r1):
    """A score for every pair (rows r0:r1 x all contigs) that orders pairs exactly as MetaBAT2's
    TNF similarity does, without evaluating it. The similarity is 1 - 1/(1 + exp(p)) with
    p = pre1 where the first model applies (pre1 > log 9) and min(pre2, log 9) elsewhere, an
    increasing function of p, so p itself ranks the pairs."""
    G = ft.tnf64[r0:r1] @ ft.tnf64.T
    G *= -2.0
    G += ft.sq[r0:r1, None]
    G += ft.sq[None, :]
    np.maximum(G, 0.0, out=G)
    d = np.sqrt(G, out=G)
    shorter = ft.logsize[r0:r1, None] <= ft.logsize[None, :]

    def coef(name):
        A, B = ft.A[name], ft.B[name]
        v = np.where(shorter, A[r0:r1, None] + B[None, :], A[None, :] + B[r0:r1, None])
        v += (ft.P15[r0:r1] * ft.W[name]) @ ft.P15.T
        return v

    pre1 = coef("c1")
    pre1 *= d
    pre1 += coef("b1")
    np.negative(pre1, out=pre1)
    pre2 = coef("c2")
    pre2 *= d
    pre2 += coef("b2")
    np.negative(pre2, out=pre2)
    np.minimum(pre2, FLOOR_PRE, out=pre2)
    return np.where(pre1 > FLOOR_PRE, pre1, pre2)


def _block(rows):
    """For a block of rows: the best TNF similarity of each row (to any other contig), and its K
    most similar is_nz contigs with their similarity. Candidates are chosen by the screening key
    and their similarities then computed exactly, as MetaBAT2 does."""
    ft, K = _FT, _K
    r0, r1 = rows
    n = ft.n
    I = np.arange(r0, r1)
    S = _screen_key(ft, r0, r1)
    S[np.arange(r1 - r0), I] = -np.inf
    # best similarity of each row: the few top candidates, recomputed exactly
    kk = min(8, n - 1)
    top = np.argpartition(-S, kk - 1, axis=1)[:, :kk]
    ii = np.repeat(I, kk)
    ex = 1.0 - tnf_dist_pairs(ft, ii, top.ravel())
    rowmax = ex.reshape(-1, kk).max(axis=1)
    # graph candidates: is_nz pairs only
    if not ft.nz.all():
        S[~(ft.nz[r0:r1, None] | ft.nz[None, :])] = -np.inf
    Kc = min(K, n - 1)
    cand = np.argpartition(-S, Kc - 1, axis=1)[:, :Kc]
    valid = np.isfinite(np.take_along_axis(S, cand, axis=1))
    ii = np.repeat(I, Kc)
    sim = (1.0 - tnf_dist_pairs(ft, ii, cand.ravel())).astype(F32).reshape(-1, Kc)
    sim[~valid] = -np.inf
    return r0, rowmax, cand.astype(np.int32), sim


def _all_pairs(ft, threads, margin=56):
    global _FT, _K
    ft.tnf64 = ft.tnf.astype(F64)
    ft.sq = (ft.tnf64 * ft.tnf64).sum(axis=1)
    _screen_tables(ft)
    _FT, _K = ft, MAX_EDGES + margin
    n = ft.n
    step = max(1, min(256, int(1_000_000 // max(n, 1))))
    blocks = [(a, min(a + step, n)) for a in range(0, n, step)]
    rowmax = np.empty(n)
    Kc = min(_K, n - 1)
    cand = np.empty((n, Kc), np.int32)
    sim = np.empty((n, Kc), F32)
    if threads > 1 and len(blocks) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(threads) as pool:
            results = pool.imap_unordered(_block, blocks, chunksize=1)
            for r0, rm, c, s in results:
                rowmax[r0:r0 + len(rm)], cand[r0:r0 + len(rm)], sim[r0:r0 + len(rm)] = rm, c, s
    else:
        for b in blocks:
            r0, rm, c, s = _block(b)
            rowmax[r0:r0 + len(rm)], cand[r0:r0 + len(rm)], sim[r0:r0 + len(rm)] = rm, c, s
    return rowmax, cand, sim


def _sample_cutoff(rowmax, n, full, rng, coverage=MAXP):
    """gen_tnf_graph_sample: the cutoff p (per mille) at which `coverage` of the sampled contigs
    have a neighbour. A contig in the sample is connected at cutoff c exactly when its best
    similarity to any other contig is >= c, so the sample's best similarities suffice."""
    m = n if full else min(n, 2500)
    idx = np.arange(n)
    left = n
    for k in range(m):                              # random_unique (partial Fisher-Yates)
        r = k + rng.rand() % left
        idx[k], idx[r] = idx[r], idx[k]
        left -= 1
    best = rowmax[idx[:m]]
    p, pp = 999, 1000
    cov = pcov = 0.0
    while p > 700:
        cutoff = p / 1000.
        cov = int((best >= cutoff).sum()) / m
        if cov >= coverage:
            if cov - coverage > coverage - pcov:
                p, cov = pp, pcov
            break
        pp, pcov = p, cov
        if p > 990:
            p -= rng.rand() % 3 + 1
        elif p > 900:
            p -= rng.rand() % 3 + 3
        else:
            p -= rng.rand() % 3 + 9
    return p


def _ptnf(rowmax, n, seed):
    rng = GlibcRand(seed)
    if n <= 25000:
        return float(_sample_cutoff(rowmax, n, True, rng))
    ptnf = 0.0
    for k in range(10):
        minp = float(_sample_cutoff(rowmax, n, False, rng))
        if minp < 701:
            minp = 700.
        ptnf += minp
        if k == 1 and ptnf / 2 < 701:
            return 700.
        if k == 9:
            ptnf /= 10.
    return ptnf


def composite(ft, threads=1, seed=1, log=print):
    """MetaBAT2's composite edge scores: (from, to, score) arrays over its TNF graph."""
    n = ft.n
    rowmax, cand, sim = _all_pairs(ft, threads)
    ptnf = _ptnf(rowmax, n, seed)
    cutoff = ptnf / 1000.
    log(f"[composite] {n} contigs, pTNF = {ptnf / 10.:.2f}")
    # gen_tnf_graph: per contig its maxEdges most similar is_nz contigs above the cutoff (ties
    # to the lower index); an edge i-j is kept from the list of the lower-numbered contig
    order = np.lexsort((cand, -sim.astype(F64)), axis=1)
    cand = np.take_along_axis(cand, order, axis=1)[:, :MAX_EDGES]
    sim = np.take_along_axis(sim, order, axis=1)[:, :MAX_EDGES]
    rows = np.repeat(np.arange(n), cand.shape[1]).reshape(cand.shape)
    keep = (sim.astype(F64) > cutoff) & (cand > rows)
    frm, to, stnf = rows[keep].astype(np.int64), cand[keep].astype(np.int64), sim[keep]
    E = len(frm)
    if E == 0:
        sys.exit("No edges were formed by TNF.")
    # abundance similarity, averaged over the samples that carry signal
    val = np.zeros(E)
    nnz = np.zeros(E, np.int64)
    for s in range(ft.nabd):
        d, nz = abd_dist_sample(ft, to, frm, s)
        val = np.where(nz, val + (1. - d), val)
        nnz += nz
    sscr = (val / nnz).astype(F32)
    abd_distr = np.sort(sscr)
    # TNF (and correlation) ranks mapped onto the abundance distribution, then geometric mean
    srt = np.sort(stnf)
    rank_tnf = (np.searchsorted(srt, stnf, "left") + 1).astype(F32)
    stnf_m = abd_distr[np.round(rank_tnf).astype(np.int64) - 1]
    wtnf = (1. / (1 + nnz)).astype(F32)
    if ft.nabd >= MIN_SAMPLE:
        scor = abd_corr(ft, to, frm).astype(F32)
        if np.isnan(scor).any():
            log(f"[composite] WARNING: {int(np.isnan(scor).sum())} edges have an undefined "
                f"abundance correlation (constant depth); MetaBAT2's ranking of those is undefined")
        srt = np.sort(scor)
        rank_cor = (np.searchsorted(srt, scor, "right")).astype(F32)
        scor_m = abd_distr[np.round(rank_cor).astype(np.int64) - 1]
        score = np.power(np.power(sscr, F32(1) - wtnf) * np.power(stnf_m, wtnf) * scor_m,
                         F32(0.5))
    else:
        score = np.power(sscr, F32(1) - wtnf) * np.power(stnf_m, wtnf)
    return frm, to, score.astype(F32), ptnf


def write_composite(ft, frm, to, score, out):
    names = ft.names
    with open(out, "w") as f:
        for a, b, s in zip(to.tolist(), frm.tolist(), score.astype(F64).tolist()):
            f.write(f'composite("{names[a]}","{names[b]}",{s:.6f}).\n')


# ---- command line ------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("features", help="read the assembly and depth table once")
    a.add_argument("--fasta", required=True)
    a.add_argument("--depth", required=True)
    a.add_argument("--minlen", type=int, default=1500)
    a.add_argument("--out", required=True)
    c = sub.add_parser("composite", help="write composite_raw.lp")
    c.add_argument("--features", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    c.add_argument("--seed", type=int, default=1)
    p = sub.add_parser("pairs", help="distances for a list of contig-name pairs")
    p.add_argument("--features", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-cv", type=float, default=MIN_CV)
    p.add_argument("--min-cv-sum", type=float, default=MIN_CV_SUM)
    args = ap.parse_args()

    if args.cmd == "features":
        if args.minlen < 1500:
            sys.exit("[Error!] Contig length < 1500 is not allowed to be used for binning.")
        build_features(args.fasta, args.depth, args.minlen, args.out)
    elif args.cmd == "composite":
        ft = Features(args.features)
        frm, to, score, _ = composite(ft, args.threads, args.seed)
        write_composite(ft, frm, to, score, args.out)
        print(f"[composite] {len(frm)} edges -> {args.out}")
    else:
        pairs = [tuple(l.split()[:2]) for l in open(args.pairs) if len(l.split()) >= 2]
        dist = pair_distances(args.features, pairs, args.min_cv, args.min_cv_sum)
        with open(args.out, "w") as f:
            for (a, b), D in dist.items():
                f.write(f'comarker_sim("{a}","{b}",{D}).\n')
        print(f"[pairs] {len(dist)} of {len(pairs)} pairs measured -> {args.out}")


if __name__ == "__main__":
    main()
