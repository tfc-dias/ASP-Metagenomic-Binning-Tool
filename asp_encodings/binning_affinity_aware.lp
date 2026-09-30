% Whole-partition metagenomic binning with affinity-aware contamination guard

% k: Number of candidate bins (overridden at runtime via -c k=...)
#const k=13.

% lengthmin: Minimum bin sequence length (200,000 bp)
#const lengthmin=200000.

% lengthmax: Maximum bin sequence length (20 Mbp), preventing degenerate mega-bins
#const lengthmax=20000000.

% threshold_completeness: Minimum distinct single-copy marker genes required per bin
#const threshold_completeness=0.

% threshold_conta: Unused in affinity-aware model (kept for argument compatibility)
#const threshold_conta=4.

% conta_simthresh: Distance threshold distinguishing intra-genome from cross-genome co-markers
#const conta_simthresh=50.

% max_penalty: Normalization constant to invert composite distance into affinity reward
#const max_penalty=100. 

bin(1..k).

% Choice rule: Assign every contig to exactly one bin
1 { att(B,C) : bin(B) } 1 :- contig(C,_).

% 1. Minimum Bin Length
:- bin(B), #sum { L,C : att(B,C), contig(C,L) } < lengthmin.

% 1b. Maximum Bin Length: forbids degenerate mega-bin formation
:- bin(B), #sum { L,C : att(B,C), contig(C,L) } > lengthmax.

% 2. Completeness (Minimum number of unique SCGs)
scg_bin(B,G) :- mg(C,G), att(B,C).
:- bin(B), #count { G : scg_bin(B,G) } < threshold_completeness.

% 3. Contamination: Affinity-aware duplicate marker constraint
dissimilar(C1,C2) :- comarker_sim(C1,C2,D), D >= conta_simthresh.
:- att(B,C1), mg(C1,G),
   att(B,C2), C1 < C2, mg(C2,G),
   dissimilar(C1,C2).

% Optimization: Maximize cumulative composite sequence affinity
affinity(C1, C2, max_penalty - D) :- composite_knn(C1, C2, D).
:~ affinity(C1, C2, W), att(B, C1), att(B, C2). [-W@1, C1, C2]

#show att/2.