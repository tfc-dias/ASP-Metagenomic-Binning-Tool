% Recruitment operator: assigns unbinned contigs to established seed bins

% threshold_conta: Duplication tolerance for recruited markers (1 = strict guard)
#const threshold_conta=1.

% conta_simthresh: Distance threshold distinguishing intra-genome from cross-genome co-markers
#const conta_simthresh=50.   % placeholder: derived per dataset and passed with -c (scripts/derive_simthresh.py)

% minaff: Minimum realized affinity required to recruit a contig (0 = disabled)
#const minaff=0.

% Pinned seed bins remain fixed
bin(B) :- pinned(B,_).
att(B,C) :- pinned(B,C).

% Assign each unbinned contig to at most one candidate bin
{ att(B,C) : bin(B) } 1 :- recruitable(C).

% Placement requires positive affinity to target bin
:- att(B,C), recruitable(C), #count { Cb : affinity(C,Cb,_), pinned(B,Cb) } = 0.
:- att(B,C), recruitable(C), #sum { W,Cb : affinity(C,Cb,W), pinned(B,Cb) } < minaff.

% Contamination guard: recruits cannot introduce duplicate dissimilar markers
dissimilar(C1,C2) :- comarker_sim(C1,C2,D), D >= conta_simthresh.
:- att(B,C1), recruitable(C1), mg(C1,G),
   att(B,C2), C1 != C2, mg(C2,G),
   dissimilar(C1,C2).

% Optimization: maximize cumulative affinity of recruited sequences
:~ affinity(C,Cb,W), att(B,C), pinned(B,Cb). [-W@1, C, Cb]

#show att/2.
