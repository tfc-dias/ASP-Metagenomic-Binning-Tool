% Merge operator: consolidates over-split bins sharing high affinity

% simthresh: Co-marker distance threshold for dissimilarity
#const simthresh=50.

% minconn: Minimum crossing affinity edge count to qualify candidate bin pairs
#const minconn=50.

bin(B)    :- seed(B,_).
binc(B,C) :- seed(B,C).

% Cross-bin connectivity (symmetrized)
cross(B1,B2,N) :- binaff(B1,B2,N).
cross(B2,B1,N) :- binaff(B1,B2,N).

% Mutual nearest-neighbor bin selection
maxcross(B,M) :- bin(B), M = #max { N : cross(B,_,N) }, M > 0.
best(B,B2)    :- cross(B,B2,N), maxcross(B,N).
mutual(B1,B2) :- best(B1,B2), best(B2,B1), B1 < B2.

% Contamination guard: block merges with conflicting dissimilar markers
dissimilar(C1,C2) :- comarker_sim(C1,C2,D), D >= simthresh.
dissimilar(C2,C1) :- comarker_sim(C1,C2,D), D >= simthresh.

blocked(B1,B2) :- mutual(B1,B2),
                  binc(B1,C1), mg(C1,G),
                  binc(B2,C2), mg(C2,G),
                  dissimilar(C1,C2).

% Merge candidate bin pairs satisfying connectivity and guard
merge(B1,B2) :- mutual(B1,B2), cross(B1,B2,N), N >= minconn, not blocked(B1,B2).

% Connected components relabeling (union-find)
link(B1,B2) :- merge(B1,B2).
link(B2,B1) :- merge(B1,B2).
reach(B,B)  :- bin(B).
reach(B,B2) :- reach(B,B1), link(B1,B2).
rep(B,M)    :- bin(B), M = #min { B2 : reach(B,B2) }.

% Output attribution
att(R,C) :- binc(B,C), rep(B,R).
#show att/2.
