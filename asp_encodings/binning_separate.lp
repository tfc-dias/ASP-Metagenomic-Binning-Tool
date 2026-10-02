% Separation operator: splits chimeric seed bins

% simthresh: Co-marker distance threshold for dissimilarity
#const simthresh=50.   % placeholder: derived per dataset and passed with -c (scripts/derive_simthresh.py)

% sepmul: Multiplier scaling cannot-link penalty against anchor affinity
#const sepmul=1.

% Sub-bin budget per seed bin
binK(B,K) :- binmax(B,K).
binK(B,1) :- seed(B,_), not binmax(B,_).

contig(C) :- seed(_,C).
1 { sub(C,S) : S=1..K } 1 :- seed(B,C), binK(B,K).

% Symmetry breaking: pin first anchor to slot 1
:- pin1(C), sub(C,S), S>1.

% Cannot-link penalty for dissimilar contig pairs sharing a marker
clpen(C1,C2,P) :- comarker_sim(C1,C2,D), D >= simthresh, P = (D - simthresh) * sepmul.
sharemarker(C1,C2) :- seed(B,C1), seed(B,C2), mg(C1,G), mg(C2,G), C1 < C2.

% Optimization: balance sequence affinity against cannot-link penalty
:~ sub(C,S), sub(A,S), affinity(C,A,W). [-W@2, C, A]
:~ sub(C1,S), sub(C2,S), sharemarker(C1,C2), clpen(C1,C2,P). [P@2, C1, C2]

% Tie-breaking: prefer slot 1
:~ sub(C,S), S>1. [1@1, C]

% Output attribution
att(B*100+S, C) :- seed(B,C), sub(C,S).
#show att/2.
