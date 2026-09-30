% Domain heuristic warm-start overlay for whole-partition binning
% Usage: clingo binning_warmstart.lp <facts...> <seed_facts.lp> --heuristic=Domain

#include "binning.lp".

% Prefer to set each seeded assignment true first (domain heuristic)
#heuristic att(B,C) : seed(B,C). [1, true]
