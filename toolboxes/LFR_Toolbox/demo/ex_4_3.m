% This example illustrates order reduction after realization (here,
% object-oriented realization) using the n-D approach.
% This example aims also at illustrating the fundamental difference
% between the 1-D (function "minlfr1") and the n-D (function "minlfr")
% approaches.

% For that purpose, let us consider an expression in which factorization
% cannot be "found" considering only one parameter at a time.
% [1/(1+d1+d2);1/(1+d1+d2)] = [1;1]*(1/(1+d1+d2))

% Realization:

   lfrs d1 d2
   S = [1/(1+d1+d2);1/(1+d1+d2)];

% The minimum forms after using the 1-D and n-D approaches are respectively
% denoted Smin1 and Smin

   Smin1 = minlfr1(S);
   Smin  = minlfr(S);

% For comparison of the results:

   size(Smin1)
   size(Smin)

% As expected Smin is minimum but Smin1 has the same complexity
% as the original LFR-object.
