% This example illustrates LFR-object realization by the tree
% decomposition, using the function "symtreed"

   syms d1 d2 d3 d4 d5
   sys_sym = [4*d1^2*d3 3*d1 0 ; d3*d5 5*d2^2*d4 d2*d4^2];

% Brute force realization

   sys1 = sym2lfr(sys_sym);
   size(sys1)

% Structured tree decomposition

   sys2 = symtreed(sys_sym);
   size(sys2)

% The size reduces from 12 to 9.

