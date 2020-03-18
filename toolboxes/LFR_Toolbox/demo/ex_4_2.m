% This example illustrates order reduction after realization (here,
% object-oriented realization) using the 1-D approach. It also shows
% that some factorizations that cannot be found during the realization
% phase might be found using the 1-D order reduction technique.

% We shall consider the right hand side member of:
% [d3*d1+d3*d2+d4*d1+d4*d2;d1+d2] = [d3+d4;1]*(d1+d2)
% (order 8 on the left, 4 on the right).

% Realization:

   lfrs d1 d2 d3 d4
   S = [d3*d1+d3*d2+d4*d1+d4*d2 ; d1+d2];

% Then, "minlfr1" is used for reducing the size of the Delta-block

   Smin = minlfr1(S);

% resulting in

   size(S)
   size(Smin)

% The order of Smin is four, that means that the factorization "was
% found".
