% The construction proposed in ex_2_11.m is treated now using the
% symbolic approach instead of the object-oriented approach.

   syms Int d1 d3
   syss = d1*[Int d1*d3]*[1 1;0 1]*[d1*Int ; d3];
   sys3 = sym2lfr(syss);

% Note that 1/s must be denoted Int in the used formula.

% Object-oriented approach (ex_2_12.m) for comparison

   lfrs Int d1 d3
   sys2 = d1*[Int d1*d3]*[1 1;0 1]*[d1*Int ; d3];

% Comparison

   distlfr(sys2,sys3)

% The distance should be zero.
