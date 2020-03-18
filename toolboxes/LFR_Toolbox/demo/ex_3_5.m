% This example illustrates the use of Maple for low order realization
% of an LFR-object.

% First, the object sys_lfr is realized without order reduction (use
% of "sym2lfr").

   syms d1 d2 d3 Int
   sys_sym1 = d1^2*Int^2 + d1*d3*Int + d1^2*d3^2;

   sys_lfr1 = sym2lfr(sys_sym1);

% Use of the Maple function "convert" for order reduction

   sys_sym2 = maple('convert',sys_sym1,'horner',d1);
   sys_lfr2 = sym2lfr(sys_sym2);

% Resulting in:

   size(sys_lfr1)
   size(sys_lfr2)

% The reduction of the number of times d1 is repeated comes from the
% Horner factorization w.r.t. d1. In this simple example we can do
% manually the same factorization.

   sys_sym3 = (d3*Int+(Int^2+d3^2)*d1)*d1
   sys_lfr3 = sym2lfr(sys_sym3);

   size(sys_lfr3)

