% This example illustrates the computation of approximation error
% bounds (function "min_max") and the use of this information for
% introducing a new uncertain parameter modelling the approximation.
% (See also ex_4_5.m).

% Let us consider an academic example in which simplifications
% can easily be guessed. Let

   lfrs d1 d2 d3 d4
   M0 = (3*d1^5+.0001*d1*d2*d3*d4)*(1-.0001*d4^4+d2^2*d3^2);
   M1 = (3*d1^5)*(1+d2^2*d3^2);

% M1 is considered as an approximation of M0. For computing the
% approximation error, we consider the difference

   DeltM = M1 - M0;
   DeltM = minlfr(DeltM);

% For computing the error bounds:

   [min_val,max_val,min_int,max_int] = min_max(DeltM);

% The results are min_val = -9.22e-04 and max_val =8.54e-04, which means
% that the approximation error varies in an interval included between these
% two values.

% It remains to build M3, and approximation of M0 better than M1. For
% that, we replace the neglected part of M1 by an uncertainty bounded
% by min_val and max_val

   lfrs err [min_val] [max_val]
   M3 = M1 + err;

   size(M3)
