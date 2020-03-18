% This example illustrates the computation of approximation error
% bounds (function "udistlfr") and the use of this information for
% introducing new uncertain parameters modelling the approximations.
% (See also ex_4_4.m).

% Approximation is performed using the function "reduclfr"

% Realization:

   lfrs a b c
   S1 = [3+0.001*a^5-b*c a^4*c;a*b*c^3+2 a^2-0.001*b*c^3+1];

% Approximation using "reduclfr" (this function handles the tolerance
% argument of order reduction algorithms in such a way that approximation
% error remains less than a given bound).

   S2 = reduclfr(S1,0.01,'a');

% Evaluation of the approximation error

   [distu,dist2,mindiff,maxdiff] = udistlfr(S1,S2);
   disp(mindiff); disp(maxdiff);

% This function returns mindiff and maxdiff that are matrices having the
% same size as S1. These matrices give term by term
% - a lower bound of the minimum value of the approximation error (mindiff)
% - an upper bound of the maximum value of the approximation error (maxdiff).

% It remains to replace the approximation errors by additional normalized
% real parameters.

   syserr = bnds2lfr('s_',mindiff,maxdiff);
   S3 = S2 + syserr;

   size(S3)

% Now S3 has a Delta-block of size 14-by-14 to be compared to the original
% size of S1 that is 23-by-23.


