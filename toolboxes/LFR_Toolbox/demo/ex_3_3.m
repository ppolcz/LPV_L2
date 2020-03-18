% This example illustrates the use of the functions "bnds2lfr"
% and "abcd2lfr". It is assumed that bounds relative to the
% elements of the matrices A,B of a are known (minA, maxA, minB
% and maxB). C and D are fixed.

% "bnds2lfr" realizes A and B as LFR-objects:

   minA = [-2 -2 -4 ; 0 -5 -5 ; 0  0 -6];
   maxA = [ 0 -2 -2 ; 0 -3 -5 ; 0  0 -6];
   disp([minA maxA])

   minB = [0.5 ; 0.5 ; 0.5];
   maxB = [1.5 ; 0.5 ; 1.5];
   disp([minB maxB])

   C=[1 1 1];
   D=1;

% The parameter names of lfrA will be "A_i_j"
% The parameter names of lfrB will be "B_i_j"

   lfrA = bnds2lfr('A_',minA,maxA);
   lfrB = bnds2lfr('B_',minB,maxB);

% Then S = [lfrA lfrB;C D] is transformed into an input/output
% LFR-object using "abcd2lfr" (3 states).

   sys = abcd2lfr([lfrA lfrB;C D],3);
   size(sys)

% Note that only parameters with non-empty interval of variation
% are created
