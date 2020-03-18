% This example uses the n-D minimization approach (function minlfr).
% It illustrates the problem due to the noncommutativity of the
% uncertain parameters.

   lfrs d1 d2 d3
   S1 = [d1*d2*d3;d2*d3*d1;d3*d1*d2];
   S2 = [d1*d2*d3;d1*d2*d3;d1*d2*d3];

% minlfr is used for reducing the size of the Delta matrices:

   S1min = minlfr(S1);
   S2min = minlfr(S2);

   size(S1min);
   size(S2min);

% In the first case Delta is 9-by-9 and in the second case 3-by-3.
% However we have "minimal objects" in both cases.

   distlfr(S1min,S2min)

% As a conclusion: commutativity must be taken into account
% before realization (i.e. using symtreed).
