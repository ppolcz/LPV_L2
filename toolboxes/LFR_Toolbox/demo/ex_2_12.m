% We solve the same problem as in ex_2_11.m using a different form
% of the considered uncertain object. This example illustrates the fact
% that the size of the block Delta is equal to the number of times
% the uncertain parameters appear in the used formula.

% Elementary LFR-objects are created

   lfrs Int d1 d3

   sys2 = d1*[Int d1*d3]*[1 1;0 1]*[d1*Int ; d3];

% Alternative form (of ex_2_11.m)

   sys1 = d1^2*Int^2 + d1*d3*Int + d1^2*d3^2;

% We can check the size

   size(sys1)
   size(sys2)

% The Delta matrix is now 5-by-5 (sys2) instead of 8-by-8 (sys1).

% Using the function distlfr we can check that we have computed
% equivalent objects:

   distlfr(sys1,sys2)

% The distance should be zero.
