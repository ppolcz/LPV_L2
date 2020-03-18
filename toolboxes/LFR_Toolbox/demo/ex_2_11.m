% Object-oriented LFR-object realization, function lfrs.

% Int is a reserved name for 1/s

   lfrs Int d1 d3
   sys1 = d1^2*Int^2 + d1*d3*Int + d1^2*d3^2;

% Result

   size(sys1)

% ex_2_12.m and ex_2_13.m proposes an alternative generation
% of the same object.

% Other example

   lfrs Int
   lfrs a b 'real' [1 1] [3 3]
   lfrs c [1] [3] [2.5]
   lfrs d 'complex' [0] [2]

   size(a*b*c*d)

