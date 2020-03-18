% In this example we illustrate the use of the function "abcd2lfr"
% that computes the input/output object corresponding to a
% quadruple (A,B,C,D). See also ex_2_7.m.

% Realization of A, B, C and D:

   lfrs x y z
   A = [1+x x*y+z;0 x^2]; B = [1;y]; C=[y+z x*y]; D=1+y^2;
   S = [A B;C D];

   sys = abcd2lfr(S,2);

% The second argument "2" specifies the size of A in S.

% In order to check this result, sys is computed in an alternative
% way:

   I2 = eye(2,2);
   sys2 = C*feedback(I2*Int,A,1)*B + D;

% Comparison

   size(sys1)
   size(sys2)
   distlfr(sys1,sys2)

% Similar results



