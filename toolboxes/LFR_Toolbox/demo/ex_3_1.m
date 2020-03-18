% This example illustrates the use of the function "rf2lfr"
% that transforms the matrix NDv = [N;D] to sys = N*D^(-1)

% Realization of NDv

   lfrs Int x y z
   NDv = [1+x*Int x*y+z;0 x^2*Int^2;1+Int 0;0 2+x*y*Int];

% Then we apply the transformation "ndv2lfr"

   sys = rf2lfr(NDv);

% In order to check the results, we generate this object in
% an alternative way

   N = [eye(2,2) zeros(2,2)]*NDv;
   D = [zeros(2,2) eye(2,2)]*NDv;
   sys2 = N*D^(-1);

% Comparison

   size(sys)
   size(sys2)
   distlfr(sys,sys2)

% The function "rf2lfr" prevents augmenting the size of the
% Delta-matrix.
