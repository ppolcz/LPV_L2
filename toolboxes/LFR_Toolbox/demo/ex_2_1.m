% Illustration of "lfr", "lfrdata", "flup"

% K_2 = I_2 (1 d2^2 + 2 d2 + 3 ), from the manual:

   D11 = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
   D12 = [0 0;0 0;1 0;0 1];
   D21 = [1 0 2 0;0 1 0 2];
   D22 = [3 0;0 3];

   blk = struct('names',{{'delta_2'}},'desc',[4;4;1;1;1;1;1;2;-1;1;0]);

   K2 = lfr(D11,D12,D21,D22,blk);

% Having such an LFR-object, it is possible to recover the matrices
% D11, D12, D21, D22 and blk by typing respectively:

   K2a = K2.a;
   K2b = K2.b;
   K2c = K2.c;
   K2d = K2.d;
   K2blk = K2.blk;

% For the same objective we can use lfrdata:

   [K2a,K2b,K2c,K2d,K2blk] = lfrdata(K2);
   size(K2)

% In order to conclude this example we can check the result by closing
% the Delta-loop. For that we shall invoke flup.

   valK2 = uplft(K2,{'delta_2'},10);
   valK2.d

% The result should be 123*eye(2,2)
