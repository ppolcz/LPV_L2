% Illustration of the function lfr/eval (see also ex_2_6.m)

% A,B,C defined as random LFR-objects

   A = rlfr(0,4,4,2,2,'a');
   B = rlfr(0,4,2,2,2,'b');
   C = rlfr(0,2,4,2,2,'c');
   D = zeros(2,2);

% Block structure corresponding to 4 states

  blk = struct('names',{{'1/s'}},'desc',[4;4;0;1;1;1;0;0;0;0;0]);

% LFR-object before evaluation

   sys = lfr(A,B,C,D,blk);
   size(sys)

% LFR-object after evaluation

   sys = eval(sys);
   size(sys)

% The uncertainties of the sub-matrices are now in the global
% matrix Delta

% Note that sys = abcd2lfr([A B;C D],4) leads to the same result,
% see ex_3_2.m

   sys2 = abcd2lfr([A B;C D],4);
   distlfr(sys,sys2)
